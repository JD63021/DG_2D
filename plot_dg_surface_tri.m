function plot_dg_surface_tri(nodeInfo, elemInfo, U, varargin)
% plot_dg_surface_tri
%   Renders a DG triangular solution as a smooth surface with only
%   element boundaries drawn (no interior sub-triangle lines).
%
%   Auto-detects P2 (6 nodes/tri) vs P1 (3 nodes/tri) from elemInfo.elements.
%
% Inputs
%   nodeInfo.x, nodeInfo.y : Nnodes×1 geometry nodes (from Gmsh)
%   elemInfo.elements      : Ne×6 (P2) or Ne×3 (P1) connectivity
%   U                      : Ndof × 1 DG vector, Ndof = Nloc*Ne
%
% Name–Value options (optional)
%   'EdgeColor' : [r g b] for element edges      (default [0 0 0])
%   'EdgeWidth' : scalar line width              (default 1.2)
%   'Colormap'  : e.g. 'turbo','parula'          (default 'turbo')
%   'Shrink'    : 0<shrink≤1 (shrink toward ctr) (default 1)
%   'View'      : 'iso' or 'top'                 (default 'iso')
%   'Subdiv'    : positive integer refinement    (default 8 for P2, 1 for P1)

% ---------- options ----------
p = inputParser;
addParameter(p,'EdgeColor',[0 0 0]);
addParameter(p,'EdgeWidth',1.2);
addParameter(p,'Colormap','turbo');
addParameter(p,'Shrink',1.0);
addParameter(p,'View','iso');
addParameter(p,'Subdiv',[]);
parse(p,varargin{:});
ec   = p.Results.EdgeColor;
ew   = p.Results.EdgeWidth;
cmap = p.Results.Colormap;
shrk = p.Results.Shrink;
viewMode = lower(p.Results.View);
subUser  = p.Results.Subdiv;

% ---------- inputs ----------
x = nodeInfo.x;  y = nodeInfo.y;
T = elemInfo.elements;           % Ne×(3 or 6)
[Ne, Nloc] = size(T);
Ndof = numel(U);
assert(Ndof == Nloc*Ne, 'U size mismatch: expected %d = %d*%d', Nloc*Ne, Nloc, Ne);

isP2 = (Nloc == 6);
if isempty(subUser)
    sub = isP2 * 8 + (~isP2) * 1;  % default: 8 for P2, 1 for P1
else
    sub = max(1, round(subUser));
end

% ---------- local reference shapes ----------
if isP2
    % P2 triangle (Gmsh node order):
    % N = [N1..N6], with ζ=1-ξ-η
    Nfun = @(xi,eta) [ ...
       (1 - xi - eta).*(2*(1 - xi - eta) - 1); ... % N1 corner
       xi.*(2*xi - 1); ...                          % N2 corner
       eta.*(2*eta - 1); ...                        % N3 corner
       4*xi.*(1 - xi - eta); ...                    % N4 edge12
       4*xi.*eta; ...                                % N5 edge23
       4*eta.*(1 - xi - eta)];                      % N6 edge31
else
    % P1 triangle: N = [λ1; λ2; λ3] with ξ=λ2, η=λ3, ζ=1-ξ-η=λ1
    Nfun = @(xi,eta) [ 1 - xi - eta; xi; eta ];
end

% ---------- build a refined reference triangulation once ----------
% Subdivide the reference triangle with a uniform barycentric grid:
% points: (i/sub, j/sub), i>=0, j>=0, i+j<=sub
pts = [];
for i=0:sub
    for j=0:(sub-i)
        pts(end+1,:) = [i/sub, j/sub]; %#ok<AGROW>
    end
end
% connectivity on the reference grid
TR = [];
idxMap = -ones(sub+1, sub+1);
k = 1;
for i=0:sub
    for j=0:(sub-i)
        idxMap(i+1,j+1) = k; k=k+1;
    end
end
for i=0:sub-1
    for j=0:(sub-1-i)
        a = idxMap(i+1,   j+1);
        b = idxMap(i+2,   j+1);
        c = idxMap(i+1, j+2);
        d = idxMap(i+2, j+2);
        TR(end+1,:) = [a b c]; %#ok<AGROW>
        if j < sub-1-i
            TR(end+1,:) = [b d c]; %#ok<AGROW>
        end
    end
end
nLocV = size(pts,1);
nLocT = size(TR,1);

% Precompute local shapes at ref vertices (nLocV×Nloc)
NlocMat = Nfun(pts(:,1), pts(:,2)).';  % (Nloc × nLocV) transpose later
% (we want nLocV×Nloc to right-multiply by Xe(=Nloc×1))
NlocMat = NlocMat.';                    % nLocV × Nloc

% ---------- assemble global refined surface (no interior lines) ----------
X = zeros(Ne*nLocV,1);
Y = zeros(Ne*nLocV,1);
Z = zeros(Ne*nLocV,1);
TT= zeros(Ne*nLocT,3);

dgDofs = reshape(1:Ndof, Nloc, Ne).';

for e=1:Ne
    ve  = T(e,:);                   % local geometry nodes
    Xe  = x(ve);  Ye = y(ve);       % geometry coords (Nloc×1)
    Ue  = U(dgDofs(e,:));           % DG dofs (Nloc×1) in local order

    % centroid for optional shrink (in XY)
    cx = mean(Xe); cy = mean(Ye);

    % evaluate at refined reference vertices
    locX = NlocMat * Xe;            % nLocV×1   (fixed)
    locY = NlocMat * Ye;            % nLocV×1
    locZ = NlocMat * Ue;            % nLocV×1

    if shrk ~= 1
        locX = cx + shrk*(locX - cx);
        locY = cy + shrk*(locY - cy);
        % do not shrink values
    end

    base = (e-1)*nLocV;
    X(base+(1:nLocV)) = locX;
    Y(base+(1:nLocV)) = locY;
    Z(base+(1:nLocV)) = locZ;

    TT((e-1)*nLocT+(1:nLocT),:) = TR + base;
end

figure;
hs = trisurf(TT, X, Y, Z);
set(hs,'FaceColor','interp','EdgeColor','none'); % hide all interior lines
colormap(cmap); colorbar; grid on; axis equal; axis tight;
shading interp; lighting gouraud; camlight headlight; material dull;

% ---------- overlay element boundary edges (curved for P2, straight for P1) ----------
hold on;
zLift = 1e-8 * (max(Z)-min(Z) + eps);  % tiny lift to avoid z-fighting

if isP2
    % P2 edges are (1--4--2), (2--5--3), (3--6--1) in local order:
    edgeTrip = [1 4 2; 2 5 3; 3 6 1];
    nSamp = max(10, 2*sub) + 1;
    t = linspace(0,1,nSamp).';           % param along edge
    % P2 edge restriction (quadratic along the edge)
    Ni = (1-t).*(1-2*t);
    Nm = 4*t.*(1-t);
    Nj = t.*(2*t-1);

    for e=1:Ne
        ve  = T(e,:);
        Xe  = x(ve);  Ye = y(ve);  Ue = U(dgDofs(e,:));
        cx = mean(Xe); cy = mean(Ye);

        for s=1:3
            loc = edgeTrip(s,:);
            Xedge = Ni*Xe(loc(1)) + Nm*Xe(loc(2)) + Nj*Xe(loc(3));
            Yedge = Ni*Ye(loc(1)) + Nm*Ye(loc(2)) + Nj*Ye(loc(3));
            Zedge = Ni*Ue(loc(1)) + Nm*Ue(loc(2)) + Nj*Ue(loc(3));
            if shrk ~= 1
                Xedge = cx + shrk*(Xedge - cx);
                Yedge = cy + shrk*(Yedge - cy);
            end
            plot3(Xedge, Yedge, Zedge + zLift, '-', 'Color', ec, 'LineWidth', ew);
        end
    end
else
    % P1: straight edges between the 3 local vertices
    edgePair = [1 2; 2 3; 3 1];
    for e=1:Ne
        ve = T(e,:);
        Xe = x(ve); Ye = y(ve); Ue = U(dgDofs(e,:));
        cx = mean(Xe); cy = mean(Ye);
        for s=1:3
            a = edgePair(s,1); b = edgePair(s,2);
            Xedge = [Xe(a); Xe(b)];
            Yedge = [Ye(a); Ye(b)];
            Zedge = [Ue(a); Ue(b)];
            if shrk ~= 1
                Xedge = cx + shrk*(Xedge - cx);
                Yedge = cy + shrk*(Yedge - cy);
            end
            plot3(Xedge, Yedge, Zedge + zLift, '-', 'Color', ec, 'LineWidth', ew);
        end
    end
end

% ---------- view ----------
switch viewMode
    case 'top'
        view(2);
    otherwise
        view(-35,25);
end
xlabel('x'); ylabel('y'); zlabel('U');
title(sprintf('DG %s surface (no interior lines)', ternary(isP2,'P2','P1')));
set(gcf,'Renderer','opengl'); % often improves line/surface compositing
end

function s = ternary(cond,a,b), if cond, s=a; else, s=b; end, end
