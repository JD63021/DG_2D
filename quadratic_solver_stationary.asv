% quadratic_solver.m  —  DG P2 scalar transport: solve, plot, broken-L2 vs 1D exact
% ---------------------------------------------------------------------------------
tic;  clear;  clc;
warning off MATLAB:NonScalarInput
warning off MATLAB:nearlySingularMatrix

%% ---------------------- User parameters ----------------------------------
penalty  = 0;       % SIPG base penalty
theta    = -1;        % +1=SIPG, 0=IIPG, -1=NIPG
shrink   = 0.90;      % visual shrink for "disconnected" look
eps_diff = 1.0;       % diffusion epsilon (used in exact solution)
beta     = [-10, 0];  % convection field beta = (bx, by)
f_coeff  = 0.0;       % no volumetric source for this test

%% ---------------------- Mesh / case --------------------------------------
test = testCase();   % user-supplied
[nodeInfo, elemInfo, boundaryInfo] = mesh5_gmsh(test.gmshVelFile);

E    = elemInfo.elements;  Ne = size(E,1);  Nloc = size(E,2);
if Nloc ~= 6, error('This driver expects P2 triangles (Nloc=6).'); end
Ndof = Ne*Nloc;

%% ---------------------- Solve (stationary) -------------------------------
% Your main should forward `beta` into the assembler (as in previous turn):
%   main1_SS(nodeInfo, elemInfo, boundaryInfo, 1, penalty, theta, beta)
[U, x, y, nodeInfo, elemInfo, boundaryInfo] = ...
    main1_SS(nodeInfo, elemInfo, boundaryInfo, 1, penalty, theta, beta);

%% ---------------------- DG dof → coords ----------------------------------
dgDofs  = reshape(1:Ndof, Nloc, Ne).';
xDG     = zeros(Ndof,1);  yDG = zeros(Ndof,1);
for e = 1:Ne
    idx = dgDofs(e,:);  nodes = E(e,:);
    xDG(idx) = nodeInfo.x(nodes);  yDG(idx) = nodeInfo.y(nodes);
end

%% ---------------------- Plot (curved P2 patches) -------------------------
figure;
plotDG_P2_curved(U, nodeInfo, elemInfo, 'shrink', shrink, 'nSub', 8, 'drawOutlines', true);
colorbar; axis equal; view(3); grid on
Pe = abs(beta(1))/max(eps_diff,eps);  % |b_x|/eps
title(sprintf('DG P2 solution – one curved patch/element  (Pe=%.3g, \\beta=(%.1f,%.1f))', ...
              Pe, beta(1), beta(2)));

%% ---------------------- Broken L2 error vs 1D exact ----------------------
[L2err, L2sq] = broken_L2_error_DG_P2_adv(U, nodeInfo, elemInfo, -beta(1), eps_diff);

% characteristic h from mean triangle area (using the 3 vertex nodes)
xy = [nodeInfo.x, nodeInfo.y];
areas = zeros(Ne,1);
for e = 1:Ne
    v = E(e,1:3);
    areas(e) = polyarea(xy(v,1), xy(v,2));
end
h = sqrt(mean(areas));

fprintf('Pe = %.6f   (|beta_x|/eps)\n', Pe);
fprintf('log10(h) = %.6f,   log10(L2err) = %.6f   (P2-tri DG)\n', ...
        log10(h), log10(L2err));
toc

% =============================== HELPERS ==================================

function [L2err, L2sq] = broken_L2_error_DG_P2_adv(U, nodeInfo, elemInfo, bx, eps_diff)
% Broken L2 for P2 triangles vs the 1-D steady advection–diffusion exact u(x):
%   -eps u'' + b u' = 0 on x∈[0,1], u(0)=1, u(1)=0 (independent of y).
T   = elemInfo.elements;              % Ne×6  [v1 v2 v3 m12 m23 m31]
Ne  = size(T,1);
Ndof= 6*Ne;
assert(numel(U)==Ndof, 'U size does not match 6*Ne for P2 DG.');

% 6-pt rule on reference triangle + P2 shapes/derivs
[Ntab, dNxi_tab, dNeta_tab, g_wt, ~] = precomputeShapeFunctionsP2_Tri();
nG = numel(g_wt);

dgDofs = reshape(1:Ndof, 6, Ne).';
L2sq = 0.0;  X = nodeInfo.x(:);  Y = nodeInfo.y(:);

b_over_eps = bx / eps_diff;
use_linear = (abs(b_over_eps) < 1e-12);
if ~use_linear, eB = exp(b_over_eps); end

for e = 1:Ne
    en = T(e,:);                 % 6 nodes
    xe = X(en);  ye = Y(en);
    Ue = U(dgDofs(e,:));         % 6×1

    for k = 1:nG
        N    = Ntab(:,k);        % 6×1
        dNxi = dNxi_tab(:,k);    % 6×1
        dNeta= dNeta_tab(:,k);   % 6×1

        % Correct P2 geometric Jacobian (this is the line that was wrong before)
        [~, ~, detJ] = q2ShapeDerivatives_AllNodes(xe, ye, dNxi, dNeta);

        xq = N.'*xe;             % physical x
        uh = N.'*Ue;

        % exact 1D solution at xq
        if use_linear
            ue = 1 - xq;
        else
            % works for both signs of b: returns u(0)=1, u(1)=0
            ue = (eB - exp(b_over_eps * xq)) / (eB - 1.0);
        end

        L2sq = L2sq + (uh - ue)^2 * detJ * g_wt(k);
    end
end
L2err = sqrt(L2sq);
end

% -------- Plotter: curved P2 patches without interior lines ---------------
function plotDG_P2_curved(U, nodeInfo, elemInfo, varargin)
p = inputParser;
p.addParameter('shrink', 0.95, @(v)isnumeric(v)&&isscalar(v)&&v>0&&v<=1);
p.addParameter('nSub',   8,    @(v)isnumeric(v)&&isscalar(v)&&v>=3);
p.addParameter('drawOutlines', true, @(v)islogical(v)&&isscalar(v));
p.parse(varargin{:});
shrink  = p.Results.shrink;
nSub    = round(p.Results.nSub);
doEdge  = p.Results.drawOutlines;

X  = nodeInfo.x(:); Y = nodeInfo.y(:);
T6 = elemInfo.elements;                 % Ne x 6
Ne = size(T6,1);
Ndof = 6*Ne; 
assert(numel(U)==Ndof, 'U must be size 6*Ne for P2 DG');

[bc, TRIloc] = baryGridAndTri(nSub);     % barycentric lattice & tris
Nb  = size(bc,1);  Nt = size(TRIloc,1);

% P2 shapes at bary points
Nsh = zeros(6, Nb);
for q = 1:Nb
    xi  = bc(q,2);   eta = bc(q,3);
    Nsh(:,q) = p2basisGmsh_ref(xi, eta);
end

totPts = Ne * Nb;   totTri = Ne * Nt;
XYZ  = zeros(totPts, 3);
TRI  = zeros(totTri, 3);
ptrP = 0;  ptrT = 0;

for e = 1:Ne
    en = T6(e,:);  xe = X(en);  ye = Y(en);
    idxDG = (e-1)*6 + (1:6);  Ue = U(idxDG);

    cx = mean(xe(1:3));  cy = mean(ye(1:3));

    Pel  = (Nsh.' * [xe(:), ye(:)]);  % Nb x 2
    Uhel = (Nsh.' * Ue(:));           % Nb x 1

    if shrink ~= 1
        Pel(:,1) = cx + shrink*(Pel(:,1) - cx);
        Pel(:,2) = cy + shrink*(Pel(:,2) - cy);
    end

    idxP = ptrP + (1:Nb);
    XYZ(idxP,1:2) = Pel;  XYZ(idxP,3) = Uhel;

    idxT = ptrT + (1:Nt);
    TRI(idxT,:) = TRIloc + ptrP;

    ptrP = ptrP + Nb;  ptrT = ptrT + Nt;
end

trisurf(TRI, XYZ(:,1), XYZ(:,2), XYZ(:,3), ...
        'FaceColor','interp','EdgeColor','none','FaceAlpha',1.0);
hold on;

if doEdge
    tEdge = linspace(0,1, max(20,2*nSub+1));
    Eset = [1 2 4; 2 3 5; 3 1 6];  % [vA vB mAB]
    for e = 1:Ne
        en = T6(e,:);  xe = X(en); ye = Y(en);
        for s = 1:3
            A = Eset(s,1); B = Eset(s,2); M = Eset(s,3);
            pA = [xe(A) ye(A)]; pB=[xe(B) ye(B)]; pM=[xe(M) ye(M)];
            L1 = (1-tEdge).*(1-2*tEdge);  L2 = 4*tEdge.*(1-tEdge);  L3 = tEdge.*(2*tEdge-1);
            xE = L1*pA(1) + L2*pM(1) + L3*pB(1);
            yE = L1*pA(2) + L2*pM(2) + L3*pB(2);
            if shrink ~= 1
                cx = mean(xe(1:3)); cy = mean(ye(1:3));
                xE = cx + shrink*(xE - cx);  yE = cy + shrink*(yE - cy);
            end
            plot3(xE, yE, zeros(size(xE)), 'Color',[0.4 0.4 0.4], 'LineWidth',0.8);
        end
    end
end
xlabel('x'); ylabel('y'); zlabel('U'); axis equal tight
hold off
end

% -------- P2 basis (Gmsh order) & barycentric grid triangulation ----------
function N = p2basisGmsh_ref(xi, eta)
zeta = 1 - xi - eta;
N1 = zeta.*(2*zeta - 1);
N2 = xi  .*(2*xi   - 1);
N3 = eta .*(2*eta  - 1);
N4 = 4*xi .* zeta;
N5 = 4*xi .* eta;
N6 = 4*eta.* zeta;
N  = [N1; N2; N3; N4; N5; N6];
end

function [bc, TRI] = baryGridAndTri(nSub)
bc   = [];  TRI = [];
ind  = nan(nSub+1, nSub+1);
ptr  = 0;
for i = 0:nSub
    for j = 0:(nSub - i)
        k = nSub - i - j;
        ptr = ptr + 1;
        ind(i+1, j+1) = ptr;
        bc(ptr,:) = [k, i, j] / nSub;   % λ1,λ2,λ3
    end
end
for i = 0:(nSub-1)
    for j = 0:(nSub-1 - i)
        a = ind(i+1,   j+1);
        b = ind(i+2,   j+1);
        c = ind(i+1, j+2);
        d = ind(i+2, j+2);
        TRI = [TRI; a b c; b d c];  %#ok<AGROW>
    end
end
end
