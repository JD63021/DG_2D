function [K,F] = build_stiffness_operator_dgs(nodeInfo,elemInfo,boundaryInfo, ...
                                                 f_coeff, pCoef, theta, g_handle)
% BUILD_STIFFNESS_OPERATOR_DG_TRI (P1 triangles, DG, strong Dirichlet)
%   -Δu = f  in Ω,  u = g on ∂Ω  (Dirichlet imposed strongly at boundary vertices)
%
% Unknowns are element-owned P1 dofs (3 per triangle).
% Interior faces use standard IP-DG:
%     K_F =  Kv + theta*Kv' + penalty,
%   with Kv = -< {∇v·n}, [u] >.
% Boundary faces: no weak (Nitsche) terms are assembled; instead we pin
% all DG dofs whose physical vertex lies on ∂Ω (strong Dirichlet).
%
% IMPORTANT FIX (triangles):
%   Edge quadrature is made ORIENTATION-AWARE. For each interior edge A→B,
%   every adjacent element uses t_local = s * t,  s ∈ {+1,-1},
%   so that both sides evaluate shape functions at the SAME physical point.
%
% Inputs
%   nodeInfo.x,y        : Nnodes×1 coordinates
%   elemInfo.elements   : Ne×3 (CCW) connectivity
%   boundaryInfo.edges  : Nbe×2 (not used here; boundaryInfo.allBoundaryNodes is)
%   f_coeff             : constant RHS coefficient (f(x,y) = f_coeff)
%   pCoef               : penalty constant on interior faces (keep constant across refinements)
%   theta               : +1=SIPG, 0=IIPG, -1=NIPG
%   g_handle(x,y)       : boundary value (default @(x,y)0)
%
% Outputs
%   K : Ndof×Ndof sparse, Ndof = 3*Ne
%   F : Ndof×1
% -------------------------------------------------------------------------

xy   = [nodeInfo.x , nodeInfo.y];
T    = elemInfo.elements;        % Ne×3
Ne   = size(T,1);
Nloc = 3;                         % P1
Ndof = 3*Ne;

if nargin<4 || isempty(f_coeff), f_coeff = 0; end
if nargin<5 || isempty(pCoef),   pCoef   = 30; end
if nargin<6 || isempty(theta),   theta   = +1; end
if nargin<7 || isempty(g_handle), g_handle = @(x,y) 0; end

f_handle = @(x,y) f_coeff;       % constant RHS

% 1D 2-point Gauss for edges (param t ∈ [-1,1])
g1D   = 1/sqrt(3);
edgeT = [-g1D, g1D];
edgeW = [ 1,    1   ];

% Local edges (node pairs) in CCW order: (1-2), (2-3), (3-1)
edgePair = [1 2; 2 3; 3 1];

% ---------- Element geometry: areas and constant P1 gradients -------------
areaK  = zeros(Ne,1);
b_all  = zeros(Ne,3);   % y_j - y_k
c_all  = zeros(Ne,3);   % x_k - x_j
for e = 1:Ne
    ve   = T(e,:);               % [i j k]
    x1 = xy(ve(1),1); y1 = xy(ve(1),2);
    x2 = xy(ve(2),1); y2 = xy(ve(2),2);
    x3 = xy(ve(3),1); y3 = xy(ve(3),2);

    areaK(e) = 0.5 * abs( (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) );

    % cyclic (i,j,k) in CCW orientation:
    b_all(e,1) = y2 - y3;   c_all(e,1) = x3 - x2;   % for N1
    b_all(e,2) = y3 - y1;   c_all(e,2) = x1 - x3;   % for N2
    b_all(e,3) = y1 - y2;   c_all(e,3) = x2 - x1;   % for N3
end

% DG dof map: e,local(1..3) -> (e-1)*3 + local
dgDofs = reshape(1:Ndof, Nloc, Ne).';

% ---------------- Triplets & RHS ------------------------------------------
maxNNZ = 9*Ne + 36*(Ne*2);       % rough upper bound
I = zeros(maxNNZ,1); J = zeros(maxNNZ,1); S = zeros(maxNNZ,1);
ptr = 1;
F   = zeros(Ndof,1);

% ----------------------------- Volume ------------------------------------
for e = 1:Ne
    Ae = areaK(e);
    be = b_all(e,:);  ce = c_all(e,:);
    % constant grads in element
    dNdx = (be ./ (2*Ae)).';       % 3×1
    dNdy = (ce ./ (2*Ae)).';       % 3×1

    % Stiffness Ke_ij = (∇Ni · ∇Nj) * |K|
    Ke = (dNdx*dNdx.' + dNdy*dNdy.') * Ae;

    % Load: f constant ⇒ ∫_K N_i dΩ = |K|/3
    Fe = (Ae/3) * f_coeff * ones(3,1);

    idx = dgDofs(e,:);
    n   = numel(idx);
    rng = ptr:ptr+n^2-1; [I(rng),J(rng)] = meshgrid(idx,idx);
    S(rng) = Ke(:); ptr = ptr + n^2;
    F(idx) = F(idx) + Fe;
end

% --------------------------- Build unique edges ---------------------------
edgeMap = containers.Map('KeyType','char','ValueType','any');
edges = struct('nodes',{},'elemM',{},'locEdgeM',{},'elemP',{},'locEdgeP',{},'len',{},'n',{});
for e = 1:Ne
    ve   = T(e,:);
    ctri = mean(xy(ve,:),1);           % centroid (for normal orientation)
    for s = 1:3
        pair = edgePair(s,:);
        % global node ids in this element's LOCAL edge order:
        globLocal = ve(pair);
        % canonical (physical) edge key uses SORTED node ids [A B] with A<B:
        nn   = sort(globLocal);
        key  = sprintf('%d_%d', nn(1), nn(2));
        if ~isKey(edgeMap, key)
            A  = xy(nn(1),:); B = xy(nn(2),:);      % physical edge A->B
            tau = B - A; L = hypot(tau(1),tau(2));
            n = [tau(2), -tau(1)]/L;               % one unit normal to A->B
            mid = 0.5*(A+B);
            if dot(n, mid - ctri) < 0, n = -n; end % orient outward from this "minus" elem
            id = numel(edges)+1;
            edges(id).nodes    = nn;               % [A B] (A<B)
            edges(id).elemM    = e;     edges(id).locEdgeM = s;
            edges(id).elemP    = 0;     edges(id).locEdgeP = 0;
            edges(id).len      = L;     edges(id).n = n;
            edgeMap(key) = id;
        else
            id = edgeMap(key);
            edges(id).elemP    = e;     edges(id).locEdgeP = s;
        end
    end
end

% ------------------------- Interior edges only ----------------------------
for id = 1:numel(edges)
    if edges(id).elemP == 0
        continue;  % boundary edge: no weak terms (strong BC later)
    end

    ed = edges(id);
    L  = ed.len;        n = ed.n.';             % 2×1
    eM = ed.elemM;      eP = ed.elemP;

    % size for penalty (keep your original |K|/|F| scaling for now)
    hM = areaK(eM)/L;      hP = areaK(eP)/L;
    invh = max(1/hM, 1/hP);
    sigma = pCoef * invh;

    % Minus/plus side constant gradients
    AeM = areaK(eM);  bM = b_all(eM,:); cM = c_all(eM,:);
    dNdxM = (bM ./ (2*AeM)).';   dNdyM = (cM ./ (2*AeM)).';
    GM    = dNdxM*n(1) + dNdyM*n(2);            % 3×1

    AeP = areaK(eP);  bP = b_all(eP,:); cP = c_all(eP,:);
    dNdxP = (bP ./ (2*AeP)).';   dNdyP = (cP ./ (2*AeP)).';
    GP    = dNdxP*n(1) + dNdyP*n(2);            % 3×1

    % Determine ORIENTATION signs sM,sP for t on each local edge
    nnAB = ed.nodes; A = nnAB(1); B = nnAB(2);   % physical edge A->B (A<B)
    % local node ids (global) for each side's edge in its local order:
    locM = T(eM, edgePair(ed.locEdgeM,:));       % [g1 g2] as stored in element eM
    locP = T(eP, edgePair(ed.locEdgeP,:));       % [h1 h2] as stored in element eP
    % sign +1 if local is [A B], -1 if local is [B A]
    sM = +1;  if locM(1)==B && locM(2)==A, sM = -1; end
    sP = +1;  if locP(1)==B && locP(2)==A, sP = -1; end

    % Edge 1D Gauss: param t in [-1,1], geom factor = L/2
    for gq = 1:2
        t   = edgeT(gq);  wg = edgeW(gq);
        sfac= (L/2) * wg;

        % P1 shapes on each side with ORIENTATION-CORRECTED t:
        tM = sM * t;   tP = sP * t;

        lpM = edgePair(ed.locEdgeM,:);   NM = zeros(3,1);
        NM(lpM(1)) = 0.5*(1 - tM);
        NM(lpM(2)) = 0.5*(1 + tM);

        lpP = edgePair(ed.locEdgeP,:);   NP = zeros(3,1);
        NP(lpP(1)) = 0.5*(1 - tP);
        NP(lpP(2)) = 0.5*(1 + tP);

        % Kv = -< {∇v·n}, [u] >
        Kv = -0.5 * [ GM*NM.' , -GM*NP.' ;  GP*NM.' , -GP*NP.' ];

        % Full interior contribution
        Kedge = (Kv + theta * Kv.') ...
                + sigma * [ NM*NM.' , -NM*NP.' ; -NP*NM.' , NP*NP.' ];
        Kedge = sfac * Kedge;

        % scatter to global
        idx = [dgDofs(eM,:) , dgDofs(eP,:)];   m = numel(idx);
        rng = ptr:ptr+m^2-1; [I(rng),J(rng)] = meshgrid(idx,idx);
        S(rng) = Kedge(:);  ptr = ptr + m^2;
    end
end

% ------------------------------ finalize K,F ------------------------------
I = I(1:ptr-1);  J = J(1:ptr-1);  S = S(1:ptr-1);
K = sparse(I,J,S,Ndof,Ndof);

% ======================== STRONG DIRICHLET at boundary vertices ========================
bVerts = unique(boundaryInfo.allBoundaryNodes(:));
if ~isempty(bVerts)
    % DG dofs whose vertex is on the boundary:
    isB   = ismember(T, bVerts.');      % Ne×3 logical
    bDofs = dgDofs(isB);  bDofs = unique(bDofs);

    % g(x,y) at those dofs:
    e_of = ceil(bDofs / Nloc);
    a_of = bDofs - (e_of-1)*Nloc;
    vIDs = T( sub2ind(size(T), e_of, a_of) );
    gx   = nodeInfo.x(vIDs);  gy = nodeInfo.y(vIDs);
    gVal = arrayfun(g_handle, gx, gy);

    % strong clamp: row=0, diag=1, RHS=g
    K(bDofs,:) = 0;
    K(sub2ind([Ndof,Ndof], bDofs, bDofs)) = 1;
    F(bDofs) = gVal;
end

end
