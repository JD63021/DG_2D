function [K,F] = build_stiffness_operator_dgcd(nodeInfo,elemInfo,boundaryInfo, ...
                                              f_coeff, pCoef, theta, g_handle, beta, ...
                                              inletVal, outletVal)
% DG(IP) for  -Δu + β·∇u = f  on P2 triangles.
% Boundary flags:
%   9 : wall (no diffusion/convective flux)
%  10 : LEFT  — enforce u = inletVal **strongly** on all edge dofs
%  11 : RIGHT — Dirichlet u = outletVal via Nitsche + upwind convection
%
% Inputs:
%   beta      : [bx,by] or @(x,y)->[bx,by] (default [0,0])
%   inletVal  : scalar for flag 10 (default 1)
%   outletVal : scalar for flag 11 (default 0)
%
% Notes:
%   • Diffusion coefficient is 1 here. If you need ε≠1, multiply the
%     diffusion volume, interior-face SIPG, and boundary Nitsche blocks by ε.

xy   = [nodeInfo.x , nodeInfo.y];
T    = elemInfo.elements;        % Ne×6 (P2, CCW)
Ne   = size(T,1);
Nloc = 6;                         % P2
Ndof = Nloc*Ne;

% defaults
if nargin<4 || isempty(f_coeff), f_coeff = 0; end
if nargin<5 || isempty(pCoef),   pCoef   = 2; end
if nargin<6 || isempty(theta),   theta   = +1; end  % SIPG
if nargin<7 || isempty(g_handle), g_handle = @(x,y) 0; end
if nargin<8 || isempty(beta),     beta     = [0,0]; end
if nargin<9 || isempty(inletVal), inletVal = 1; end      % flag 10
if nargin<10|| isempty(outletVal),outletVal= 0; end      % flag 11

% beta as function handle
if isa(beta,'function_handle')
    beta_fun = beta;
else
    bconst = beta(:).';
    beta_fun = @(x,y) bconst;
end

pDeg = 2;                          % P2
Cdeg = (pDeg+1)^2;
f_handle = @(x,y) f_coeff;

% --- volume quadrature + P2 tables (you already have these helpers) ---
[Nvol, dNxi_vol, dNeta_vol, g_wt, g_pt] = precomputeShapeFunctionsP2_Tri();
nG = numel(g_wt);

% --- 1D edge Gauss (3-pt) ---
g1 = sqrt(3/5);
edgeT = [-g1, 0, g1];
edgeW = [ 5/9, 8/9, 5/9 ];

% local vertex pairs per edge
edgePair = [1 2; 2 3; 3 1];

% element areas (vertex-only) for hF
areaK = zeros(Ne,1);
for e = 1:Ne
    ve = T(e,1:3);
    x1=xy(ve(1),1); y1=xy(ve(1),2);
    x2=xy(ve(2),1); y2=xy(ve(2),2);
    x3=xy(ve(3),1); y3=xy(ve(3),2);
    areaK(e) = 0.5 * abs( (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) );
end

% -------- map boundary edge (n1,n2) -> flag from boundaryInfo -----------
edgeFlagMap = containers.Map('KeyType','char','ValueType','int32');
if isfield(boundaryInfo,'edges') && ~isempty(boundaryInfo.edges)
    E = boundaryInfo.edges;
    Ff= boundaryInfo.edgeFlags(:);
    for k = 1:size(E,1)
        a = E(k,1); b = E(k,2);
        key = sprintf('%d_%d', min(a,b), max(a,b));
        edgeFlagMap(key) = Ff(k);
    end
else
    warning('boundaryInfo.edges/edgeFlags missing; flags unknown.');
end

% DG dof map
dgDofs = reshape(1:Ndof, Nloc, Ne).';

% Triplets + RHS
maxNNZ = Ne * 2200;
I = zeros(maxNNZ,1); J = zeros(maxNNZ,1); S = zeros(maxNNZ,1);
ptr = 1;
F   = zeros(Ndof,1);

% ---------- collectors for strong Dirichlet on flag 10 ----------
strongRows = [];   % global dofs on flag 10
strongVals = [];   % their values (all = inletVal)

% ============================== VOLUME ====================================
for e = 1:Ne
    nodes = T(e,:);
    xe = xy(nodes,1);  ye = xy(nodes,2);
    Ke = zeros(Nloc,Nloc);
    Fe = zeros(Nloc,1);

    for q = 1:nG
        Nv    = Nvol(:,q);
        dNxi  = dNxi_vol(:,q);
        dNeta = dNeta_vol(:,q);

        [dNx, dNy, detJ] = q2ShapeDerivatives_AllNodes(xe, ye, dNxi, dNeta);

        % diffusion volume:  ∫ ∇φ_i·∇φ_j
        Ke = Ke + (dNx*dNx.' + dNy*dNy.') * (detJ * g_wt(q));

        % convection volume (strong form):  -∫ u (β·∇v)
        xq = Nv.'*xe;  yq = Nv.'*ye;
        b  = beta_fun(xq,yq);  bx = b(1);  by = b(2);
        Kconv_vol = -(bx*dNx + by*dNy) * (Nv.') * (detJ * g_wt(q));
        Ke = Ke + Kconv_vol;

        % RHS
        Fe = Fe + Nv * f_handle(xq,yq) * (detJ * g_wt(q));
    end

    idx = dgDofs(e,:); n = numel(idx);
    rng = ptr:ptr+n^2-1; [I(rng),J(rng)] = meshgrid(idx,idx);
    S(rng) = Ke(:); ptr = ptr + n^2;
    F(idx) = F(idx) + Fe;
end

% ======================== UNIQUE EDGE LIST ================================
edgeMap = containers.Map('KeyType','char','ValueType','any');
edges = struct('nodes',{},'elemM',{},'locEdgeM',{},'elemP',{},'locEdgeP',{},'len',{},'n',{});
for e = 1:Ne
    ve3 = T(e,1:3);
    ctri = mean(xy(ve3,:),1);
    for s = 1:3
        pairLoc  = edgePair(s,:);
        globPair = ve3(pairLoc);
        nn = sort(globPair);
        key = sprintf('%d_%d', nn(1), nn(2));
        if ~isKey(edgeMap, key)
            A  = xy(nn(1),:); B = xy(nn(2),:);
            tau = B - A; L = hypot(tau(1),tau(2));
            n = [tau(2), -tau(1)]/L;
            mid = 0.5*(A+B);
            if dot(n, mid - ctri) < 0, n = -n; end
            id = numel(edges)+1;
            edges(id).nodes = nn;
            edges(id).elemM = e;              edges(id).locEdgeM = s;
            edges(id).elemP = 0;              edges(id).locEdgeP = 0;
            edges(id).len   = L;              edges(id).n = n;
            edgeMap(key) = id;
        else
            id = edgeMap(key);
            edges(id).elemP = e;              edges(id).locEdgeP = s;
        end
    end
end

% ============================= EDGE TERMS =================================
for id = 1:numel(edges)
    ed = edges(id);
    L  = ed.len;     n = ed.n.';           % 2×1 column
    eM = ed.elemM;

    % orientation wrt physical A→B
    nnAB = ed.nodes;  A = nnAB(1);  B = nnAB(2);
    locM_v = T(eM, edgePair(ed.locEdgeM,:));
    sM = +1;  if locM_v(1)==B && locM_v(2)==A, sM = -1; end
    nodesM = T(eM,:);  xeM = xy(nodesM,1); yeM = xy(nodesM,2);

    if ed.elemP ~= 0
        % ---------------- INTERIOR EDGE ----------------
        eP = ed.elemP;
        locP_v = T(eP, edgePair(ed.locEdgeP,:));
        sP = +1;  if locP_v(1)==B && locP_v(2)==A, sP = -1; end
        nodesP = T(eP,:);  xeP = xy(nodesP,1); yeP = xy(nodesP,2);

        hM = 2*areaK(eM)/L;  hP = 2*areaK(eP)/L;
        invh = max(1/hM, 1/hP);
        sigma = pCoef * Cdeg * invh;

        for gq = 1:numel(edgeT)
            t  = edgeT(gq);   wg = edgeW(gq);
            sf = (L/2) * wg;

            % minus side eval (M)
            tM = sM * t;
            iM = edgePair(ed.locEdgeM,1);
            jM = edgePair(ed.locEdgeM,2);
            cM = setdiff(1:3, [iM jM]);
            lamM = zeros(3,1);
            lamM(iM) = 0.5*(1 - tM);
            lamM(jM) = 0.5*(1 + tM);
            lamM(cM) = 0;
            xiM  = lamM(2);  etaM = lamM(3);
            [NM, dNxiM, dNetaM] = p2basis_ref(xiM, etaM);
            [dNxM, dNyM, ~] = q2ShapeDerivatives_AllNodes(xeM, yeM, dNxiM, dNetaM);
            GM = dNxM*n(1) + dNyM*n(2);
            xqM = NM.'*xeM;  yqM = NM.'*yeM;

            % plus side eval (P)
            tP = sP * t;
            iP = edgePair(ed.locEdgeP,1);
            jP = edgePair(ed.locEdgeP,2);
            cP = setdiff(1:3, [iP jP]);
            lamP = zeros(3,1);
            lamP(iP) = 0.5*(1 - tP);
            lamP(jP) = 0.5*(1 + tP);
            lamP(cP) = 0;
            xiP  = lamP(2);  etaP = lamP(3);
            [NP, dNxiP, dNetaP] = p2basis_ref(xiP, etaP);
            [dNxP, dNyP, ~] = q2ShapeDerivatives_AllNodes(xeP, yeP, dNxiP, dNetaP);
            GP = dNxP*n(1) + dNyP*n(2);
            xqP = NP.'*xeP;  yqP = NP.'*yeP;

            % diffusion (SIPG)
            Kv = -0.5 * [ GM*NM.' , -GM*NP.' ;  GP*NM.' , -GP*NP.' ];
            Kedge = (Kv + theta*Kv.') ...
                    + sigma * [ NM*NM.' , -NM*NP.' ; -NP*NM.' , NP*NP.' ];

            % convection (upwind interior)
            bxby = beta_fun(0.5*(xqM+xqP), 0.5*(yqM+yqP));
            bn = bxby(1)*n(1) + bxby(2)*n(2);
            bn_plus  = max(bn,0);
            bn_minus = min(bn,0);

            Kedge_conv = ...
                bn_plus  * [ NM*NM.' ,  zeros(Nloc) ; -NP*NM.' , zeros(Nloc) ] + ...
                bn_minus * [ zeros(Nloc) ,  NM*NP.' ; zeros(Nloc) , -NP*NP.' ];

            Kedge = sf * (Kedge + Kedge_conv);

            idx = [dgDofs(eM,:) , dgDofs(eP,:)];   m = numel(idx);
            rng = ptr:ptr+m^2-1; [I(rng),J(rng)] = meshgrid(idx,idx);
            S(rng) = Kedge(:);  ptr = ptr + m^2;
        end

    else
        % ---------------- BOUNDARY EDGE ----------------
        key = sprintf('%d_%d', min(A,B), max(A,B));
        if isKey(edgeFlagMap, key)
            flag = edgeFlagMap(key);
        else
            flag = -999;
        end

        hM = 2*areaK(eM)/L;
        sigma = pCoef * Cdeg * (1/hM);

        % ---- STRONG DIRICHLET on flag 10: collect dofs and skip assembly
        if flag == 10
            locEdge  = ed.locEdgeM;               % 1..3
            locNodes = [ edgePair(locEdge,:) , 3+locEdge ];  % [i j mid]
            strongRows = [ strongRows ; dgDofs(eM,locNodes).' ];
            strongVals = [ strongVals ; inletVal*ones(numel(locNodes),1) ];
            continue
        end

        % For other flags, proceed with boundary quadrature
        for gq = 1:numel(edgeT)
            t  = edgeT(gq);   wg = edgeW(gq);
            sf = (L/2) * wg;

            tM = sM * t;
            iM = edgePair(ed.locEdgeM,1);
            jM = edgePair(ed.locEdgeM,2);
            cM = setdiff(1:3, [iM jM]);
            lamM = zeros(3,1);
            lamM(iM) = 0.5*(1 - tM);
            lamM(jM) = 0.5*(1 + tM);
            lamM(cM) = 0;

            xiM  = lamM(2);  etaM = lamM(3);
            [NM, dNxiM, dNetaM] = p2basis_ref(xiM, etaM);
            [dNxM, dNyM, ~] = q2ShapeDerivatives_AllNodes(xeM, yeM, dNxiM, dNetaM);
            GM = dNxM*n(1) + dNyM*n(2);

            xq = NM.'*xeM;  yq = NM.'*yeM;   gqVal = g_handle(xq,yq);
            bxy = beta_fun(xq,yq);  bn = bxy(1)*n(1) + bxy(2)*n(2);

            switch flag
                case 9
                    % wall: no diffusive or convective flux -> nothing
                    Kb = zeros(Nloc); Fb = zeros(Nloc,1);
                case 11
                    % RIGHT Dirichlet = outletVal via Nitsche
                    gqVal = outletVal;
                    Kb = ( - (GM*NM.')  - theta*(NM*GM.')  + sigma*(NM*NM.') );
                    Fb = ( - GM*gqVal   + sigma*NM*gqVal );
                    % Convection upwind at boundary
                    if bn >= 0
                        % outflow: natural
                        Kb_conv = bn * (NM*NM.');  Fb_conv = 0;
                    else
                        % inflow: inject boundary value (0 usually)
                        Kb_conv = zeros(Nloc);     Fb_conv = bn * (NM * gqVal);
                    end
                    Kb = sf * (Kb + Kb_conv);  Fb = sf * (Fb + Fb_conv);
                otherwise
                    % fallback: treat as Dirichlet with g_handle
                    Kb = ( - (GM*NM.')  - theta*(NM*GM.')  + sigma*(NM*NM.') );
                    Fb = ( - GM*gqVal   + sigma*NM*gqVal );
                    if bn >= 0
                        Kb_conv = bn * (NM*NM.');  Fb_conv = 0;
                    else
                        Kb_conv = zeros(Nloc);     Fb_conv = bn * (NM * gqVal);
                    end
                    Kb = sf * (Kb + Kb_conv);  Fb = sf * (Fb + Fb_conv);
            end

            if any(Kb(:)) || any(Fb)
                idx = dgDofs(eM,:);  m = numel(idx);
                rng = ptr:ptr+m^2-1; [I(rng),J(rng)] = meshgrid(idx,idx);
                S(rng) = Kb(:);  ptr = ptr + m^2;
                F(idx) = F(idx) + Fb;
            end
        end
    end
end

% ------------------------------ finalize ----------------------------------
I = I(1:ptr-1);  J = J(1:ptr-1);  S = S(1:ptr-1);
K = sparse(I,J,S,Ndof,Ndof);

% ----------------- enforce strong Dirichlet on flag 10 --------------------
if ~isempty(strongRows)
    [strongRows, ia] = unique(strongRows, 'stable');
    strongVals = strongVals(ia);

    mask = false(Ndof,1);  mask(strongRows) = true;

    % move known values to RHS, zero columns
    F = F - K(:,mask) * strongVals;
    K(:,mask) = 0;

    % set rows to identity, RHS to values
    K(mask,:) = 0;
    K = K + sparse(strongRows,strongRows,1,Ndof,Ndof);
    F(mask) = strongVals;
end

end % ================= end function ========================================

% -------------------- local helper: P2 basis on reference tri --------------
function [N, dNxi, dNeta] = p2basis_ref(xi, eta)
zeta = 1 - xi - eta;   % λ1
N1 = zeta*(2*zeta - 1);  N2 = xi*(2*xi - 1);  N3 = eta*(2*eta - 1);
N4 = 4*xi*zeta;          N5 = 4*xi*eta;       N6 = 4*eta*zeta;
N  = [N1; N2; N3; N4; N5; N6];
dN1_dxi  = 1 - 4*zeta;   dN1_deta = 1 - 4*zeta;
dN2_dxi  = (2*xi - 1) + 2*xi;      dN2_deta = 0;
dN3_dxi  = 0;                       dN3_deta = (2*eta - 1) + 2*eta;
dN4_dxi  = 4*(zeta - xi);          dN4_deta = -4*xi;
dN5_dxi  = 4*eta;                  dN5_deta = 4*xi;
dN6_dxi  = -4*eta;                 dN6_deta = 4*(zeta - eta);
dNxi  = [dN1_dxi;  dN2_dxi;  dN3_dxi;  dN4_dxi;  dN5_dxi;  dN6_dxi];
dNeta = [dN1_deta; dN2_deta; dN3_deta; dN4_deta; dN5_deta; dN6_deta];
end
