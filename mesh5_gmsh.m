function [nodeInfo, elemInfo, boundaryInfo] = mesh5_gmsh(gmshFile)
% mesh5_gmsh (P2 triangles aware)
%   Reads a 2D Gmsh-exported mesh struct 'msh' from a MATLAB script or MAT
%   file (the argument is the filename of the script that defines 'msh').
%
%   Supports:
%     - P2 triangles:   msh.TRIANGLES6  (Ne × 7)  [n1 n2 n3 n4 n5 n6 physID]
%     - P1 triangles:   msh.TRIANGLES   (Ne × 4)  [n1 n2 n3 physID]   (fallback)
%     - Quadratic lines:msh.LINES3      (Nb × 4)  [n1 n2 nmid physID]
%     - Linear lines:   msh.LINES       (Nb × 3)  [n1 n2 physID]      (fallback)
%
% Outputs
%   nodeInfo:
%     .x, .y        : Nnodes×1 coordinates
%     .id           : 1…Nnodes
%
%   elemInfo:
%     .elements     : Ne×Nloc connectivity (Nloc=6 if TRIANGLES6, else 3)
%     .physIDs      : Ne×1 physical region IDs
%     .order        : 2 for P2, 1 for P1
%
%   boundaryInfo:
%     .edgesP2      : Nb×3 quadratic boundary edges [n1 n2 nmid]   (if LINES3)
%     .edges        : Nb×2 linear endpoints [n1 n2]                (always set)
%     .edgeFlags    : Nb×1 physical IDs (from LINES/LINES3)
%     .allBoundaryNodes : unique node IDs on the boundary (includes mids if P2)
%     .flag_<id>    : (optional) list of node IDs that belong to boundary group <id>
%
% Notes
%   • All triangles are re‑oriented to CCW using vertex coordinates.
%     For P2, if a triangle is CW we swap vertex nodes 2↔3 **and** mid‑edge
%     nodes 4↔6 so the (12,23,31) mid‑edge convention remains correct.
%   • This keeps your P2 DG assembly consistent with the Gmsh node ordering.

% ----------------- 1) load mesh struct -----------------
if ischar(gmshFile)
    run(gmshFile);   % should define variable 'msh'
else
    error('gmshFile must be a filename string that defines variable ''msh''.');
end
if ~exist('msh','var')
    error('No variable ''msh'' found after running %s.', gmshFile);
end

% ----------------- 2) node coordinates -----------------
coords    = msh.POS;            % Nnodes×2 or ×3
numNodes  = size(coords,1);
nodeInfo.x  = coords(:,1);
nodeInfo.y  = coords(:,2);
nodeInfo.id = (1:numNodes).';

% For convenience in orientation checks
x = nodeInfo.x;  y = nodeInfo.y;

% ----------------- 3) elements (prefer P2 triangles) ----
isP2 = false;
if isfield(msh,'TRIANGLES6') && ~isempty(msh.TRIANGLES6)
    tri6 = msh.TRIANGLES6;     % Ne × 7
    if size(tri6,2) < 7
        error('msh.TRIANGLES6 must be Ne×7: [n1 n2 n3 n4 n5 n6 physID].');
    end
    elems = tri6(:,1:6);       % Ne×6
    phys  = tri6(:,7);         % Ne×1
    isP2  = true;
elseif isfield(msh,'TRIANGLES') && ~isempty(msh.TRIANGLES)
    tri3 = msh.TRIANGLES;      % Ne × 4
    if size(tri3,2) < 4
        error('msh.TRIANGLES must be Ne×4: [n1 n2 n3 physID].');
    end
    elems = tri3(:,1:3);       % Ne×3
    phys  = tri3(:,4);         % Ne×1
    isP2  = false;
else
    error('No TRIANGLES6 (P2) or TRIANGLES (P1) found in msh.');
end
Ne = size(elems,1);

% --------- 3a) Enforce CCW orientation for every triangle ----------
% signed area A2 = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);  CCW if A2>0
if isP2
    v1 = elems(:,1); v2 = elems(:,2); v3 = elems(:,3);
    A2 = (x(v2)-x(v1)).*(y(v3)-y(v1)) - (x(v3)-x(v1)).*(y(v2)-y(v1));
    cw = (A2 < 0);   % clockwise triangles
    if any(cw)
        % swap vertex 2 and 3
        tmp = elems(cw,2);              elems(cw,2) = elems(cw,3);  elems(cw,3) = tmp;
        % swap mid-edge nodes 4 and 6 (keep mid 5 as is)
        tmp = elems(cw,4);              elems(cw,4) = elems(cw,6);  elems(cw,6) = tmp;
    end
else
    v1 = elems(:,1); v2 = elems(:,2); v3 = elems(:,3);
    A2 = (x(v2)-x(v1)).*(y(v3)-y(v1)) - (x(v3)-x(v1)).*(y(v2)-y(v1));
    cw = (A2 < 0);
    if any(cw)
        % swap vertex 2 and 3 to make CCW
        tmp = elems(cw,2);              elems(cw,2) = elems(cw,3);  elems(cw,3) = tmp;
    end
end

elemInfo.elements = elems;       % Ne×6 (P2) or Ne×3 (P1)
elemInfo.physIDs  = phys;        % Ne×1
elemInfo.order    = 2*isP2 + 1*(~isP2);   % 2 for P2, 1 for P1

% ----------------- 4) boundary edges -------------------
boundaryInfo = struct();

edgesP2 = [];  edgesLin = [];  eFlags = [];

if isfield(msh,'LINES3') && ~isempty(msh.LINES3)
    % Quadratic boundary lines: [n1 n2 nmid physID]
    L3 = msh.LINES3;
    if size(L3,2) < 4
        error('msh.LINES3 must be Nb×4: [n1 n2 nmid physID].');
    end
    edgesP2 = L3(:,1:3);           % Nb × 3
    edgesLin= L3(:,1:2);           % endpoints only
    eFlags  = L3(:,4);             % phys IDs
elseif isfield(msh,'LINES') && ~isempty(msh.LINES)
    % Linear boundary lines: [n1 n2 physID]
    L1 = msh.LINES;
    if size(L1,2) < 3
        error('msh.LINES must be Nb×3: [n1 n2 physID].');
    end
    edgesLin = L1(:,1:2);
    eFlags   = L1(:,3);
    % No mid-nodes available
else
    warning('No boundary lines (LINES3/LINES) found → boundaryInfo will be empty.');
end

boundaryInfo.edgesP2    = edgesP2;              % may be empty
boundaryInfo.edges      = edgesLin;             % endpoints (always set, may be empty)
boundaryInfo.edgeFlags  = eFlags(:);            % phys IDs (may be empty)

% unique nodes on boundary (include mids if present)
if ~isempty(edgesP2)
    boundaryInfo.allBoundaryNodes = unique([edgesP2(:,1); edgesP2(:,2); edgesP2(:,3)]);
elseif ~isempty(edgesLin)
    boundaryInfo.allBoundaryNodes = unique([edgesLin(:,1); edgesLin(:,2)]);
else
    boundaryInfo.allBoundaryNodes = [];
end

% group by physical ID (optional convenience fields)
if ~isempty(eFlags)
    flags = unique(eFlags);
    for k = 1:numel(flags)
        f = flags(k);
        mask = (eFlags == f);
        if ~isempty(edgesP2)
            nodes = unique([edgesP2(mask,1); edgesP2(mask,2); edgesP2(mask,3)]);
        else
            nodes = unique([edgesLin(mask,1); edgesLin(mask,2)]);
        end
        boundaryInfo.(sprintf('flag_%d',f)) = nodes;
    end
end
end
