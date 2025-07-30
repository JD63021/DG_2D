function uv_new = update_bc(uv, boundaryInfo, nodeInfo, Nxy, t, bcFlags, inletProfile)
% update_bc:
%   Applies (Dirichlet) boundary conditions to the solution vector uv.
%   Originally for 2D velocity, but simplified to allow 1D or 2D usage.
%
% Inputs:
%   uv           : solution vector. In a 2D velocity problem, this might be
%                  [U_x; U_y; possibly pressure], but in 1D conduction it might
%                  just be [T].
%   boundaryInfo : structure with boundary sets, e.g. boundaryInfo.flag_1, etc.
%   nodeInfo     : structure with node coordinates (nodeInfo.velocity.x, etc.)
%   Nxy          : number of DOFs for the first component (e.g. # of nodes if single DOF).
%   t            : current time (if your inletProfile depends on it).
%   bcFlags      : struct with fields (for example) .inlet and .wall,
%                  each containing boundary flag IDs.
%   inletProfile : either a numeric constant or a function handle for the BC.
%                  If it's a function handle, we call:  inletProfile(t, xVal)
%
% Output:
%   uv_new       : solution vector with boundary conditions applied.

uv_new = uv;

%% Example usage:
% - If you are in a 1D conduction scenario, you probably have just one DOF per node.
%   Then "globalRow" can just return nodeID itself if you like,
%   or you can adapt it accordingly.
%
% - If you are in 2D Navier–Stokes, you might have 2 velocity DOFs plus pressure, etc.

% ---------------------------------------------------------------------
%  Helper function: get the global row index for the 'x' (or 'y') DOF
    function row = globalRow_1D(vNode)
        % For a single-DOF 1D case, the global DOF = vNode.
        % If you have a 2D velocity approach, you’d do:
        %   row = vNode  for the 'x' component
        %   row = vNode + Nxy for 'y'
        row = vNode;
    end

% ---------------------------------------------------------------------
% 1) Apply "inlet" BC (Dirichlet).
if isfield(bcFlags, 'inlet')
    inFlagList = bcFlags.inlet;
    if ~iscell(inFlagList), inFlagList = num2cell(inFlagList); end
    % (If bcFlags.inlet is a single integer, that’s fine too.)

    for iF = 1:numel(inFlagList)
        flagVal = inFlagList{iF};
        % If bcFlags.inlet is numeric, remove the curly braces usage or adapt as needed.
        if isnumeric(flagVal), flagVal = flagVal(1); end
        flagStr = ['flag_' num2str(flagVal)];
        if isfield(boundaryInfo, flagStr)
            inletNodes = boundaryInfo.(flagStr)(:).';  % row vector
        else
            inletNodes = [];
        end

        for nodeID = inletNodes
            % For a 1D case, we only have one DOF -> row = nodeID:
            rowX = globalRow_1D(nodeID);

            % Obtain the node’s x-coordinate (or y, if relevant). For 1D we just do x:
            xVal = nodeInfo.x(nodeID);

            % If inletProfile is a function handle, call it; else assume numeric.
            if isa(inletProfile, 'function_handle')
                Uin = inletProfile(t, xVal);
            else
                Uin = inletProfile;  % numeric
            end

            uv_new(rowX) = Uin;
        end
    end
end

% ---------------------------------------------------------------------
% 2) Apply "wall" BC (e.g. zero velocity or zero temperature).
if isfield(bcFlags, 'wall')
    wFlagList = bcFlags.wall;
    if ~iscell(wFlagList), wFlagList = num2cell(wFlagList); end

    for iF = 1:numel(wFlagList)
        flagVal = wFlagList{iF};
        if isnumeric(flagVal), flagVal = flagVal(1); end
        flagStr = ['flag_' num2str(flagVal)];
        if ~isfield(boundaryInfo, flagStr), continue; end
        sideNodes = boundaryInfo.(flagStr)(:);
        sideNodes = unique(sideNodes);

        for nodeID = sideNodes'
            rowX = globalRow_1D(nodeID);
            uv_new(rowX) = 0;
        end
    end
end

end

