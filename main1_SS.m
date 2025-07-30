function [U, x,y, nodeInfo, elemInfo, boundaryInfo] = ...
    main1_SS(  ...
             nodeInfo, elemInfo, boundaryInfo, time,penalty,epsilon,beta)
% main1_SS - Stationary solver for the Navier–Stokes equations.
%
% Inputs:
%   U1, U2, U3       : initial guesses for x-velocity, y-velocity, and pressure.
%   D, Re, o, gamma, mu, rho : physical parameters.
%   nt               : number of iterations (for steady state, typically 1).
%   nodeInfo, elemInfo, boundaryInfo : mesh and boundary data.
%   time             : simulation time (used for evaluating BCs).
%   bcFlags          : structure with boundary flags and an optional inletProfile.
%   inletProfile     : function handle for the inlet velocity profile.
%
% Outputs:
%   U                : final solution vector.
%   x, y             : velocity node coordinates.
%   x1, y1           : pressure node coordinates.
%   z1, z2, z3, z    : solution components for visualization.
%   T                : final "time" (remains the same for a steady state).
%   FF, J, ...       : additional outputs from the Newton iteration.

%% 1) Precompute shape functions and related quantities.
% [quad1_s3, quad2_s3, dq1_s3, dq2_s3, lin1_s3, lin2_s3, ...
%  g_wt_s3, g_pt_s3] = precomputeShapeFunctions2();

Nxy = length(nodeInfo.x);
% Npr = max(elemInfo.presElements(:));

%% 2) Initialize the solution.
% U = [U1; U2; U3];  % initial solution vector
currentTime = time;  % For a stationary problem, time is arbitrary.




   f_coeff = 10;
   f_handle = ones(numel(nodeInfo.x),1)*10;

   % [K, F] = build_stiffness_operator_dg( nodeInfo, elemInfo, boundaryInfo, 10,3000);

   [K, F] = build_stiffness_operator_dgcd( nodeInfo, elemInfo, boundaryInfo, 0,penalty,epsilon, @(x,y) 0, beta);

   % [K, F] = build_stiffness_operator_s3( nodeInfo, elemInfo, boundaryInfo, f_handle);

   
    % After you assemble K,F (or inside main1_SS right after [K,F]=...):
xy  = [nodeInfo.x, nodeInfo.y];
E   = elemInfo.elements;
areas = arrayfun(@(e) polyarea(xy(E(e,:),1), xy(E(e,:),2)), (1:size(E,1))');
domainArea = sum(areas);          % should be ≈ 1 for your unit square
fprintf('RHS check: sum(F)=%.6f,  expected ≈ f_coeff*Area=%.6f\n', ...
        sum(F), f_coeff*domainArea);

  


  
    % Compute the Newton update.
    U = K \ F;
  




% Set node coordinates.
x = nodeInfo.x; y = nodeInfo.y;


end

