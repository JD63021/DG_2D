function testCase = testCase()
% fpc - Test case configuration for the flow past cylinder problem.
%
% This file specifies:
%   (i) the gmsh files to run (now two files: one for velocity and one for pressure),
%   (ii) any velocity boundary flags to exclude (if desired),
%   (iii) pressure boundary flags and values (if desired), and 
%   (iv) the boundary flags to use for enforcing Dirichlet conditions.
%
% For fpc, we choose:
%   - Top boundary flag = 2.
%   - Typical side and bottom boundary flags = [1, 4, 5, 6, 7, 8].
%
% Note: fpc2q (velocity grid) and fpc2l (pressure grid) replace the older fpc2r.
% FPC9 GIVES BEST FIT AS OF YET!!!

    testCase.gmshVelFile  = 'squarecd1.m';  % gmsh file for velocity grid (Q2)
    % testCase.gmshPresFile = 'quadp1.m';  % gmsh file for pressure grid (Q1)
    
    testCase.excludedVelFlags = [];             % velocity boundary flags to exclude
    testCase.pBoundFlags      = [];                % (no pressure BCs in this case)
    testCase.pBoundVals       = [];
    testCase.boundaryFlags.inlet = [];              % top boundary flag (for inlet)
    testCase.boundaryFlags.wall  = [9,10,11]; % side and bottom boundaries
    testCase.order = 1;
    testCase.constantsource = 0;

    % Specify the inlet velocity profile (function handle).
    % This function should accept: (t, y, H) and return the x-velocity.
    testCase.inletProfile   = @(x) 1;
    testCase.inletProfileSS = 1;%4*0.3*(y*(H-y))/(H^2);
    % testCase.corner = 2;%113;
end

