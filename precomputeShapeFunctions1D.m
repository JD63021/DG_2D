function [N, dN, gp, w] = precomputeShapeFunctions1D(order, nGP)
% precomputeShapeFunctions1D:
%   Returns shape function values (N) and derivatives (dN) for a 1D Lagrange
%   element of specified 'order' (1, 2, or 3), along with a Gauss rule on [0,1].
%
%   The integration is performed using a Gauss rule on [0,1]. By default,
%   nGP is set to 3 for orders 1 and 2 and to 5 for order 3. (This ensures that
%   when using cubic elements the error norm is integrated accurately.)
%
%   For quadratic and cubic elements the output ordering is chosen to mimic the
%   GMsh convention. For a quadratic element (3 nodes) the ordering is:
%     N1 (node at 0), N3 (node at 1), N2 (mid node).
%   For a cubic element (4 nodes) we assume the standard nodal positions:
%     ξ₁ = 0,  ξ₂ = 1/3,  ξ₃ = 2/3,  ξ₄ = 1,
%   but GMsh reorders the nodes as:
%     N1, N4, N2, N3.
%
% Inputs:
%   order : 1 => Linear (2 nodes)
%           2 => Quadratic (3 nodes)
%           3 => Cubic    (4 nodes)
%   nGP   : (optional) number of Gauss points to use on [0,1]. 
%           Default: nGP = 3 for order 1 & 2; nGP = 5 for order 3.
%
% Outputs:
%   N  : shape function values, size = (numLocalNodes) x (numGauss)
%   dN : shape function derivatives wrt ξ, same size as N.
%   gp : Gauss points in [0,1]
%   w  : Gauss weights for [0,1]

    %% 1) Determine number of Gauss points if not provided.
    if nargin < 2
        if order == 3
            nGP = 5;
        else
            nGP = 3;
        end
    end

    %% 2) Define Gauss integration rule on [0,1]
    % For nGP==3, use the original 3-point rule.
    if nGP == 3
        gp = [0.112701665379258, 0.5, 0.887298334620742];
        w  = [0.277777777777778, 0.444444444444444, 0.277777777777778];
    elseif nGP == 5
        % 5-point Gauss-Legendre rule on [0,1]:
        gp = [0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992];
        w  = [0.11846344, 0.23931434, 0.28444444, 0.23931434, 0.11846344];
    else
        error('nGP must be either 3 or 5.');
    end
    numGauss = length(gp);

    %% 3) Set number of local nodes
    switch order
        case 1
            numNodes = 2;  % linear
        case 2
            numNodes = 3;  % quadratic
        case 3
            numNodes = 4;  % cubic
        otherwise
            error('Order must be 1, 2, or 3.');
    end

    %% 4) Allocate output arrays
    N  = zeros(numNodes, numGauss);
    dN = zeros(numNodes, numGauss);

    %% 5) Define shape functions and derivatives for each order.
    switch order
        case 1  % Linear elements
            for ig = 1:numGauss
                xi = gp(ig);
                % Standard linear Lagrange on [0,1]:
                % L1 = 1 - xi, L2 = xi.
                N(1,ig) = 1 - xi;
                N(2,ig) = xi;
                
                dN(1,ig) = -1;
                dN(2,ig) =  1;
            end

        case 2  % Quadratic elements (using GMsh ordering: N1, N3, N2)
            for ig = 1:numGauss
                xi = gp(ig);
                % Polynomials chosen so that:
                % At xi=0: N1=1, N2=0, N3=0;
                % At xi=0.5: N1=0, N2=0, N3=1;
                % At xi=1: N1=0, N2=1, N3=0.
                N(1,ig) = 1 - 3*xi + 2*xi^2;   % node at 0
                N(3,ig) = 4*xi - 4*xi^2;         % mid node (at 0.5)
                N(2,ig) = -xi + 2*xi^2;          % node at 1
                
                dN(1,ig) = -3 + 4*xi;
                dN(3,ig) = 4 - 8*xi;
                dN(2,ig) = -1 + 4*xi;
            end

        case 3  % Cubic elements (using standard Lagrange on [0,1] then GMsh ordering)
            % Standard nodal positions for cubic element:
            nodes_std = [0, 1/3, 2/3, 1];
            % Precompute denominators for each basis function.
            denom = zeros(1,4);
            for i = 1:4
                prod_val = 1;
                for j = 1:4
                    if j ~= i
                        prod_val = prod_val * (nodes_std(i) - nodes_std(j));
                    end
                end
                denom(i) = prod_val;
            end
            
            % Allocate temporary arrays for standard ordering (Lagrange basis)
            L  = zeros(4, numGauss);
            dL = zeros(4, numGauss);
            
            % Compute the standard Lagrange basis functions and their derivatives.
            for ig = 1:numGauss
                xi_val = gp(ig);
                for i = 1:4
                    % Basis function L_i(xi):
                    prod_val = 1;
                    for j = 1:4
                        if j ~= i
                            prod_val = prod_val * (xi_val - nodes_std(j));
                        end
                    end
                    L(i,ig) = prod_val / denom(i);
                    
                    % Derivative dL_i(xi) using the product rule:
                    sum_val = 0;
                    for k = 1:4
                        if k ~= i
                            prod_der = 1;
                            for j = 1:4
                                if j ~= i && j ~= k
                                    prod_der = prod_der * (xi_val - nodes_std(j));
                                end
                            end
                            sum_val = sum_val + prod_der;
                        end
                    end
                    dL(i,ig) = sum_val / denom(i);
                end
            end
            
            % Reorder the standard basis functions to match GMsh ordering.
            % Standard ordering is: [L1, L2, L3, L4] corresponding to nodes [0, 1/3, 2/3, 1].
            % GMsh ordering is: [N1, N4, N2, N3] i.e. first corner, second corner, then the interior nodes.
            N(1,:) = L(1,:); % node at 0
            N(2,:) = L(4,:); % node at 1
            N(3,:) = L(2,:); % node at 1/3
            N(4,:) = L(3,:); % node at 2/3
            
            dN(1,:) = dL(1,:);
            dN(2,:) = dL(4,:);
            dN(3,:) = dL(2,:);
            dN(4,:) = dL(3,:);
    end
end
