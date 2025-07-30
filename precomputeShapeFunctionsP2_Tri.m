function [N, dNxi, dNeta, g_wt, g_pt] = precomputeShapeFunctionsP2_Tri()
% precomputeShapeFunctionsP2_Tri:
%  Returns shape functions (N), derivatives (dNxi, dNeta),
%  and 6-point Gauss integration rule (g_pt, g_wt)
%  for a P2 triangle in the reference domain 0<=xi, 0<=eta, xi+eta<=1.
%
%  Gmsh node ordering for P2 (anticlockwise):
%    Node1 -> corner1
%    Node2 -> corner2
%    Node3 -> corner3
%    Node4 -> mid-edge(corner1-corner2)
%    Node5 -> mid-edge(corner2-corner3)
%    Node6 -> mid-edge(corner3-corner1)

    % Get 6-point Gauss rule over the reference triangle
    [g_pt, g_wt] = triGaussPoints6();
    numGauss = length(g_wt);

    numNodes = 6;  % 6 local nodes for P2
    N      = zeros(numNodes, numGauss);
    dNxi   = zeros(numNodes, numGauss);
    dNeta  = zeros(numNodes, numGauss);

    for k = 1:numGauss
        xi  = g_pt(k,1);
        eta = g_pt(k,2);

        % P2 shape functions matching Gmsh anticlockwise node ordering
        [Ni, dNxi_k, dNeta_k] = p2basisGmsh(xi, eta);

        N(:,k)     = Ni;
        dNxi(:,k)  = dNxi_k;
        dNeta(:,k) = dNeta_k;
    end
end

%--------------------------------------------------------------------
function [N, dNxi, dNeta] = p2basisGmsh(xi, eta)
% p2basisGmsh:
%   Returns the 6 shape functions and their partial derivatives
%   wrt (xi, eta) for a quadratic triangle, with node ordering:
%     (1) corner1, (2) corner2, (3) corner3,
%     (4) edge12,  (5) edge23,  (6) edge31.
%
% Let zeta = 1 - xi - eta.

    zeta = 1 - xi - eta;  % corner1 (anticlockwise from node1)

    % shape functions
    N1 = zeta*(2*zeta - 1);  % node1
    N2 = xi*(2*xi - 1);      % node2
    N3 = eta*(2*eta - 1);    % node3
    N4 = 4*xi*zeta;          % mid-edge(1-2)
    N5 = 4*xi*eta;           % mid-edge(2-3)
    N6 = 4*eta*zeta;         % mid-edge(3-1)

    N = [N1; N2; N3; N4; N5; N6];

    % partial derivatives wrt xi, eta
    % zeta = 1 - xi - eta => dzeta/dxi = -1, dzeta/deta = -1

    dN1_dxi = (2*(zeta) - 1)*(-1) + zeta*(2*(-1));
        % expanded: dN1_dxi = - (2zeta - 1) - 2*zeta
        %           = -2zeta + 1 - 2zeta = 1 - 4zeta
    dN1_deta = 1 - 4*zeta; % same pattern

    dN2_dxi  = (2*xi - 1) + xi*(2); % chain rule for xi*(2xi-1)
    dN2_deta = 0; % no eta in N2

    dN3_dxi  = 0; % no xi in N3
    dN3_deta = (2*eta - 1) + eta*(2);

    dN4_dxi  = 4*zeta + 4*xi*(-1);  % 4(zeta dxi + xi dzeta/dxi= -xi)
               % = 4*(zeta - xi)
    dN4_deta = 4*( xi * (-1) );     % = -4 xi

    dN5_dxi  = 4*eta;
    dN5_deta = 4*xi;

    dN6_dxi  = 4*( eta * (-1) );  % = -4 eta
    dN6_deta = 4*( zeta + eta*(-1) );
               % = 4*(zeta - eta)

    dNxi  = [ dN1_dxi;  dN2_dxi;  dN3_dxi;  dN4_dxi;  dN5_dxi;  dN6_dxi ];
    dNeta = [ dN1_deta; dN2_deta; dN3_deta; dN4_deta; dN5_deta; dN6_deta ];
end

%--------------------------------------------------------------------
function [g_pt, g_wt] = triGaussPoints6()
% triGaussPoints6: returns a 6-point integration rule on the
%  reference triangle (xi,eta>=0, xi+eta<=1).
% Each row of g_pt is (xi, eta), and the corresponding weight is in g_wt.
% This rule integrates polynomials exactly up to degree 4.
%
% Weights sum to area = 1/2 for the standard reference triangle,
% so each weight is the sub-area contribution.
%
% Source of these points is standard in many FE references.

    g_pt = [ ...
      0.4459484909, 0.4459484909;
      0.4459484909, 0.1081030182;
      0.1081030182, 0.4459484909;
      0.0915762135, 0.0915762135;
      0.0915762135, 0.8168475730;
      0.8168475730, 0.0915762135 ];

    g_wt = 0.5*[ ...
      0.2233815897;
      0.2233815897;
      0.2233815897;
      0.1099517437;
      0.1099517437;
      0.1099517437 ];
end
