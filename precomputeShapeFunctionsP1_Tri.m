function [N, dNxi, dNeta, g_wt, g_pt] = precomputeShapeFunctionsP1_Tri()
% precomputeShapeFunctionsP1_Tri:
%  Same pattern, but for a linear (3-node) triangle in anticlockwise order:
%    Node1: corner1
%    Node2: corner2
%    Node3: corner3
%
% shape: N1=1 - xi - eta, N2=xi, N3=eta
% partial derivatives are constant.

    [g_pt, g_wt] = triGaussPoints6(); % we can reuse the same 6-pt rule
    numGauss = length(g_wt);

    numNodes = 3;
    N      = zeros(numNodes, numGauss);
    dNxi   = zeros(numNodes, numGauss);
    dNeta  = zeros(numNodes, numGauss);

    for k = 1:numGauss
        xi  = g_pt(k,1);
        eta = g_pt(k,2);

        % shape
        N1 = 1 - xi - eta;
        N2 = xi;
        N3 = eta;

        % derivatives
        dN1_dxi = -1;  dN1_deta = -1;
        dN2_dxi =  1;  dN2_deta =  0;
        dN3_dxi =  0;  dN3_deta =  1;

        N(:,k)      = [N1; N2; N3];
        dNxi(:,k)   = [dN1_dxi; dN2_dxi; dN3_dxi];
        dNeta(:,k)  = [dN1_deta; dN2_deta; dN3_deta];
    end
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
