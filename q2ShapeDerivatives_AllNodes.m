
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated function for shape derivatives at a single Gauss point
% using precomputed full shape functions and derivative multipliers.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dNxAll, dNyAll, detJ] = q2ShapeDerivatives_AllNodes( xcoords, ycoords, dNxi, dNeta )
    % Compute the Jacobian components from the derivative multipliers
    dX_dxi   = sum( xcoords .* dNxi );
    dX_deta  = sum( xcoords .* dNeta );
    dY_dxi   = sum( ycoords .* dNxi );
    dY_deta  = sum( ycoords .* dNeta );
    
    J = [ dX_dxi,  dY_dxi ;
          dX_deta, dY_deta ];
    
    detJ = dX_dxi*dY_deta - dX_deta*dY_dxi;
    invJ = inv(J);
    
    % Transform derivatives to physical coordinates
    dNxAll = invJ(1,1)*dNxi + invJ(1,2)*dNeta;
    dNyAll = invJ(2,1)*dNxi + invJ(2,2)*dNeta;
end


