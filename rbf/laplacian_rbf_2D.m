function [ value ] = laplacian_rbf_2D( c, r )
% Laplace rbf of NMQ
% c is the shape parameter or a list of shape parameters when usimg variable shape
% r is the distance matrix
import rbf.*;
% value = 2 * c.^2 .* (rbf(c, r) .^ 2 + 1) ./rbf(c, r).^3;

value = ( 2 * (c.^ 2)./ rbf(c, r) ) - ( (c .^ 4) .* ( r.^2)./(rbf(c, r).^3) );
end