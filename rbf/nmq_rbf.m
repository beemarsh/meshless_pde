function [ value ] = nmq_rbf( c, r )
% Calculate RBF MQ = sqrt (1 + ξ^2 * r^2 )

value = sqrt(1+(c.*r).^2);

end