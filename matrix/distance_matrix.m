% This function D calculates the distance between each points in the matrix x and y.
% Here, x is a n*2 matrix and y is a m*2 matrix.
% Each row denotes a set of x,y collocation points
% Using bsxfun, we calculate the difference between each points.
% Using hypot, we calculate the Euclidean distance.
% Essentially, all its doing is: sqrt( (x_1 - x_2)^2 + (y_1 - y_2)^2 )

function D = distance_matrix(x, y)
% Calculates the Euclidean distance between points in matrices x and y
%
% Inputs:
%   x: n×2 matrix where each row is a point [x,y]
%   y: m×2 matrix where each row is a point [x,y]
%
% Output:
%   D: n×m matrix of Euclidean distances between points in x and y

D = hypot(bsxfun(@minus,x(:,1),y(:,1)'), bsxfun(@minus,x(:,2),y(:,2)'));
end