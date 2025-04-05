
function [x,y] = uniform_circle(N, x_center, y_center, R)
% Generate N quasi-uniform points inside a circle using corrected polar coordinates
% INPUTS:
%   x_center, y_center: Center coordinates of the circle
%   R: Radius of the circle
%   N: Number of points to generate
% OUTPUT:
%   points: N x 2 matrix of [x, y] coordinates

if nargin < 2, x_center = 0; end
if nargin < 3, y_center = 0; end
if nargin < 4, R = 1; end

bases = [2, 3]; % Prime bases for dimensions
halton = sobolset(2); % Use Sobol for better performance
points_unit = net(halton, N);

% Map to polar coordinates (corrected radial distribution)
r = sqrt(points_unit(:,1)); % r ∈ [0,1]
theta = 2 * pi * points_unit(:,2); % θ ∈ [0,2π)

% Convert to Cartesian and scale/center
x = x_center + R * r .* cos(theta);
y = y_center + R * r .* sin(theta);


% x =x'; y= y';
end