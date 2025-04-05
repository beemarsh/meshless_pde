function coordinates = Cassini(N)
% Function to generate N points for the Cassini domain
% INPUT: N - Number of points
% OUTPUT: x, y - Cartesian coordinates of the points

theta = linspace(0, 2*pi, N); % Generate N values for theta
r = (cos(3*theta) + sqrt(2 - sin(3*theta).^2)).^(1/3); % Compute r(theta)

% Convert polar coordinates to Cartesian
x = r .* cos(theta);
y = r .* sin(theta);

coordinates = [x',y'];
end
