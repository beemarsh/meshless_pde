%Lshape_Domain = LDomain([-1,1,-1,1],[-1,0,-1,0]);
%Lshape_Domain = Lshape_Domain.generateLShape(100);
%Lshape_Domain.scatterPlot('L-Shaped Domain [-1,1]\[-1,0]');

domain = [-1,1];


% Define the size of the grid
n = 100; % Grid size
n_boundary = 10;

% Orientation: Specify the orientation of the L-shape
% Options: 'top-left', 'top-right', 'bottom-left', 'bottom-right'
orientation = 'bottom-left';

% Generate the L-shaped grid using numgrid
L_grid = numgrid('L', n);

% Rotate the L-shaped grid based on the specified orientation
switch orientation
    case 'top-left'
        % No rotation needed
    case 'top-right'
        L_grid = rot90(L_grid, 1); % Rotate 90 degrees clockwise
    case 'bottom-left'
        L_grid = rot90(L_grid, -1); % Rotate 90 degrees counterclockwise
    case 'bottom-right'
        L_grid = rot90(L_grid, 2); % Rotate 180 degrees
    otherwise
        error('Invalid orientation. Choose from ''top-left'', ''top-right'', ''bottom-left'', or ''bottom-right''.');
end

% Define the domain for the meshgrid
x_domain = linspace(domain(1), domain(2), n);
%y_domain = linspace(-1, 1, n);

% Generate the meshgrid
[X, Y] = meshgrid(x_domain, x_domain);

% Extract points corresponding to the L-shaped domain
L_mask = L_grid > 0; % Logical mask where L_grid > 0

% Apply the mask to get points inside the L-shaped domain
X_L = X(L_mask);
Y_L = Y(L_mask);

% Generate boundary points
%Only for bottom left

top_points = [linspace(-1,1,n_boundary);ones(1,n_boundary)];
top_right_points = [ones(1,n_boundary);linspace(-1,1,n_boundary)];

%Boundary points where they are reduced to half
top_left = [-1 * ones(1,n_boundary/2);linspace(0,1,n_boundary/2)];
bottom_L = [linspace(-1,0,n_boundary/2); zeros(1,n_boundary/2)];
right_L = [zeros(1,n_boundary/2);linspace(-1,0,n_boundary/2)];
bottom_bottom = [linspace(0,1,n_boundary/2); -1 * ones(1,n_boundary/2)];

boundary_points = [top_points, top_right_points,top_left, bottom_L, right_L, bottom_bottom];


% Plot the points in the L-shaped domain
figure;
scatter(X_L, Y_L, 40, 'filled');
hold on;
scatter(boundary_points(1,:),boundary_points(2,:), 40, 'filled','r');
grid on;
axis equal;
xlabel('X');
ylabel('Y');
title(['L-Shaped Domain Points - Orientation: ', strrep(orientation, '-', ' ')]);
hold off;
