% Poison equation
% comparason between ghost point and Kansa method
% Kansas method on a square domain.
clear;
warning off all;

addpath('./functions/');
addpath('./domain')
addpath('./rbf/');
addpath('./matrix/')

u = @(x, y) exp(x + y);
laplace_u = @(x, y) 2*exp(x + y);

%% define number of points
Nb = 100;
N = 400;
Ni = N - Nb;

theta = linspace(0, 2*pi, Nb+1);
theta(end)=[];


r_parametric = @(theta, scale) scale.*(cos(3*theta) + (2 - sin(3*theta).^2).^0.5).^(1/3);
r = @(theta) r_parametric(theta, 1);

square_domain = Square([0,1]);
square_domain = square_domain.generateGrid(19);

boundary_pts = square_domain.BoundaryPoints;


% boundary_pts = [r(theta).*cos(theta); r(theta).*sin(theta)]';
interior_pts = square_domain.InteriorPoints;

% [interior_pts, ~, domain_center] = within_polygon(boundary_pts, Ni);
all_pts = [interior_pts; boundary_pts];

%% test points
% t = linspace(0, 2*pi, 50);
% test_boundary_pts = [r(t).*cos(t); r(t).*sin(t)]';

x_c = 0.5;
y_c=0.5; W = 1; H = 1; NUM_TEST = 50;

% Compute min and max coordinates
x_min = x_c - W/2; x_max = x_c + W/2;
y_min = y_c - H/2; y_max = y_c + H/2;

% Generate parameter t (from 0 to 1)
t = linspace(0, 1, NUM_TEST+1); % Extra point for closure
t(end) = []; % Remove duplicate last point

% Compute x and y coordinates
x = zeros(size(t));
y = zeros(size(t));

% Assign points based on t range (moving counterclockwise)
for i = 1:NUM_TEST
    if t(i) < 1/4
        x(i) = x_min + 4 * t(i) * W;
        y(i) = y_min;
    elseif t(i) < 1/2
        x(i) = x_max;
        y(i) = y_min + 4 * (t(i) - 1/4) * H;
    elseif t(i) < 3/4
        x(i) = x_max - 4 * (t(i) - 1/2) * W;
        y(i) = y_max;
    else
        x(i) = x_min;
        y(i) = y_max - 4 * (t(i) - 3/4) * H;
    end
end

test_boundary_pts = [x;y]';


[test_interior_pts, ~, ~] = within_polygon(test_boundary_pts, 100); %test points
test_points = [test_interior_pts;test_boundary_pts];

%% exact
u_exact = u(test_points(:, 1), test_points(:, 2));

rhs = [laplace_u(interior_pts(:, 1),interior_pts(:, 2)); u(boundary_pts(:,1), boundary_pts(:,2))];


for R = linspace(0.5, 50, 100)

    %% ghost points
    [x_ghost, y_ghost] = fabric_pattern(N,1);
    ghost_pts =[R*x_ghost, R*y_ghost];

    %% Matrix
    dm_test_points = distance_matrix(test_points, ghost_pts);
    dm_interior_pts = distance_matrix(interior_pts, ghost_pts);
    dm_boundary_pts = distance_matrix(boundary_pts, ghost_pts);

    %% Shape parameter (the same length as basis)
    shape_franke = franke(N, 2*R);
    [shapes, shape_min, shape_max] = variable_shape(length(ghost_pts), shape_franke, 0.3, 0.5);

    %% computation
    basis_test_points = nmq_rbf(shapes, dm_test_points);
    basis_interior_pts = laplacian_rbf_2D(shapes, dm_interior_pts);
    basis_boundary_pts = nmq_rbf(shapes, dm_boundary_pts);
    basis_matrix = [basis_interior_pts; basis_boundary_pts];
    %% solve
    % First find the coefficient vector
    % α = A^−1 * b
    alpha = basis_matrix \ rhs;

    % Find the solution at the test points
    u_approx = basis_test_points * alpha;
    diff = u_exact - u_approx;
    error = norm(diff,Inf);

    fprintf('R =%4.1f,  Maxerr = %10.3e,  c_i=%5.3f,  c_range=(%5.3f, %5.3f)\n', R, max(error), shape_franke, shape_min, shape_max );
end