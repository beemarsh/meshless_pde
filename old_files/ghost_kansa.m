addpath('./domain/');
addpath('./functions/');
addpath('./plot/')
addpath('./rbf/');
addpath('./matrix/')

m=[1 2 1 2 3 1 3 2 4 ];
l=[1 1 2 2 1 3 2 3 1 ];


radii = linspace(1, 10, 500);


%count the time of the main part of the code
t0=clock;

Nb = 100;
N = 400;
Ni = N- Nb;

% Center of the Square
x_c= 0.5; y_c= 0.5;
% length of square
L = 1;

[x_b , y_b] = rectangle_boundary(x_c,y_c, L,L, Nb);
% boundary_pts = Cassini(Nb);
boundary_pts = [x_b; y_b]';

[interior_pts, ~, domain_center] = within_polygon(boundary_pts, Ni);
all_pts = [interior_pts; boundary_pts];


% Store the result for eigenvalues. Each column stores the approximate eigenvalues for each radius.

result_eigenvalues = zeros(length(m), length(radii));

error_eigenvalues = zeros(length(m), length(radii));

t0=clock;

% for R = linspace(0.5, 50, 100)
i=0;
for R= radii
    i=i+1;
    % ghost points
    [x_ghost, y_ghost] = fabric_pattern(N, 1, x_c, y_c,R);
    % ghost_pts =[R*x_ghost, R*y_ghost];
    ghost_pts =[x_ghost, y_ghost];

    DM_interior= distance_matrix(interior_pts, ghost_pts);
    DM_boundary= distance_matrix(boundary_pts, ghost_pts);
    DM_all= [DM_interior; DM_boundary];

    % Shape parameter (the same length as basis)
    shape_franke = franke(N, 2*R);

    shapes = 2.69;
    % [shapes, shape_min, shape_max] = variable_shape(length(ghost_pts), shape_franke, 0.5, 0.5);

    % computation
    % basis_test_points = nmq_rbf(shapes, dm_test_points);
    basis_interior_pts = laplacian_rbf_2D(shapes, DM_interior);
    basis_boundary_pts = nmq_rbf(shapes, DM_boundary);

    % We have a linear system: -A α = λ B α
    A = [basis_interior_pts; basis_boundary_pts];

    B = [rbf(shapes, DM_interior); zeros(Nb, N)];

    % Now let: V = A * inverse(B)
    V = A \ B;

    % Now, we calculate the eigenvalues.
    % Note that we are using negative because our equation has it.
    [alpha, lambda] = eigs(-V, length(m), 0);

    % extract the first length(m) eigenvalues;
    % We need to extract them from the diagonal.

    approximate_eigenvalues = real(diag(lambda));

    alpha = real(alpha);

    % Compute the exact eigenvalues and then calculate the relative error
    exact_eigenvalues = pi^2*(m.^2.+l.^2)';
    relative_error = (abs(approximate_eigenvalues' - pi^2*(m.^2.+l.^2)))';

    result_eigenvalues(:,i) = approximate_eigenvalues;

    error_eigenvalues(:,i) = relative_error;
end

cpu=etime(clock,t0); % count the time of the main part of the code

fprintf('Run time : %6.2f\n',cpu);

eigenvalue_Max_errors = zeros(1, length(radii));
for i=1:length(radii)
    eigenvalue_Max_errors(i) = max(error_eigenvalues(:,i));
end

figure;
plot(radii, eigenvalue_Max_errors, 'o-');
xlabel('Shape Parameter');
ylabel('Max Error');
title('Max Error of USING GHOST POINTS(FIXED SHAPE PARAMETER)');


function [x, y] = rectangle_boundary(x_c, y_c, W, H, N)
% Generates N boundary points of a rectangle centered at (x_c, y_c)
% INPUTS:
%   x_c, y_c - Center of the rectangle
%   W, H - Width and Height of the rectangle
%   N - Total number of boundary points
% OUTPUTS:
%   x, y - Cartesian coordinates of boundary points

% Compute min and max coordinates
x_min = x_c - W/2; x_max = x_c + W/2;
y_min = y_c - H/2; y_max = y_c + H/2;

% Generate parameter t (from 0 to 1)
t = linspace(0, 1, N+1); % Extra point for closure
t(end) = []; % Remove duplicate last point

% Compute x and y coordinates
x = zeros(size(t));
y = zeros(size(t));

% Assign points based on t range (moving counterclockwise)
for i = 1:N
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

end
