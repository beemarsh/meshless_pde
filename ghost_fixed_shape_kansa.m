clear all;

addpath('./functions/');
addpath('./plot/')
addpath('./rbf/');
addpath('./matrix/')
addpath('./domain/');

%% Define the parameters
m=[1 2 1 2 3 1 3 2 4 ];
l=[1 1 2 2 1 3 2 3 1 ];

radii = [1]

shapes = [9];

N = 441; % Total number of points
Nb = 80; % Number of boundary points
Ni = N- Nb; % Number of interior points

bottom_left = [0, 0]; % Bottom left corner of the rectangle
L = 1; % Length of the rectangle
W = 1; % Width of the rectangle

%% Define the domain

rectangle = Rectangle(bottom_left, L, W, Nb, Ni);

interior_pts = rectangle.InteriorPoints;
boundary_pts = rectangle.BoundaryPoints;
all_pts = rectangle.Coordinates;

center = rectangle.Center;

% rectangle.plot("Rectangle Domain");

%% Construct the varaibles to store the results

result_eigenvalues = zeros(length(m), length(shapes), length(radii));

error_eigenvalues = zeros(length(m), length(shapes), length(radii));

%% Start the clock and run the loop

t0=clock;
i=0;

for R= radii
    i=i+1;
    % Generate Ghost points
    % [x_ghost, y_ghost] = fabric_pattern(N, 1, center(1), center(2), R);
    % [x_ghost, y_ghost] = uniform_circle(N, center(1), center(2), R);
    % ghost_pts =[x_ghost, y_ghost];
    ghost_pts = all_pts;

    DM_interior= distance_matrix(interior_pts, ghost_pts);
    DM_boundary= distance_matrix(boundary_pts, ghost_pts);
    DM_all= [DM_interior; DM_boundary];

    %% Loop for each shape parameter
    j=0;
    for shape=shapes
        j=j+1;
        basis_interior_pts = laplacian_rbf_2D(shape, DM_interior);
        basis_boundary_pts = nmq_rbf(shape, DM_boundary);

        % We have a linear system: -A α = λ B α
        A = [basis_interior_pts; basis_boundary_pts];

        B = [rbf(shape, DM_interior); zeros(Nb, N)];

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
        relative_error = (abs(approximate_eigenvalues - exact_eigenvalues))';

        result_eigenvalues(:,j, i) = approximate_eigenvalues;

        error_eigenvalues(:,j,i ) = relative_error;

        format short
        fprintf('The first %d exact eigenvalues/numerical eigenvalues/Relative errors are:\n',length(m))
        [pi^2*(m.^2.+l.^2)' approximate_eigenvalues (abs(approximate_eigenvalues'-pi^2*(m.^2.+l.^2)))']
    end
    %%


end

%% Print the time taken
cpu=etime(clock,t0); % count the time of the main part of the code

fprintf('Run time : %6.2f\n',cpu);

%% Now plot the maximum error

% Plot the errors of eigenvalues
% plot_eigenvalue_errors(nodes,shapes,error_eigenvalues,l,m);

% eigenvalue_Max_errors = zeros(length(shapes), length(radii));
% for i=1:length(radii)
%     for j=1:length(shapes)
%         % Find the maximum error for each shape parameter
%         eigenvalue_Max_errors(j, i) = max(error_eigenvalues(:,j,i));
%     end
% end

% for i=1:length(radii)
%     % Find the maximum error for each shape parameter

%     [min_error, min_idx] = min(eigenvalue_Max_errors(:,i));
%     best_shape = shapes(min_idx);

%     figure;
%     plot(shapes, eigenvalue_Max_errors(:,i)', 'o-');
%     hold on;

%     % Highlight the minimum error point with a distinct marker
%     plot(best_shape, min_error, 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Red circle

%     % Add a text label with an arrow
%     text(best_shape, min_error, ...
%         sprintf('\\leftarrow Minimum Error: %.2e\n(Shape = %.2f)', min_error, best_shape), ...
%         'VerticalAlignment', 'middle', ...
%         'HorizontalAlignment', 'right', ...
%         'FontSize', 10);
%     xlabel('Shape Parameter');
%     ylabel('Max Error');
%     title("Max Error of RADII = " + num2str(radii(i)));




%     % eigenvalue_Max_errors(:, i) = max(error_eigenvalues(:, :, i), [], 1);
% end


