% This script approximates the eigenvalues and eigenmodes of the Laplace operator on a unit square using the Kansa method.

clear all;
warning off;
addpath('domain/');
addpath('matrix/');
addpath('rbf/');
addpath("plot/")
addpath("shapes/")

%% Define the parameters

shapes= [2.19]; % Shape parameter using in RBF.
% shapes = linspace(1,5,500);

nodes=[20]; % Grid sizes for the domain.

m=[1 2 1 2 3 1 3 2 4 ];
l=[1 1 2 2 1 3 2 3 1 ];

%% Store the results

result_eigenvalues = zeros(length(m), length(shapes), length(nodes)); % Approximate eigenvalues

error_eigenvalues = zeros(length(m), length(shapes), length(nodes)); % Absolute error in eigenvalues

max_err_eigenmodes = zeros(length(m), length(shapes), length(nodes)); % Maximum error in eigenmodes

relative_err_eigenmodes = zeros(length(m), length(shapes), length(nodes)); % Relative error in eigenmodes

rms_err_eigenmode = zeros(length(m), length(shapes), length(nodes)); % RMS error in eigenmodes

%% Start the timer
t0=clock;

%% Loop through the nodes

i=0; % Initialize the node counter
for nn=nodes
    i=i+1; % Increment the node counter
    j=0; % Initialize the shape parameter counter

    %% Define the domain
    square_domain = Square([0, 1]); % Create a square domain
    square_domain = square_domain.generateGrid(nn); % Generate a grid of size (N+1)^2

    coordinates = [square_domain.Coordinates, zeros(length(square_domain.Coordinates),1)]; % Add a third column for eigenmodes

    % square_domain.scatterPlot("Square Domain", true, true); % Plot the square domain

    % Store the interior points and their indices
    interior_pts = square_domain.InteriorPoints;
    interior_idx = square_domain.InteriorIndices;

    % Store the boundary points and their indices
    boundary_pts = square_domain.BoundaryPoints;
    boundary_idx = square_domain.BoundaryIndices;

    % Calculate the size of the set of points
    num_interior_pts = length(interior_pts(:,1));
    num_boundary_pts = length(boundary_pts(:,1));

    num_total_pts = num_interior_pts + num_boundary_pts;

    X = square_domain.Grid.X;
    Y = square_domain.Grid.Y;

    %% Loop through the shape parameters
    for shape = shapes
        j=j+1;

        %% Main computation
        DM = distance_matrix(coordinates(:,1:2), coordinates(:,1:2)); % Find the distance matrix

        % We have a linear system: -A α = λ B α
        A_i = laplacian_rbf_2D(shape, DM(interior_idx,:)); % Laplacian RBF for interior points

        B = nmq_rbf(shape, DM); % RBF for all points

        V = A_i / B; % V = A inv(B) Thus, we get the system: -V α = λ α

        V = V(:, interior_idx); % Make V Square. Extract the entries corresponding to interior collocation points

        % Note that we are using negative because our equation has it.
        [alpha, lambda] = eigs(-V, length(m), 0); % Calculate the eigenvalues.

        approximate_eigenvalues = real(diag(lambda)); % Extract the eigenvalues from the diagonal matrix

        alpha = real(alpha); % Extract the coefficient vector

        %% Calculate the eigenvalues
        exact_eigenvalues = pi^2*(m.^2.+l.^2)'; % Calculate the exact eigenvalues
        eigenvalue_abs_err = (abs(approximate_eigenvalues' - pi^2*(m.^2.+l.^2)))'; %

        % The first layer stores the ith node. Each column stores for each shape parameter.
        result_eigenvalues(:,j,i) = approximate_eigenvalues; % Store the approximate eigenvalues

        error_eigenvalues(:,j,i) = eigenvalue_abs_err; % Store error in eigenvalues

        format short
        fprintf('The first %d exact eigenvalues/numerical eigenvalues/Relative errors are:\n',length(m))
        [pi^2*(m.^2.+l.^2)' approximate_eigenvalues (abs(approximate_eigenvalues'-pi^2*(m.^2.+l.^2)))']


        %% Calculate Eigenmodes
        exact_eigenmode = sin(pi*coordinates(:,1)*m).*sin(pi*coordinates(:,2)*l);
        normf = sqrt(sum(exact_eigenmode.^2));

        % % Calculate relative, RMS, and absolute errors in eigenmodes
        for k=1:length(m)
            firsteigmode = exact_eigenmode(interior_idx,k); % kth eigenmode
            error_eigenmode_k = abs(abs(firsteigmode)-normf(k)*abs(alpha(1:num_interior_pts,k))); % Error in kth eigenmode

            max_err_eigenmodes(k, j, i) = max(error_eigenmode_k); % Store max error in kth eigenmode

            relative_err_eigenmodes(k, j, i) = max(error_eigenmode_k/max(firsteigmode)); % Store relative error in kth eigenmode

            rms_err_eigenmode(k,j,i) = sqrt(sum(error_eigenmode_k.^2)/num_interior_pts); % Store RMS error in kth eigenmode
        end

        %% Plot the exact, numerical, and absolute errors for eigenmodes.

        % plot_exact_eigenmodes(length(m), exact_eigenmode, X, Y, shape); % Plot the exact eigenmodes

        % plot_numeric_eigenmodes(length(m), coordinates, X, Y, interior_idx, num_interior_pts, normf, alpha, shape); % Plot the numerical eigenmodes

        % plot_abs_error(length(m), coordinates, X, Y, interior_idx, num_interior_pts, normf, alpha, exact_eigenmode, shape); % Plot the absolute errors of the eigenmodes

    end
end

%% Stop the timer
cpu=etime(clock,t0); % count the time of the main part of the code

fprintf('Run time : %6.2f\n',cpu);

%% Find the ideal shape parameter for eigenvalue.
[best_shapes_eigvalue, best_errors_eigvalue] = best(error_eigenvalues, "Shapes", shapes, "Nodes", nodes, true);

% [best_shapes_eigmodes, best_errors_eigmodes] = best(max_err_eigenmodes, "Shapes", shapes, "Nodes", nodes, true);


%% Plot the Errors

% plot_eigenvalue_errors(nodes,shapes,error_eigenvalues,l,m); %Plot the eigenvalues errors

% plot_eigenmode_abs_error(nodes,shapes,max_err_eigenmodes,l,m); % Plot the absolute error of eigenmodes

% plot_eigenmode_rel_error(nodes,shapes,relative_err_eigenmodes,l,m); % Plot the relative error of eigenmodes

% plot_eigenmode_rms_error(nodes,shapes,rms_err_eigenmode,l,m); % Plot the RMS error of eigenmodes

%%