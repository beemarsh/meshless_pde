% Kansas method on a square domain.
addpath('./domain/');
addpath('./functions/');


% We use different shape parameters
% shapes = [3 2 1 0.1 0.01 0.001 0.0001 0.00001 0.000001,0.0000001];

% We use different number of nodes
nodes=[19];

m=[1 2 1 2 3 1 3 2 4 ];
l=[1 1 2 2 1 3 2 3 1 ];

delta = 0.3;
epsilon = 0.5;

% Radius of the spiral ghost points.
radii = [0.5, 1 , 1.5, 2, 2.5, 3];

% Store the result for eigenvalues in 3 columns:
% Each layer stores the approximate eigenvalues for each size of the node.
% Each column stores the approximate eigenvalues for each shape parameter.
% The first layer stores the ith node. Each column stores for each shape parameter.
result_eigenvalues = zeros(length(m), length(radii), length(nodes));

% Store the error in eigenvalues.
error_eigenvalues = zeros(length(m), length(radii), length(nodes));

% Store results for eigenmode.
max_err_eigenmodes = zeros(length(m), length(radii), length(nodes)); % Store maximum error

relative_err_eigenmodes = zeros(length(m), length(radii), length(nodes)); % Store relative error

rms_err_eigenmode = zeros(length(m), length(radii), length(nodes)); % Store RMS error


%count the time of the main part of the code
t0=clock;


%Keep track of nodes
i=0;
for nn=nodes
    i=i+1;
    j=0;
    %First define the domain
    % Unit square domain
    square_domain = Square([0, 1]);
    % Generate a grid of size N^2
    square_domain = square_domain.generateGrid(nn);

    % Store the coordinates.
    % Here coordinates have two columns, X and Y values.
    % Now add a third column for storing the eigenmodes associated with the point. It is filled with 1s.
    coordinates = [square_domain.Coordinates, zeros(length(square_domain.Coordinates),1)];

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

    % Diameter and Center of the smallest circle enclosing the interpolation points
    [ center, r] = smallest_circle(boundary_pts);

    % Calculate the shape parameter using Modified Franke
    c_f = franke(num_total_pts, 2*r);

    % Now get the variable shapes.
    shapes = variable_shape(num_total_pts,c_f, delta, epsilon);

    for R = radii
        j=j+1; % Keep track of the index of the shape parameter

        square_domain = square_domain.generateGhostPoints(1, 0.5,0.5, R);

        % square_domain.scatterPlot("Square Domain", true, true);

        % Store ghost points
        ghost_pts = square_domain.GhostPoints;

        % We find the distance matrix
        DM = D(coordinates(:,1:2), ghost_pts(:,1:2));
        % shapes = reshape(var_shapes, size(DM));

        % Calculate RBF MQ = sqrt (1 + ξ^2 * r^2 )
        % Only for collocation points
        MQ = rbf(shapes, DM(interior_idx,:));

        % Now, RBF for all points
        mq = rbf(shapes, DM);

        % Now we take the laplacian of MQ RBF.
        L_MQ = ( 2 * (shapes.^2)./ MQ ) - ( (shapes .^ 4) .* ( DM(interior_idx,:).^2)./(MQ.^3) );

        % We have a linear system: -A α = λ B α
        % We recast this form into standard eigenvalue problem.
        % Let V = A inv(B)
        % Thus, we get the system: -V α = λ α
        V = L_MQ / mq;

        % To calculate eigenvalues, first we make the matrix V square.
        % For that, we extract the entries corresponding to interior
        % collocation points
        V = V(:, interior_idx);

        % Now, we calculate the eigenvalues.
        % Note that we are using negative because our equation has it.
        [alpha, lambda] = eigs(-V, length(m), 0);

        % extract the first length(m) eigenvalues;
        % We need to extract them from the diagonal.

        approximate_eigenvalues = real(diag(lambda));

        alpha = real(alpha);

        % Compute the exact eigenvalues and then calculate the relative error
        exact_eigenvalues = pi^2*(m.^2.+l.^2)';
        relative_error = (abs(approximate_eigenvalues' - pi^2*(m.^2.+l.^2))./(pi^2*(m.^2.+l.^2)))';

        % Store approximate eigenvalues, exact eigenvalues and relative error in eigenvalues.
        % The first layer stores the ith node. Each column stores for each shape parameter.
        result_eigenvalues(:,j,i) = approximate_eigenvalues;

        % Same way store the errors
        error_eigenvalues(:,j,i) = relative_error;

        % Compute the exact eigenmodes
        exact_eigenmode = sin(pi*coordinates(:,1)*m).*sin(pi*coordinates(:,2)*l);
        normf = sqrt(sum(exact_eigenmode.^2));
        % After computing eigenmodes, we will calculate errors: relative, max,
        % and RMS.
        for k=1:length(m)
            firsteigmode = exact_eigenmode(interior_idx,k);
            error_eigenmode_k = abs(abs(firsteigmode)-normf(k)*abs(alpha(1:num_interior_pts,k)));

            max_err_eigenmodes(k, j, i) = max(error_eigenmode_k); % Just like the eigenvalues, we store the maximum error in a 3D matrix.

            relative_err_eigenmodes(k, j, i) = max(error_eigenmode_k/max(firsteigmode));

            rms_err_eigenmode(k,j,i) = sqrt(sum(error_eigenmode_k.^2)/num_interior_pts);
        end


        % Plot the exact eigenmodes
        % plot_exact_eigenmodes(length(m), exact_eigenmode, X, Y, R);

        % Plot the numerical eigenmodes
        % plot_numeric_eigenmodes(length(m), coordinates, X, Y, interior_idx, num_interior_pts, normf, alpha, R);

        % Plot the absolute errors of the eigenmodes
        % plot_abs_error(length(m), coordinates, X, Y, interior_idx, num_interior_pts, normf, alpha, exact_eigenmode, shape);

    end

end

cpu=etime(clock,t0); % count the time of the main part of the code

fprintf('Run time : %6.2f\n',cpu);

% Plot the errors of eigenvalues
plot_eigenvalue_errors(nodes,radii,error_eigenvalues,l,m);

% Plot the absolute error of eigenmodes
% plot_eigenmode_abs_error(nodes,shapes,max_err_eigenmodes,l,m);

% Plot the relative error of eigenmodes
% plot_eigenmode_rel_error(nodes,shapes,relative_err_eigenmodes,l,m);

% Plot the RMS error of eigenmodes
% plot_eigenmode_rms_error(nodes,shapes,rms_err_eigenmode,l,m);

% save_matrix(log10(error_eigenvalues), nodes, shapes, 'error_eigenvalues_kansa_ghost');