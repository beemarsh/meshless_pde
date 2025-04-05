% Kansas method on a square domain.
addpath('./domain/');
addpath('./functions/');
addpath('./plot/')

% We use different shape parameters
% shapes = [ 11 10 9 8];
% shapes = linspace(0.1 , 50, 1000);
shapes= [9];
% We use different number of nodes
nodes=[19];

m=[1 2 1 2 3 1 3 2 4 ];
l=[1 1 2 2 1 3 2 3 1 ];

% This function D calculates the distance between each points in the matrix x and y.
% Here, x is a n*2 matrix and y is a m*2 matrix.
% Each row denotes a set of x,y collocation points
% Using bsxfun, we calculate the difference between each points.
% Using hypot, we calculate the Euclidean distance.
% Essentially, all its doing is: sqrt( (x_1 - x_2)^2 + (y_1 - y_2)^2 )
D = @(x,y) hypot(bsxfun(@minus,x(:,1),y(:,1)'),bsxfun(@minus,x(:,2),y(:,2)'));

% Store the result for eigenvalues in 3 columns:
% Each layer stores the approximate eigenvalues for each size of the node.
% Each column stores the approximate eigenvalues for each shape parameter.
% The first layer stores the ith node. Each column stores for each shape parameter.
result_eigenvalues = zeros(length(m), length(shapes), length(nodes));

% Store the error in eigenvalues.
error_eigenvalues = zeros(length(m), length(shapes), length(nodes));

% Store results for eigenmode.
max_err_eigenmodes = zeros(length(m), length(shapes), length(nodes)); % Store maximum error

relative_err_eigenmodes = zeros(length(m), length(shapes), length(nodes)); % Store relative error

rms_err_eigenmode = zeros(length(m), length(shapes), length(nodes)); % Store RMS error

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

        % square_domain.scatterPlot("Square Domain", true, true);

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

        for shape = shapes

                j=j+1;

                % We find the distance matrix
                DM = D(coordinates(:,1:2), coordinates(:,1:2));

                % Calculate RBF MQ = sqrt (1 + ξ^2 * r^2 )
                % Only for collocation points
                MQ = sqrt( ( shape^2 * DM(interior_idx,:).^2 ) + 1 );

                % Now, RBF for all points
                mq = sqrt( ( shape^2 * DM.^2 ) + 1 );

                % Now we take the laplacian of MQ RBF.
                L_MQ = ( 2 * (shape ^ 2)./ MQ ) - ( (shape ^ 4) * ( DM(interior_idx,:).^2)./(MQ.^3) );

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
                relative_error = (abs(approximate_eigenvalues' - pi^2*(m.^2.+l.^2)))';

                format short
                fprintf('The first %d exact eigenvalues/numerical eigenvalues/Relative errors are:\n',length(m))
                [pi^2*(m.^2.+l.^2)' approximate_eigenvalues (abs(approximate_eigenvalues'-pi^2*(m.^2.+l.^2)))']

                % Store approximate eigenvalues, exact eigenvalues and relative error in eigenvalues.
                % The first layer stores the ith node. Each column stores for each shape parameter.
                result_eigenvalues(:,j,i) = approximate_eigenvalues;

                % Same way store the errors
                error_eigenvalues(:,j,i) = relative_error;

                % Compute the exact eigenmodes
                % exact_eigenmode = sin(pi*coordinates(:,1)*m).*sin(pi*coordinates(:,2)*l);
                % normf = sqrt(sum(exact_eigenmode.^2));

                % % After computing eigenmodes, we will calculate errors: relative, max,
                % % and RMS.
                % for k=1:length(m)
                %         firsteigmode = exact_eigenmode(interior_idx,k);
                %         error_eigenmode_k = abs(abs(firsteigmode)-normf(k)*abs(alpha(1:num_interior_pts,k)));

                %         max_err_eigenmodes(k, j, i) = max(error_eigenmode_k); % Just like the eigenvalues, we store the maximum error in a 3D matrix.

                %         relative_err_eigenmodes(k, j, i) = max(error_eigenmode_k/max(firsteigmode));

                %         rms_err_eigenmode(k,j,i) = sqrt(sum(error_eigenmode_k.^2)/num_interior_pts);
                % end


                % Plot the exact eigenmodes
                % plot_exact_eigenmodes(9, exact_eigenmode, X, Y, shape);

                % Plot the numerical eigenmodes
                % plot_numeric_eigenmodes(9, coordinates, X, Y, interior_idx, num_interior_pts, normf, alpha, shape);

                % Plot the absolute errors of the eigenmodes
                % plot_abs_error(9, coordinates, X, Y, interior_idx, num_interior_pts, normf, alpha, exact_eigenmode, shape);

        end
end

cpu=etime(clock,t0); % count the time of the main part of the code

fprintf('Run time : %6.2f\n',cpu);

%Plot the eigenvalues errors
% plot_eigenvalue_errors(nodes,shapes,error_eigenvalues,l,m);


% eigenvalue_Max_errors = zeros(1, length(shapes));
% for i=1:length(shapes)
%     eigenvalue_Max_errors(i) = max(error_eigenvalues(:,i,1));
% end

% [min_error, min_idx] = min(eigenvalue_Max_errors);
% best_shape = shapes(min_idx);

% figure;
% plot(shapes, eigenvalue_Max_errors, 'o-');
% hold on;

% % Highlight the minimum error point with a distinct marker
% plot(best_shape, min_error, 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Red circle

% % Add a text label with an arrow
% text(best_shape, min_error, ...
%     sprintf('\\leftarrow Minimum Error: %.2e\n(Shape = %.2f)', min_error, best_shape), ...
%     'VerticalAlignment', 'middle', ...
%     'HorizontalAlignment', 'right', ...
%     'FontSize', 10);
% xlabel('Shape Parameter');
% ylabel('Max Error');
% title('Max Error of USING GHOST POINTS(FIXED SHAPE PARAMETER)');

% Plot the absolute error of eigenmodes
% plot_eigenmode_abs_error(nodes,shapes,max_err_eigenmodes,l,m);

% Plot the relative error of eigenmodes
% plot_eigenmode_rel_error(nodes,shapes,relative_err_eigenmodes,l,m);

% Plot the RMS error of eigenmodes
% plot_eigenmode_rms_error(nodes,shapes,rms_err_eigenmode,l,m);

% save_matrix(log10(error_eigenvalues), nodes, shapes, 'error_eigenvalues_kansa');