% Kansas method on a square domain.
addpath('./domain/');
addpath('./functions/');

%First define the domain
% Unit square domain
square_domain = Square([0, 1]);

% We use different shape parameters
shapes = [9 10 11];

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
% First: approximate eigenvalues
% Second: exact eigenvalues
% Third:  error in eigenvalues
result_eigenvalues = zeros([length(m), 3]);

% Store results for eigenmode.
max_err_eigenmodes = zeros(length(shapes), length(m)); % Store maximum error
relative_err_eigenmodes = zeros(length(shapes), length(m)); % Store relative error
rms_err_eigenmode = zeros(length(shapes), length(m)); % Store RMS error

% Generate a grid of size N^2
square_domain = square_domain.generateGrid(20);
square_domain.scatterPlot("Square Domain", true, false);

% Store the coordinates.
coordinates = square_domain.Coordinates; % Here coordinates have two columns, X and Y values.
% Now add a third column for storing the eigenmodes associated with the point.
coordinates = [coordinates, zeros(length(coordinates),1)];

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


% Keeps count when looping
count = 0;
% For each different shape parameter, we approximate values.

%shape = shapes(1);
for shape = shapes

count = count+1;

% We find the distance matrix
DM = D(coordinates(:,1:2), coordinates(:,1:2));

% Calculate RBF MQ = sqrt (1 + ξ^2 * r^2 )
% Only for collocation points
MQ =  (1/(9 * (shape ^ 2))) .* ( ( (4 + (shape * DM(interior_idx,:)).^2  ) .* (sqrt(1 + (shape * DM(interior_idx,:)).^2 ))) - 3 * log(sqrt(1 + (shape * DM(interior_idx,:)).^2 ) + 1));


% Now, RBF for all points
mq =  (1/(9 * (shape ^ 2))) * ( ( (4 + (shape * DM).^2  ) * (sqrt(1 + (shape * DM).^2 ))) - 3 * log(sqrt(1 + (shape * DM).^2 ) + 1));


% Now we take the laplacian of MQ RBF.
L_MQ = sqrt( ( shape^2 * DM(interior_idx,:).^2 ) + 1 );

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

% Store approximate eigenvalues, exact eigenvalues and relative error
% in eigenvalues.
result_eigenvalues(:,1:3)=[approximate_eigenvalues, exact_eigenvalues, relative_error];

% Compute the exact eigenmodes
exact_eigenmode = sin(pi*coordinates(:,1)*m).*sin(pi*coordinates(:,2)*l);
normf = sqrt(sum(exact_eigenmode.^2));

% After computing eigenmodes, we will calculate errors: relative, max,
% and RMS.
for k=1:length(m)
        firsteigmode = exact_eigenmode(interior_idx,k);
        error_eigenmode_k = abs(abs(firsteigmode)-normf(k)*abs(alpha(1:num_interior_pts,k)));

        max_err_eigenmodes(count, k) = max(error_eigenmode_k);

        relative_err_eigenmodes(count, k) = max(error_eigenmode_k/max(firsteigmode));

        rms_err_eigenmode(count, k) = sqrt(sum(error_eigenmode_k.^2)/num_interior_pts);
end


% Plot the exact eigenmodes
plot_exact_eigenmodes(9, exact_eigenmode, X, Y, shape);

% Plot the numerical eigenmodes
plot_numeric_eigenmodes(9, coordinates, X, Y, interior_idx, num_interior_pts, normf, alpha, shape);

% Plot the absolute errors of the eigenmodes
plot_abs_error(9, coordinates, X, Y, interior_idx, num_interior_pts, normf, alpha, exact_eigenmode, shape);

end