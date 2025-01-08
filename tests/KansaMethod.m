% Generate a Lshaped Domain

new_L = LDomain([-3,3],[-3,0],[-3,0], 'bottom-left');

new_L = new_L.generateLShape(40);

new_L = new_L.generateBoundaryPoints(40);

new_L = new_L.generateGhostPoints(10, 2, 0.1);

new_L.scatterPlot("L-Shaped Domain",true, false);

%new_L.scatterPlot("L-Shaped Domain",true, true);

X_interior = new_L.LPoints(1,:); % Interior Collocation Points
Y_interior = new_L.LPoints(2,:);

X_interior = X_interior(:); Y_interior = Y_interior(:); % Convert to column vectors

X_boundary = new_L.BoundaryPoints(1,:); %Boundary collocation points
Y_boundary = new_L.BoundaryPoints(2,:);
X_boundary = X_boundary(:); Y_boundary = Y_boundary(:); % Convert to column vectors

% Shape parameter
c = 0.81;

function [lambda, alpha] = kansa_eigenvalues(Xint, Yint, Xbd, Ybd, c)
    % Combine all points into a single list
    X = [Xint; Xbd];  
    Y = [Yint; Ybd];
    N = length(X);  % Total number of collocation points
    Ni = length(Xint);  % Number of interior points
    Nb = length(Xbd);  % Number of boundary points

    % Compute pairwise distances between points
    [X1, X2] = meshgrid(X, X);
    [Y1, Y2] = meshgrid(Y, Y);
    R = sqrt((X1 - X2).^2 + (Y1 - Y2).^2);

    % Define the Multiquadric (MQ) RBF and its Laplacian
    Phi = sqrt(R.^2 + c^2);
    Phi_lap = (2 ./ sqrt(R.^2 + c^2)) - (R.^2 ./ (R.^2 + c^2).^(3/2));

    % Construct A and B matrices
    A = zeros(N, N);
    B = zeros(N, N);

    % Apply Laplacian operator to the first Ni interior points
    A(1:Ni, :) = Phi_lap(1:Ni, :);
    B(1:Ni, :) = Phi(1:Ni, :);

    % Apply Dirichlet boundary conditions (u = 0 at boundary)
    A(Ni+1:end, :) = Phi(Ni+1:end, :);
    B(Ni+1:end, :) = zeros(Nb, N);

    % Solve the generalized eigenvalue problem
    [V, D] = eig(A, B);
    lambda = diag(D); % Eigenvalues
    alpha = V; %Corresponding eigenvectors

end

% lambda = eigenvalue, alpha = eigen funciton
[lambda, alpha] = kansa_eigenvalues(X_interior, Y_interior, X_boundary, Y_boundary, c);

% 
% % Print eigenvalues
% disp('Eigenvalues:');
disp(lambda);
% 
% % Print eigenfunctions
% disp('Eigenfunctions (First Few Entries for Each):');
% num_eigenfunctions = size(alpha, 2); % Number of eigenfunctions
% for i = 1:num_eigenfunctions
%     disp(['Eigenfunction ', num2str(i), ':']);
%     disp(alpha(1:10, i)); % Print the first 10 entries of each eigenfunction
% end


% Generate grid for interpolation
[xq, yq] = meshgrid(linspace(-3, 3, 100), linspace(-3, 3, 100));

L_shape_mask = (xq >= -3 & xq <= 3 & yq >= -3 & yq <= 3) & ...
               ~((xq <= 0 & yq <= 0)); % Exclude bottom-left quadrant


% Number of eigenfunctions to plot
num_eigenfunctions = 8;

figure;
for i = 1:num_eigenfunctions
    % Extract real part of eigenfunction
    eigenfunction_values = real(alpha(1:length(X_interior), i)); 

    % Interpolate eigenfunction values
    eigenfunction_grid = griddata(X_interior, Y_interior, eigenfunction_values, xq, yq, 'cubic');

    % Apply the mask
    eigenfunction_grid(~L_shape_mask) = NaN;

    % Plot contour
    subplot(4, 2, i); % Arrange in a 4x2 grid
    contourf(xq, yq, eigenfunction_grid, 200, 'LineColor', 'none'); % 200 levels
    colorbar;
    title(['Eigenfunction ', num2str(i), ' (\lambda = ', num2str(real(lambda(i))), ')']);
    xlabel('X');
    ylabel('Y');
    axis equal;
    grid on;
end
