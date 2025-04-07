function plot_abs_error(n , coordinates, X, Y, interior_idx, num_interior_pts, normf, alpha, exact_eigenmode, shape, description)
% PLOT_ABS_ERROR Plots the absolute contour errors of eigenmodes.
%
% Inputs:
% n: The number of eigenmodes to plot.
% All the inputs are the varialbles used in kansa method.
if nargin < 11
    description = "Numerical Eigenmodes"; % Default description if not provided
end
% Create a new figure for the plots
figure;

% Loop through each eigenmode and plot the absolute error
for j=1:n
    % Calculate the absolute difference between the exact eigenmode and the computed eigenmode at the interior points.
    coordinates(interior_idx, 3) = abs(abs(exact_eigenmode(interior_idx,j))-abs(normf(j)*alpha(1:num_interior_pts,j)));

    % Reshape the error values to match the grid dimensions for plotting.
    Z=reshape(coordinates(:,3),size(X));

    % Create a subplot within the figure for the current eigenmode.
    subplot(3,3,j)

    % Generate a surface plot of the absolute error.
    surf(X,Y,Z);

    % Hold the current plot to allow for potential overlays.
    hold on;

    % Set the axis to be equal.
    axis equal;

    % Set the axis limits to [0, 1] for both x and y.
    axis([0,1,0,1]);

    % Set the title of the subplot to display the shape parameter used.
    title("Shape = " + shape);
end

sgtitle(description, 'FontSize', 14, 'FontWeight', 'bold');
end