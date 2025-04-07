function plot_numeric_eigenmodes(n, coordinates, X, Y, interior_idx, num_interior_pts, normf, alpha, shape, description)
%PLOT_NUMERIC_EIGENMODES Plots the computed numerical eigenmodes as surface plots.
%
%   This function takes the computed eigenmodes and visualizes them as
%   surface plots. It reshapes the eigenmode values onto a grid defined by X and Y
%   and then uses the `surf` function to create the plots.
%
% Inputs:
%   n: The number of eigenmodes to plot.
% All the inputs are the varialbles used in kansa method.
% Create a figure for the plots

if nargin < 10
    description = "Numerical Eigenmodes"; % Default description if not provided
end

figure;

% Loop through each eigenmode and plot it
for j=1:n
    % Store the j-th numerical eigenmode in the 3rd column of the coordinates matrix,
    % but only for the interior points.
    coordinates(interior_idx , 3) = normf(j)*alpha(1:num_interior_pts,j);

    % Reshape the values to match the grid for plotting
    Z = reshape(coordinates(:,3),size(X));

    % Create a subplot for the current eigenmode
    subplot(3,3,j)

    % Generate the surface plot
    surf(X,Y,Z);

    % Add plot configurations
    hold on;
    axis equal;
    axis([0,1,0,1]);
    title("Shape = " + shape);
end

sgtitle(description, 'FontSize', 14, 'FontWeight', 'bold');

end