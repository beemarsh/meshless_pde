function plot_exact_eigenmodes(n, exact_eigenmode, X, Y, shape, description)
%PLOT_EXACT_EIGENMODES Plots the exact eigenmodes as surface plots.
%
%   This function visualizes the exact eigenmodes by plotting them as surface plots
%   on a grid defined by X and Y.
%
%   % Inputs:
% n: The number of eigenmodes to plot.
% All the inputs are the varialbles used in kansa method.

% Create a new figure for the plots

if nargin < 6
    description = "Exact Eigenmodes"; % Default description if not provided
end

figure

% Loop through the first 9 eigenmodes and plot them
for j=1:n
    % Get the jth eigenmode
    j_th_eigenmode = exact_eigenmode(:,j);

    % Reshape the matrix to get the Z values for the surface plot
    Z = reshape(j_th_eigenmode,size(X));

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