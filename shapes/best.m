function [best_values, min_errors] = best(errors, p1_name, p1, p2_name, p2, showPlot)
% errors is a 3D matrix. Each row contains the eigenvalues or eigenmodes. Each column contains the value for each p1. p1 can be shape or radius or whatever. Each layer contains the values of each layer. p2 can be number of nodes or radius.
% For example: p1 can be shapes and p2 can be nodes

if nargin < 6
    showPlot = false;
end

best_values = zeros(1, length(p2));
min_errors = zeros(1, length(p2));

for i=1:length(p2)

    max_errors_i = zeros(1, length(p1));
    for j=1:length(p1)
        max_errors_i(j) = max(errors(:,j,i));
    end

    [min_error, min_idx] = min(max_errors_i);
    best_value_i = p1(min_idx);

    best_values(i) = best_value_i;
    min_errors(i) = min_error;

    if showPlot
        figure;
        plot(p1, max_errors_i, 'o-');
        hold on;

        % Highlight the minimum error point with a distinct marker
        plot(best_value_i, min_error, 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Red circle

        % Add a text label with an arrow
        text(best_value_i, min_error,    sprintf('Minimum Error: %.2e\n(Shape = %.2f)', min_error, best_value_i), 'FontSize', 10);
        xlabel(p1_name);
        ylabel('Max Error');
        title("Max Error (" + p2_name + " = " + p2(i) + " )");
    end
end

end