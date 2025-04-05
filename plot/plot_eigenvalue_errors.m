function plot_eigenvalue_errors(nodes,shapes,error_eigenvalues,l,m)
    color=['r','b','m','g','m','c','y','k','k'];

    % Plot the error of eigenvalues w.r.t. sp and nodes
    figure
    nrows = floor(length(shapes)/2) + mod(length(shapes),2);
    for j = 1:length(shapes)
            shape = shapes(j);
            subplot(nrows,2,j);
    
            for i = 1:length(nodes)
    
                    scatter3((nodes(i)+1)^2*ones(1,length(m)), 1:length(m), (squeeze(error_eigenvalues(:,j,i))), 'filled', color(i));
                    hold on;
                    % Use plot3 to connect the points with lines
                    plot3((nodes(i)+1)^2*ones(1,length(m)), 1:length(m), (squeeze(error_eigenvalues(:,j,i))), '-o', 'LineWidth', 1.5);
                    hold on;
            end
            % Labeling the axes
            xlabel('Number of nodes');
            ylabel('Indices of eigenvalues');
            zlabel('(Errors of eigenvalues)');
            title("Error of eigenvalues vs. the Num of nodes (shape = " + num2str(shape) + ")");
            grid on;
    
    end
    hold off;
end