function plot_eigenmode_abs_error(nodes,shapes,max_err_eigenmodes,l,m)
    color=['r','b','m','g','m','c','y','k','k'];

    % Plot the error of eigenvalues w.r.t. sp and nodes
    figure
    nrows = floor(length(shapes)/2) + mod(length(shapes),2);
    for j = 1:length(shapes)
            shape = shapes(j);
            subplot(nrows,2,j);
    
            for i = 1:length(nodes)
    
                    scatter3((nodes(i)+1)^2*ones(1,length(m)), 1:length(m), log10(squeeze(max_err_eigenmodes(:,j,i))), 'filled', color(i));
                    hold on;
                    % Use plot3 to connect the points with lines
                    plot3((nodes(i)+1)^2*ones(1,length(m)), 1:length(m), log10(squeeze(max_err_eigenmodes(:,j,i))), '-o', 'LineWidth', 1.5);
                    hold on;
            end
            % Labeling the axes
            xlabel('Number of nodes');
            ylabel('Indices of eigenvalues');
            zlabel('log10(Max Error of eigenmodes)');
            title("Max error of eigenmodes vs. the number of nodes (shape = " + num2str(shape) + ")");
            grid on;
    
    end
    hold off;
end