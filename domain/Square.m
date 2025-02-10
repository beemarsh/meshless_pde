classdef Square
    properties        
        SPoints        % (2 x N) Interior (collocation) points in the square
        BoundaryPoints % (2 x M) Points on the boundary of the square
        GhostPoints    % (2 x K) "Ghost" points outside the square
        Grid           % Structure containing grid information for surface plotting
    end

    properties(SetAccess = private)
        Domain % [xmin, xmax, ymin, ymax] specifying the square domain
    end

    methods
        % Constructor
        function obj = Square(domain)
            % domain should be a 1x4 vector: [xmin, xmax, ymin, ymax]
            if nargin > 0
                validateattributes(domain, {'numeric'}, {'vector', 'numel', 4});
                obj.Domain = domain;
            end
        end
        
        % Generate interior (collocation) points in the square
                % Generate interior (collocation) points in the square
                function obj = generateGrid(obj, N)
            % N: number of intervals in each coordinate direction,
            % so there will be N+1 grid points per direction.
            xmin = obj.Domain(1); xmax = obj.Domain(2);
            ymin = obj.Domain(3); ymax = obj.Domain(4);
            
            % Create a grid with N+1 points along each axis.
            x_vals = linspace(xmin, xmax, N+1);
            y_vals = linspace(ymin, ymax, N+1);
            [X, Y] = meshgrid(x_vals, y_vals);
            
            % Save the full grid for plotting
            obj.Grid.Full.X = X;
            obj.Grid.Full.Y = Y;
            
            % Combine the grid points into a 2 x ((N+1)^2) array.
            % Each column is a point [x; y].
            allPoints = [X(:)'; Y(:)'];
            
            % Identify the boundary points.
            % A point is on the boundary if its x value equals xmin or xmax,
            % or its y value equals ymin or ymax (allowing for numerical tolerance).
            tol = 1e-12;
            isBoundary = (abs(allPoints(1,:) - xmin) < tol) | (abs(allPoints(1,:) - xmax) < tol) | ...
                         (abs(allPoints(2,:) - ymin) < tol) | (abs(allPoints(2,:) - ymax) < tol);
            
            % Separate the points.
            obj.BoundaryPoints = allPoints(:, isBoundary)';
            obj.SPoints = allPoints(:, ~isBoundary)';
        end

        
        
        % Generate ghost points outside the square domain
        function obj = generateGhostPoints(obj, numPoints, width, spread)
            % numPoints: number of ghost points per edge per layer
            % width: number of ghost layers (lines offset from the boundary)
            % spread: the offset distance added per ghost layer
            
            xmin = obj.Domain(1); xmax = obj.Domain(2);
            ymin = obj.Domain(3); ymax = obj.Domain(4);
            ghostArray = [];
            
            for i = 1:width
                offset = i * spread;
                % Top ghost layer: offset upward from top boundary
                x_top = linspace(xmin - offset, xmax + offset, numPoints);
                y_top = (ymax + offset) * ones(1, numPoints);
                
                % Bottom ghost layer: offset downward from bottom boundary
                x_bottom = linspace(xmin - offset, xmax + offset, numPoints);
                y_bottom = (ymin - offset) * ones(1, numPoints);
                
                % Left ghost layer: offset leftward from left boundary
                y_left = linspace(ymin - offset, ymax + offset, numPoints);
                x_left = (xmin - offset) * ones(1, numPoints);
                
                % Right ghost layer: offset rightward from right boundary
                y_right = linspace(ymin - offset, ymax + offset, numPoints);
                x_right = (xmax + offset) * ones(1, numPoints);
                
                % Concatenate ghost points for this layer
                ghostLayer = [ [x_top, x_bottom, x_left, x_right]', ...
                               [y_top, y_bottom, y_left, y_right]' ];
                ghostArray = [ghostArray; ghostLayer];
            end
            
            % Remove duplicate rows, if any
            ghostArray = unique(ghostArray, 'rows', 'stable');
            % Store as 2 x K array
            obj.GhostPoints = ghostArray';
        end
        
        % Scatter plot of interior, boundary, and ghost points
        function scatterPlot(obj, plot_title, showBoundary, showGhost)
            figure;
            % Plot interior (collocation) points in blue
            scatter(obj.SPoints(:,1), obj.SPoints(:,2), 10, 'filled', 'b');
            hold on;
            % Plot boundary points in red (if requested)
            if nargin > 2 && showBoundary && ~isempty(obj.BoundaryPoints)
                scatter(obj.BoundaryPoints(:,1), obj.BoundaryPoints(:,2), 15, 'filled', 'r');
            end
            % Plot ghost points in green (if requested)
            if nargin > 3 && showGhost && ~isempty(obj.GhostPoints)
                scatter(obj.GhostPoints(:,1), obj.GhostPoints(:,2), 15, 'filled', 'g');
            end
            hold off;
            axis equal;
            xlabel('X');
            ylabel('Y');
            title(plot_title);
            legendEntries = {'Interior Points'};
            if nargin > 2 && showBoundary
                legendEntries{end+1} = 'Boundary Points';
            end
            if nargin > 3 && showGhost
                legendEntries{end+1} = 'Ghost Points';
            end
            legend(legendEntries, 'Location', 'best');
        end
    end
end
