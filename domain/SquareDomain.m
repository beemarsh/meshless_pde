classdef SquareDomain
    properties        
        SPoints        % (2 x N) Interior (collocation) points in the square
        BoundaryPoints % (2 x M) Points on the boundary of the square
        GhostPoints    % (2 x K) "Ghost" points outside the square
    end

    properties(SetAccess = private)
        Domain % [xmin, xmax, ymin, ymax] specifying the square domain
    end

    methods
        % Constructor
        function obj = SquareDomain(domain)
            % domain should be a 1x4 vector: [xmin, xmax, ymin, ymax]
            if nargin > 0
                validateattributes(domain, {'numeric'}, {'vector', 'numel', 4});
                obj.Domain = domain;
            end
        end
        
        % Generate interior (collocation) points in the square
        function obj = generateSquare(obj, numPoints)
            % numPoints: number of points in each direction
            xmin = obj.Domain(1); xmax = obj.Domain(2);
            ymin = obj.Domain(3); ymax = obj.Domain(4);
            % Create a uniform grid over the square
            x_domain = linspace(xmin, xmax, numPoints);
            y_domain = linspace(ymin, ymax, numPoints);
            [X, Y] = meshgrid(x_domain, y_domain);
            % Store as a 2 x N array
            obj.SPoints = [X(:)'; Y(:)'];
        end
        
        % Generate boundary points along the four edges of the square
        function obj = generateBoundaryPoints(obj, numPoints)
            % numPoints: number of points per edge
            xmin = obj.Domain(1); xmax = obj.Domain(2);
            ymin = obj.Domain(3); ymax = obj.Domain(4);
            
            % Top edge (y = ymax)
            x_top = linspace(xmin, xmax, numPoints);
            y_top = ymax * ones(1, numPoints);
            
            % Bottom edge (y = ymin)
            x_bottom = linspace(xmin, xmax, numPoints);
            y_bottom = ymin * ones(1, numPoints);
            
            % Left edge (x = xmin)
            y_left = linspace(ymin, ymax, numPoints);
            x_left = xmin * ones(1, numPoints);
            
            % Right edge (x = xmax)
            y_right = linspace(ymin, ymax, numPoints);
            x_right = xmax * ones(1, numPoints);
            
            % Combine all boundary points
            allPoints = [x_top, x_bottom, x_left, x_right;
                         y_top, y_bottom, y_left, y_right];
            
            % Remove duplicate points (corners appear twice)
            unique_points = unique(allPoints', 'rows', 'stable')';
            obj.BoundaryPoints = unique_points;
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
            scatter(obj.SPoints(1,:), obj.SPoints(2,:), 10, 'filled', 'b');
            hold on;
            % Plot boundary points in red (if requested)
            if nargin > 2 && showBoundary && ~isempty(obj.BoundaryPoints)
                scatter(obj.BoundaryPoints(1,:), obj.BoundaryPoints(2,:), 15, 'filled', 'r');
            end
            % Plot ghost points in green (if requested)
            if nargin > 3 && showGhost && ~isempty(obj.GhostPoints)
                scatter(obj.GhostPoints(1,:), obj.GhostPoints(2,:), 15, 'filled', 'g');
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
