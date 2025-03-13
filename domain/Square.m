classdef Square
    properties        
        SPoints        % (2 x N) Interior (collocation) points in the square
        BoundaryPoints % (2 x M) Points on the boundary of the square
        GhostPoints    % (2 x K) "Ghost" points outside the square
        Grid           % Structure containing grid information for surface plotting
        GhostGrid
        Size            % Size of the grid.
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
              obj.Size = N;     %assign
            % N: number of intervals in each coordinate direction,
            % so there will be N+1 grid points per direction.
            xmin = obj.Domain(1); xmax = obj.Domain(2);
            ymin = obj.Domain(3); ymax = obj.Domain(4);
            
            % Create a grid with N+1 points along each axis.
            x_vals = linspace(xmin, xmax, N+1);
            y_vals = linspace(ymin, ymax, N+1);
            [X, Y] = meshgrid(x_vals, y_vals);
            
            % Save the full grid for plotting
            obj.Grid.X = X;
            obj.Grid.Y = Y;
            
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
        function obj = generateGhostPoints(obj, width, spread)
    % Extract domain boundaries
    xmin = obj.Domain(1);
    xmax = obj.Domain(2);
    ymin = obj.Domain(3);
    ymax = obj.Domain(4);
    
    % Generate original grid points (from generateGrid)
    N = obj.Size;
    x_vals = linspace(xmin, xmax, N+1);
    y_vals = linspace(ymin, ymax, N+1);
    
    % Generate ghost layers for x-direction
    left_x = xmin - spread * (1:width);
    right_x = xmax + spread * (1:width);
    x_extended = [left_x, x_vals, right_x];
    
    % Generate ghost layers for y-direction
    bottom_y = ymin - spread * (1:width);
    top_y = ymax + spread * (1:width);
    y_extended = [bottom_y, y_vals, top_y];
    
    % Create meshgrid for the extended domain
    [X, Y] = meshgrid(x_extended, y_extended);
    
    % Save the ghost grid for plotting
    obj.GhostGrid.X = X;
    obj.GhostGrid.Y = Y;
    
    % Combine all extended points into a 2xK array
    allGhostPoints = [X(:)'; Y(:)'];
    
    % Identify points outside the original domain (with tolerance)
    tol = 1e-12;
    isGhost = (allGhostPoints(1,:) < xmin - tol) | (allGhostPoints(1,:) > xmax + tol) | ...
              (allGhostPoints(2,:) < ymin - tol) | (allGhostPoints(2,:) > ymax + tol);
    
    % Assign ghost points
    obj.GhostPoints = allGhostPoints(:, isGhost)';
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
