classdef CircularDomain
    %CDomain A class for creating and visualizing a 2D circular domain
    %   Example usage:
    %       circleObj = CDomain([0, 0], 1);       % Circle centered at (0,0) with radius 1
    %       circleObj = circleObj.generateCDomain(50);     % Generate ~50x50 interior grid
    %       circleObj = circleObj.generateBoundaryPoints(50); % 50 boundary points
    %       circleObj.scatterPlot('My Circle', true);      % Plot interior+boundary

    properties
        CPoints        % (N x 2) The (x,y) points that lie inside the circular domain
        BoundaryPoints % (2 x M) The (x, y) points that lie on the boundary of the circle
        GhostPoints    % Not used here, but shown for consistency
    end

    properties(SetAccess = private)
        Center % [cx, cy] center of the circle
        Radius % Scalar radius
    end

    methods
        % Constructor
        function obj = CircularDomain(center, radius)
            if nargin > 0
                obj.Center = center;
                obj.Radius = radius;
            end
        end

        % Generate interior points by discretizing a bounding box
        % and keeping those inside the circle
        function obj = generateCDomain(obj, numPoints)
            % 1) Compute bounding box around the circle
            x_min = obj.Center(1) - obj.Radius;
            x_max = obj.Center(1) + obj.Radius;
            y_min = obj.Center(2) - obj.Radius;
            y_max = obj.Center(2) + obj.Radius;

            % 2) Create a uniform grid in the bounding box
            x_domain = linspace(x_min, x_max, numPoints);
            y_domain = linspace(y_min, y_max, numPoints);
            [X, Y] = meshgrid(x_domain, y_domain);

            % 3) Keep only points inside the circle: (x-cx)^2 + (y-cy)^2 <= r^2
            dist_sq = (X - obj.Center(1)).^2 + (Y - obj.Center(2)).^2;
            mask = dist_sq <= obj.Radius^2;

            % 4) Extract (x, y) coordinates for those points
            X_in = X(mask);
            Y_in = Y(mask);

            % 5) Store them in CPoints as an N x 2 array
            obj.CPoints = [X_in, Y_in];
        end

        % Generate boundary points around the circle
        function obj = generateBoundaryPoints(obj, numPoints)
            % 1) Parametric angle from 0 to 2*pi
            theta = linspace(0, 2*pi, numPoints);

            % 2) Evaluate the circle parametric equations
            x_b = obj.Center(1) + obj.Radius * cos(theta);
            y_b = obj.Center(2) + obj.Radius * sin(theta);

            % 3) Store as a 2 x M matrix for consistency with your LDomain style
            obj.BoundaryPoints = [x_b; y_b];
        end

        % Scatter plot of the interior and (optionally) the boundary
        function scatterPlot(obj, plot_title, showBoundary)
            figure;
            % 1) Scatter the interior points (if they exist)
            if ~isempty(obj.CPoints)
                scatter(obj.CPoints(:,1), obj.CPoints(:,2), 10, 'filled','b');
                hold on;
            end

            % 2) Scatter the boundary points in red
            if nargin > 2 && showBoundary && ~isempty(obj.BoundaryPoints)
                scatter(obj.BoundaryPoints(1,:), obj.BoundaryPoints(2,:), 15, 'filled', 'r');
            end

            axis equal;
            xlabel('X');
            ylabel('Y');
            title(plot_title);

            if nargin > 2 && showBoundary
                legend({'Collocation Points', 'Boundary Points'}, 'Location', 'best');
            else
                legend({'Collocation Points'}, 'Location', 'best');
            end
        end
    end
end

