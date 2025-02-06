classdef Circular
   %CIRCULARDOMAIN A class for creating and visualizing a 2D circular domain
   %   Example usage:
   %       circleObj = CircularDomain([0, 0], 1);          % Circle centered at (0,0) with radius 1
   %       circleObj = circleObj.generateCDomain(50);      % Generate ~50x50 interior grid
   %       circleObj = circleObj.generateBoundaryPoints(50); % 50 boundary points
   %       circleObj = circleObj.generateGhostPoints(5, 0.05, 50); % 5 offset rings, each 0.05 away, 50 points each
   %       circleObj.scatterPlot('My Circle', true);       % Plot interior+boundary (and if you modify scatterPlot to show ghost points, that too)
   properties
       CPoints        % (N x 2) The (x,y) points that lie inside the circular domain
       BoundaryPoints % (2 x M) The (x, y) points that lie on the boundary of the circle
       GhostPoints    % (2 x K) The (x, y) "ghost" points outside the circle
   end
   properties(SetAccess = private)
       Center % [cx, cy] center of the circle
       Radius % Scalar radius
   end
   methods
       % Constructor
       function obj = Circular(center, radius)
           if nargin > 0
               obj.Center = center;
               obj.Radius = radius;
           end
       end
       % Generate interior points
       function obj = generateCDomain(obj, numPoints)
           % 1) Compute bounding box
           x_min = obj.Center(1) - obj.Radius;
           x_max = obj.Center(1) + obj.Radius;
           y_min = obj.Center(2) - obj.Radius;
           y_max = obj.Center(2) + obj.Radius;
           % 2) Create uniform grid
           x_domain = linspace(x_min, x_max, numPoints);
           y_domain = linspace(y_min, y_max, numPoints);
           [X, Y] = meshgrid(x_domain, y_domain);
           % 3) Keep only points inside the circle
           dist_sq = (X - obj.Center(1)).^2 + (Y - obj.Center(2)).^2;
           mask = dist_sq <= obj.Radius^2;
           % 4) Extract (x, y)
           X_in = X(mask);
           Y_in = Y(mask);
           % 5) Store
           obj.CPoints = [X_in, Y_in];
       end
       % Generate boundary points
       function obj = generateBoundaryPoints(obj, numPoints)
           theta = linspace(0, 2*pi, numPoints);
           x_b = obj.Center(1) + obj.Radius * cos(theta);
           y_b = obj.Center(2) + obj.Radius * sin(theta);
           % Store as 2 x M
           obj.BoundaryPoints = [x_b; y_b];
       end
       % Generate ghost points outside the circular boundary
       function obj = generateGhostPoints(obj, width, spread, numPoints)
           % width: number of offset "rings" outside the circle
           % spread: radial distance between each ring
           % numPoints: number of points per ring
           %
           % Example:
           %   width = 5, spread = 0.05, numPoints = 50
           %   => 5 concentric offset circles, each 0.05 bigger in radius
           %      than the last, each discretized into 50 points.
           % Preallocate or collect in a local array
           ghostArray = [];
           for i = 1:width
               % Radius of the i-th offset circle
               offsetRadius = obj.Radius + i * spread;
               % Discretize in angle
               theta = linspace(0, 2*pi, numPoints);
               % (x, y) for this ring
               x_ghost = obj.Center(1) + offsetRadius * cos(theta);
               y_ghost = obj.Center(2) + offsetRadius * sin(theta);
               % Append to ghostArray as row (x, y)
               ghostArray = [ ghostArray; x_ghost(:), y_ghost(:) ];
           end
           % Remove duplicates (if any) and store as 2xK
           ghostArray = unique(ghostArray, 'rows', 'stable');
           obj.GhostPoints = ghostArray';
       end
       % Scatter plot
       function scatterPlot(obj, plot_title, showBoundary, showGhost)
           figure;
           % 1) Scatter the interior points
           if ~isempty(obj.CPoints)
               scatter(obj.CPoints(:,1), obj.CPoints(:,2), 10, 'filled','b');
               hold on;
           end
           % 2) Scatter the boundary points in red
           if nargin > 2 && showBoundary && ~isempty(obj.BoundaryPoints)
               scatter(obj.BoundaryPoints(1,:), obj.BoundaryPoints(2,:), 15, 'filled', 'r');
           end
           % 3) Scatter the ghost points in green
           if nargin > 3 && showGhost && ~isempty(obj.GhostPoints)
               scatter(obj.GhostPoints(1,:), obj.GhostPoints(2,:), 15, 'filled', 'g');
           end
           axis equal;
           xlabel('X');
           ylabel('Y');
           title(plot_title);
           if nargin > 3 && showGhost
               legend({'Collocation Points', 'Boundary Points', 'Ghost Points'}, 'Location', 'best');
           elseif nargin > 2 && showBoundary
               legend({'Collocation Points', 'Boundary Points'}, 'Location', 'best');
           else
               legend({'Collocation Points'}, 'Location', 'best');
           end
       end
   end
end
