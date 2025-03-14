% It generates a square grid of collocation points and boundary points.
% It takes Domain as input which is a 1x2 vector: [min, max] where min and max are the boundaries of the square.
% It generates a grid of size N+1, where N is the number of collocation points.

classdef Square
    properties
        Coordinates           % (2 x K) All points in the square.

        InteriorPoints        % (2 x N) Interior (collocation) points in the square
        InteriorIndices       % Indices of interior points

        BoundaryPoints        % (2 x M) Points on the boundary of the square
        BoundaryIndices       % Indices of boundary points


        GhostPoints             % (2 x K) "Ghost" points outside the square
        Grid                    % Structure containing grid information for surface plotting
        GhostGrid
        Size                    % Size of the grid.
    end

    properties(SetAccess = private)
        Domain % [min, max] specifying the square domain
    end

    methods
        % Constructor
        function obj = Square(domain)
            % domain should be a 1x2 vector: [min, max]
            if nargin > 0
                validateattributes(domain, {'numeric'}, {'vector', 'numel', 2});
                obj.Domain = domain;

            end
        end

        % Generate interior (collocation) points in the square
        % Generate interior (collocation) points in the square
        function obj = generateGrid(obj, N)

            N = N+1; % Adding one for boundary points

            obj.Size = N;

            points = linspace(obj.Domain(1), obj.Domain(2),N);

            % Create a grid of size N.
            [X, Y] = meshgrid(points, points);

            % Save the full grid for plotting
            obj.Grid.X = X;
            obj.Grid.Y = Y;

            % Flatten the grid
            x_values = reshape(X,[],1); % This forms a column vector of x-coordinates of size N^2
            y_values = reshape(Y,[],1); % This forms a column vector of y-coordinates of size N^2

            % The following is a 3-column matrix
            % First two columns are x and y values of nodes.
            % Third column is index for either interior or boundary points.
            % 1 indicates boundary points, 0 indicates interior points
            coordinates = [x_values, y_values, zeros(length(x_values),1)];

            % Separate boundary and interior points.
            % Each value stores the index of the boundary points.
            bottom_boundary_indices = find(y_values==0 & x_values~=1); % Bottom boundary such that y=0
            right_boundary_indices=find(x_values==1 & y_values~=1); % Right boundary such that x=1
            top_boundary_indices=find(y_values==1 & x_values~=0); % Top boundary such that y=1
            left_boundary_indices=find(x_values==0 & y_values~=0); % Left boundary such that x=0

            % Assign whether the point is a boundary point or not.
            % For every index in the boundary, assign the value 1.
            coordinates(bottom_boundary_indices,3)=ones(length(bottom_boundary_indices),1); % Bottom boundary
            coordinates(right_boundary_indices,3)=ones(length(right_boundary_indices),1); % Right boundary
            coordinates(top_boundary_indices,3)=ones(length(top_boundary_indices),1); % Top boundary
            coordinates(left_boundary_indices,3)=ones(length(left_boundary_indices),1); % Left boundary

            % Store the indices of interior points
            interior_indices = find(coordinates(:,3)==0);

            % Store the actual interior collocation points
            interior_points = coordinates(interior_indices,1:2);% interior points

            % Store the indices of boundary points
            boundary_indices = find(coordinates(:,3)>0);

            % Store the actual boundary points
            boundary_points = coordinates(boundary_indices,1:2);% boundary points

            % Save all the required information in the object
            obj.InteriorPoints = interior_points;
            obj.InteriorIndices = interior_indices;

            obj.BoundaryPoints = boundary_points;
            obj.BoundaryIndices = boundary_indices;

            obj.Coordinates = coordinates(:,1:2); % All points in the square
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
            scatter(obj.InteriorPoints(:,1), obj.InteriorPoints(:,2), 10, 'filled', 'b');
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
