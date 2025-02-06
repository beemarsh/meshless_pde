classdef IrregularDomain
    % IrregularDomain A class for creating and visualizing an irregular 2D domain
    %   Example usage:
    %       vertices = [0, 1, 1, 0; 0, 0, 1, 1]; % Square
    %       hole = {[0.4, 0.6, 0.6, 0.4; 0.4, 0.4, 0.6, 0.6]}; % Square hole
    %       domainObj = IrregularDomain(vertices, hole);
    %       domainObj = domainObj.generateInteriorPoints(100);
    %       domainObj = domainObj.generateBoundaryPoints(20);
    %       domainObj.scatterPlot('Irregular Domain with Hole', true);

    properties
        InteriorPoints   % (N x 2) Interior (x, y) points
        BoundaryPoints   % (2 x M) Boundary (x, y) points
        GhostPoints      % Not used, but kept for consistency
    end

    properties(SetAccess = private)
        Vertices       % (2 x P) Matrix of polygon vertices [x; y]
        HoleVertices   % Cell array of (2 x Q_i) matrices for each hole
    end

    methods
        % Constructor
        function obj = IrregularDomain(vertices, holeVertices)
            if nargin > 0
                obj.Vertices = vertices;
                if nargin > 1
                    obj.HoleVertices = holeVertices;
                else
                    obj.HoleVertices = {};
                end
            end
        end

        % Generate interior points
        function obj = generateInteriorPoints(obj, numPointsPerSide)
            % 1) Compute bounding box
            x_min = min(obj.Vertices(1, :));
            x_max = max(obj.Vertices(1, :));
            y_min = min(obj.Vertices(2, :));
            y_max = max(obj.Vertices(2, :));

            % 2) Create a uniform grid
            x_domain = linspace(x_min, x_max, numPointsPerSide);
            y_domain = linspace(y_min, y_max, numPointsPerSide);
            [X, Y] = meshgrid(x_domain, y_domain);

            % 3) Point-in-polygon test for the main polygon
            in = inpolygon(X, Y, obj.Vertices(1, :), obj.Vertices(2, :));

            % 4) Exclude points inside any holes
            if ~isempty(obj.HoleVertices)
                for i = 1:length(obj.HoleVertices)
                    hole = obj.HoleVertices{i};
                    in_hole = inpolygon(X, Y, hole(1, :), hole(2, :));
                    in = in & ~in_hole;
                end
            end

            % 5) Extract interior points
            X_in = X(in);
            Y_in = Y(in);

            % 6) Store interior points
            obj.InteriorPoints = [X_in, Y_in];
        end

        % Generate boundary points
        function obj = generateBoundaryPoints(obj, numPointsPerEdge)
            % Initialize empty matrices for boundary points
            boundaryX = [];
            boundaryY = [];
            
            % Main polygon boundary
            numVertices = size(obj.Vertices, 2);
            for i = 1:numVertices
                % Current and next vertex (with wrap-around)
                v1 = obj.Vertices(:, i);
                v2 = obj.Vertices(:, mod(i, numVertices) + 1);
                
                % Generate points along the edge
                x_edge = linspace(v1(1), v2(1), numPointsPerEdge);
                y_edge = linspace(v1(2), v2(2), numPointsPerEdge);
                
                % Append to boundary arrays
                boundaryX = [boundaryX, x_edge];
                boundaryY = [boundaryY, y_edge];
            end
            
            % Handle holes if any
            if ~isempty(obj.HoleVertices)
                for h = 1:length(obj.HoleVertices)
                    hole = obj.HoleVertices{h};
                    numHoleVertices = size(hole, 2);
                    for i = 1:numHoleVertices
                        v1 = hole(:, i);
                        v2 = hole(:, mod(i, numHoleVertices) + 1);
                        
                        % Generate points along the hole edge
                        x_edge = linspace(v1(1), v2(1), numPointsPerEdge);
                        y_edge = linspace(v1(2), v2(2), numPointsPerEdge);
                        
                        % Append to boundary arrays
                        boundaryX = [boundaryX, x_edge];
                        boundaryY = [boundaryY, y_edge];
                    end
                end
            end

            % Combine and remove duplicates
            boundaryPoints = [boundaryX; boundaryY]';
            [uniqueRows, ~] = unique(boundaryPoints, 'rows', 'stable');
            obj.BoundaryPoints = uniqueRows';
        end

        % Scatter plot of the interior and boundary points
        function scatterPlot(obj, plot_title, showBoundary)
            figure;
            hold on;
            
            % Plot interior points
            if ~isempty(obj.InteriorPoints)
                scatter(obj.InteriorPoints(:,1), obj.InteriorPoints(:,2), 10, 'filled', 'b');
            end
            
            % Plot boundary points
            if showBoundary && ~isempty(obj.BoundaryPoints)
                scatter(obj.BoundaryPoints(1,:), obj.BoundaryPoints(2,:), 15, 'filled', 'r');
            end
            
            % Plot the polygon edges
            plot([obj.Vertices(1,:), obj.Vertices(1,1)], [obj.Vertices(2,:), obj.Vertices(2,1)], 'k-', 'LineWidth', 1.5);
            
            % Plot holes if any
            if ~isempty(obj.HoleVertices)
                for h = 1:length(obj.HoleVertices)
                    hole = obj.HoleVertices{h};
                    plot([hole(1,:), hole(1,1)], [hole(2,:), hole(2,1)], 'k--', 'LineWidth', 1.5);
                end
            end
            
            hold off;
            axis equal;
            xlabel('X');
            ylabel('Y');
            title(plot_title);
            legendEntries = {'Interior Points'};
            if showBoundary
                legendEntries{end+1} = 'Boundary Points';
            end
            legend(legendEntries, 'Location', 'best');
        end
    end
end
