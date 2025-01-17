classdef LDomain
    properties        
        LPoints % The Generate 2 Dimensional LShapeDomain (It is generated only after running the generate function
        GhostPoints
        BoundaryPoints
        Orientation
    end

    properties(SetAccess = private)
        Domain %   domain Square A vector specifying the [X,Y]^2 of the full domain
        XBounds %   A vector specifying the [xmin, xmax] of the excluded region
        YBounds %    A vector specifying the [xmin, xmax] of the excluded region
    end


    methods
        function obj = LDomain(domain, x_bounds, y_bounds, orientation)
            if nargin > 0
                % Validate input arguments
                %validateattributes(domain, {'numeric'}, {'size', [1, 2,3,4]});
               % validateattributes(bounds, {'numeric'}, {'size', [1, 2,3,4]});
               
                % Set the properties
                obj.Domain=domain;
                obj.XBounds=x_bounds;
                obj.YBounds=y_bounds;
                obj.Orientation=orientation;
            end
        end

        function obj = generateLShape(obj, numPoints)

                    % Generate the L-shaped grid using numgrid
                    L_grid = numgrid('L', numPoints);   

                    % Rotate the L-shaped grid based on the specified orientation
                    switch obj.Orientation
                        case 'top-left'
                            % No rotation needed
                        case 'top-right'
                            L_grid = rot90(L_grid, 1); % Rotate 90 degrees clockwise
                        case 'bottom-left'
                            L_grid = rot90(L_grid, -1); % Rotate 90 degrees counterclockwise
                        case 'bottom-right'
                            L_grid = rot90(L_grid, 2); % Rotate 180 degrees
                        otherwise
                            error('Invalid orientation. Choose from ''top-left'', ''top-right'', ''bottom-left'', or ''bottom-right''.');
                    end

                    % Define the domain for the meshgrid
                    x_domain = linspace(obj.Domain(1), obj.Domain(2), numPoints);
                    y_domain = linspace(obj.Domain(1), obj.Domain(2), numPoints);
                    %y_domain = linspace(-1, 1, n);
                    
                    % Generate the meshgrid
                    [X, Y] = meshgrid(x_domain, y_domain);

                    % Extract points corresponding to the L-shaped domain
                    L_mask = L_grid > 0; % Logical mask where L_grid > 0
                    
                    % Apply the mask to get points inside the L-shaped domain
                    X_L = X(L_mask);
                    Y_L = Y(L_mask);

                    % Assign the value to the object
                    obj.LPoints = [X_L, Y_L];
              
        end

        function obj = generateBoundaryPoints(obj, numPoints)
            % Generate boundary points
            %Only for bottom left

            switch obj.Orientation
                        case 'top-left'
                            a = [linspace(obj.XBounds(2),obj.Domain(2),ceil(numPoints/2)); obj.Domain(2)*ones(1,ceil(numPoints/2))]; %Top right half                           
                            b = [obj.XBounds(2)*ones(1,ceil(numPoints/2)); linspace(obj.YBounds(1), obj.YBounds(2), ceil(numPoints/2))]; %Top mid left half 
                            %b=[;];
                            c = [linspace(obj.XBounds(1), obj.XBounds(2), ceil(numPoints/2)); obj.YBounds(1) * ones(1,ceil(numPoints/2))]; %Left mid
                            d = [obj.Domain(1) * ones(1, ceil(numPoints/2)); linspace(obj.Domain(1), obj.YBounds(1), ceil(numPoints/2))]; %Left Bottom Half
                            e = [linspace(obj.Domain(1),obj.Domain(2),numPoints); obj.Domain(1) * ones(1,numPoints)]; %Bottom
                            f = [obj.Domain(2)*ones(1,numPoints);linspace(obj.Domain(1),obj.Domain(2),numPoints)]; %Right
                            
                        case 'top-right'
                            a = [linspace(obj.Domain(1),obj.XBounds(1),ceil(numPoints/2)); obj.Domain(2)*ones(1,ceil(numPoints/2))]; %Top left half                           
                            b = [obj.Domain(1) * ones(1,numPoints);linspace(obj.Domain(1),obj.Domain(2),numPoints)]; %Left
                            c = [linspace(obj.Domain(1),obj.Domain(2),numPoints);obj.Domain(1) * ones(1,numPoints)]; %Bottom
                            d = [obj.Domain(2) * ones(1, ceil(numPoints/2)); linspace(obj.Domain(1), obj.YBounds(1), ceil(numPoints/2))]; %Right Bottom Half
                            e = [linspace(obj.XBounds(1),obj.XBounds(2),ceil(numPoints/2));obj.YBounds(1)*ones(1,ceil(numPoints/2))]; %Right Mid
                            f = [obj.XBounds(1) * ones(1,ceil(numPoints/2));linspace(obj.YBounds(1),obj.YBounds(2),ceil(numPoints/2))]; %Top Mid

                        case 'bottom-left'
                            a = [linspace(obj.Domain(1),obj.Domain(2),numPoints); obj.Domain(2)*ones(1,numPoints)]; %Top
                            b = [obj.Domain(1) * ones(1,ceil(numPoints/2));linspace(obj.YBounds(2),obj.Domain(2),ceil(numPoints/2))]; %Left mid half
                            c = [linspace(obj.XBounds(1),obj.XBounds(2),ceil(numPoints/2)); obj.YBounds(2)*ones(1,ceil(numPoints/2))]; %Mid Left
                            d = [obj.XBounds(2)*ones(1,ceil(numPoints/2));linspace(obj.YBounds(1),obj.YBounds(2),ceil(numPoints/2))]; %Left Bottom Mid half
                            e = [linspace(obj.XBounds(2),obj.Domain(2),ceil(numPoints/2)); obj.Domain(1) * ones(1,ceil(numPoints/2))]; %Bottom right half
                            f = [obj.Domain(2)*ones(1,numPoints);linspace(obj.Domain(1),obj.Domain(2),numPoints)]; %Right

                            
                        case 'bottom-right'
                            a = [linspace(obj.Domain(1),obj.Domain(2),numPoints); obj.Domain(2)*ones(1,numPoints)]; %Top                           
                            b = [obj.Domain(1) * ones(1,numPoints);linspace(obj.Domain(1),obj.Domain(2),numPoints)]; %Left
                            c = [linspace(obj.Domain(1),obj.XBounds(1),ceil(numPoints/2));obj.Domain(1) * ones(1,ceil(numPoints/2))]; %Bottom Left Half
                            d = [obj.XBounds(1) * ones(1, ceil(numPoints/2)); linspace(obj.YBounds(1), obj.YBounds(2), ceil(numPoints/2))]; %Right Mid
                            e = [linspace(obj.XBounds(1),obj.XBounds(2),ceil(numPoints/2));obj.YBounds(2)*ones(1,ceil(numPoints/2))]; %Right Mid
                            f = [obj.Domain(2) * ones(1,ceil(numPoints/2));linspace(obj.YBounds(2),obj.Domain(2),ceil(numPoints/2))]; %Top Right

                        otherwise
                            error('Invalid orientation. Choose from ''top-left'', ''top-right'', ''bottom-left'', or ''bottom-right''.');
            end
            
            %The corner poitns are duplicated. Now remove duplicate points

            % Transpose the matrix to work with rows as individual points
            points_transposed = [a,b,c,d,e,f]';

            % Use the `unique` function to find unique rows
            [unique_rows, ~] = unique(points_transposed, 'rows', 'stable');
            %Duplication remove completed

            % Transpose back to maintain the 2xM format and store it
            obj.BoundaryPoints = unique_rows';
        end

        function scatterPlot(obj, plot_title, boundary)
                    figure;
                    scatter(obj.LPoints(:,1), obj.LPoints(:,2), 10, 'filled','b');
                    hold on;
                    if boundary
                        scatter(obj.BoundaryPoints(1,:),obj.BoundaryPoints(2,:),15,'filled','r');
                    end
                    axis equal;
                    xlabel('X');
                    ylabel('Y');
                    title(plot_title);
                    legend({'Collocation Points', 'Boundary Points'}, 'Location', 'best');

            end
    end

end






     