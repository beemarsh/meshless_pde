classdef LDomain
    properties        
        LPoints % The Generate N Dimensional LShapeDomain (It is generated only after running the generate function
    end

    properties(SetAccess = private)
        Domain %The Entire Domain [Square Domain]
        Dimension=2
        Bounds % To create an L shaped domain from a square domain A, we do A/B which is the difference of set A and B
    end


    methods
        function obj = LDomain(dimension, domain, bounds)
            if nargin > 0
                % Validate input arguments
                validateattributes(dimension, {'numeric'}, {'nonempty'});
                validateattributes(domain, {'numeric'}, {'size', [1, 2], 'increasing'});
                validateattributes(bounds, {'numeric'}, {'size', [1, 2], 'increasing'});
               
                % Set the properties
                obj.Domain=domain;
                obj.Dimension=dimension;
                obj.Bounds=bounds;
            end
        end

        function obj = generateLShape(obj, numPoints)
                % Validate the number of points
                validateattributes(numPoints, {'numeric'}, {'scalar', 'positive', 'integer'}, 'DomainProperties', 'numPoints');

                x = linspace(obj.Domain(1), obj.Domain(2), numPoints);
                [X, Y] = meshgrid(x,x);

                % Create the L-shaped domain mask
                L_mask = ~(X >= obj.Bounds(1) & X <= obj.Bounds(2) & Y >= obj.Bounds(1) & Y <= obj.Bounds(2));

                % Extract points within the L-shaped domain
                X_L = X(L_mask);
                Y_L = Y(L_mask);

                % Combine into a single matrix of points
                obj.LPoints = [X_L, Y_L];      
            end

            function scatterPlot(obj, plot_title)
                    figure;
                    scatter(obj.LPoints(:,1), obj.LPoints(:,2), 10, 'filled');
                    axis equal;
                    xlabel('X');
                    ylabel('Y');
                    title(plot_title);
            end
    end

end






     