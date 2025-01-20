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
                    obj.LPoints = [X_L, Y_L]';
              
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

        function obj = generateGhostPoints(obj, numPoints, width, spread)
                    %Width: Number of layers (parallel offset lines)
                    % For example: Width = 4 means 4 concentric offsets

                    %numPoints: Points along each line

                    %spread:  Spacing along each line

                   
                    switch obj.Orientation
                        case 'top-left'
                            %Top Right Half Part
                            a=zeros(ceil(numPoints/2) * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=ceil(numPoints/2);
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints/2 number
                                %of points. 
                                % X spans from XBounds(2) to Domain(2)
                                % Y is fixed at YBounds(2)

                                a(start:e_nd,1:2) = [linspace(obj.XBounds(2) - i * spread, obj.Domain(2) + i*spread, ceil(numPoints/2)); (obj.Domain(2) + i*spread )*ones(1,ceil(numPoints/2))]';
                                start = ceil(numPoints/2)*i + 1;
                                e_nd = e_nd + ceil(numPoints/2);
                            end

                             %Middle Top L
                            b=zeros(ceil(numPoints/2) * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=ceil(numPoints/2);
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints number
                                %of points. 
                                % X constant at XBounds(2)
                                % Y from YBounds(1) to YBounds(2)
                                
                                b(start:e_nd,1:2) = [(obj.XBounds(2) - i*spread) * ones(1,ceil(numPoints/2)); linspace(obj.YBounds(1) + width * spread, obj.YBounds(2) + i*spread, ceil(numPoints/2))]';
                                start = ceil(numPoints/2)*i + 1;
                                e_nd = e_nd + ceil(numPoints/2);
                            end

                            %Middle Left Half Part
                            c=zeros(ceil(numPoints/2) * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=ceil(numPoints/2);
                            for i=1:width % We generate `width` number of lines
                                %   In each line we generate numPoints/2 number
                                %   of points. 
                                % X from XBounds(1) to XBounds(2)
                                % Y fixed at YBounds(1)
                                c(start:e_nd,1:2) = [linspace(obj.XBounds(1) - i * spread, obj.XBounds(2) -  spread, ceil(numPoints/2)); (obj.YBounds(1) + i * spread) * ones(1,ceil(numPoints/2))]';
                                start = ceil(numPoints/2) + 1;
                                e_nd = e_nd + ceil(numPoints/2);
                            end

                            %Right Bottom Half
                            d=zeros(ceil(numPoints/2) * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=ceil(numPoints/2);
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints/2 number
                                %of points. 
                                % X remains same - Domain(1)
                                % Y from Domain(1) to YBounds(1)
                                d(start:e_nd,1:2) = [(obj.Domain(1) - i*spread) * ones(1,ceil(numPoints/2)); linspace(obj.Domain(1) - i * spread, obj.YBounds(1) + i * spread, ceil(numPoints/2))]';
                                start = ceil(numPoints/2)*i + 1;
                                e_nd = e_nd + ceil(numPoints/2);
                            end

                            %Bottom Part
                            e=zeros(numPoints * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=numPoints;
                            for i=1:width % We generate `width` number of lines
                                %   In each line we generate numPoints/2 number
                                %   of points. 
                                % Y axis remain same. Only increase Y axis by `spread` for each line produced.
                                % Generate a linearly spaced X numPoints/2 points from
                                % domain(1) to bounds(1) each increased by
                                % spread. On the right hand side, dont
                                % encroach the X boundary.

                                e(start:e_nd,1:2) = [linspace(obj.Domain(1) - i * spread, obj.Domain(2) + i * spread, numPoints); (obj.Domain(1) - i * spread) * ones(1,numPoints)]';
                                start = numPoints + 1;
                                e_nd = e_nd + numPoints;
                            end

                            % Right
                            f=zeros(numPoints * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=numPoints;
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints/2 number
                                %of points. 
                                % X axis remain same. Only change x values
                                % by `spread` for each line produced.

                                % Generate a linearly spaced Y values numPoints points from
                                % domain(1) to domain(2) each increased by
                                % spread.
                                
                                f(start:e_nd,1:2) = [(obj.Domain(2) + i*spread) * ones(1,numPoints); linspace(obj.Domain(1) - i * spread, obj.Domain(2) + i * spread, numPoints) ]';
                                start = numPoints*i + 1;
                                e_nd = e_nd + numPoints;
                            end                             
                        %---------------------------------------Top Left Ends--------------------------------------------------------------------                            
                        case 'top-right'
                            %Top Left Half Part
                            a=zeros(ceil(numPoints/2) * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=ceil(numPoints/2);
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints/2 number
                                %of points. 
                                % X- axis from Domain(1) to XBounds(1) with
                                % spread added on both sides
                                % Y axis fixed at Domain(2) with added
                                % spread

                                a(start:e_nd,1:2) = [linspace(obj.Domain(1) - i * spread, obj.XBounds(1) + i*spread, ceil(numPoints/2)); (obj.Domain(2) + i*spread )*ones(1,ceil(numPoints/2))]';
                                start = ceil(numPoints/2)*i + 1;
                                e_nd = e_nd + ceil(numPoints/2);
                            end

                             %Left Part
                            b=zeros(numPoints * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=numPoints;
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints number
                                %of points. 
                                % X remain constant at Domain(1) - i *
                                % spread
                                % Y goes from Domain(1) to Domain(2) with
                                % spread on both sides
                                

                                b(start:e_nd,1:2) = [(obj.Domain(1) - i*spread) * ones(1,numPoints); linspace(obj.Domain(1) - i * spread, obj.Domain(2) + i*spread, numPoints)]';
                                start = numPoints*i + 1;
                                e_nd = e_nd + numPoints;
                            end

                            %Bottom Part
                            c=zeros(numPoints * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=numPoints;
                            for i=1:width % We generate `width` number of lines
                                %   In each line we generate numPoints/2 number
                                %   of points. 
                                % Y axis remain same. Only increase Y axis by `spread` for each line produced.
                                % Generate a linearly spaced X numPoints/2 points from
                                % domain(1) to bounds(1) each increased by
                                % spread. On the right hand side, dont
                                % encroach the X boundary.

                                c(start:e_nd,1:2) = [linspace(obj.Domain(1) - i * spread, obj.Domain(2) + i * spread, numPoints); (obj.Domain(1) - i * spread) * ones(1,numPoints)]';
                                start = numPoints + 1;
                                e_nd = e_nd + numPoints;
                            end

                            %Right Half
                            d=zeros(ceil(numPoints/2) * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=ceil(numPoints/2);
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints/2 number
                                %of points. 
                                % X remains same - Domain(2)
                                % Y from Domain(1) to Domain to YBounds(1)
                                d(start:e_nd,1:2) = [(obj.Domain(2) + i*spread) * ones(1,ceil(numPoints/2)); linspace(obj.Domain(1) - i * spread, obj.YBounds(1) + i * spread, ceil(numPoints/2))]';
                                start = ceil(numPoints/2)*i + 1;
                                e_nd = e_nd + ceil(numPoints/2);
                            end

                            %Middle Right Half L Shape
                            e=zeros(ceil(numPoints/2) * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=ceil(numPoints/2);
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints/2 number
                                %of points. 
                                % X from XBounds(1) to XBounds(2)
                                % But remember not to encroach boundary on
                                % left side
                                % Y remains the same at YBounds(1)
                                
                                e(start:e_nd,1:2) = [linspace(obj.XBounds(1) + spread, obj.XBounds(2) + i * spread, ceil(numPoints/2)); (obj.YBounds(1) + i*spread) * ones(1,ceil(numPoints/2)) ]';
                                start = ceil(numPoints/2)*i + 1;
                                e_nd = e_nd + ceil(numPoints/2);
                            end

                            % Top Middle Half
                            f=zeros(ceil(numPoints/2) * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=ceil(numPoints/2);
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints/2 number
                                %of points. 
                                % X constant at XBounds(1)

                                % Y from YBounds(1) to YBounds(2) with
                                % spread
                                
                                f(start:e_nd,1:2) = [(obj.XBounds(1) + i*spread) * ones(1,ceil(numPoints/2)); linspace(obj.YBounds(1) + width*spread, obj.YBounds(2) + i * spread, ceil(numPoints/2)) ]';
                                start = ceil(numPoints/2)*i + 1;
                                e_nd = e_nd + ceil(numPoints/2);
                            end                     
                        %---------------------------------------Top Right Ends--------------------------------------------------------------------

                        case 'bottom-left'
                            %Top Part
                            a=zeros(numPoints * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=numPoints;
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints number
                                %of points. 
                                % Increase x axis by `spread` on both ends.
                                % Increase y axis by `spread` for each line

                                a(start:e_nd,1:2) = [linspace(obj.Domain(1) - i * spread, obj.Domain(2) + i*spread, numPoints); (obj.Domain(2) + i*spread )*ones(1,numPoints)]';
                                start = numPoints*i + 1;
                                e_nd = e_nd + numPoints;
                            end

                             %Top Left Half Part
                            b=zeros(ceil(numPoints/2) * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=ceil(numPoints/2);
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints/2 number
                                %of points. 
                                % X axis remain same. Only decrease x axis by `spread` for each line produced.
                                % Generate a linearly spaced points for Y axis from
                                % ybounds(2) to domain(2) each increased by
                                % spread.
                                

                                b(start:e_nd,1:2) = [(obj.Domain(1) - i*spread) * ones(1,ceil(numPoints/2)); linspace(obj.YBounds(2) - i * spread, obj.Domain(2) + i*spread, ceil(numPoints/2))]';
                                start = ceil(numPoints/2)*i + 1;
                                e_nd = e_nd + ceil(numPoints/2);
                            end

                            %Middle Right to Left
                            c=zeros(ceil(numPoints/2) * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=ceil(numPoints/2);
                            for i=1:width % We generate `width` number of lines
                                %   In each line we generate numPoints/2 number
                                %   of points. 
                                % Y axis remain same. Only increase Y axis by `spread` for each line produced.
                                % Generate a linearly spaced X numPoints/2 points from
                                % domain(1) to bounds(1) each increased by
                                % spread. On the right hand side, dont
                                % encroach the X boundary.

                                c(start:e_nd,1:2) = [linspace(obj.XBounds(1) - i * spread, obj.XBounds(2) - spread, ceil(numPoints/2)); (obj.YBounds(2) - i * spread) * ones(1,ceil(numPoints/2))]';
                                start = ceil(numPoints/2)*i + 1;
                                e_nd = e_nd + ceil(numPoints/2);
                            end

                            %Middle Bottom L Shape
                            d=zeros(ceil(numPoints/2) * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=ceil(numPoints/2);
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints/2 number
                                %of points. 
                                % X axis remain same. Only increase X axis
                                % by `spread` for each line produced. Here
                                % we need to be careful because if we
                                % increase X as usual by spread on the top
                                % end we will encroach the boundary.
                                % Threrfore, at the top end, we only
                                % subtract by spread.

                                d(start:e_nd,1:2) = [(obj.XBounds(2) - i*spread) * ones(1,ceil(numPoints/2)); linspace(obj.YBounds(1) - i * spread, obj.YBounds(2) - width * spread, ceil(numPoints/2))]';
                                start = ceil(numPoints/2)*i + 1;
                                e_nd = e_nd + ceil(numPoints/2);
                            end

                            %Bottom Right Half Shape
                            e=zeros(ceil(numPoints/2) * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=ceil(numPoints/2);
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints/2 number
                                %of points. 
                                % Y axis remain same. Only decrease Y axis
                                % by `spread` for each line produced.

                                % Generate a linearly spaced X values numPoints/2 points from
                                % xbounds(2) to domain(2) each increased by
                                % spread.
                                
                                e(start:e_nd,1:2) = [linspace(obj.XBounds(2) - i * spread, obj.Domain(2) + i * spread, ceil(numPoints/2)); (obj.YBounds(1) - i*spread) * ones(1,ceil(numPoints/2)) ]';
                                start = ceil(numPoints/2)*i + 1;
                                e_nd = e_nd + ceil(numPoints/2);
                            end

                            % Right
                            f=zeros(numPoints * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=numPoints;
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints/2 number
                                %of points. 
                                % X axis remain same. Only change x values
                                % by `spread` for each line produced.

                                % Generate a linearly spaced Y values numPoints points from
                                % domain(1) to domain(2) each increased by
                                % spread.
                                
                                f(start:e_nd,1:2) = [(obj.Domain(2) + i*spread) * ones(1,numPoints); linspace(obj.Domain(1) - i * spread, obj.Domain(2) + i * spread, numPoints) ]';
                                start = numPoints*i + 1;
                                e_nd = e_nd + numPoints;
                            end                     
                        %---------------------------------------Bottom Left Ends--------------------------------------------------------------------

                        case 'bottom-right'
                            %Top Part
                            a=zeros(numPoints * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=numPoints;
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints number
                                %of points. 
                                % Increase x axis by `spread` on both ends.
                                % Increase y axis by `spread` for each line

                                a(start:e_nd,1:2) = [linspace(obj.Domain(1) - i * spread, obj.Domain(2) + i*spread, numPoints); (obj.Domain(2) + i*spread )*ones(1,numPoints)]';
                                start = numPoints*i + 1;
                                e_nd = e_nd + numPoints;
                            end

                            %Left Part
                            b=zeros(numPoints * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=numPoints;
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints number
                                %of points. 
                                % X axis remain same. Only increase x axis by `spread` for each line produced.
                                % Generate a linearly spaced points from
                                % domain(1) to domain(2) each increased by
                                % spread.
                                

                                b(start:e_nd,1:2) = [(obj.Domain(1) - i*spread) * ones(1,numPoints); linspace(obj.Domain(1) - i * spread, obj.Domain(2) + i*spread, numPoints)]';
                                start = numPoints*i + 1;
                                e_nd = e_nd + numPoints;
                            end

                            %Bottom Left Half Part
                            c=zeros(ceil(numPoints/2) * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=ceil(numPoints/2);
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints/2 number
                                %of points. 
                                % Y axis remain same. Only increase Y axis by `spread` for each line produced.
                                % Generate a linearly spaced X numPoints/2 points from
                                % domain(1) to bounds(1) each increased by
                                % spread.
                                

                                c(start:e_nd,1:2) = [linspace(obj.Domain(1) - i * spread, obj.XBounds(1) + i*spread, ceil(numPoints/2)); (obj.Domain(1) - i*spread) * ones(1,ceil(numPoints/2))]';
                                start = ceil(numPoints/2)*i + 1;
                                e_nd = e_nd + ceil(numPoints/2);
                            end

                            %Middle Bottom L Shape
                            d=zeros(ceil(numPoints/2) * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=ceil(numPoints/2);
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints/2 number
                                %of points. 
                                % X axis remain same. Only increase X axis
                                % by `spread` for each line produced. Here
                                % we need to be careful because if we
                                % increase X as usual by spread on the top
                                % end we will encroach the boundary.
                                % Threrfore, at the top end, we only
                                % subtract by spread.
                                % Generate a linearly spaced X numPoints/2 points from
                                % domain(1) to bounds(1) each increased by
                                % spread.
                                

                                d(start:e_nd,1:2) = [(obj.XBounds(1) + i*spread) * ones(1,ceil(numPoints/2)); linspace(obj.YBounds(1) - i * spread, obj.YBounds(2) - spread, ceil(numPoints/2))]';
                                start = ceil(numPoints/2)*i + 1;
                                e_nd = e_nd + ceil(numPoints/2);
                            end

                            %Middle Right L Shape
                            e=zeros(ceil(numPoints/2) * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=ceil(numPoints/2);
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints/2 number
                                %of points. 
                                % Y axis remain same. Only decrease Y axis
                                % by `spread` for each line produced.

                                % Generate a linearly spaced X values numPoints/2 points from
                                % xbounds(1) to xbounds(1) each increased by
                                % spread. Similarly, dont encroach the
                                % boundary for the left end point.
                                

                                e(start:e_nd,1:2) = [linspace(obj.XBounds(1) + (width) * spread, obj.XBounds(2) + i * spread, ceil(numPoints/2)); (obj.YBounds(2) - i*spread) * ones(1,ceil(numPoints/2)) ]';
                                start = ceil(numPoints/2)*i + 1;
                                e_nd = e_nd + ceil(numPoints/2);
                            end

                            %Top Right Half
                            f=zeros(ceil(numPoints/2) * width, 2);
                            start=1; % Start from 1 to numPoints
                            e_nd=ceil(numPoints/2);
                            for i=1:width % We generate `width` number of lines
                                %In each line we generate numPoints/2 number
                                %of points. 
                                % X axis remain same. Only change x values
                                % by `spread` for each line produced.

                                % Generate a linearly spaced Y values numPoints/2 points from
                                % ybounds(2) to domain(2) each increased by
                                % spread.
                                

                                f(start:e_nd,1:2) = [(obj.XBounds(2) + i*spread) * ones(1,ceil(numPoints/2)); linspace(obj.YBounds(2) - i * spread, obj.Domain(2) + i * spread, ceil(numPoints/2)) ]';
                                start = ceil(numPoints/2)*i + 1;
                                e_nd = e_nd + ceil(numPoints/2);
                            end

                        otherwise
                            error('Invalid orientation. Choose from ''top-left'', ''top-right'', ''bottom-left'', or ''bottom-right''.');
                    end

                    % Combine the rows together                  
                    points = [a;b;c;d;e;f];

                    % Use the `unique` function to find unique rows
                    [unique_rows, ~] = unique(points, 'rows', 'stable');
                    %Duplication remove completed


                    % Transpose back to maintain the 2xM format and store it
                    obj.GhostPoints = unique_rows';
        end

        function scatterPlot(obj, plot_title, boundary, ghost)
                    figure;
                    scatter(obj.LPoints(1,:), obj.LPoints(2,:), 10, 'filled','b');
                    hold on;
                    if boundary
                        scatter(obj.BoundaryPoints(1,:),obj.BoundaryPoints(2,:),15,'filled','r');
                    end

                    if ghost
                       scatter(obj.GhostPoints(1,:),obj.GhostPoints(2,:),15,'filled','g');
                    end
                    hold off;

                    axis equal;
                    xlabel('X');
                    ylabel('Y');
                    title(plot_title);
                    legend({'Collocation Points', 'Boundary Points', 'GhostPoints'}, 'Location', 'best');

            end
    end

    methods(Access=private)
      
    end
            
end






     