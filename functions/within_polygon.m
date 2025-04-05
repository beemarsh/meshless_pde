function [ points, D, center ] = within_polygon( polygon_points, number_of_points, uniform )
    %WITHIN_POLYGON_POINT_GENERATOR Generates points within a polygon
    %   polygon_points: polygon described with points in [x1, y1; x2, y2;....; xn, yn]
    %   number_of_points: number of points to generate
    %   uniform: (optional) if true, generates uniformly distributed points
    %           default is false
    
    if nargin < 3
        uniform = false;
    end
    
    center = [mean(polygon_points(:, 1)), mean(polygon_points(:, 2))];
    
    x_d = abs(max(polygon_points(:, 1)) - min(polygon_points(:, 1)));
    y_d = abs(max(polygon_points(:, 2)) - min(polygon_points(:, 2)));
    
    D = norm([x_d, y_d]);
    
    if uniform
        % Create uniform grid with inset from boundary
        x_min = min(polygon_points(:, 1));
        x_max = max(polygon_points(:, 1));
        y_min = min(polygon_points(:, 2));
        y_max = max(polygon_points(:, 2));
        
        % Add small inset (1% of domain size) to avoid boundary
        inset_x = 0.01 * x_d;
        inset_y = 0.01 * y_d;
        
        x_min = x_min + inset_x;
        x_max = x_max - inset_x;
        y_min = y_min + inset_y;
        y_max = y_max - inset_y;
        
        n = ceil(sqrt(number_of_points * x_d/y_d));
        m = ceil(number_of_points/n);
        
        x = linspace(x_min, x_max, n);
        y = linspace(y_min, y_max, m);
        [X, Y] = meshgrid(x, y);
        samples = [X(:) Y(:)];
    else
        % Create circle domain (original random distribution)
        samples = disk_2d(D, center, 10000);
    end
    
    %% generate points
    is_within_domain = inpolygon(samples(:,1), samples(:,2), polygon_points(:,1), polygon_points(:,2));
    samples = samples(is_within_domain, :);
    
    if uniform
        % For uniform distribution, take first number_of_points that are inside
        if size(samples, 1) < number_of_points
            warning('Not enough points inside polygon. Returning all available points.');
            points = samples;
        else
            points = samples(1:number_of_points, :);
        end
    else
        % For random distribution, take random subset
        if size(samples, 1) < number_of_points
            warning('Not enough points inside polygon. Returning all available points.');
            points = samples;
        else
            idx = randperm(size(samples, 1), number_of_points);
            points = samples(idx, :);
        end
    end
    
end