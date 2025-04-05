function [center, radius] = smallest_circle(points)
    % Finds the smallest circle that encloses a set of points
    %
    % Inputs:
    %   points: nx2 matrix where each row is a point [x,y]
    %
    % Outputs:
    %   center: [x,y] coordinates of circle center
    %   radius: radius of the enclosing circle
    
    % Randomize points to improve average running time
    points = points(randperm(size(points,1)), :);
    
    % Find minimum enclosing circle
    [center, radius] = welzl_algorithm(points, [], size(points,1));
end

function [center, radius] = welzl_algorithm(points, R, n)
    % Welzl's algorithm implementation
    % R is the set of points on the circle boundary
    
    % Base cases
    if n == 0 || size(R,1) == 3
        [center, radius] = min_circle_from_points(R);
        return;
    end
    
    % Pick a point
    p = points(n,:);
    
    % Recursively compute minimum circle without p
    [center, radius] = welzl_algorithm(points, R, n-1);
    
    % If p is inside current circle, return current solution
    if norm(p - center) <= radius
        return;
    end
    
    % Otherwise, p must be on the boundary of the new circle
    [center, radius] = welzl_algorithm(points, [R; p], n-1);
end

function [center, radius] = min_circle_from_points(R)
    % Compute minimum circle from 0-3 points
    n = size(R,1);
    
    if n == 0
        center = [0 0];
        radius = 0;
    elseif n == 1
        center = R(1,:);
        radius = 0;
    elseif n == 2
        center = 0.5 * (R(1,:) + R(2,:));
        radius = 0.5 * norm(R(1,:) - R(2,:));
    else
        % Three points define a unique circle (if not collinear)
        A = R(2,:) - R(1,:);
        B = R(3,:) - R(1,:);
        D = 2 * (A(1)*B(2) - A(2)*B(1));
        
        if abs(D) < 1e-10  % Points are collinear
            % Find the two most distant points
            dists = pdist2(R, R);
            [~, idx] = max(dists(:));
            [i, j] = ind2sub(size(dists), idx);
            center = 0.5 * (R(i,:) + R(j,:));
            radius = 0.5 * norm(R(i,:) - R(j,:));
        else
            % Compute circle center using perpendicular bisectors
            U = [A(1)^2 + A(2)^2, B(1)^2 + B(2)^2];
            center = R(1,:) + [U(1)*B(2) - U(2)*A(2), A(1)*U(2) - B(1)*U(1)] / D;
            radius = norm(R(1,:) - center);
        end
    end
end