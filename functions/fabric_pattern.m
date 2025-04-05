function [x, y] = fabric_pattern(n, alpha, x_center, y_center, scale)
    % Default values if not provided
    if nargin < 3, x_center = 0; end
    if nargin < 4, y_center = 0; end
    if nargin < 5, scale = 1; end  % Scale factor to control pattern size

    b = round(alpha * sqrt(n));
    phi = (sqrt(5) + 1)/2;  % Golden ratio
    x = zeros(n, 1); 
    y = x;
    
    for k = 1:n
        r = radius(k, n, b);
        theta = 2 * pi * k / phi^2;
        % Apply scaling and shifting
        x(k) = x_center + scale * r * sin(theta); 
        y(k) = y_center + scale * r * cos(theta);
    end
end

function r = radius(k, n, b)
    if k > n - b
        r = 1;
    else
        r = sqrt(k - 1/2) / sqrt(n - (b + 1)/2);
    end
end