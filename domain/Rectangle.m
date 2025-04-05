classdef Rectangle
    properties
        Coordinates           % (2 x K) All points in the rectangle.

        InteriorPoints        % (2 x N) Interior (collocation) points in the rectangle

        BoundaryPoints        % (2 x M) Points on the boundary of the rectangle

        Nb                % Number of boundary points
        Ni                % Number of interior points
        N               % Total number of points

        Center              % Center of the rectangle
    end


    methods
        % Constructor
        function obj = Rectangle(bottom_left, L, W, Nb, Ni)
            % bottom_left should be a 1x2 vector: [x, y]
            % L and W are the length and width of the rectangle
            % Nb is the number of boundary points
            % Ni is the number of interior points
            if nargin > 0
                validateattributes(bottom_left, {'numeric'}, {'vector', 'numel', 2});
                validateattributes(L, {'numeric'}, {'scalar', 'positive'});
                validateattributes(W, {'numeric'}, {'scalar', 'positive'});
                validateattributes(Nb, {'numeric'}, {'scalar', 'integer', 'positive'});
                validateattributes(Ni, {'numeric'}, {'scalar', 'integer', 'positive'});

                obj.Nb = Nb;
                obj.Ni = Ni;
                obj.N = Nb + Ni;

                % Generate boundary points
                [x_b, y_b] = rectangle_boundary_from_corner(bottom_left(1), bottom_left(2), L, W, Nb);
                boundary_pts = [x_b; y_b]';
                obj.BoundaryPoints = boundary_pts;

                % Generate interior points
                [interior_pts, ~, center] = within_polygon(boundary_pts, Ni,true);

                obj.InteriorPoints = interior_pts;
                obj.Center = center;

                obj.Coordinates = [interior_pts; boundary_pts];

            end

        end

        function plot(obj, title_str)
            % Plot the rectangle and its points
            figure;
            hold on;
            plot(obj.BoundaryPoints(:, 1), obj.BoundaryPoints(:, 2), 'ro', 'MarkerSize', 5, 'DisplayName', 'Boundary Points');
            plot(obj.InteriorPoints(:, 1), obj.InteriorPoints(:, 2), 'bo', 'MarkerSize', 5, 'DisplayName', 'Interior Points');
            plot(obj.Center(1), obj.Center(2), 'go', 'MarkerSize', 10, 'DisplayName', 'Center');
            axis equal;
            grid on;
            title(title_str);
            xlabel('X-axis');
            ylabel('Y-axis');
            legend show;
        end


    end
end

function [x, y] = rectangle_boundary_from_corner(x_bl, y_bl, W, H, N)
% Generates N boundary points of a rectangle starting from bottom left corner
% INPUTS:
%   x_bl, y_bl - Bottom left corner coordinates
%   W, H - Width and Height of the rectangle
%   N - Total number of boundary points
% OUTPUTS:
%   x, y - Arrays of x and y coordinates of boundary points

% Generate parameter t (from 0 to 1)
t = linspace(0, 1, N+1);
t(end) = []; % Remove last point to avoid duplication

% Initialize arrays
x = zeros(size(t));
y = zeros(size(t));

% Assign points based on t range (moving counterclockwise)
for i = 1:N
    if t(i) < 1/4
        % Bottom edge
        x(i) = x_bl + 4 * t(i) * W;
        y(i) = y_bl;
    elseif t(i) < 1/2
        % Right edge
        x(i) = x_bl + W;
        y(i) = y_bl + 4 * (t(i) - 1/4) * H;
    elseif t(i) < 3/4
        % Top edge
        x(i) = x_bl + W - 4 * (t(i) - 1/2) * W;
        y(i) = y_bl + H;
    else
        % Left edge
        x(i) = x_bl;
        y(i) = y_bl + H - 4 * (t(i) - 3/4) * H;
    end
end
end
