function demo_multiple_offsets
    % Define an L-shaped boundary in a 2 x M array:
    %   B(1,:) = x-coordinates
    %   B(2,:) = y-coordinates
    B = [
        0   0   0   0   1   2   3   3   3  ;  % x
        0   1   2   3   0   0   0   1   2   % y
    ];
    
    % Number of layers (parallel offset lines)
    nLines     = 2;       % e.g. 4 concentric offsets
    % Points along each line
    nPts       = 100;      
    % Spacing between each line
    spreadDist = 0.2;     
    
    % Compute offsets
    [xAll, yAll] = offsetBoundaryMultiple(B, nLines, nPts, spreadDist);
    
    % Plot
    figure('Color','white'); 
    hold on; axis equal; grid on;
    
    % Plot original boundary
    plot(B(1,:), B(2,:), 'b-o','LineWidth',1.5, ...
         'DisplayName','Original L Boundary');
     
    % Plot each offset line
    colors = lines(nLines);  % just for distinct colors
    for i = 1:nLines
        scatter(xAll(i,:), yAll(i,:), 15, 'r');
    end
    
    legend('Location','best');
    title('Multiple Offset Lines around an L-Shape');
end


function [xAll, yAll] = offsetBoundaryMultiple(B, nLines, nPts, spreadDist)
%OFFSETBOUNDARYMULTIPLE Generate multiple parallel offset lines
%   around a given polygonal boundary B.
%
%   INPUTS:
%       B          : 2 x M array of boundary coordinates
%                    B(1,:) = x, B(2,:) = y
%       nLines     : Number of offset lines to produce
%       nPts       : Number of points per offset line
%       spreadDist : Distance between successive offset lines
%
%   OUTPUTS:
%       xAll, yAll : Each is nLines x nPts, where row i
%                    corresponds to the i-th offset line.

    % Pre-allocate output arrays
    xAll = zeros(nLines, nPts);
    yAll = zeros(nLines, nPts);
    
    % For each line i, offset distance is i * spreadDist
    for i = 1:nLines
        currentOffset = i * spreadDist;
        [xOff, yOff] = offsetBoundarySingle(B, nPts, currentOffset);
        xAll(i,:) = xOff;
        yAll(i,:) = yOff;
    end
end

function [xOffset, yOffset] = offsetBoundarySingle(B, nPts, offsetDist)
%OFFSETBOUNDARYSINGLE Generate nPts equally spaced points along B,
%   then shift them outward by offsetDist.
%
%   B          : 2 x M array [x; y] for boundary
%   nPts       : number of sample points along boundary
%   offsetDist : outward shift distance
%
%   xOffset, yOffset : row vectors of length nPts

    % Extract X and Y
    x = B(1,:);
    y = B(2,:);
    
    % Ensure boundary is closed
    if x(end) ~= x(1) || y(end) ~= y(1)
        x(end+1) = x(1);
        y(end+1) = y(1);
    end
    
    % Segment lengths, cumulative length
    segLengths = sqrt(diff(x).^2 + diff(y).^2);  
    cumLen = [0, cumsum(segLengths)];
    totalLen = cumLen(end);
    
    % Sample points at these perimeter locations
    t = linspace(0, totalLen, nPts);
    xInterp = interp1(cumLen, x, t);
    yInterp = interp1(cumLen, y, t);
    
    % Allocate
    xOffset = zeros(1, nPts);
    yOffset = zeros(1, nPts);
    
    % Compute outward normal for each point
    for i = 1:nPts
        % indices for backward/forward difference
        if i == 1
            iPrev = nPts;
            iNext = 2;
        elseif i == nPts
            iPrev = nPts-1;
            iNext = 1;
        else
            iPrev = i - 1;
            iNext = i + 1;
        end
        
        % Tangent vector
        tx = xInterp(iNext) - xInterp(iPrev);
        ty = yInterp(iNext) - yInterp(iPrev);
        tLen = sqrt(tx^2 + ty^2);
        
        % Outward normal (assuming B is CCW)
        nx = -ty / tLen;
        ny =  tx / tLen;
        
        % Offset
        xOffset(i) = xInterp(i) + offsetDist * nx;
        yOffset(i) = yInterp(i) + offsetDist * ny;
    end
end

demo_multiple_offsets();