addpath('domain');

function plotLShapeDomain()
    % L-shaped Domain
    % First instantiate a class
    new_L = Lshape([-3,3],[-3,0],[-3,0], 'bottom-left');
    
    % Then generate the collocation points
    % This will only generate interior points
    new_L = new_L.generateLShape(40);
    
    % Then, generate boundary points
    new_L = new_L.generateBoundaryPoints(40);
    
    %Then, generate Ghost points
    new_L = new_L.generateGhostPoints(10, 2, 0.1);
    
    %Then, scatter plot the domain.
    new_L.scatterPlot("L-Shaped Domain",true, true);
% End L Shaped Domain
% -------------------------------------------------------------------------
end


% Circular Domain
function plotCircularDomain()
    
    % First instantiate a class
    circleObj = Circular([0, 0], 1); % Center at 0,0 and radius = 1 unit.
    
    % Generate 50x50 interior grid.
    circleObj = circleObj.generateCDomain(50);
    
    % Generate 50 boundary points
    circleObj = circleObj.generateBoundaryPoints(50);
    
    %Generate ghost points:
    circleObj = circleObj.generateGhostPoints(5, 0.05, 50); %
    
    % Scatter plot the points
    circleObj.scatterPlot('Circle with Ghost Points', true, true);

% End Circular Domain
% -------------------------------------------------------------------------
end

function plotIrregularDomain()
    % Irregular Domain
    
    % Define main polygon vertices (e.g., a pentagon)
    mainVertices = [
        0, 1, 2, 1.5, 0.5;
        0, 0, 1, 2, 1.5
    ];
    
    % Define a hole (e.g., a smaller triangle)
    hole1 = [
        0.8, 1.2, 1;
        0.8, 0.8, 1
    ];
    
    % Create the IrregularDomain object
    domainObj = Irregular(mainVertices, {hole1});
    
    % Generate interior points
    domainObj = domainObj.generateInteriorPoints(100);
    
    % Generate boundary points
    domainObj = domainObj.generateBoundaryPoints(20);
    
    % Plot the domain
    domainObj.scatterPlot('Irregular Domain with a Hole', true);
% End Irregular Domain
% -------------------------------------------------------------------------
end

function plotSquareDomain()

    % Create a square domain object (domain: [xmin, xmax, ymin, ymax])
    squareObj = Square([0, 1, 0, 1]);
    
    % Generate interior points (a 21x21 grid)
    squareObj = squareObj.generateSquare(21);
    
    % Generate boundary points (50 points per edge)
    squareObj = squareObj.generateBoundaryPoints(50);
    
    % Generate ghost points (3 layers, each with 50 points per edge, spread of 0.05)
    squareObj = squareObj.generateGhostPoints(50, 3, 0.05);
    
    % Plot the square domain with interior, boundary, and ghost points
    squareObj.scatterPlot('Square Domain', true, true);

end


% Function Calls for plots here
plotLShapeDomain();

plotCircularDomain();

plotIrregularDomain();

plotSquareDomain();