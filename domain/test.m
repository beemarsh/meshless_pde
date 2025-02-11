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
