% Create a circle at (0,0) with radius 1
circleObj = CDomain([0, 0], 1);

% Generate ~50x50 grid points inside the circle
circleObj = circleObj.generateCDomain(50);

% Generate 50 boundary points
circleObj = circleObj.generateBoundaryPoints(50);

% Plot both interior (blue) and boundary (red) points
circleObj.scatterPlot('My Circular Domain', true);