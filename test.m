new_L = LDomain([-3,3],[-3,0],[-3,0], 'bottom-left');
new_L = new_L.generateLShape(100);
new_L = new_L.generateBoundaryPoints(100);


new_L = new_L.generateGhostPoints(100, 2, 0.1);
new_L.scatterPlot("L-Shaped Domain",true, true);

% Create a circle at (0,0) with radius 1
%circleObj = CircularDomain([0, 0], 1);

% Generate ~50x50 grid points inside the circle
%circleObj = circleObj.generateCDomain(50);

% Generate 50 boundary points
%circleObj = circleObj.generateBoundaryPoints(50);

% Plot both interior (blue) and boundary (red) points
%circleObj.scatterPlot('My Circular Domain', true);

