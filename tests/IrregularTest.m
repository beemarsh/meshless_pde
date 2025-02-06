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
domainObj = IrregularDomain(mainVertices, {hole1});

% Generate interior points
domainObj = domainObj.generateInteriorPoints(100);

% Generate boundary points
domainObj = domainObj.generateBoundaryPoints(20);

% Plot the domain
domainObj.scatterPlot('Irregular Domain with a Hole', true);
