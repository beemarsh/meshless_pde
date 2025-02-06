new_L = LDomain([-3,3],[-3,0],[-3,0], 'bottom-left');
new_L = new_L.generateLShape(40);
new_L = new_L.generateBoundaryPoints(40);

new_L = new_L.generateGhostPoints(10, 2, 0.1);
new_L.scatterPlot("L-Shaped Domain",true, true);

% % Create a circular domain of radius 1 centered at (0,0)
% circleObj = CircularDomain([0, 0], 1);
% % Generate ~50x50 interior grid
% circleObj = circleObj.generateCDomain(50);
% % Generate 50 boundary points
% circleObj = circleObj.generateBoundaryPoints(50);
% % Generate ghost points:
% %   width = 5 offset circles
% %   spread = 0.05 away each
% %   50 points per circle
% circleObj = circleObj.generateGhostPoints(5, 0.05, 50);
% % Visualize everything
% % pass true for boundary, and true for ghost
% circleObj.scatterPlot('Circle with Ghost Points', true, true);

