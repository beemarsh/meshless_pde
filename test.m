new_L = LDomain([-3,3],[0,3],[-3,0], 'bottom-right');
new_L = new_L.generateLShape(100);
new_L = new_L.generateBoundaryPoints(100);
new_L.scatterPlot("L-Shaped Domain",true);