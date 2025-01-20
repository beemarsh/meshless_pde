# Welcome!

This is a part of my research with Dr. Huiqing Zhu. It is about  using meshless methods to approximate PDEs. 

In this repo, we have programmed only for the need of the research. We are aware that the code isn't efficient and can be optimized for better performance.

We have used various domains: L-shaped, Circular, Spherical, and Irregular. Different classes are implemented in the code. It is easy to import the class and use at any time. The short implementation details for the classes are listed below.

### LDomain
This class is used to create a square shaped L-shape domain.
- Arguments: `(Domain, XBounds, YBounds, orientation)`
	- `Domain`: This is an array of two numbers. It is the complete square domain, eg: `[-2,2] ` forms a square domain from `x=-2, y=-2` to `x=2,y=2`
	- `XBounds`: This is an array of two numbers. The range provided here will be used as a reference to remove the square that forms an L-shape. eg: `[0,2]` means the points in the range 0 to 2 will be omitted that forms an L-shape.
	- `YBounds` : This is an array of two numbers. The range provided here will be used as a reference to remove the square that forms an L-shape. eg: `[0,2]` means the points in the range 0 to 2 will be omitted that forms an L-shape.
	- `orientation` : It takes only four values namely: `bottom-right`,`bottom-left`,`top-right`, and `top-left`. The orientation decides where the shape of L will be formed.
		-	```matlab
			new_L = LDomain([-3,3],[-3,0],[-3,0], 'bottom-left');
			```
- Methods:
	- `generateLShape(numPoints)`: This function will generate collocation points inside the given domain. It first forms a grid of size `numPoints*numPoints`. For `numPoints=10`, it will form a 10x10. The grid is numbered with numbers>=0. The code will populate the locations in the grid that are not equal to zero with values linearly spaced inside the domain. This function will only generate collocation points inside the boundary but not on the boundary iteself.
		-	```matlab
			new_L = LDomain([-3,3],[-3,0],[-3,0], 'bottom-left');
			new_L = new_L.generateLShape(100); % A 100*100 grid
			```

	-	`generateBoundaryPoints(numPoints)`: This function will generate points only on the boundary. It will generate `numPoints` points on each side of the boundary. For example: for `numPoints=10`, there will be 10 points on each side of the square. However, there will be `numPoints/2` points where the domain forms an L-shape. The points are linearly spaced just like collocation points.
		-	```matlab
			new_L = LDomain([-3,3],[-3,0],[-3,0], 'bottom-left');
			new_L = new_L.generateBoundaryPoints(100);
			```

	-	`generateGhostPoints(numPoints, width, spread)`: This function will generate points just outside the boundary. The `width` parameter defines how many lines of points to be made outside the boundary. `spread` parameter defines how close the lines are to each other or near the boundary for the first line.
		-	```matlab
			new_L = LDomain([-3,3],[-3,0],[-3,0], 'bottom-left');
			new_L = new_L.generateGhostPoints(100, 2, 0.1);
			```

	-	`scatterPlot(plot_title, boundary, ghost)`: This function scatter plot the collocation points. If `boundary` is set to `true`, it will plot the boundary points as well. Similarly, if `ghost` is set to `true`, it will plot the ghost points.
		-	```matlab
			new_L = LDomain([-3,3],[-3,0],[-3,0], 'bottom-left');

			new_L = new_L.generateLShape(100);

			new_L = new_L.generateBoundaryPoints(100);

			new_L = new_L.generateGhostPoints(100, 2, 0.1);

			new_L.scatterPlot("L-Shaped Domain",true, true);
			```

![L-shaped Domain Example Graph](https://github.com/beemarsh/meshless_pde/blob/main/images/L_shape_example.png?raw=true)

