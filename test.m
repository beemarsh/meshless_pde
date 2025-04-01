% Generate some sample data
x = rand(10,1);  % 10 random x points
y = rand(10,1);  % 10 random y points
z = rand(10,1);  % 10 random z points

% Create scatter plot
figure;
scatter3(x, y, z, 'filled');  % 'filled' makes solid markers
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
grid on;