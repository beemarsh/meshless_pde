% Franke's Function
f = @(x) 0.75*exp(-1/4*(9*x(1)-2)^2 - 1/4*(9*x(2)-2)^2) + ...
    0.75*exp(-(9*x(1)+1)^2/49-(9*x(2)+1)/10) + ...
    0.5*exp(-(9*x(1)-7)^2/4-(9*x(2)-3)^2/4) - ...
    0.2*exp(-(9*x(1)-4)^2 - (9*x(2)-7)^2);

% Points from 0 to 1 with a 0.01 steps
t_evaluate = 0:0.01:1;

% Creating a grid
[X,Y]=meshgrid(t_evaluate);

% Evaluating the Z Value usinf Franke's Function
Z = zeros(size(X))
for i = 1:length(t_evaluate)
    for j = 1:length(t_evaluate)
        Z(i,j) = f([X(i,j); Y(i,j)]);
    end
end 

% Creating a surface plot of the actual value
figure
surf(X,Y,Z);
title("Plot of Actual Values of the Franke's Function")
hold on;

% Generating a random 2 dimensional data using Halton Sequence
num_samples = 100;
halton_seq = haltonset(2);
random_points = net(halton_seq, num_samples)'; % A column vector representing random x and y values from halton sequence

Z_samples = zeros(num_samples,1);

for i=1:num_samples
    Z_samples(i) = f(random_points(:,i));
end

% Scatter Plotting the sampled points in the 3D surface
plot3(random_points(1,:), random_points(2,:), Z_samples, 'o', 'LineWidth',6);


% Kernel Functions and Shape Parameter
mu = 2;
%K = @(x,y) exp(1/mu*x'*y);
K = @(x,y) exp(-1/mu*norm(x-y)^2);

% Gram Matrix
GramMatrix = zeros(num_samples);

for i = 1:num_samples
    for j = 1:num_samples
        GramMatrix(i,j) = K(random_points(:,i),random_points(:,j));
    end
end

% Weights

Weights = pinv(GramMatrix)*Z_samples;

% Evaluating the Approximation
Z_approximation = zeros(size(X));

for i = 1:length(t_evaluate)
    for j = 1:length(t_evaluate)
        sum = 0;
        for ii = 1:num_samples
            sum = sum + Weights(ii)*K([X(i,j);Y(i,j)],random_points(:,ii));
        end
        Z_approximation(i,j) = sum;
    end
end

% Plot

figure

surf(X,Y,Z_approximation);
title('RBF Approximation of Frankes Function')

figure

surf(X,Y,abs(Z-Z_approximation));
title('Error between true and approximated values')