addpath('./domain/')
clear all;
clear;


% Kansas method on a square domain.

%First define the domain
% Unit square domain
square_domain = Square([0, 1, 0, 1]);
square_domain = square_domain.generateGrid(20);
%square_domain.scatterPlot("Square Domain", true, false);

% This contains two columns, first X values and second Y values
collocation_points = square_domain.SPoints;
boundary_points = square_domain.BoundaryPoints;
all_points = [collocation_points; boundary_points];

num_collocation_points = length(collocation_points(:,1));
num_boundary_points = length(boundary_points(:,1));
num_total_points = num_collocation_points + num_boundary_points;


X = square_domain.Grid.X;
Y = square_domain.Grid.Y;

%This function D calculates the distance between each points in the matrix x and y.
%Here, x is a n*2 matrix and y is a m*2 matrix.
% Each row denotes a set of x,y collocation points
% Using bsxfun, we calculate the difference between each points.
% Using hypot, we calculate the Euclidean distance.
% Essentially, all its doing is: sqrt( (x_1 - x_2)^2 + (y_1 - y_2)^2 )
D = @(x,y) hypot(bsxfun(@minus,x(:,1),y(:,1)'),bsxfun(@minus,x(:,2),y(:,2)'));

% We try for different shape parameters
shapes = [9 10 11];

m=[1 2 1 2 3 1 3 2 4 ];
l=[1 1 2 2 1 3 2 3 1 ];

% Here, we have store the result for eigenvalues in 3 columns:
% approximate eigenvalues
% exact eigenvalues
% error_eigenvalues
result_eigenvalues = zeros([length(m), 3]);

% The following variables hold results for eigenmode
% Stores the errors in eigenmodes.
max_err_eigenmodes = zeros(length(shapes), length(m));
relative_err_eigenmodes = zeros(length(shapes), length(m));
rms_err_eigenmode = zeros(length(shapes), length(m));

% Keeps count when looping
count = 0;
% For each different shape parameter, we approximate values.

shape = shapes(1);
%for shape = shapes

    count = count+1;

    % We find the distance matrix
    DM = D(all_points,all_points);

    % Calculate RBF MQ = sqrt (1 + ξ^2 * r^2 )
    % Only for collocation points
    MQ = sqrt( ( shape^2 * DM(1:num_collocation_points,:).^2 ) + 1 );

    % Now, RBF for all points
    mq = sqrt( ( shape^2 * DM.^2 ) + 1 );

    % Now we take the laplacian of MQ RBF.
    L_MQ = ( 2 * (shape ^ 2)./ MQ ) - ( (shape ^ 4) * ( DM(1:num_collocation_points,:).^2)./(MQ.^3) );

    % We have a linear system: -A α = λ B α
    % We recast this form into standard eigenvalue problem.
    % Let V = A inv(B)
    % Thus, we get the system: -V α = λ α
    V = L_MQ / mq;
    
    % To calculate eigenvalues, first we make the matrix V square.
    % For that, we extract the entries corresponding to interior
    % collocation points
    V = V(:,1:num_collocation_points);



    % Now, we calculate the eigenvalues.
    % Note that we are using negative because our equation has it.
    [alpha, lambda] = eigs(-V, length(m), 0);

    % extract the first length(m) eigenvalues;
    % We need to extract them from the diagonal.

    approximate_eigenvalues = real(diag(lambda));

    alpha = real(alpha);

    % Compute the exact eigenvalues and then calculate the relative error
    exact_eigenvalues = pi^2*(m.^2.+l.^2)';
    relative_error = (abs(approximate_eigenvalues' - pi^2*(m.^2.+l.^2))./(pi^2*(m.^2.+l.^2)))';

    % Store approximate eigenvalues, exact eigenvalues and relative error
    % in eigenvalues.
    result_eigenvalues(:,1:3)=[approximate_eigenvalues, exact_eigenvalues, relative_error];

    % Compute the exact eigenmodes
    exact_eigenmode = sin(pi*all_points(:,1)*m).*sin(pi*all_points(:,2)*l);
    normf = sqrt(sum(exact_eigenmode.^2));

    % After computing eigenmodes, we will calculate errors: relative, max,
    % and RMS.
    for k=1:length(m)
            firsteigmode = exact_eigenmode(1:num_collocation_points,k);
            error_eigenmode_k = abs(abs(firsteigmode)-normf(k)*abs(alpha(1:num_collocation_points,k)));

            max_err_eigenmodes(count, k) = max(error_eigenmode_k);

            relative_err_eigenmodes(count, k) = max(error_eigenmode_k/max(firsteigmode));

            rms_err_eigenmode(count, k) = sqrt(sum(error_eigenmode_k.^2)/num_collocation_points);
    end


    
   %plot the contours of the exact eigenmodes
   figure 
   for j=1:9
    j_th_eigenmode = exact_eigenmode(:,j);
    Z = reshape(j_th_eigenmode,size(X));
    subplot(3,3,j)
    %contourf(X,Y,Z,20);
    surf(X,Y,Z);
    hold on;
    axis equal;
    axis([0,1,0,1]);
    title("Exact Eigenmodes for Shape Parameter : " + shape);
   end
   
   % Numerical Eigenmodes
figure
   for j=1:9
    num_eigenmodes = zeros(num_total_points,1);
    num_eigenmodes(1:num_collocation_points,:) = normf(j)*alpha(1:num_collocation_points,j);
    num_eigenmodes(num_collocation_points+1:num_total_points) = exact_eigenmode(num_collocation_points+1, j);
    Z = reshape(num_eigenmodes , size(X));
    subplot(3,3,j)
    %contourf(X,Y,Z,20);
    surf(X,Y,Z);
    hold on;
    axis equal;
    axis([0,1,0,1]);
    title("Numerical Eigenmodes for Shape Parameter : " + shape);
   end

    % Plot the absolute contour errors of the eigenmodes
   figure 
   for j=1:9
    errors = zeros(num_total_points,1);
    errors(1:num_collocation_points,:) = abs(abs(exact_eigenmode(1:num_collocation_points,j))-abs(normf(j)*alpha(1:num_collocation_points,j)));
    
    Z=reshape(errors,size(X));
    subplot(3,3,j)
    %contourf(X,Y,Z,20);
    surf(X,Y,Z);
    hold on;
    axis equal;
    axis([0,1,0,1]);
    title("Absolute Errors for Shape Parameter : " + shape);
   end

%end


% We have saved eigenvalues and their errors. We can always plot them here.