clear all;
warning off;


addpath("functions/")
%PDE: laplacian u=-lambda u;
%B.C.: u=0;
%analytical eigenvalue parameters of the PDE (Note: the order of the pairs may not consistent with the order obtained from the numerical method)
m=[1 2 1 2 3 1 3 2 4 ]; %u(x,y)=sin(m*pi*x)*sin(l*pi*y) with lambda=(m^2+l^2)*pi^2
l=[1 1 2 2 1 3 2 3 1 ];
%m=[1 2 1 2 3 1 3 2 4 1 3 4 2 5 1]; %15
%l=[1 1 2 2 1 3 2 3 1 4 3 2 4 1 5];
% number of centers/collocation points(nodes)
nodes=[20];
%shape parameter values
% sp= [2.5]
sp = linspace(2, 5,1000);
% sp = [2.727]


%define a function that can calculate distances from x to y (x and y could be vectors containing nodes'coordinates )
D = @(x,y) hypot(bsxfun(@minus,x(:,1),y(:,1)'),bsxfun(@minus,x(:,2),y(:,2)'));


%define vectors to collect numerial results (eigenvalues, errors of numerical eigenvalues and eigenmodes)
smallest_eigenvalues=zeros(length(nodes),length(sp),length(m));
Err_eigenvalues=zeros(length(nodes),length(sp),length(m));
Err_eigenmodes_max=zeros(length(nodes),length(sp),length(m));
Err_eigenmodes_relative=zeros(length(nodes),length(sp),length(m));
Err_eigenmodes_rms=zeros(length(nodes),length(sp),length(m));


%count the time of the main part of the code
t0=clock;

%% ---------------------------------------------------------

% For each shape parameter, and for each nodes, compute the global coefficient matrix and its first several eigenvalues/eigenmodes
ii=0;
%per nodes' size
for nn=nodes
  ii=ii+1;jj=0;


  %define centers/collocation points in coor (here centers and collocation points are the same set)
  t=linspace(0,1,nn+1);
  [X,Y]=meshgrid(t,t);
  x1=reshape(X,[],1);y1=reshape(Y,[],1);
  % coor has three columes:col 1 and 2: x and y values of nodes; col 3: index for interior or boundary pts
  coor=[x1 y1 zeros(length(x1),1)];
  p1=find(y1==0 & x1~=1);p2=find(x1==1 & y1~=1);q1=find(y1==1 & x1~=0);q2=find(x1==0 & y1~=0);
  coor(p1,3)=ones(length(p1),1);coor(p2,3)=2*ones(length(p2),1);
  coor(q1,3)=3*ones(length(q1),1);coor(q2,3)=4*ones(length(q2),1);
  %define interior and boundary nodes' indices
  int_ind=find(coor(:,3)==0);intnode=coor(int_ind,1:2);% interior points
  int_ind_b=find(coor(:,3)>0);
  bdpt=coor(int_ind_b,1:2);% boundary points
  n=length(coor(:,1));ni=length(intnode(:,1));nb=n-ni;% # of interior, boundary, and total points

  [x_ghost, y_ghost] = fabric_pattern(n, 2.6, 0.5, 0.5, 0.8);

  ghost_pts = [x_ghost, y_ghost];
  %%plot all nodes
  % figure('Units', 'inches', 'Position', [1 1 6 4]); % 6 inches wide, 4 inches high

  % plot(coor(:,1),coor(:,2), 'bo','MarkerSize', 3);
  % hold on;

  % plot(x_ghost, y_ghost, 'ro', 'MarkerSize', 3);
  % title("Domain and Ghost Points")
  % xlabel("X");
  % ylabel("Y");

  % axis equal;
  % print('images/ghost_domain_r_8.png', '-dpng', '-r300');




  %per shape parameter
  for c=sp
    jj=jj+1;
    % Global RBFCM/MAPS algorithm
    cc=c;% rescaled shape parameter
    DM= D(coor(:,1:2), ghost_pts); %Calculating the distance matrix
    MQ = sqrt(cc^2*DM(int_ind,:).^2 + 1);  %Calculating the RBF \phi from the distance vector
    mq=sqrt(cc^2*DM.^2+1); %Calculating the RBF \phi from the distance matrix
    L_MQ=2*cc^2./MQ-cc^4*DM(int_ind,:).^2./(MQ.^3);%Calculating the L\phi from the distance vector
    A=L_MQ/mq; %Calculating the local weights \alpha

    %extract entries corresponding to interior nodes
    A=A(:,int_ind);


    % compute numerical eigencalues/eigenmodes:lambda stores the first length(m)
    % eigenvalues; alpha stores the corresponding eigenmode function values at interior nodes
    [alpha, lambda] = eigs(-A,length(m),0); % solving the eigenvalue problem p. 8 Driscol
    lam= real(diag(lambda)');% extract the first length(m) eigenvalues;
    alpha=real(alpha);
    %End of localized RBFCM/LMAPS


    %  fprintf('The first %d exact eigenvalues are:\n %s.\n', length(m),num2str(pi^2*(m.^2.+l.^2)))
    %  fprintf('The first %d numerical eigenvalues are:\n %s.\n', length(m),num2str(lam))
    %  fprintf('Relative errors of the first %d exact eigenvalues are:\n %s.\n', length(m),num2str(abs(lam-pi^2*(m.^2.+l.^2))./(pi^2*(m.^2.+l.^2))))
    %  format short
    %  fprintf('The first %d exact eigenvalues/numerical eigenvalues/Relative errors are:\n',length(m))
    %  [pi^2*(m.^2.+l.^2)' lam' (abs(lam-pi^2*(m.^2.+l.^2))./(pi^2*(m.^2.+l.^2)))']


    %compute the exact/analytical eigenmodes and compare them with the numerical eigenmodes
    firsteigmode=sin(pi*coor(:,1)*m).*sin(pi*coor(:,2)*l);
    normf=sqrt(sum(firsteigmode.^2)); %compute the maginititudes of eigenmodes because the numerical eigenmodes were defaultly normalized


    %plot exact/numerical eigenmodes and their differences when the first node
    %number and the first sp are used
    %  if ii==1&&jj==1
    %    k=2;%plot the second exact/numerical eigenmodes
    %    % plot the second exact eigenmodes
    %    figure
    %    scatter3(coor(:,1),coor(:,2),firsteigmode(:,k),3);
    %    axis([0,1,0,1,-1,1]);
    %    % plot the second exact eigenmodes
    %    figure
    %    coor(:,4)=zeros(n,1);%Define a new column for coor to store this numerical eigenmode
    %    coor(int_ind,4)=alpha(1:ni,k);%Assign this numerical eigenmode values at interior modes (it remains zero at boundary nodes due to the homogeneous boundary conditions)
    %    scatter3(coor(:,1),coor(:,2),normf(k)*coor(:,4),3);
    %    axis([0,1,0,1,-1,1]);
    %    %plot  their differences (errors) at all nodes
    %    figure
    %    scatter3(coor(:,1),coor(:,2), abs(abs(firsteigmode(:,k))-normf(k)*abs(coor(:,4))),3);%toc
    %  end


    %collecting data of eigenvalues
    smallest_eigenvalues(ii,jj,:)=lam;
    Err_eigenvalues(ii,jj,:)=abs(lam-pi^2*(m.^2.+l.^2))./(pi^2*(m.^2.+l.^2));

    %  collecting data of eigenmodes
    for k=1:length(m)
      firsteigmode_interior=firsteigmode(int_ind,k);
      error_eigenmode_k=abs(abs(firsteigmode_interior)-normf(k)*abs(alpha(1:ni,k)));
      Err_eigenmodes_max(ii,jj,k)=max(error_eigenmode_k);
      Err_eigenmodes_relative(ii,jj,k)=max(error_eigenmode_k/max(firsteigmode_interior));
      Err_eigenmodes_rms(ii,jj,k)=sqrt(sum(error_eigenmode_k.^2)/ni);
    end




  end % end of for c=sp
end %end of for nn=nodes





cpu=etime(clock,t0);
fprintf('cputime= %6.2f\n',cpu);

% eigenmode_Max_errors = zeros(1, length(sp));
% for i=1:length(sp)
%     eigenmode_Max_errors(i) = max(Err_eigenmodes_max(1,i,:));
% end

% [min_error, min_idx] = min(eigenmode_Max_errors);
% best_shape = sp(min_idx);

% figure;
% plot(sp, eigenmode_Max_errors, 'o-');
% hold on;

% % Highlight the minimum error point with a distinct marker
% plot(best_shape, min_error, 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Red circle

% % Add a text label with an arrow
% text(best_shape, min_error, ...
%     sprintf('\\leftarrow Minimum Error: %.2e\n(RADIUS = %.2f)', min_error, best_shape), ...
%     'VerticalAlignment', 'middle', ...
%     'HorizontalAlignment', 'right', ...
%     'FontSize', 10);
% xlabel('Shapes');
% ylabel('Max Error');
% title('Max Error Eigenmodes)');

%end counting time
%--------------------------------------------------------------------

color=['r','b','m','g','m','c','y','k','k'];
% I.----------------------------------------------------




% II.----------------------------------------------------
%plot the absolution max error of eigenmodes w.r.t. sp and nodes
% figure
% subplot(1,2,1) %the first plot in this figure
%   for ii=1:length(nodes)
%     jj=1;%plot the results when the first sp is used
%     scatter3((nodes(ii)+1)^2*ones(1,length(m)), 1:length(m), log10(squeeze(Err_eigenmodes_max(ii,jj,:))), 'filled', color(ii));
%     hold on;
%     % Use plot3 to connect the points with lines
%     plot3((nodes(ii)+1)^2*ones(1,length(m)), 1:length(m), log10(squeeze(Err_eigenmodes_max(ii,jj,:))), '-o', 'LineWidth', 1.5);
%     hold on;
%   end;
%   hold off;
% % Labeling the axes
% xlabel('Number of nodes');
% ylabel('Indices of eigenvalues');
% zlabel('log10(Max Error of eigenmodes)');
% title('3D Scatter Plot of max error of eigenmodes vs. the number of nodes');
% grid on;


% subplot(1,2,2)
%   for jj=1:length(sp)
%     ii=1;%plot the results when the first #of nodes is used
%     scatter3(sp(jj)*ones(1,length(m)), 1:length(m), log10(squeeze(Err_eigenmodes_max(ii,jj,:))), 'filled', color(jj));
%     hold on;
%     % Use plot3 to connect the points with lines
%     plot3(sp(jj)*ones(1,length(m)), 1:length(m), log10(squeeze(Err_eigenmodes_max(ii,jj,:))), '-o', 'LineWidth', 1.5);
%     hold on;
%   end;
%   hold off;
% % Labeling the axes
% xlabel('Values of the RBF shape parameter');
% ylabel('Indices of eigenvalues');
% zlabel('log10(Max Error of eigenmodes)');
% title('3D Scatter Plot of max error of eigenmodes vs. the shape parameter');
% grid on;






% plot the error of eigenvalues w.r.t. sp and nodes
% figure
% subplot(1,2,1) %the first plot in this figure
%   for ii=1:length(nodes)
%     jj=1;%plot the results when the first sp is used
%     scatter3((nodes(ii)+1)^2*ones(1,length(m)), 1:length(m), (squeeze(Err_eigenvalues(ii,jj,:))), 'filled', color(ii));
%     hold on;
%    % Use plot3 to connect the points with lines
%     plot3((nodes(ii)+1)^2*ones(1,length(m)), 1:length(m), (squeeze(Err_eigenvalues(ii,jj,:))), '-o', 'LineWidth', 1.5);
%     hold on;
%   end;
%   hold off;
% % Labeling the axes
% xlabel('Number of nodes');
% ylabel('Indices of eigenvalues');
% zlabel('log10(Errors of eigenvalues)');
% title('3D Scatter Plot of error of eigenvalues vs. the number of nodes');
% grid on;


% subplot(1,2,2)
%   for jj=1:length(sp)
%     ii=1;%plot the results when the first #of nodes is used
%     scatter3(sp(jj)*ones(1,length(m)), 1:length(m), log10(squeeze(Err_eigenvalues(ii,jj,:))), 'filled', color(jj));
%     hold on;
%     % Use plot3 to connect the points with lines
%     plot3(sp(jj)*ones(1,length(m)), 1:length(m), log10(squeeze(Err_eigenvalues(ii,jj,:))), '-o', 'LineWidth', 1.5);
%     hold on;
%   end;
%   hold off;
% % Labeling the axes
% xlabel('Values of the RBF shape parameter');
% ylabel('Indices of eigenvalues');
% zlabel('log10(Errors of eigenvalues)');
% title('3D Scatter Plot of error of eigenvalues vs. the shape parameter');
% grid on;

Err_eigenvalues = (Err_eigenmodes_max);

eigenvalue_Max_errors = zeros(1, length(sp));
for i=1:length(sp)
  eigenvalue_Max_errors(i) = max(Err_eigenvalues(1, i, :));
end

[min_error, min_idx] = min(eigenvalue_Max_errors);
best_shape = sp(min_idx);

figure('Units', 'inches', 'Position', [1 1 6 4]); % 6 inches wide, 4 inches high
plot(sp, (eigenvalue_Max_errors), 'r-', 'LineWidth', 1.8);
hold on;

% Highlight the minimum error point with a distinct marker
plot(best_shape, min_error, 'bo', 'MarkerSize', 8, 'LineWidth', 5); % Red circle

% Add a text label with an arrow
text(best_shape, min_error, ...
  sprintf('Max Err: %.2e\n(c = %.6f)', min_error, best_shape), 'FontSize', 10);
xlabel('Shape Parameter');
ylabel('(Max  Absolute Error)');

title('Ghost Method (Max Absolute Error vs Shape Parameter)');
print('images/ghost_eigenmode_best_shape_r_0.8.png', '-dpng', '-r300');

% display(max(log10(Err_eigenvalues(:))));