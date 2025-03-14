m=[1 2 1 2 3 1 3 2 4 ];
l=[1 1 2 2 1 3 2 3 1 ];

%m=[1 2 1 2 3 1 3 2 4 1 3 4 2 5 1]; %15
%l=[1 1 2 2 1 3 2 3 1 4 3 2 4 1 5];


% number of centers/collocation points(nodes)
nodes=[20];


%shape parameter values
sp=[ 9];%[0.1:0.1:1.01 1:0.5:100];


%This function D calculates the distance between each points in the matrix
%x and y.
%Here, x is a n*2 matrix and y is a m*2 matrix.
% Each row denotes a set of x,y collocation points
% Using bsxfun, we calculate the difference between each points.
%Using hypot, we calculate the Euclidean distance.
%Essentially, all its doing is: sqrt( (x_1 - x_2)^2 + (y_1 - y_2)^2 )

D = @(x,y) hypot(bsxfun(@minus,x(:,1),y(:,1)'),bsxfun(@minus,x(:,2),y(:,2)'));


%define vectors to collect numerical results (eigenvalues, errors of numerical eigenvalues and eigenmodes)
% These are the variables where the final calculations will be stored.
smallest_eigenvalues=zeros(length(nodes),length(sp),length(m));
Err_eigenvalues=zeros(length(nodes),length(sp),length(m));
Err_eigenmodes_max=zeros(length(nodes),length(sp),length(m));
Err_eigenmodes_relative=zeros(length(nodes),length(sp),length(m));
Err_eigenmodes_rms=zeros(length(nodes),length(sp),length(m));


%count the time of the main part of the code
t0=clock;


%---------------------------------------------------------
% For each shape parameter, and for each nodes, compute the global coefficient matrix and its first several eigenvalues/eigenmodes
ii=0;

% For each node size, we loop
for nn=nodes

    ii=ii+1;
    jj=0;

    % All it does is just based on the number of nodes (nn), it will generate
    %points. For eg nn=20, will make a square grid of size 20*20.

    %define centers/collocation points in coor (here centers and collocation points are the same set)
    t=linspace(0,1,nn+1);
    [X,Y]=meshgrid(t,t);

    x1=reshape(X,[],1);
    y1=reshape(Y,[],1);

    % Here, we have 3 columns in coor.
    % First two columns are x and y values of nodes.
    % Third column is index for interior or boundary pts
    % Therefore, coor stores x values, y values, and index for interior or boundary pts

    % Don't worry about the code down here.
    % It is just an algorithm that is separating boundary points and interior
    % points

    % int_ind is the index of interior points in coor.
    % intnode are the actual interior collocation points.

    % int_ind_b is the index of boundary points in coor
    % bdpt are the actual boundary points.

    coor=[x1 y1 zeros(length(x1),1)];


    p1=find(y1==0 & x1~=1);
    p2=find(x1==1 & y1~=1);
    q1=find(y1==1 & x1~=0);
    q2=find(x1==0 & y1~=0);

    coor(p1,3)=ones(length(p1),1);
    coor(p2,3)=2*ones(length(p2),1);
    coor(q1,3)=3*ones(length(q1),1);coor(q2,3)=4*ones(length(q2),1);

    


    %define interior and boundary nodes' indices
    int_ind=find(coor(:,3)==0);
    intnode=coor(int_ind,1:2);% interior points
    int_ind_b=find(coor(:,3)>0);
    bdpt=coor(int_ind_b,1:2);% boundary points
    n=length(coor(:,1));ni=length(intnode(:,1));nb=n-ni;% # of interior, boundary, and total points

    %figure;
    %scatter(intnode(:,1),intnode(:,2), 15, 'filled', 'g');
    %hold on;
    %scatter(bdpt(:,1),bdpt(:,2), 15, 'filled', 'b');

    
    % Interior and Boundary Points Figuring Out Complete
    % ***************************************************

    
    %per shape parameter
    for c=sp
        jj=jj+1;
        % Global RBFCM/MAPS algorithm

        % We don't need to worry about the algorithm here as well.
        % As described in the paper, we are just finding matrix A and L.
        % The solution system is A*α = λ*B*α
        % The matrix A = [ L; P].
        % L is the for interior points and P is for boundary points

        % here, DM = sqrt( (x1-x2)^2 + (y1-y2)^2 ). Basically, euclidean norm.

        % MQ is the RBF. MQ = sqrt (1 + ξ^2 * r^2 ). Here, r^2 is basically DM.

        % For matrix L, we only do it for interior points. We, also need to take
        % the laplacian which transforms phi. It is reflected in L_MQ.

        cc=c;% re-scaled shape parameter
        DM= D(coor(:,1:2), coor(:,1:2)); %Calculating the distance matrix
        

        %Calculating the RBF (Only for interior points)
        MQ = sqrt(cc^2*DM(int_ind,:).^2 + 1);

        %Calculating the RBF \phi. (For all points)
        mq = sqrt(cc^2*DM.^2+1);

        % We take the laplacian of φ. This gives the matrix L. Note: This is
        % only for interior points
        L_MQ = 2*cc^2./MQ-cc^4*DM(int_ind,:).^2./(MQ.^3);

        

        % From the solution system, we get α values by α = B * inverse(A).
        % Here, B is L_MQ and A is mq. / operator finds the inverse.
        % Don't be confused because we are taking laplacian of B. If it was
        % Kansas method, we would take laplacian of A. But for MAPS we will
        % take laplacian of B.

       

        A = L_MQ/mq; %Calculating the local weights \alpha
      
        %extract entries corresponding to interior nodes
        A = A(:,int_ind);

        % Now, need to compute the eigenvalues
        % Here, alpha stores the corresponding eigen functions at interior
        % nodes.
        % lambda stores first M or length(m) eigenvalues.

        [alpha, lambda] = eigs(-A, length(m), 0);

        % extract the first length(m) eigenvalues;
        % We need to extract them from the diagonal.
        lam = real(diag(lambda)');

        alpha = real(alpha);

        % Now we are done approximating eigenmodes and eigenvalues
        % ******************************************************************

        exact_eigenvalues = pi^2*(m.^2.+l.^2)';
        approximate_eigenvalues = lam';

        relative_error = (abs(lam- pi^2*(m.^2.+l.^2))./(pi^2*(m.^2.+l.^2)))';

        %compute the exact/analytical eigenmodes and compare them with the numerical eigenmodes
        exact_eigenmode = sin(pi*coor(:,1)*m).*sin(pi*coor(:,2)*l);
        normf = sqrt(sum(exact_eigenmode.^2)); %compute the magnitudes of eigenmodes because the numerical eigenmodes were default normalized

        % Collecting data of eigenvalues
        % All this code is doing below is storing the values in the above
        % declared vector. We are also calculating errors and storing them.
        smallest_eigenvalues(ii,jj,:)=lam;
        Err_eigenvalues(ii,jj,:)=abs(lam-pi^2*(m.^2.+l.^2))./(pi^2*(m.^2.+l.^2));

        %collecting data of eigenmodes
        for k=1:length(m)
            firsteigmode_interior=exact_eigenmode(int_ind,k);
            error_eigenmode_k=abs(abs(firsteigmode_interior)-normf(k)*abs(alpha(1:ni,k)));
            
            Err_eigenmodes_max(ii,jj,k)=max(error_eigenmode_k);
            Err_eigenmodes_relative(ii,jj,k)=max(error_eigenmode_k/max(firsteigmode_interior));
            Err_eigenmodes_rms(ii,jj,k)=sqrt(sum(error_eigenmode_k.^2)/ni);
        end
    end
end