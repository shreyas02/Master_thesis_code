clc
clear

%% Definition of the domain 

% fine domain - 
% Length of the domain
L = 1 ;
% Number of elements
nElem = 1000 ;
% Number of nodes per element
nne = 2 ;
% Boundary nodes
BCNodes = [1; nElem+1];
BCValues = [0; 0];
% Coordinates of the nodes
crd = [0:L/nElem:L]';
ndof = size(crd,1) ;

%% defining subdomains 

num  = ndof - 2;

%number of subdomains 
n = input("Number of Subdomains -");
overlap = 2;

%approx size of each subdomain 
sizef =  num./n ; 

%connectivity of each subdomain
con = zeros(n,2);
con(1,1) = 1; % starting value 
con(1,2) = 1 + sizef;
for i = 2:n
        con(i,1) = con(i-1,1) + sizef;
        con(i,2) = con(i,1) + sizef;
        if(con(i,2) > num)
            con(i,2) = num;
        end
end
con = floor(con);
for i = 2:n-1
    con(i,2) = con(i,2) + overlap;
    con(i,1) = con(i,1) - overlap;
end

%restriction operators 

for i = 1:n
        before = con(i,1)-1;
        after = num - con(i,2);
        sizef = con(i,2) - con(i,1);
        temp1 = zeros(sizef+1,before);
        temp2 = eye(sizef+1);
        temp3 = zeros(sizef+1,after);
        temp = [temp1,temp2,temp3];
        restriction{i} = temp;
end
clear temp temp1 temp2 temp3

%% Definition of the coarse space
% Number of elements in coarse space -
ndofC = n + 1;
step = ndof/ndofC; % approx step size 

conC = int32(linspace(1,ndof,ndofC));
crdC = crd(conC);

%% Construction of the coarse restiction operators 

r0 = zeros(ndofC-2,ndof-2);
for i = 2:ndofC-1
    xn = crdC(i);
    xnm1 = crdC(i-1);
    xnp1 = crdC(i+1);
    for j = 2:ndof-1
        x = crd(j);
        if x > xnm1 && x <= xn
            test = (x - xnm1)/(xn - xnm1);
            r0(i-1,j-1) = test;
            continue
        end
        if x > xn && x <= xnp1 
            test = (xnp1 - x)/(xnp1 - xn);
            r0(i-1,j-1)= test;
            continue
        end
    end
end

%% Finite element matrix construction

% Source term 
f = 1;

% Connectivity matrix of the elements
conn = zeros(nElem,nne); 
nn = 1;
for i=1:nElem
    nn = nn-1;
    for j=1:nne
        nn = nn+1; 
        conn(i,j) = nn;
    end
end

% Definition of Gauss Quadrature and Shape Function Space

% Location of Gauss points
gP = [-1/sqrt(3), 1/sqrt(3)] ;
% Weights of Gauss points
gW = [1,  1] ;
% Number of Gauss points
nQuad = length(gW) ;
% Shape functions
N(1,:) = 0.5.*(1-gP(1,:)) ;
N(2,:) = 0.5.*(1+gP(1,:)) ;
% Gradient of shape functions
Nx(1,:) = -0.5.*ones(1,2) ;
Nx(2,:) =  0.5.*ones(1,2) ;

% Formation of Local to Global DOF Mapping

ii = zeros(nne^2*nElem,1); 
jj = zeros(nne^2*nElem,1);
index = 0;
for i = 1:nne
   for j = 1:nne
      ii(index+1:index+nElem) = double(conn(:,i)); 
      jj(index+1:index+nElem) = double(conn(:,j));  
      index = index + nElem;
   end
end

%Initialization of Solution Vector 
ndof = size(crd,1) ;
Sol.u = zeros(ndof,1);

% Satisfy boundary conditions
Sol.u(BCNodes) = BCValues ;

% Initialization for Element-level Calculations 

xx = zeros(size(conn));
yy = zeros(size(conn));

for i=1:nne
   xx(:,i) = crd(conn(:,i),1);
end

%%%%% Element-level Evaluation of Matrices %%%%%conn

% Gauss quadrature loop
for p = 1:nQuad  
    
    % Jacobian evaluation and its absolute value
    J = [xx*Nx(:,p)];
    absJac = abs(J);
    
    % Evaluation of gradients of shape functions
    DNDx = (1./absJac)*Nx(:,p)' ;  
            
    index = 0;
    for i = 1:nne
        for j = 1:nne
            
            % Galerkin diffusion term
            Kij = gW(p)*(DNDx(:,i).*DNDx(:,j));
            Kij = Kij.*absJac;
            sK(index+1:index+nElem,p) = Kij; 
                    
            % Galerkin source term
            Fij = gW(p)*(N(i,p)*N(j,p));
            Fij = Fij.*absJac ;
            sF(index+1:index+nElem,p) = Fij;

            index = index + nElem;
        end
    end
end

% Summation of All Quadrature Data for Numerical Integration
sK = sum(sK,2);
sF = sum(sF,2);

%Assembly of Local Matrices to Global Matrix 
K= sparse(ii,jj,sK,ndof,ndof); 
F_global = sparse(ii,jj,sF,ndof,ndof); 

%Defining force matrix 
F = F_global*(f.*ones(ndof,1)) ;

%%%%% Incorporating Boundary Conditions to the RHS %%%%%
F = F - K*Sol.u ;

% Selecting unknown degrees of freedom for solution
unKnowns = setdiff([1:ndof]',BCNodes) ;

K = K(unKnowns,unKnowns);
F = F(unKnowns);

%% SOlving the system 

%calculating the preconditioner
precond = sparse(zeros(num,num));

parfor i = 1:n
    tempRes = cell2mat(restriction(i));
    temp = sparse(tempRes'*inverse(tempRes*K*tempRes')*tempRes);  %inverse function using LU factorization 
    precond = precond + temp;
end

temp = sparse(r0'*inverse(r0*K*r0')*r0);
precond = precond + temp;

clear temp tempRes

% conjugate gradient
[Sol.u(unKnowns),cgnorm,iter] = cg(K,F,precond);

%% post processing 

% Exact solution
x_Refined = [0:L/1000:L]';
u_exact = -(f/2).*x_Refined.^2 + (f*L/2).*x_Refined ;

figure(2);
plot(crd,Sol.u,':sk',x_Refined,u_exact,'-r','LineWidth',3,'MarkerSize',15);
set(gca,'TickLabelInterpreter','latex','FontSize',30);
xlabel('$x$','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
legend('Numerical','Exact','Interpreter','latex');

H = mean(con(:,2)-con(:,1));
O = overlap;
condition = condest(precond*K);
prop = 1 + H/O;
fprintf('Condition Number - : %e, \n', condition);
% fprintf("H = %e, Overlap = %e,\n ",H,O);
fprintf('Proportionality constant - : %e, \n', prop);

save 2x300.mat n iter

