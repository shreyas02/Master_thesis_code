clc
clear all

%Definition of Discretized Domain
% for i=1:size(BCTop,1)
%     bb = ismember(BCNodes, BCTop(i,1));
%     BCValues(bb~=0) = 1 ;
% end
% Length of the domain
Lx = 1 ;
Ly = 1 ;
% Number of elements
nElemX = 20 ;
nElemY = 20;
nElem = nElemX*nElemY ;
nCrds = (nElemX+1)*(nElemY+1) ;
% Number of nodes per element
nne = 4 ;
% Source term
f =-1 ;
% Boundary nodes
BCBottom = [1:1:nElemX+1]';
BCLeft = [1:nElemX+1:nCrds]';
BCRight = [nElemX+1:nElemX+1:nCrds]';
BCTop = [nCrds-nElemX:1:nCrds]';
BCNodes = unique([BCBottom; BCLeft; BCRight; BCTop]);
BCValues = zeros(size(BCNodes,1),1);
% for i=1:size(BCTop,1)
%     bb = ismember(BCNodes, BCTop(i,1));
%     BCValues(bb~=0) = 1 ;
% end

% Coordinates of the nodes
x = linspace(0,Lx,nElemX+1);
y = linspace(0,Ly,nElemY+1);
[X,Y] = meshgrid(x,y);
crd = [reshape(X',nCrds,1) reshape(Y',nCrds,1)];

% Connectivity matrix of the elements
conn = zeros(nElem,nne); 
nn = 1;
a0 = 0 ;
a1 = 0 ;
for i=1:nElem
    conn(i,1) = i+a0 ;
    conn(i,2) = conn(i,1)+1 ;
    conn(i,3) = i+nElemX+2+a1 ;
    conn(i,4) = conn(i,3)-1 ;
    if (mod(i,nElemX)==0)
        a0 = a0+1 ;
        a1 = a1+1 ;
    end
end

%%%%% Definition of Gauss Quadrature and Shape Function Space %%%%%

% Location of Gauss points
gP = [-1/sqrt(3),  1/sqrt(3), 1/sqrt(3), -1/sqrt(3);
      -1/sqrt(3), -1/sqrt(3), 1/sqrt(3),  1/sqrt(3)];
% Weights of Gauss points
gW = [1,  1,  1,  1] ;
% Number of Gauss points
nQuad = length(gW) ;
% Shape functions
N(1,:) = 0.25.*(1-gP(1,:)).*(1-gP(2,:)) ;
N(2,:) = 0.25.*(1+gP(1,:)).*(1-gP(2,:)) ;
N(3,:) = 0.25.*(1+gP(1,:)).*(1+gP(2,:)) ;
N(4,:) = 0.25.*(1-gP(1,:)).*(1+gP(2,:)) ;
% Gradient of shape functions
Nx(1,:) = -0.25.*(1-gP(2,:)) ;
Nx(2,:) =  0.25.*(1-gP(2,:)) ;
Nx(3,:) =  0.25.*(1+gP(2,:)) ;
Nx(4,:) = -0.25.*(1+gP(2,:)) ;
Ny(1,:) = -0.25.*(1-gP(1,:)) ;
Ny(2,:) = -0.25.*(1+gP(1,:)) ;
Ny(3,:) =  0.25.*(1+gP(1,:)) ;
Ny(4,:) =  0.25.*(1-gP(1,:)) ;

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

%Initialization for Element-level Calculations 

xx = zeros(size(conn));
yy = zeros(size(conn));

for i=1:nne
   xx(:,i) = crd(conn(:,i),1);
   yy(:,i) = crd(conn(:,i),2);
end

%Element-level Evaluation of Matrices 

% Gauss quadrature loop
for p = 1:nQuad  
    
    % Jacobian evaluation and its absolute value
    J = [xx*[Nx(:,p)], yy*[Nx(:,p)],...
         xx*[Ny(:,p)], yy*[Ny(:,p)]];
    absJac =abs( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
    
    % Evaluation of gradients of shape functions 
    DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,2))*Ny(:,p)')./repmat(absJac,1,nne);
    DNDy = ((-J(:,3))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(absJac,1,nne);
            
    index = 0;
    for i = 1:nne
        for j = 1:nne
            
            % Galerkin diffusion term
            Kij_1 = -gW(p)*(DNDx(:,i).*DNDx(:,j));
            Kij_2 = -gW(p)*(DNDy(:,i).*DNDy(:,j));
            Kij_1 = Kij_1.*absJac;
            Kij_2 = Kij_2.*absJac;
            sK_1(index+1:index+nElem,p) = Kij_1;
            sK_2(index+1:index+nElem,p) = Kij_2;
                    
            % Galerkin source term
            Fij = gW(p)*(N(p,i)*N(p,j));
            Fij = Fij.*absJac ;
            sF(index+1:index+nElem,p) = Fij;

            index = index + nElem;
        end
    end
end

% Summation of All Quadrature Data for Numerical Integration
sK_1 = sum(sK_1,2);
sK_2 = sum(sK_2,2);
sF = sum(sF,2);

% Assembly of Local Matrices to Global Matrix 
K1_global = sparse(ii,jj,sK_1,ndof,ndof); 
K2_global = sparse(ii,jj,sK_2,ndof,ndof); 
F_global = sparse(ii,jj,sF,ndof,ndof); 

% Defining of K Matrix and RHS Vector 
K = K1_global + K2_global ;
RHS = F_global*(f.*ones(ndof,1)) ;

%  Incorporating Boundary Conditions to the RHS
RHS = RHS - K*Sol.u ;

% Selecting unknown degrees of freedom for solution
unKnowns = setdiff([1:ndof]',BCNodes) ;

%% subdomains 
temp = 0;
count = 1;
for i = 1:nElem
    for j = 1:nne
        if ismember(conn(i,j),BCNodes)
            temp = 1;
        end
    end
    if temp == 1
        temp = 0;
        continue
    end
    ids(count) = i;
    count = count + 1;
end

for i = 1 : size(ids,2)
    temp = floor(ids(i)/nElemX);
    tempy(i) = temp + 1;
    tempx(i) = ids(i) - temp*nElemX;
end

ids = [ids',tempx',tempy'];
clear tempx tempy

% size of the subdomain 
sizeX = 2;
sizeY = 2;

%overlap region 
overX = 1;
overY = 1;
nextX = sizeX - overX;
nextY = sizeY - overY;
subdLoc = [];
subd = [];

% all the element numbers for each subdomain  
for i = 1:size(ids,1)
    x = ids(i,2);
    y = ids(i,3);
    if (ids(i,2) + sizeX) > nElemX || (ids(i,3) + sizeY) > nElemY
            continue
    end
    subdLoc = [subdLoc,ids(i,1)];
    for j = 1:size(ids,1) 
        
        if (ids(i,2) + sizeX) <= ids(j,2) ||  (ids(i,3) + sizeY) <= ids(j,3)
            continue
        end
        if ids(j,2)>=x && ids(j,3) >=y && i ~= j  
            %write here
            subdLoc = [subdLoc,ids(j,1)];
        end
    end
    subd = [subd;subdLoc];
    subdLoc = [];
end

% node numbers in each subdomain 

subdnode = [];

for i = 1:size(subd,1)
    temp = [];
    for j = 1: size(subd,2)
        element = subd(i,j);
        temp = [temp,conn(element,:)];
    end
    temp = unique(temp);  
    subdnode = [subdnode;temp];
    clear temp
end

% construction of restriction operators 

for i = 1 : size(subd,1)
    temp = zeros(size(subdnode,2) ,ndof);
    for j = 1:size(subdnode,2)
        t = subdnode(i,j);
        temp(j,t) = 1;
    end
    restriction(:,:,i) = temp;
    clear temp
end

%number of subdomains 
n = size(subd,1);

% calculating the preconditioner 
precond = sparse(zeros(ndof,ndof));
for i = 1:n
    temp = sparse(restriction(:,:,i)'*inverse(restriction(:,:,i)*K*restriction(:,:,i)')*restriction(:,:,i));  %inverse function using LU factorization 
    precond = precond + temp;
end
clear temp

precond = precond(unKnowns,unKnowns);

% taking the unkowns 
K = K(unKnowns,unKnowns);
RHS = RHS(unKnowns);

% conjugate gradient
[Sol.u(unKnowns),cgnorm,iter] = cg(K,RHS,precond);

%%  Post-processing 
figure(1);
mesh(X',Y',reshape(Sol.u,(nElemX+1),(nElemY+1)));
set(gca,'TickLabelInterpreter','latex','FontSize',30);
g1 = colorbar;
colormap('jet');
set(g1,'TickLabelInterpreter','latex','FontSize',30);
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$u$','Interpreter','latex');

condition = condest(precond*K);
size1 = sqrt(sizeY*sizeX);
overlap = sqrt(overY*overX);
Prop = 1 + 1/(size1*overlap);
fprintf('Cndition Number - : %e, \n', condition);
fprintf('Proportionality constant - : %e, \n', condition);

% save xx6.mat Prop condition
save xy8.mat n iter















