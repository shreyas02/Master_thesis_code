include("scicomp.jl")
using Plots
using SparseArrays
using LinearAlgebra

#Length of the Domain 
Lx = 5;
Ly = 5;
#Number of elements 
nElemX = 10;
nElemY = 10;
nElem = nElemX*nElemY;
nCrds = (nElemX+1).*(nElemY+1);
#Number of nodes per element 
nne = 4;
#Source Term 
f = 1; 
#Boundary Nodes 
BCBottom = linArr(1,nElemX + 1,1);
BCLeft = linArr(1,nCrds,nElemX+1);
BCRight = linArr(nElemX+1, nCrds, nElemX+1);
BCTop = linArr(nCrds-nElemX,nCrds,1);
BCNodes = sort(unique([BCBottom;BCLeft;BCRight;BCTop]));
BCValues = zeros(size(BCNodes,1),1);
#Coordinates of the nodes
x = LinRange(0,Lx,nElemX + 1);
y = LinRange(0,Ly,nElemY + 1);
X,Y = meshgrid(x,y);
crd = [reshape(X',nCrds,1) reshape(Y',nCrds,1)];
#Connectivity matrix of the elements 
conn = rand(nElem,nne).*0;
nn = 1;
let a0 = 0;
let a1 = 0; 
for i = 1:nElem
    conn[i,1] = i + a0;
    conn[i,2] = conn[i,1] + 1;
    conn[i,3] = i + nElemX + 2 + a1;
    conn[i,4] = conn[i,3] - 1;
    if mod(i,nElemX) == 0
        a0 = a0+1;
        a1 = a1+1;
    end
end
end
end
#Location of Gauss points 
gP = [ -1/sqrt(3)  1/sqrt(3) 1/sqrt(3) -1/sqrt(3) ; -1/sqrt(3)  -1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];
#Weights of gauss points 
gW = [1 1 1 1];
#Number of Gauss points 
nQuad = length(gW);
#Shape Functions 
N = rand(4,4).*0; #initialization
N[1,:] = 0.25.*(1 .- gP[1,:]).*(1 .- gP[2,:]);
N[2,:] = 0.25.*(1 .+ gP[1,:]).*(1 .- gP[2,:]);
N[3,:] = 0.25.*(1 .+ gP[1,:]).*(1 .+ gP[2,:]);
N[4,:] = 0.25.*(1 .- gP[1,:]).*(1 .+ gP[2,:]);
#Gradient of shape functions 
Nx = rand(4,4).*0; #initialization
Nx[1,:] = -0.25.*(1 .- gP[2,:]);
Nx[2,:] = 0.25.*(1 .- gP[2,:]);
Nx[3,:] = 0.25.*(1 .+ gP[2,:]);
Nx[4,:] = -0.25.*(1 .+ gP[2,:]);
Ny = rand(4,4).*0; #initialization
Ny[1,:] = -0.25.*(1 .- gP[1,:]);
Ny[2,:] = -0.25.*(1 .+ gP[1,:]);
Ny[3,:] = 0.25.*(1 .+ gP[1,:]);
Ny[4,:] = 0.25.*(1 .- gP[1,:]);
#Formation of local to global DOF mapping 
ii = zeros(nne^2*nElem,1);
jj = zeros(nne^2*nElem,1);
let index = 0;
for i = 1:nne
    for j = 1:nne
        ii[index+1:index+nElem] = conn[:,i];
        jj[index+1:index+nElem] = conn[:,j];
        index = index+nElem;
    end
end
end
#Initialization of the soultion vector 
ndof = size(crd,1);
Solu = zeros(ndof,1);
#All Nodes for this consructed mesh
nodes = rand(ndof,1).*0;
let test = 1;
for i = 1:ndof
    nodes[i] = test;
    test = test + 1;
end
end
#Satisfy boundary conditions 
Solu[BCNodes] = BCValues;
#Initialization for Element 
xx = zeros(size(conn));
yy = zeros(size(conn));
for i = 1:nElem
    for j = 1:nne
        xx[i,j] = crd[Int(conn[i,j]),1];
        yy[i,j] = crd[Int(conn[i,j]),2];
    end
end
#Element level evaluation of matrices 
#Initialization of the Jacobian Matrix 
J = rand(nElem,nne).*0;
absJac = rand(nElem,1).*0;
#Initializations of vector form matrices 
sK_1 = rand(nne^2*nElem,1).*0;
sK_2 = rand(nne^2*nElem,1).*0;
sF = rand(nne^2*nElem,1).*0;
#Gauss quadrature loop 
for p = 1:nQuad   
    #Jacobian evaluation and its absolute value
    global J = [xx*Nx[:,p]  yy*Nx[:,p] xx*Ny[:,p]  yy*Ny[:,p]];
    global absJac = J[:,1].*J[:,4] .- J[:,2].*J[:,3];
    #Evaluation of gradients of the shape function 
    global DNDx = (J[:,4]*Nx[:,p]' .- J[:,2]*Ny[:,p]')./(absJac.*ones(nElem,nne));
    global DNDy = (-J[:,3]*Nx[:,p]' .+ J[:,1]*Ny[:,p]')./(absJac.*ones(nElem,nne));

    let index = 0
    for i = 1:nne 
        for j = 1:nne
            #Galerkin diffusion term 
            Kij_1 = -gW[p]*(DNDx[:,i].*DNDx[:,j]);
            Kij_2 = -gW[p]*(DNDy[:,i].*DNDy[:,j]);
            Kij_1 = Kij_1.*absJac;
            Kij_2 = Kij_2.*absJac;
            sK_1[index+ 1:index + nElem] .+= Kij_1;
            sK_2[index+ 1:index + nElem] .+= Kij_2;
            #Galerkin source term 
            Fij = gW[p]*(N[i,p].*N[j,p]);
            Fij = Fij.*absJac;
            sF[index + 1:index+nElem] .+= Fij;
            #updating index
            index = index + nElem;
        end
    end
    end #End for index = 0 
end
#Assembly of local Matrices to Global Matrices 
K = sprand(ndof,ndof,0.2).*0;
F = sprand(ndof,ndof,0.2).*0;
sK = sK_1 + sK_2;
for i = 1:nne^2*nElem
    K[Int(ii[i]),Int(jj[i])] = K[Int(ii[i]),Int(jj[i])] + sK[i];
    F[Int(ii[i]),Int(jj[i])] = F[Int(ii[i]),Int(jj[i])] + sF[i];
end
#Source Vector 
F = F*(rand(ndof,1).^0 .*f);
#Incorporatring the Dirichlet boundary conditions
F = F - K*Solu;
unKnowns = Int.(setdiff(nodes,BCNodes));
K = submat(K,unKnowns);
F = subvec(F,unKnowns);
#Solu[unKnowns] = K\F;
Solu[unKnowns] ,nrm, iter = cg(K,F,I);
#Postprocessing 
s = surface(X,Y,reshape(Solu,(nElemX+1, nElemY+1)), title="2D Poisson Problem", xlabel="x axis", ylabel="y axis", zlabel = "Range")
display(s)
