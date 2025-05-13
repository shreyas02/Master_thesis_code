include("scicomp.jl")
using Plots
using SparseArrays
using LinearAlgebra
#Length of the Domain
L = 1; 
#Number of the elements 
nElem = 100;
#Number of nodes per element 
nne = 2; 
#Source Term 
f = 1; 
#Boundary Nodes 
BCNodes = [1;nElem+1];
#Boundary Values
BCValues = [0;0];
#Number of nodes 
ndof = Int(nElem + 1);
#All Nodes for this consructed mesh
nodes = rand(ndof,1).*0;
let test = 1;
for i = 1:ndof
    nodes[i] = test;
    test = test + 1;
end
end
#Coordinates array 
crd = rand(ndof,1).*0;
dx = L/nElem;
let temp = 0;
for i in 1:ndof
    crd[i] = temp;
    temp = temp + dx;
end
end
#Connectivity matrix construction 
cnn = rand(nElem,nne).*0;
let nn = 1;
for i in 1:nElem
    nn = nn - 1;
    for j = 1:nne
        nn = nn + 1;
        cnn[i,j] = nn;
    end
end
end
#Definition of Gauss Quadrature and shape function space 
#Location of Gauss Points 
gP = [-1/sqrt(3),1/sqrt(3)];
#Weights of Gauss Points 
gW = [1,1];
#Number of Gauss Points 
nQuad = length(gW);
#Shape Functions 
N = rand(2,2).*0; #initialization 
N[1,:] = 0.5.*(1 .- gP);
N[2,:] = 0.5.*(1 .+ gP);
#Gradient of shape Functions
Nx = rand(2,2).*0; #initialization
Nx[1,:] = [-0.5,-0.5];
Nx[2,:] = [0.5,0.5];
#Formation of local to global DOF Mapping 
ii = rand(nne^2*nElem,1).*0;
jj = rand(nne^2*nElem,1).*0;
let index = 0;
for i = 1:nne
    for j = 1:nne
        ii[index+1:index+nElem] = cnn[:,i];
        jj[index+1:index+nElem] = cnn[:,j];
        index = index + nElem;
    end
end
end
#Initialization of the Solution Vector 
Solu = rand(ndof,1).*0;
Solu[BCNodes] = BCValues;
#initialization of element level calculation 
i = size(cnn,1);
j = size(cnn,2);
xx = rand(i,j).*0;
for i = 1:nne
    for j = 1:size(cnn,1)
        xx[j,i] = crd[Int(cnn[j,i])];
    end
end
#Element level Evaluation of Matrices 
#Initializations of vector form matrices 
sK = rand(nne^2*nElem,nQuad).*0;
sF = rand(nne^2*nElem,nQuad).*0;
#Gauss quadrature loop 
for p = 1:nQuad
    #Jacobian evaluation and its absolute value
    global J = xx*Nx[:,p];
    global absJac = abs.(J);
    #Evaluation of gradients of shape functions 
    global DNDx = (1 ./absJac)*Nx[:,p]';
    let index = 0
    for i = 1:nne
        for j = 1:nne
            #Galerkin Diffusion Term 
            Kij = gW[p]*(DNDx[:,i].*DNDx[:,j]);
            Kij = Kij.*absJac;
            sK[index + 1:index+nElem,p] = Kij;
            #Galerkin Source Term 
            Fij = gW[p]*(N[i,p]*N[j,p]);
            Fij = Fij.*absJac;
            sF[index + 1:index+nElem,p] = Fij;
            #updating index
            index = index + nElem;
        end
    end
    end #end for let index variable declaration 
end
#Summation of quadrature data for numerical integration 
#Initializations of vector form matrices global 
global gK = rand(nne^2*nElem,1).*0;
global gF = rand(nne^2*nElem,1).*0;
for i = 1:nne^2*nElem
    for j = 1:nQuad
        gK[i] = gK[i] + sK[i,j];
        gF[i] = sF[i] + sF[i,j];
    end
end
#Assembly of local Matrices to Global Matrices 
K = sprand(ndof,ndof,0.2).*0;
F = sprand(ndof,ndof,0.2).*0;
for i = 1:nne^2*nElem
    K[Int(ii[i]),Int(jj[i])] = K[Int(ii[i]),Int(jj[i])] + gK[i];
    F[Int(ii[i]),Int(jj[i])] = F[Int(ii[i]),Int(jj[i])] + gF[i];
end
#Source Vector 
F = F*(rand(ndof,1).^0 .*f);
#Incorporatring the Dirichlet boundary conditions
F = F - K*Solu;
#Solving for unknowns 
unKnowns = Int.(setdiff(nodes,BCNodes));
K = submat(K,unKnowns);
F = subvec(F,unKnowns);
Solu[unKnowns] ,nrm, iter = cg(K,F,I);
plot(crd, Solu, title="1D Poisson Problem", xlabel="Domain", ylabel="Range")


