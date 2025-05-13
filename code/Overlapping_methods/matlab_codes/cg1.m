clc
clear 

a = [4,1;1,3];
b = [1;2];

x = a\b;
a = sparse(a);


%%conjugate gradient method 
n = size(a,1);
%initial guess
xk = ones(n,1);
rk = b - a*xk;
pk = rk;
k = 0;
toll = 1e-3;

%iteration starts 
for i = 1 : 1000
    alphaNum = rk'*rk;
    alphaDen = pk'*a*pk;
    alpha = alphaNum/alphaDen;
    xk1 = xk + alpha.*pk;
    rk1 = rk - alpha.*a*pk;
    beta = rk1'*rk1;
    beta = beta/alphaNum;
    pk1 = rk1 + beta.*pk;
    pk = pk1;
    rk = rk1;
    xk = xk1;
    if(norm(rk1) < toll)
        break
    end
end

[xn,nrm] = cg(a,b,inv(a));
[L,U] = LUfactor(a);

xn1 = forward(L,b);
xn2 = backward(U,b);
Ain = inverse(a);


