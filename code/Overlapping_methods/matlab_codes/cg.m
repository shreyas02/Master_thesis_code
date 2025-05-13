function [xk,nrm,iter] = cg(a,b,precond)

%%conjugate gradient method 
n = size(a,1);
%initial guess
xk = ones(n,1);
rk = b - a*xk;
zk = precond*rk;
pk = zk;
k = 0;
toll = 1e-6;

%iteration starts 
for i = 1 : 1000
    alphaNum = rk'*zk;
    alphaDen = pk'*a*pk;
    alpha = alphaNum/alphaDen;
    xk1 = xk + alpha.*pk;
    rk1 = rk - alpha.*a*pk;
    zk1 = precond*rk1;
    beta = zk1'*rk1;
    betaDen = zk'*rk;
    beta = beta/betaDen;
    pk1 = zk1 + beta.*pk;
    pk = pk1;
    rk = rk1;
    xk = xk1;
    zk = zk1;
    if(norm(rk1) < toll)
        break
    end
end
nrm = norm(rk1);
fprintf('Iterations - : %e, \n', i);
fprintf('Condition Number - : %e, \n', cond(precond*a));
iter = i;
end

