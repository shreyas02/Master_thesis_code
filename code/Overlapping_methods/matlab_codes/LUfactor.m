function [L,U] = LUfactor(A)

%LU factorisation 
n = size(A,1);
L = eye(n);
U = A;
for k = 1:n-1
    for l = k+1:n
        m = U(l,k)/U(k,k);
        U(l,k)=0;
        for c = k+1:n
            U(l,c) = U(l,c) - m*U(k,c);
        end
        L(l,k) = m;
    end
end

end


