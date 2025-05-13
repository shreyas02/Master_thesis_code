function Ain = inverse(A)

n = size(A,1);
[L,U] = LUfactor(A);
Ain = zeros(n,n);
I = eye(n);
for i = 1:n
    b = I(:,i);
    c = forward(L,b);
    Ain(:,i) = backward(U,c);
end

end