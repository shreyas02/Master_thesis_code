%backward solve

function y = backward(U,b)
n = size(U,1);
y = b;
y(n) = y(n)/U(n,n);
for k = n-1:-1:1
    for c = n:-1:k+1
        y(k) = y(k)-U(k,c)*y(c);
    end
    y(k) = y(k)/U(k,k);
end
end
