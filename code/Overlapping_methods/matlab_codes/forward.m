%forward solve

function y = forward(L,b)
n = size(L,1);
y = b;
for k = 2:n
    for c = 1:k-1
        y(k) = y(k)-L(k,c)*y(c);
    end
end
end