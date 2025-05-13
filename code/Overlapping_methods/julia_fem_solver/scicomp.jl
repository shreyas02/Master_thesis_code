function submat(mat,unKnowns)
    NumNodes = size(unKnowns,1);
    matF = sprand(NumNodes,NumNodes,0.2).*0;
    for i = 1:NumNodes
        for j = 1:NumNodes
            x = Int(unKnowns[i]);
            y = Int(unKnowns[j]);
            matF[i,j] = mat[x,y];
        end
    end
    return matF
end

function subvec(vec,unKnowns)
    NumNodes = size(unKnowns,1);
    vecF = zeros(NumNodes,1);
    for i = 1:NumNodes
            x = Int(unKnowns[i]);
            vecF[i] = vec[x]
    end
    return vecF
end

function linArr(start,N,step) #N is the max value 
    arr = [start];
    for i = 2:N
        if(arr[i-1]+step<=N)
            arr = vcat(arr,arr[i-1]+step)
        else
            break
        end
    end
    return arr
end

function meshgrid(xin,yin)
    nx=length(xin)
    ny=length(yin)
    xout=zeros(ny,nx)
    yout=zeros(ny,nx)
    for jx=1:nx
        for ix=1:ny
            xout[ix,jx]=xin[jx]
            yout[ix,jx]=yin[ix]
        end
    end
    return (x=xout, y=yout)
end

function cg(a,b,precond)
    #conjugate gradient method 
    n = size(a,1);
    #Initial guess
    xk = ones(n,1);
    rk = b - a*xk;
    zk = precond*rk;
    pk = zk;
    k = 0;
    toll = 1e-6;
    #iteration process starts 
    for i = 1:1000
        alphaNum = rk'*zk;
        alphaDen = pk'*a*pk;
        alpha = alphaNum/alphaDen;
        global xk1 = xk .+ alpha.*pk;
        global rk1 = rk .- alpha.*a*pk;
        global zk1 = precond*rk1;
        beta = zk1'*rk1;
        betaDen = zk'*rk;
        beta = beta/betaDen;
        pk1 = zk1 + beta.*pk;
        pk = pk1;
        rk = rk1;
        xk = xk1;
        zk = zk1;
        if norm(rk1) < toll
            break
        end
    end
    nrm = norm(rk1);
    iter = i;
    return (xk , nrm, iter)
end

function backward(U,b)
    n = size(U,1);
    y = b;
    y[n] = y[n]/U[n,n];
    for k = n-1:-1:1
        for c = n:-1:k+1
            y[k] = y[k] - U[k,c]*y[c];
        end
        y[k] = y[k]/U[k,k]; 
    end
    return y 
end

function forward(L,b)
    n = size(L,1);
    y = b;
    for k = 2:n
        for c = 1:k-1
            y[k] = y[k] - L[k,c]*y[c];
        end
    end
    return y
end

function LUfactor(A)
    n = size(A,1);
    L = I(n).*ones(n,n);
    U = A;
    for k = 1:n-1
        for l = k+1:n
            m = U[l,k]/U[k,k];
            U[l,k] = 0;
            for c = k+1:n
                U[l,c] = U[l,c] - m*U[k,c];
            end
            L[l,k] = m;
        end
    end 
    return (L,U)
end

function inverse(A)
    n = size(A,1);
    L,U = LUfactor(A);
    Ain = zeros(n,n);
    eye = I(n).*ones(n,n);
    for i = 1:n
        b = eye[:,i];
        c = forward(L,b);
        Ain[:,i] = backward(U,c);
    end
    return Ain
end