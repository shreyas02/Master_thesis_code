using Plots
include("ThinFSIFunc.jl")

# Robin Robin Testing 2
# Syntax for the Test Function is 
# MaxIter, Iter, Res = ThinFSI(α = Relaxation Parameter ,αf = Robin Parameter ,αs = Robin Parameter ,ρs = Solid Density)
# This test is carried out only for the first time step 

# Dirichlet Neumann Tests 
# ρs = 45, α = varied from 0 10 1 with 0.1 increments 

function RN2(Dir,ρs)
    relaxation = collect(0:0.01:1)
    let arr = []
    for i = 1:101
        println("Test Step is ",i)
        MaxIter, Iter, Res = ThinFSI(relaxation[i] , Dir, Dir , ρs)
        arr = [arr ; MaxIter]
    end
    plot!(relaxation,arr,label = "ρs = "*string(ρs)*" αf,αs= "*string(Dir))
    end
end

function DN1(ρs)
    relaxation = collect(0:0.01:1)
    let arr = []
    for i = 1:101
        println("Test Step is ",i)
        MaxIter, Iter, Res = ThinFSI(relaxation[i] , 1e10, 0 , ρs)
        arr = [arr ; MaxIter]
    end
    plot!(relaxation,arr,label = "ρs = "*string(ρs)*" Dirichlet-Neumann")
    end
end

p = plot(xlabel= "Relaxation Parameter", ylabel = "Iterations", title = "Implicit Robin Neumann")

RN2(2100,20)
RN2(1050,20)
RN2(525,20)
RN2(262.5,20)


savefig("RobinNeumann3.png")
