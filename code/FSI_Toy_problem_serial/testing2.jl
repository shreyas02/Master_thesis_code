using Plots
include("ThinFSIFunc.jl")

# Robin Neumann Testing 1
# Syntax for the Test Function is 
# MaxIter, Iter, Res = ThinFSI(α = Relaxation Parameter ,αf = Robin Parameter ,αs = Robin Parameter ,ρs = Solid Density)
# This test is carried out only for the first time step 

# Dirichlet Neumann Tests 
# ρs = 45, α = varied from 0 10 1 with 0.1 increments 

function RN1(Dir,ρs)
    relaxation = collect(0:0.01:1)
    let arr = []
    for i = 1:101
        println("Test Step is ",i)
        MaxIter, Iter, Res = ThinFSI(relaxation[i] , Dir, 0 , ρs)
        arr = [arr ; MaxIter]
    end
    plot!(relaxation,arr,label = "ρs = "*string(ρs)*" αf = "*string(Dir))
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

RN1(262.5,20)
DN1(20)

RN1(565.5,45)
DN1(45)
plot!(ylim = (0,50))
plot!(title = "")
savefig("RobinNeumann.png")
