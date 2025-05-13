using Plots
include("ThinFSIFunc.jl")

# Dirichlet Neumann Testing
# Syntax for the Test Function is 
# MaxIter, Iter, Res = ThinFSI(α = Relaxation Parameter ,αf = Robin Parameter ,αs = Robin Parameter ,ρs = Solid Density)
# This test is carried out only for the first time step 

# Dirichlet Neumann Tests 
# ρs = 45, α = varied from 0 10 1 with 0.1 increments 

function DN1(ρs)
    relaxation = collect(0:0.01:1)
    let arr = []
    for i = 1:101
        println("Test Step is ",i)
        MaxIter, Iter, Res = ThinFSI(relaxation[i] , 1e10, 0 , ρs)
        arr = [arr ; MaxIter]
    end
    plot!(relaxation,arr,label = "ρs = "*string(ρs))
    end
end

p = plot(xlabel= "Relaxation Parameter", ylabel = "Iterations", title = "Implicit Dirichlet Neumann")

DN1(30)
DN1(35)
DN1(40)
DN1(45)
DN1(50)
plot!(ylim = (0,50))
savefig("DirichletNeumann.png")
