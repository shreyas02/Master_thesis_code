using LinearAlgebra
using FillArrays, BlockArrays

using Gridap
using Gridap.ReferenceFEs, Gridap.Algebra, Gridap.Geometry, Gridap.FESpaces
using Gridap.CellData, Gridap.MultiField, Gridap.Algebra
using PartitionedArrays

using GridapSolvers
using GridapSolvers.LinearSolvers, GridapSolvers.MultilevelTools
using GridapSolvers.BlockSolvers: LinearSystemBlock, BiformBlock, BlockTriangularSolver
using MAT, Plots

include("modrr1")

function dirichletNeumann4(a1,a2,interface)
    relaxation = collect(0:0.01:1)
    let arr = []
    for i = 1:101
        temp = robin(relaxation[i],interface,a1,a2,Inf,0)
        arr = [arr ; temp]
    end
    plot!(relaxation,arr,ylims=(0,50),label = "interface = "*string(interface))
    end
end

p = plot(xlabel= "Relaxation Parameter", ylabel = "Iterations", title = "a₁ = 1 and a₂ = 10")

dirichletNeumann4(1,10,1)
dirichletNeumann4(1,10,2)
dirichletNeumann4(1,10,2.5)
dirichletNeumann4(1,10,3)
dirichletNeumann4(1,10,4)

savefig("DN4")  
