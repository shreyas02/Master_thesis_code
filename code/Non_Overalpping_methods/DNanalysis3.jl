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

function dirichletNeumann3(a1,a2)
    interface = 2.5
    relaxation = collect(0:0.01:1)
    let arr = []
    for i = 1:101
        temp = robin(relaxation[i],interface,a1,a2,Inf,0)
        arr = [arr ; temp]
    end
    plot!(relaxation,arr,ylims=(0,50),label = "a₁ = "*string(a1)*" a₂ = "*string(a2))
    end
end

p = plot(xlabel= "Relaxation Parameter", ylabel = "Iterations")
 
dirichletNeumann3(1,1)
dirichletNeumann3(10,1)
dirichletNeumann3(1,10)
dirichletNeumann3(100,1)
dirichletNeumann3(1,100)

savefig("DN3") 

