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

function dirichletRobin7(a1,a2,interface,neu)
    relaxation = collect(0:0.01:1)
    let arr = []
    for i = 1:101
        temp = robin(relaxation[i],interface,a1,a2,Inf,neu)
        arr = [arr ; temp]
    end
    plot!(relaxation,arr,ylims=(0,500),label = "αs = "*string(neu))
    end
end

p = plot(xlabel= "Relaxation Parameter", ylabel = "Iterations", title = "a₁ = 10, a₂ = 1, interface = 4")

dirichletRobin7(10,1,4,0)
dirichletRobin7(10,1,4,0.25)
dirichletRobin7(10,1,4,0.5)
dirichletRobin7(10,1,4,1)

savefig("DR7")
