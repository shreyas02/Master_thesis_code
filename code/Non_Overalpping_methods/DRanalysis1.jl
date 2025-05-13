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

function dirichletRobin1(a1,a2,interface,neu)
    relaxation = collect(0:0.01:1)
    let arr = []
    for i = 1:101
        temp = robin(relaxation[i],interface,a1,a2,Inf,neu)
        arr = [arr ; temp]
    end
    plot!(relaxation,arr,ylims=(0,50),label = "αs = "*string(neu))
    end
end

p = plot(xlabel= "Relaxation Parameter", ylabel = "Iterations")

dirichletRobin1(1,1,2.5,0)
dirichletRobin1(1,1,2.5,0.25)
dirichletRobin1(1,1,2.5,0.5)
dirichletRobin1(1,1,2.5,1)
dirichletRobin1(1,1,2.5,1.5)
dirichletRobin1(1,1,2.5,2)

savefig("DR1")
