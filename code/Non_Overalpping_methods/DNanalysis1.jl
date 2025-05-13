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

function dirichletNeumann1(interface)
    relaxation = collect(0:0.01:1)
    let arr = []
    for i = 1:101
        temp = robin(relaxation[i],interface,1,1,Inf,0)
        arr = [arr ; temp]
    end
    plot!(relaxation,arr,ylims=(0,50),label = "interface = "*string(interface))
    end
end

plot(xlabel= "Relaxation Parameter", ylabel = "Iterations")

dirichletNeumann1(1)
dirichletNeumann1(2)
dirichletNeumann1(2.5)
dirichletNeumann1(3)
dirichletNeumann1(4)

savefig("DN1") 
