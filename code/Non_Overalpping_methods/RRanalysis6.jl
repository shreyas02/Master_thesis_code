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

function RobinRobin6(a1,a2,interface,rob)
    relaxation = collect(0:0.01:1)
    let arr = []
    for i = 1:101
        temp = robin(relaxation[i],interface,a1,a2,rob,rob)
        arr = [arr ; temp]
    end
    plot!(relaxation,arr,ylims=(0,100),label = "αf,αs = "*string(rob))
    end
end

p = plot(xlabel= "Relaxation Parameter", ylabel = "Iterations", title = "a₁ = 10, a₂ = 1, interface = 1")

RobinRobin6(10,1,1,0.25)
RobinRobin6(10,1,1,0.5)
RobinRobin6(10,1,1,1)
RobinRobin6(10,1,1,1.5)
RobinRobin6(10,1,1,2)

savefig("RR6")
