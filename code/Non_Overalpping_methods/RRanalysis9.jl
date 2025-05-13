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

function RobinRobin9(a1,a2,interface,rob)
    relaxation = collect(0:0.01:1)
    let arr = []
    for i = 1:101
        temp = robin(relaxation[i],interface,a1,a2,rob,rob)
        arr = [arr ; temp]
    end
    plot!(relaxation,arr,ylims=(0,50),label = "αf,αs = "*string(rob))
    end
end

p = plot(xlabel= "Relaxation Parameter", ylabel = "Iterations", title = "a₁ = 1, a₂ = 10, interface = 4")

RobinRobin9(1,10,4,0.25)
RobinRobin9(1,10,4,0.5)
RobinRobin9(1,10,4,1)
RobinRobin9(1,10,4,1.5)
RobinRobin9(1,10,4,2)

savefig("RR9")
