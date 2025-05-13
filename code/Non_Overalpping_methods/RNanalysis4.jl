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

function robinNeumann4(a1,a2,interface,dir)
    relaxation = collect(0:0.01:1)
    let arr = []
    for i = 1:101
        temp = robin(relaxation[i],interface,a1,a2,dir,0)
        arr = [arr ; temp]
    end
    plot!(relaxation,arr,ylims=(0,200),label = "αf = "*string(dir))
    end
end

p = plot(xlabel= "Relaxation Parameter", ylabel = "Iterations",title = "a₁ = 10, a₂ = 1")

robinNeumann4(10,1,2.5,Inf)
robinNeumann4(10,1,2.5,10)
robinNeumann4(10,1,2.5,50)
robinNeumann4(10,1,2.5,100)

savefig("RN4") 
