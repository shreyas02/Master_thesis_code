using LinearAlgebra
using FillArrays, BlockArrays

using Revise
using Gridap
using Gridap.ReferenceFEs, Gridap.Algebra, Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.MultiField, Gridap.Algebra

using GridapSolvers
using GridapSolvers.LinearSolvers, GridapSolvers.MultilevelTools
using GridapSolvers.BlockSolvers: LinearSystemBlock, BiformBlock, BlockTriangularSolver

using WriteVTK

# Definitions 
# MatrixModifier setup 

struct MatrixModifier <: LinearSolver
    params::Dict
    ls::LinearSolver
end

struct MatrixModifierSymbolicSetup <: SymbolicSetup
    ls::MatrixModifier
    ss::SymbolicSetup
end

mutable struct MatrixModifierNumericalSetup{T<:AbstractMatrix} <: NumericalSetup
    ss::MatrixModifierSymbolicSetup
    A::T
    ns::NumericalSetup
    BlockStructure::Tuple
    Q::AbstractMatrix
end

function Gridap.Algebra.symbolic_setup(ls::MatrixModifier, mat::AbstractMatrix) 
    MatrixModifierSymbolicSetup(ls,Gridap.Algebra.symbolic_setup(ls.ls,mat))
end

function Gridap.Algebra.numerical_setup(ss::MatrixModifierSymbolicSetup,A::AbstractMatrix,ws::Tuple{Vararg{Real}}) 
    UI = ss.ls.params["UI"]
    UΓ = ss.ls.params["UΓ"]
    PI = ss.ls.params["PI"]
    SI = ss.ls.params["SI"]
    SΓ = ss.ls.params["SΓ"]
    αf = ss.ls.params["αf"]
    αs = ss.ls.params["αs"]
    # Modifying the Jacobian matrix
    A[SΓ,UI] = A[UΓ,UI]
    A[SΓ,PI] = A[UΓ,PI]
    A[SΓ,UΓ] = A[UΓ,UΓ]
    A[UΓ,PI] = 0*A[UΓ,PI]
    A[UΓ,UI] = 0*A[UΓ,UI]
    A[UΓ,UΓ] = ws[1]*I(size(UΓ)[1])
    A[UΓ,SΓ] = -ws[2]*I(size(UΓ)[1])
    # Creating the Block Structure in the Jacobain Matrix 
    tempf = size([UΓ;UI;PI])[1]
    temps = size([SΓ;SI])[1] + tempf
    rowSize = [tempf , temps - tempf]  
    columnSize = [tempf , temps - tempf] 
    A = BlockArray(A, rowSize, columnSize)
    # Creation of the Robin-Robin mapping matrix 
    Q = BlockArray(zeros(size(A)), blocksizes(A,1), blocksizes(A,2))
    Q[UI,UI] = 1.0*I(size(UI)[1])
    Q[PI,PI] = 1.0*I(size(PI)[1])
    Q[UΓ,UΓ] = (αf)*I(size(UΓ)[1])
    Q[SΓ,SΓ] = 1.0*I(size(SΓ)[1]) 
    Q[UΓ,SΓ] = 1.0*I(size(UΓ)[1])
    Q[SΓ,UΓ] = -(αs)*I(size(SΓ)[1])
    if αf > 1e5
        Q[UΓ,UΓ] = (1.0)*I(size(UΓ)[1])
        Q[UΓ,SΓ] = (1/αf)*I(size(SΓ)[1])
    end
    if αs > 1e5
        Q[SΓ,UΓ] = -(1)*I(size(UΓ)[1])
        Q[SΓ,SΓ] = (1/αs)*I(size(SΓ)[1])
    end
    BlockStructure = (rowSize,columnSize)
    # Modifying the final matrix to include robin conditions 
    A = Q*A
    MatrixModifierNumericalSetup(ss,A,Gridap.Algebra.numerical_setup(ss.ss,A),BlockStructure,Q)
end

function Gridap.Algebra.numerical_setup!(ns::MatrixModifierNumericalSetup,A::AbstractMatrix,ws::Tuple{Vararg{Real}})
    UI = ns.ss.ls.params["UI"]
    UΓ = ns.ss.ls.params["UΓ"]
    PI = ns.ss.ls.params["PI"]
    SI = ns.ss.ls.params["SI"]
    SΓ = ns.ss.ls.params["SΓ"]
    # Modifying the Jacobian matrix
    A[SΓ,UI] = A[UΓ,UI]
    A[SΓ,PI] = A[UΓ,PI]
    A[SΓ,UΓ] = A[UΓ,UΓ]
    A[UΓ,PI] = 0.0*A[UΓ,PI]
    A[UΓ,UI] = 0.0*A[UΓ,UI]
    A[UΓ,UΓ] = ws[1]*I(size(UΓ)[1])
    A[UΓ,SΓ] = -ws[2]*I(size(UΓ)[1])
    # Creating the Blcok Structure
    A = BlockArray(A, ns.BlockStructure[1], ns.BlockStructure[2])
    # Modifying the final matrix to include robin conditions 
    A = ns.Q*A
    ns.A = A
    Gridap.Algebra.numerical_setup!(ns.ns,A)
end

function Gridap.Algebra.solve!(x::AbstractVector,ns::MatrixModifierNumericalSetup,b::AbstractVector)
    solve!(x,ns.ns,b)
end

function Gridap.Algebra.solve!(
x::AbstractVector,
ls::MatrixModifier, lop::Gridap.ODEs.LinearStageOperator,
ns::Nothing
)

J = lop.J
ss = Gridap.Algebra.symbolic_setup(ls, J)
ns = Gridap.Algebra.numerical_setup(ss, J, lop.ws)

r = lop.r
rmul!(r, -1)
# Modifying the residual # JX = r
r[ls.params["SΓ"]] = r[ls.params["SΓ"]] + r[ls.params["UΓ"]]
r[ls.params["UΓ"]] = (lop.usx[2][ls.params["SΓ"]] - lop.usx[1][ls.params["UΓ"]])
# Creating the Block Structre for the final vector 
r = BlockVector(r, ns.BlockStructure[2])
# Modifying the final vector to include robin conditions
r = ns.Q*r
# Solving
solve!(x, ns, r)
ns
end

function Gridap.Algebra.solve!(
x::AbstractVector,
ls::MatrixModifier, lop::Gridap.ODEs.LinearStageOperator,
ns
)
if !lop.reuse
    J = lop.J
    Gridap.Algebra.numerical_setup!(ns, J, lop.ws) 
end
r = lop.r
rmul!(r, -1)
# Modifying the residual # JX = r
r[ls.params["SΓ"]] = r[ls.params["SΓ"]] + r[ls.params["UΓ"]]
r[ls.params["UΓ"]] = (lop.usx[2][ls.params["SΓ"]] - lop.usx[1][ls.params["UΓ"]])
# Creating the Block Structre for the final vector 
r = BlockVector(r, ns.BlockStructure[2])
# Modifying the final vector to include robin conditions
r = ns.Q*r
# Solving
solve!(x, ns, r)
ns
end

# Test Function

function ThinFSI(α,αf,αs,ρs)
    
    # FSI Problem Definition starts here 

    model = DiscreteModelFromFile("mesh/square.json")

    # Define reference FE (Q2/P1(disc) pair)
    order = 2
    reffeᵤ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
    reffeₚ = ReferenceFE(lagrangian,Float64,order-1)
    reffeₛ = ReferenceFE(lagrangian, Float64, order )

    # Define triangulation and integration measure
    degree = 2*order + 1
    Ωf = Interior(model)
    dΩf = Measure(Ωf,degree)
    Σ = BoundaryTriangulation(model,tags="Top")
    dΣ = Measure(Σ,degree)

    # Define Test FE Functions 
    V = TestFESpace(Ωf,reffeᵤ,dirichlet_tags = ["Bottom","BCTop"], 
                        dirichlet_masks =[(false,true),(false,true)] ,conformity=:H1)
    Q = TestFESpace(Ωf,reffeₚ,dirichlet_tags = ["Left","Right"], conformity=:H1)
    S = TestFESpace(Σ,reffeₛ,dirichlet_tags = ["BCTop"],conformity=:H1)
    Y = MultiFieldFESpace([V,Q,S])

    # Define trial FESpaces from Dirichlet values
    gu(t) = x -> VectorValue(0.0,0.0) # Velocity BCs
    gd(t) = x -> 0.0 # Displacement BCs
    gpr(t) = x -> 0.0 # Right Pressure BCs
    function gpl(t) # Left Pressure BCs
        function inner_function(x)
            if t > 0.005
                return 0.0
            else
                return 20000
            end
        end
        return inner_function
    end
    U = TransientTrialFESpace(V,[gu,gu])
    P = TransientTrialFESpace(Q,[gpl,gpr])
    D = TransientTrialFESpace(S, [gd])
    X = TransientMultiFieldFESpace([U,P,D])

    # Initial conditions function
    h_u(t) = x -> VectorValue(0.0,0.0) # Velocity
    h_d(t) = x -> 0.0 # Displacement
    function h_p(t) # Left Pressure BCs
        function inner_function(x)
            if t > 0.005
                return 0.0
            else
                if x[1] == 0
                    return 20000
                else
                    return 0
                end
            end
        end
        return inner_function
    end

    # Initial Velocity and Acelleration terms 
    v_p(t) = x -> 0.0 # Pressure
    v_u(t) = x -> VectorValue(0.0,0.0) # Velocity
    v_d(t) = x -> 0.0 # Displacement

    # Define bilinear and linear form
    ρf = 1;
    # ρs already taken from user
    ϵ = 0.1;
    L = 1e5; # L = ϵ*E/(R²(1-ν²))
    jac(t,(du,dp,dd),(v,q,s)) = ∫((∇⋅v)*dp + (∇⋅du)*q)dΩf + ∫(L*dd*s)dΣ
    jac_t(t,(dtu,dtp,dtd),(v,q,s)) = ∫(ρf*v⋅dtu)dΩf
    jac_tt(t,(dttu,dttp,dttd),(v,q,s)) = ∫(ρs*ϵ*s*dttd)dΣ 
    l(t,(v,q,s)) = ∫(v⋅VectorValue(0.0, 0.0))dΩf + ∫(s*0.0)dΣ

    # Build affine FE operator
    op = TransientLinearFEOperator((jac, jac_t, jac_tt), l, X, Y)

    # Time numerics 
    t0, tF = 0.0, 0.001 # 0 to 3 seconds
    dt = 0.001

    ##### Getting the common degrees of freedom
    # Function to get the range for each field
    num_fields = 3
    sizes = [Gridap.FESpaces.num_free_dofs(X[i]) for i in 1:num_fields]
    function get_range(field)
        start = sum(sizes[1:field-1]) + 1
        stop = sum(sizes[1:field])
        return start:stop
    end

    IdCord2D = t -> x -> VectorValue(0.0,1e10*x[1]+x[2])
    IdCord1D = t -> x -> Float64(1e10*x[1]+x[2])

    cord = interpolate_everywhere([IdCord2D(t0) , IdCord1D(t0), IdCord1D(t0)], X(t0))
    cord = get_free_dof_values(cord)
    cord1 = interpolate_everywhere(IdCord2D(t0), U(t0))
    cord1 = get_free_dof_values(cord1)
    cord2 = interpolate_everywhere(IdCord1D(t0), D(t0))
    cord2 = get_free_dof_values(cord2)

    common_cord = intersect(cord1,cord2)
    UΓ = vcat([findall(x -> x == elem, cord1) for elem in common_cord]...)
    SΓ = vcat([findall(x -> x == elem, cord2) for elem in common_cord]...)
    SolidRange = get_range(3)
    SΓ = SolidRange[SΓ] #SΓ DOFs
    UI = collect(get_range(1))
    PI = collect(get_range(2))
    SI = collect(get_range(3))
    SI = [x for x in SI if !(x in SΓ)]
    UI = [x for x in UI if !(x in UΓ)]

    # Loading all the parameters 
    params = Dict("UI"=>UI,"PI"=>PI,"SI"=>SI,"UΓ"=>UΓ,"SΓ"=>SΓ,"αf"=>αf,"αs"=>αs)

    # Creation of the relaxation parameter
    # α relaxation parameter already taken from user  
    ω = ones(num_free_dofs(X))
    ω[SΓ] = α*ω[SΓ]

    # System Solver Definition 
    solver_blocks = LUSolver()
    bblocks = [LinearSystemBlock() LinearSystemBlock();
                LinearSystemBlock() LinearSystemBlock()]
    coeffs = [1.0 0.0;
                1.0 1.0]
    P = BlockTriangularSolver(bblocks,[solver_blocks,solver_blocks],coeffs,:lower)
    #SystemSolver = GMRESSolver(20;Pl=P,restart=false,m_add=1,maxiter=100,atol=1e-12,rtol=1.e-6,verbose=true)
    SystemSolver = RichardsonLinearSolver(ω,51;Pl=P,rtol=1e-8,atol=1e-6,verbose=true)

    SolverStokes = MatrixModifier(params,SystemSolver)

    # ODE Solver 
    solver_ode = Gridap.ODEs.GeneralizedAlpha2(SolverStokes,dt,0.5)

    # Imposing Intitial conditions
    x_initial = interpolate_everywhere([h_u(t0) , h_p(t0), h_d(t0)], X(t0))
    v_initial = interpolate_everywhere([v_u(t0) , v_p(t0), v_d(t0)], X(t0))
    a_initial = interpolate_everywhere([v_u(t0) , v_p(t0), v_d(t0)], X(t0))

    xₜ = solve(solver_ode, op, t0, tF, (x_initial,v_initial,a_initial))

    if !isdir("tmp")
        mkdir("tmp")
    end

    pvd_Ωf = paraview_collection("results/TestThinFSIFluids")
    pvd_Σ = paraview_collection("results/TestThinFSISolids")

    pvd_Ωf[0] = createvtk(Ωf, "tmp/resultsFluid_0" * ".vtu", cellfields=["u" => x_initial[1] , "p" => x_initial[2]])
    pvd_Σ[0] = createvtk(Σ, "tmp/resultsSolid_0" * ".vtu", cellfields=["d" => x_initial[3]])

    for (tn, (uhn,phn,dhn)) in xₜ
        pvd_Ωf[tn] = createvtk(Ωf, "tmp/resultsFluid_$tn" * ".vtu", cellfields=["u" => uhn, "p" => phn])
        pvd_Σ[tn] = createvtk(Σ, "tmp/resultsSolid_$tn" * ".vtu", cellfields=["d" => dhn])
        println(tn)
    end
    vtk_save(pvd_Ωf)
    vtk_save(pvd_Σ)
    
    MaxIter = SystemSolver.log.num_iters
    Iter = collect(1:MaxIter+1)
    Res = SystemSolver.log.residuals
    Res = Res[Iter]
   
    # Returning data
    return MaxIter, Iter, Res
end
