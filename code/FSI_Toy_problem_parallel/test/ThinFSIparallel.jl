using Gridap, GridapGmsh
using Gridap.ReferenceFEs, Gridap.Algebra, Gridap.Geometry, Gridap.FESpaces, Gridap.MultiField
using Gridap.CellData, Gridap.MultiField, Gridap.Algebra

using GridapSolvers
using GridapSolvers.LinearSolvers, GridapSolvers.MultilevelTools
using GridapSolvers.BlockSolvers: LinearSystemBlock, BiformBlock, BlockTriangularSolver

using LinearAlgebra, SparseArrays, FillArrays, BlockArrays, WriteVTK

using GridapDistributed, PartitionedArrays, MPI
using Revise, Logging, Plots

include("../shared_lib/mpicpp.jl")

#######################
## Begin Trilinos Solve 
struct TrilinosSolve <: Gridap.Algebra.LinearSolver
end

struct TrilinosSolveSymbolicSetup <: Gridap.Algebra.SymbolicSetup
    solver
end

mutable struct TrilinosSolveNumericalSetup{T<:AbstractMatrix} <: NumericalSetup 
    solver
    A::T
end

function Gridap.Algebra.symbolic_setup(solver::TrilinosSolve, mat::AbstractMatrix) 
    return TrilinosSolveSymbolicSetup(solver)
end

function Gridap.Algebra.numerical_setup(ss::TrilinosSolveSymbolicSetup,A::AbstractMatrix) 
    return TrilinosSolveNumericalSetup(ss.solver,A)
end

function Gridap.Algebra.numerical_setup(ss::TrilinosSolveSymbolicSetup,A::AbstractMatrix, x::AbstractVector) 
    return TrilinosSolveNumericalSetup(ss.solver,A)
end

function Gridap.Algebra.numerical_setup!(ns::TrilinosSolveNumericalSetup,A::AbstractMatrix)
    return ns
end

function Gridap.Algebra.solve!(x::AbstractVector,ns::TrilinosSolveNumericalSetup,b::AbstractVector)
    map(ns.A.matrix_partition,ns.A.row_partition,ns.A.col_partition, b.vector_partition, 
        x.vector_partition, x.index_partition) do LocMat, RowMap, ColMap, LocRhs, LocSoln, xmap
        OwnRowMap = RowMap[own_to_local(RowMap)];
        mpicpp.TrilinosParallel(LocMat.nzval, LocMat.rowval .- 1, 
            LocMat.colptr .- 1, LocRhs, LocSoln, OwnRowMap .- 1, ColMap .- 1,
            size(ns.A)[1], size(OwnRowMap)[1], size(LocMat)[2],
            own_to_local(xmap) .-1, own_to_local(RowMap) .-1)
    end
end
## End Trilinos Solve 
#####################

################################
## Begin convert_interfaceDOFmap
function convert_interfaceDOFmap(dofmap, dommap)
    conv_dofmap = Int[]
    for i = 1:size(dofmap)[1]
        if(dofmap[i] in dommap)
            index = findfirst(x -> x == dofmap[i], dommap)
            push!(conv_dofmap,index)
        end
    end
    return conv_dofmap
end
## End convert_interfaceDOFmap
##############################

##################################
## Begin convert_interfaceDOFmap_2
function convert_interfaceDOFmap_2(dofmap, dommap)
    conv_dofmap = Int[]
    for i = 1:size(dofmap)[1]
        if(dofmap[i] in dommap)
            push!(conv_dofmap,dofmap[i])
        end
    end
    return conv_dofmap
end
## End convert_interfaceDOFmap_2
################################

###########################
## Begin domain_range_map()
function domain_range_map(DomainMap, RangeMap)
    dom_ran_map =  Int[]
    for i = 1:size(RangeMap)[1]
        if (RangeMap[i] in DomainMap)
            index = findfirst(x -> x == RangeMap[i], DomainMap)
            push!(dom_ran_map,index)
        end
    end
    return dom_ran_map
end
## End domain_range_map()
#########################

#######################
# Begin get_FESpace_map
function get_FESpace_map(V::FESpace) # Change to Distributed FEspace
    test = get_free_dof_values(zero(V))
    FESpace_map = []
    map(test.index_partition) do locPartition
        FESpace_map = locPartition
    end
    return FESpace_map
end
# End get_FESpace_map
#####################

#####################################
## Begin get_background_cell_to_faces
function get_background_cell_to_faces(BoundaryTrian::Triangulation{}, reffe_bg::ReferenceFE{}, face_to_cell_glue::Gridap.Geometry.FaceToCellGlue{})
    boundary_dim = Gridap.Geometry.num_cell_dims(BoundaryTrian)
    face_to_bgface = face_to_cell_glue.face_to_bgface
    face_to_cell = face_to_cell_glue.face_to_cell
    background_model = Gridap.Geometry.get_background_model(BoundaryTrian)
    background_dim = Gridap.Geometry.num_cell_dims(background_model)
    @assert background_dim - 1 == boundary_dim "Inconsistant dimensions of the triangulations"
    grid_topology_bg = Gridap.Geometry.get_grid_topology(background_model)
    # Get parent cell_to_faces
    cell_to_faces_bg = Gridap.Geometry.get_cell_faces(grid_topology_bg)
    # Get the offset 
    d_to_offset = Gridap.Geometry.get_offsets(grid_topology_bg)
    # Offset face_to_bgface
    face_to_bgface = face_to_bgface .+ d_to_offset[background_dim]
    boundary_bg_cell_to_faces = Int32[] 
    node_permutations = Gridap.Geometry.get_face_nodes(reffe_bg,boundary_dim) 
    face_nodes = Gridap.Geometry.get_face_own_nodes(reffe_bg,boundary_dim)
    ptrs = Int32[1]
    for elem in face_to_cell
        bg_cell = cell_to_faces_bg[elem]
        for i = 1:size(face_nodes)[1]
            if(bg_cell[face_nodes[i][1]] in face_to_bgface)
                for j = 1:size(node_permutations[i])[1]
                    push!(boundary_bg_cell_to_faces,bg_cell[node_permutations[i][j]]);
                end
                push!(ptrs,ptrs[end] + size(node_permutations[i])[1]); continue;
            end
        end
    end
    return Gridap.Arrays.Table(boundary_bg_cell_to_faces,ptrs)
end
## End get_background_cell_to_faces
###################################

#########################
## Begin get_node_dof_ids
function get_node_dof_ids(reffe::ReferenceFE{}, gmsh_tag::Vector{String}, V::FESpace; vec_comp::Int = 1)
    # Get cell dofs 
    cell_dofs = Gridap.FESpaces.get_cell_dof_ids(V)

    # Get the triangulation 
    trian = Gridap.FESpaces.get_triangulation(V)

    # Get the Discrete Model from the triangulation 
    discrete_model = Gridap.Geometry.get_active_model(trian) 

    # Get dimensions
    D = Gridap.Geometry.num_cell_dims(trian)

    # Get the Grid Topology 
    grid_topology = Gridap.Geometry.get_grid_topology(discrete_model)

    # Get the cell to d_to_dface mapping along with the offset 
    d_to_offset = Gridap.Geometry.get_offsets(grid_topology)
    cell_to_faces = Gridap.Geometry.get_cell_faces(grid_topology) 

    # Exit the function if the size of cell_to_faces is zero or the triangulation is not there in the local subdomain
    if (size(cell_to_faces)[1] == 0) 
        return [],[] 
    end

    # Get the facelabeling of the model 
    facelabeling = Gridap.Geometry.get_face_labeling(discrete_model)

    # Get the bool mask vectors for the tag given
    d_to_dface_to_tag = [ Gridap.Geometry.get_face_tag_index(facelabeling,gmsh_tag,d)  for d in 0:D] 

    # Get the number of cells 
    numcells = size(cell_to_faces)[1]
    @assert size(cell_to_faces)[1] == size(cell_dofs)[1] "Inconsistant number of cells; quitting program"
    numfaces = size(cell_to_faces[1])[1]

    # Understand the dof arrangement in reffe 
    dof_to_node = Gridap.ReferenceFEs.get_dof_to_node(reffe)
    numnodes = maximum(dof_to_node)
    numelem = size(dof_to_node)[1]
    noderept = numelem/numnodes # To know number of components of the vector 
    @assert mod(numelem,numnodes) == 0 "Inconsistent number of components in the vector"

    # Create a cellfield of bools to know if a element has a tag
    bool_cell_field = SparseArrays.spzeros(numcells,numnodes)
    # Looping over all cells 
    for cell = 1:numcells 
        for node = 1:numfaces
            node_id = cell_to_faces[cell][node]
            for j = 1:size(d_to_offset)[1]
                if (node_id > d_to_offset[j] && j == size(d_to_offset)[1]) # Last face node check 
                    node_id = node_id - d_to_offset[j] # offset node_id
                    bool_cell_field[cell,node] = d_to_dface_to_tag[j][node_id]; continue;
                end
                if (node_id > d_to_offset[j] && node_id <= d_to_offset[j+1])
                    node_id = node_id - d_to_offset[j] # offset node_id
                    bool_cell_field[cell,node] = d_to_dface_to_tag[j][node_id]; continue;
                end
            end
        end
    end
    # New Comparing the bool matrix to the cell dof matrix and add the intface DOFs
    face_own_nodes = vcat([Gridap.ReferenceFEs.get_face_own_nodes(reffe,d) for d in 0:D]...)
    @assert size(face_own_nodes)[1] == numfaces "Inconsistent number of faces"
    VΓ = Int[]
    VΓ_dface = Int[]
    for cell = 1:numcells
        for node = 1:numfaces
            if(bool_cell_field[cell,node] == 1)
                dof_node = face_own_nodes[node] .+ numnodes*(vec_comp - 1)
                dface = cell_to_faces[cell][node].*(ones(size(dof_node)[1]))
                loc_dof = cell_dofs[cell][dof_node]
                append!(VΓ,loc_dof)
                append!(VΓ_dface,dface)
            end
        end
    end
    indices = [findfirst(==(v), VΓ) for v in sort(unique(VΓ))]
    VΓ = VΓ[indices]
    VΓ_dface = VΓ_dface[indices]
    return VΓ, VΓ_dface
end
## End get_node_dof_ids
#######################

######################
## Begin get_dface_map
function get_dface_map(S::FESpace, reffe_bg::ReferenceFE{}, face_to_cell_glue::Gridap.Geometry.FaceToCellGlue{})
    # Get the triangulation 
    trian = Gridap.FESpaces.get_triangulation(S)
    # Get the active model 
    active_model = Gridap.Geometry.get_active_model(trian) 
    # Get  Grid Topology 
    grid_topology = Gridap.Geometry.get_grid_topology(active_model)
    # Get cell_to_faces 
    cell_to_faces_active = vcat(Gridap.Geometry.get_cell_faces(grid_topology)...)
    # Get corresponding background cell_to_faces
    cell_to_faces_background = vcat(get_background_cell_to_faces(trian,reffe_bg,face_to_cell_glue)...) # Calling new function 
    # Remove repetitions in the dface ids 
    indices = [findfirst(==(v), cell_to_faces_active) for v in unique(cell_to_faces_active)]
    # Preparing ordered lists
    cell_to_faces_active = cell_to_faces_active[indices]
    cell_to_faces_background = cell_to_faces_background[indices]
    return [cell_to_faces_active cell_to_faces_background]
end
## End get_dface_map
####################

####################################
## Begin get_corresponding_dface_ids
function get_corresponding_dface_ids(SΓ_dface::AbstractVector, S_to_V_map::AbstractMatrix, pos::Int)
    SVΓ_dface = Int[]
    indices = [findfirst(==(v), S_to_V_map[:,pos]) for v in SΓ_dface]
    SVΓ_dface = S_to_V_map[indices,2]
    return SVΓ_dface
end
## End get_corresponding_dface_ids
##################################

#####################
## Begin order_dof_id
function order_dof_ids(VΓ::AbstractVector, VΓ_dface::AbstractVector, SΓ::AbstractVector, SVΓ_dface::AbstractVector, SΓ_dface::AbstractVector)
    indices_VΓ = [findfirst(==(v), VΓ_dface) for v in sort(VΓ_dface)]
    indices_SΓ = [findfirst(==(v), SVΓ_dface) for v in sort(SVΓ_dface)]
    VΓ_dface = VΓ_dface[indices_VΓ]
    SΓ_dface = SΓ_dface[indices_SΓ]
    SVΓ_dface = SVΓ_dface[indices_SΓ]
    @assert VΓ_dface == SVΓ_dface "Inconsistent DOF order"
    VΓ = VΓ[indices_VΓ]
    SΓ = SΓ[indices_SΓ]
    return 
end
## End order_dof_id
###################

##############################
## Begin get_interface_dof_ids
function get_interface_dof_ids(V::FESpace,S::FESpace,reffeᵤ::ReferenceFE{},reffeₛ::ReferenceFE{},mesh_tag::Vector{String}, face_to_cell_glue::Gridap.Geometry.FaceToCellGlue{}; vec_comp1::Int = 1, vec_comp2::Int = 1)
    VΓ, VΓ_dface =  get_node_dof_ids(reffeᵤ, mesh_tag, V; vec_comp = vec_comp1)
    SΓ, SΓ_dface =  get_node_dof_ids(reffeₛ, mesh_tag, S; vec_comp = vec_comp2)
    S_to_V_map = get_dface_map(S,reffeᵤ,face_to_cell_glue)
    SVΓ_dface = get_corresponding_dface_ids(SΓ_dface,S_to_V_map,1)
    order_dof_ids(VΓ, VΓ_dface, SΓ, SVΓ_dface, SΓ_dface)
    @assert size(VΓ) == size(SΓ) "Size of Interface DOFs along the FESpaces does not match"
    return VΓ, VΓ_dface, SΓ, SΓ_dface, SVΓ_dface
end
## End get_interface_dof_ids
############################

############################
# Begin MatrixModifier setup 
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
    x_new::AbstractVector
end

function Gridap.Algebra.symbolic_setup(ls::MatrixModifier, mat::AbstractMatrix) 
    # Getting the Interface DOFs for Matrix map
    if (!haskey(ls.params,"VΓRow"))
        VΓ = ls.params["VΓ"]
        SΓ = ls.params["SΓ"]
        dommap_V = ls.params["dommap_V"]
        dommap_S = ls.params["dommap_S"]
        VΓRow, VΓDom, SΓRow, SΓDom = [],[],[],[] # Initialisation
        map(blocks(mat)[2,2].row_partition, blocks(mat)[1,1].row_partition) do rowmap_V, rowmap_S
            DomainRowMap_V = domain_range_map(dommap_V,rowmap_V)
            OwnDomainRowMap_V = DomainRowMap_V[own_to_local(rowmap_V)]
            VΓRow = convert_interfaceDOFmap(VΓ,OwnDomainRowMap_V)
            VΓDom = convert_interfaceDOFmap_2(VΓ,OwnDomainRowMap_V)
            DomainRowMap_S = domain_range_map(dommap_S,rowmap_S)
            OwnDomainRowMap_S = DomainRowMap_S[own_to_local(rowmap_S)]
            SΓRow = convert_interfaceDOFmap(SΓ,OwnDomainRowMap_S)
            SΓDom = convert_interfaceDOFmap_2(SΓ,OwnDomainRowMap_S)
        end
        new_entries = Dict("VΓRow" => VΓRow, "VΓDom" => VΓDom, "SΓRow" => SΓRow, "SΓDom" => SΓDom)
        merge!(ls.params, new_entries)
        delete!(ls.params, "VΓ") # Delete old parameters 
        delete!(ls.params, "SΓ") # Delete old parameters 
        delete!(ls.params, "dommap_S") # Delete old parameters 
        delete!(ls.params, "dommap_V") # Delete old parameters 
    end
    MatrixModifierSymbolicSetup(ls,Gridap.Algebra.symbolic_setup(ls.ls,mat))
end

function Gridap.Algebra.numerical_setup(ss::MatrixModifierSymbolicSetup,A::AbstractMatrix,ws::Tuple{Vararg{Real}})  
    # Block Solve Allocation
    x_new = allocate_in_domain(A); fill!(x_new,0.0)
    # Matrix modification 
    map(blocks(A)[2,2].matrix_partition,blocks(A)[2,1].matrix_partition,blocks(A)[1,2].matrix_partition, 
        blocks(A)[1,1].matrix_partition) do A22,A21,A12,A11
        # Robin parameters
        αf_dir = ss.ls.params["αf"]; αf_neu = 1.0;
        αs_dir = ss.ls.params["αs"]; αs_neu = 1.0;
        if(αf_dir > 1e-4)
            αf_neu = 1/αf_dir;
            αf_dir = 1.0;
        end
        if(αs_dir > 1e-4)
            αs_neu = 1/αs_dir;
            αs_dir = 1.0;
        end
        # Nuemann interface condition
        velBlock = A22[ss.ls.params["VΓRow"],:]; 
        dispBlock = A11[ss.ls.params["SΓRow"],:];

        A12[ss.ls.params["SΓRow"],:] = αs_neu.*velBlock; 
        A11[ss.ls.params["SΓRow"],:] = αs_neu.*dispBlock;

        A22[ss.ls.params["VΓRow"],:] = αf_neu.*velBlock;
        A21[ss.ls.params["VΓRow"],:] = αf_neu.*dispBlock;

        # Dirichlet Interface condition 
        A22[ss.ls.params["VΓRow"],ss.ls.params["VΓRow"]] = A22[ss.ls.params["VΓRow"],ss.ls.params["VΓRow"]] + αf_dir.*ws[1]*sparse(I(size(ss.ls.params["VΓRow"])[1])); 
        A21[ss.ls.params["VΓRow"],ss.ls.params["SΓRow"]] = A21[ss.ls.params["VΓRow"],ss.ls.params["SΓRow"]] + -αf_dir.*ws[2]*sparse(I(size(ss.ls.params["VΓRow"])[1])); 

        A12[ss.ls.params["SΓRow"],ss.ls.params["VΓRow"]] = A12[ss.ls.params["SΓRow"],ss.ls.params["VΓRow"]] + -αs_dir.*ws[1]*sparse(I(size(ss.ls.params["SΓRow"])[1])); 
        A11[ss.ls.params["SΓRow"],ss.ls.params["SΓRow"]] = A11[ss.ls.params["SΓRow"],ss.ls.params["SΓRow"]] + αs_dir.*ws[2]*sparse(I(size(ss.ls.params["SΓRow"])[1]));
    end
    MatrixModifierNumericalSetup(ss,A,Gridap.Algebra.numerical_setup(ss.ss,A),x_new)
end

function Gridap.Algebra.numerical_setup!(ns::MatrixModifierNumericalSetup,A::AbstractMatrix,ws::Tuple{Vararg{Real}})
    ns.x_new = 0*ns.x_new; # reinitializing solution vector
    map(blocks(A)[2,2].matrix_partition,blocks(A)[2,1].matrix_partition,blocks(A)[1,2].matrix_partition, 
        blocks(A)[1,1].matrix_partition) do A22,A21,A12,A11
        # Robin parameters
        αf_dir = ns.ss.ls.params["αf"]; αf_neu = 1.0;
        αs_dir = ns.ss.ls.params["αs"]; αs_neu = 1.0;
        if(αf_dir > 1e-4)
            αf_neu = 1/αf_dir;
            αf_dir = 1.0;
        end
        if(αs_dir > 1e-4)
            αs_neu = 1/αs_dir;
            αs_dir = 1.0;
        end
        # Nuemann interface condition
        velBlock = A22[ns.ss.ls.params["VΓRow"],:]; 
        dispBlock = A11[ns.ss.ls.params["SΓRow"],:];

        A12[ns.ss.ls.params["SΓRow"],:] = αs_neu.*velBlock; 
        A11[ns.ss.ls.params["SΓRow"],:] = αs_neu.*dispBlock;

        A22[ns.ss.ls.params["VΓRow"],:] = αf_neu.*velBlock;
        A21[ns.ss.ls.params["VΓRow"],:] = αf_neu.*dispBlock;

        # Dirichlet Interface condition 
        A22[ns.ss.ls.params["VΓRow"],ns.ss.ls.params["VΓRow"]] = A22[ns.ss.ls.params["VΓRow"],ns.ss.ls.params["VΓRow"]] + αf_dir.*ws[1]*sparse(I(size(ns.ss.ls.params["VΓRow"])[1])); 
        A21[ns.ss.ls.params["VΓRow"],ns.ss.ls.params["SΓRow"]] = A21[ns.ss.ls.params["VΓRow"],ns.ss.ls.params["SΓRow"]] + -αf_dir.*ws[2]*sparse(I(size(ns.ss.ls.params["VΓRow"])[1])); 

        A12[ns.ss.ls.params["SΓRow"],ns.ss.ls.params["VΓRow"]] = A12[ns.ss.ls.params["SΓRow"],ns.ss.ls.params["VΓRow"]] + -αs_dir.*ws[1]*sparse(I(size(ns.ss.ls.params["SΓRow"])[1])); 
        A11[ns.ss.ls.params["SΓRow"],ns.ss.ls.params["SΓRow"]] = A11[ns.ss.ls.params["SΓRow"],ns.ss.ls.params["SΓRow"]] + αs_dir.*ws[2]*sparse(I(size(ns.ss.ls.params["SΓRow"])[1]));    
    end
    # ns.A = A
    Gridap.Algebra.numerical_setup!(ns.ns,A)
end

function Gridap.Algebra.solve!(x::AbstractVector,ns::MatrixModifierNumericalSetup,b::AbstractVector)
    solve!(ns.x_new,ns.ns,b)
    copy!(x, ns.x_new)
end

# Solver for GridapODE
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
    # Modifying the residual vector Jx = r
    map(blocks(r)[1].vector_partition,blocks(r)[2].vector_partition, 
        blocks(lop.usx[2])[1].vector_partition, blocks(lop.usx[1])[2].vector_partition) do r1,r2,usx2,usx1
        αf_dir = ls.params["αf"]; αf_neu = 1.0;
        αs_dir = ls.params["αs"]; αs_neu = 1.0;
        if(αf_dir > 1e-4)
            αf_neu = 1/αf_dir;
            αf_dir = 1.0;
        end
        if(αs_dir > 1e-4)
            αs_neu = 1/αs_dir;
            αs_dir = 1.0;
        end
        nuemann_block = r1[ls.params["SΓRow"]] + r2[ls.params["VΓRow"]];
        dirichlet_block = usx2[ls.params["SΓDom"]] - usx1[ls.params["VΓDom"]];
        r1[ls.params["SΓRow"]] = αs_neu.*nuemann_block - αs_dir.*dirichlet_block;
        r2[ls.params["VΓRow"]] = αf_neu.*nuemann_block + αf_dir.*dirichlet_block;
    end
    # Solving
    solve!(x, ns, r)
    ns
end

# Solver for GridapODE
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
    # Modifying the residual vector Jx = r
    map(blocks(r)[1].vector_partition,blocks(r)[2].vector_partition, 
        blocks(lop.usx[2])[1].vector_partition, blocks(lop.usx[1])[2].vector_partition) do r1,r2,usx2,usx1
        αf_dir = ls.params["αf"]; αf_neu = 1.0;
        αs_dir = ls.params["αs"]; αs_neu = 1.0;
        if(αf_dir > 1e-4)
            αf_neu = 1/αf_dir;
            αf_dir = 1.0;
        end
        if(αs_dir > 1e-4)
            αs_neu = 1/αs_dir;
            αs_dir = 1.0;
        end
        nuemann_block = r1[ls.params["SΓRow"]] + r2[ls.params["VΓRow"]];
        dirichlet_block = usx2[ls.params["SΓDom"]] - usx1[ls.params["VΓDom"]];
        r1[ls.params["SΓRow"]] = αs_neu.*nuemann_block - αs_dir.*dirichlet_block;
        r2[ls.params["VΓRow"]] = αf_neu.*nuemann_block + αf_dir.*dirichlet_block;
    end
    # Solving
    solve!(x, ns, r)
    ns
end
# End MatrixModifier setup 
##########################
  

#################################
# Starting FSI Problem Definition 
function main(distribute , parts)
    # ranks and MPI related 
    ranks = distribute(LinearIndices((prod(parts),)))
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(MPI.COMM_WORLD)

    length = 6;

    # Defining the  model 
    domain = (0,length,0,1)
    partition = (length*7,5)
    model = CartesianDiscreteModel(ranks, parts, domain, partition)

    # Define tags in the model
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels,"Top",[6,])
    add_tag_from_tags!(labels,"Left",[1,7,3])
    add_tag_from_tags!(labels,"Right",[2,4,8])
    add_tag_from_tags!(labels,"Bottom",[1,5,2])
    add_tag_from_tags!(labels,"BCTop",[3,4])
 
    # Define reference FE (Q2/P1(disc) pair)
    order = 2
    reffeᵤ = Gridap.ReferenceFEs.LagrangianRefFE(VectorValue{2,Float64},QUAD, order)
    reffeₚ = ReferenceFE(lagrangian,Float64,order-1;space=:P)
    reffeₛ = Gridap.ReferenceFEs.LagrangianRefFE(Float64,SEGMENT, order)

    # Define triangulation and integration measure
    degree = 2*order + 1
    Ωf = Triangulation(model)
    dΩf = Measure(Ωf,degree)
    Σ = BoundaryTriangulation(model,tags=["Top"])
    dΣ = Measure(Σ,degree)
    Σleft = BoundaryTriangulation(model,tags=["Left"])
    dΣleft = Measure(Σleft,degree)
    Σright = BoundaryTriangulation(model,tags=["Right"])
    dΣright = Measure(Σright,degree)

    # Define Test FE Functions 
    mfs = BlockMultiFieldStyle(2,(1,2),(1,2,3))
    V = TestFESpace(Ωf,reffeᵤ,dirichlet_tags = ["Bottom","BCTop"], 
                        dirichlet_masks =[(false,true),(false,true)])
    Q = TestFESpace(Ωf,reffeₚ;conformity=:L2,constraint=:zeromean)
    S = TestFESpace(Σ,reffeₛ,dirichlet_tags = ["BCTop"])
    Yfluid = MultiFieldFESpace([V,Q])
    Y = MultiFieldFESpace([S,V,Q]; style = mfs)
   
    # Define trial FESpaces from Dirichlet values
    gu(t) = x -> VectorValue(0.0,0.0) # Velocity BCs
    gd(t) = x -> 0.0 # Displacement BCs
    gpr(t) = x -> 0.0 # Right Pressure BCs
    function gpl(t) # Left Pressure BCs
        function inner_function(x)
            if t > 0.005
                return 0.0
            else
                return 20000.0
            end
        end
        return inner_function
    end
    U = TransientTrialFESpace(V,[gu,gu])
    P = TransientTrialFESpace(Q)
    D = TransientTrialFESpace(S, [gd])
    Xfluid = TransientMultiFieldFESpace([U,P])
    X = TransientMultiFieldFESpace([D,U,P]; style = mfs)

    # Initial conditions function
    h_u(t) = x -> VectorValue(0.0,0.0) # Velocity
    h_d(t) = x -> 0.0 # Displacement
    function h_p(t) # Pressure 
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

    # Robin-Robin parameters 
    αs = 0.0

    # Define bilinear and linear form
    ρf = 1.0;
    ρs = 45;
    ϵ = 0.1;Inf
    L = 1.0e5; # L = ϵ*E/(R²(1-ν²))

    αf = ((ρs*ϵ/0.001) + L*0.001)*(0.1)

    jac(t,(dd,du,dp),(s,v,q)) = ∫(-(∇⋅v)*dp + (∇⋅du)*q)dΩf + ∫(L*dd*s)dΣ
    jac_t(t,(dtd,dtu,dtp),(s,v,q)) = ∫(ρf*v⋅dtu)dΩf
    jac_tt(t,(dttd,dttu,dttp),(s,v,q)) = ∫(ρs*ϵ*s*dttd)dΣ
    l(t,(s,v,q)) = ∫(v⋅VectorValue(0.0, 0.0) + (q*0.0))dΩf + ∫(s*0.0)dΣ - ∫(gpr(t)*v⋅VectorValue(-1,0))dΣright - ∫(gpl(t)*v⋅VectorValue(-1,0))dΣleft

    # Build affine FE operator
    op = TransientLinearFEOperator((jac, jac_t, jac_tt), l, X, Y)

    # Extracting the interface DOFs
    VΓ ,VΓ_dface ,SΓ , SΓ_dface ,SVΓ_dface = [] ,[] ,[] ,[] ,[]
    map(local_views(V), local_views(S), local_views(Σ)) do locV, locS, locΣ
        #### Major limitation of this method is that the triangulation of S should be exactly Σ
        #### It can be overcomed by using ```get_glue(get_triangulation(S),Val(dim)).tface_to_mface]``` which gives face_to_bgface but no cell_to_faces
        VΓ,VΓ_dface,SΓ, SΓ_dface,SVΓ_dface  =  get_interface_dof_ids(locV ,locS ,
            reffeᵤ ,reffeₛ ,["Top"],locΣ.parent.glue ;vec_comp1 = 2, vec_comp2 = 1)
    end

    # Get the Gridap initialised solution vector partition map
    dommap_V = get_FESpace_map(Xfluid(0.0))
    dommap_S = get_FESpace_map(X(0.0)[1])
    params = Dict("VΓ" =>VΓ, "SΓ" => SΓ, "dommap_V" =>dommap_V, "dommap_S" =>dommap_S, "αf" =>αf, "αs" =>αs)
    
    # Time numerics 
    t0, tF = 0.0, 3.0
    dt = 0.001

    # System Solver Definition 
    
    # Defining the preconditioner
    solverBlock = TrilinosSolve()
    
    coeffs = [1.0 1.0;
              0.0 1.0;]

    bblocks = [LinearSystemBlock() LinearSystemBlock();
            LinearSystemBlock() LinearSystemBlock();]

    P = BlockTriangularSolver(bblocks,[solverBlock,solverBlock],coeffs,:upper)

    solver = FGMRESSolver(20,P;Pl=nothing,restart=false,m_add=1,maxiter=1000,atol=1e-8,rtol=1.e-6,verbose=i_am_main(ranks))

    SysSolver = MatrixModifier(params,solver)

    # ODE Solver 
    solver_ode = Gridap.ODEs.GeneralizedAlpha2(SysSolver,dt,0.5)

    # Imposing Intitial conditions
    x_initial = interpolate_everywhere([h_d(t0) , h_u(t0), h_p(t0)], X(t0))
    v_initial = interpolate_everywhere([v_d(t0) , v_u(t0), v_p(t0)], X(t0))
    a_initial = interpolate_everywhere([v_d(t0) , v_u(t0), v_p(t0)], X(t0))

    # Initializing the solving process
    xₜ = solve(solver_ode, op, t0, tF, (x_initial,v_initial,a_initial))

    # Post processing 
    if !isdir("postprocess/tmp")  
        mkdir("postprocess/tmp")
    end

    # Creating pvd files
    pvd_Ωf = createpvd(ranks,"postprocess/DNThinFSIFluids") 
    pvd_Σ = createpvd(ranks,"postprocess/DNThinFSISolids")

    # Storing initial conditions in the pvd file 
    pvd_Ωf[0] = createvtk(Ωf, "postprocess/tmp/resultsFluid_0", cellfields=["u" => x_initial[2] , "pressure" => x_initial[3]])
    pvd_Σ[0] = createvtk(Σ, "postprocess/tmp/resultsSolid_0", cellfields=["displacement" => x_initial[1]])

    # Ploting the residual vs iteration plot
    function plot_residual(solver)
        itermax = solver.log.num_iters
        res = solver.log.residuals[2:itermax+1]
        iter = collect(2:itermax+1) .- 1
        p = plot(iter,res, yscale=:log10,label="Robin-Neumann",color=:black,linewidth = 2, marker = true)
        p = plot!(xlabel = "Iteration",ylabel = "Residual")
        savefig("residual_plot.png")
        println("residuals - ",res)
        println("iterations - ", iter)
    end

    timer1 = MPI.Wtime() # Adding a timer 
    timestep = 0; 

    # Solving and storing results in the pvd file 
    mpicpp.KokkosInitialize()
    with_logger(SimpleLogger(stderr, Logging.Error)) do # Warnings in the vtu files. Need to fix later 
        for (tn, (dhn,uhn,phn)) in xₜ
            pvd_Ωf[tn] = createvtk(Ωf, "postprocess/tmp/resultsFluid_$tn", cellfields=["u" => uhn, "pressure" => phn])
            pvd_Σ[tn] = createvtk(Σ, "postprocess/tmp/resultsSolid_$tn", cellfields=["displacement" => dhn])
            timer2 = MPI.Wtime() # Adding a timer 
            SolvingTime = timer2 - timer1
            timestep += 1
            if(rank == 0) println("Solve time at timestep $(timestep) is ",SolvingTime) end
            timer1 = MPI.Wtime() # Adding a timer
            # Plotting residual plot
            if(rank == 0 && tn == tF)
                plot_residual(solver)
            end
        end
    end

    mpicpp.KokkosFinalize()

    # Sving the pvd files 
    MPI.Barrier(comm)
    if(rank == 0)
        savepvd(pvd_Ωf)
        savepvd(pvd_Σ)
        println("Results Saved!")
    end
end

# Ending FSI Problem Definition 
###############################

# Calling the main function 
with_mpi() do distribute
    rank_partition = (MPI.Comm_size(MPI.COMM_WORLD),1)
    main(distribute, rank_partition)
end
