using Gridap

# Adding the model to the script 
domain = (0,5,0,5)
partition = (10,10)
model = CartesianDiscreteModel(domain,partition)
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"sides",[1,2,3,4,5,6,7,8])

# Defining the finite element spaces 
order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
V0 = TestFESpace(model,reffe;conformity=:H1,dirichlet_tags=["sides"])

g(x) = 0.0
Ug = TrialFESpace(V0,g)

# Setting up numerical integration 
degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

# Weak form 

function α(x)
    temp1 = x[1]
    if temp1<=2.5 
        return 1
    end
    if temp1>2.5
        return 1
    end
end

f(x) = 1

a(u,v) = ∫(α*∇(v)⋅∇(u) )*dΩ
b(v) = ∫( v * f )*dΩ 

# Building the Finite element problem 

op = AffineFEOperator(a,b,Ug,V0)

# Choosing the type of solver 

ls = LUSolver()
solver = LinearFESolver(ls)

# Solving the Linear system of equations

uh = solve(solver,op)

# Post processing 
writevtk(Ω,"square",cellfields=["uh"=>uh])