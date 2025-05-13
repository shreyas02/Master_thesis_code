using Gridap
using Gridap.Io
using GridapGmsh

model = GmshDiscreteModel("square.msh")

writevtk(model,"square")

fn = "square.json"
to_json_file(model,fn) 
