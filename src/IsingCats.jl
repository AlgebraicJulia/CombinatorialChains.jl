module IsingCats

using Catlab, Catlab.Theories, Catlab.CategoricalAlgebra, Catlab.CategoricalAlgebra.CSets
using Catlab.CategoricalAlgebra.DPO, Catlab.Graphs, Catlab.Present, Catlab.Graphics


export IsingModel, SchemaIsingModel, calculate_hamiltonian

function calculate_hamiltonian(J::Number, μ::Number, ising_model::CSet)
  return J * length(ising_model.tables.E) - μ * (length(ising_model.tables.V1) - length(ising_model.tables.V2))
end

# Shema for the (two state) Ising model 
@present SchemaIsingModel(FreeSchema) begin
  V1::Ob
  V2::Ob
  E::Ob
  L1::Ob
  L2::Ob
  src1::Hom(L1,V1)
  tgt1::Hom(L1,V1)
  src2::Hom(L2,V2)
  tgt2::Hom(L2,V2)
  p::Hom(E, V1)
  q::Hom(E, V2)
  
end

# Make it into a type??
const IsingModel = ACSetType(SchemaIsingModel, index=[:src1,:tgt1,:src2,:tgt2]) 

end