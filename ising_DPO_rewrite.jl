using Catlab, Catlab.Theories, Catlab.CategoricalAlgebra, Catlab.CategoricalAlgebra.CSets
using Catlab.CategoricalAlgebra.DPO, Catlab.Graphs, Catlab.Present, Catlab.Graphics


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


# Create an instance of that schema
i = IsingModel()
add_parts!(i, :V1, 3)
add_parts!(i, :V2, 1)
add_parts!(i, :E, 2, p=[1,3], q=[1,1])
add_parts!(i, :L1, 2, src1=[1,2], tgt1=[2,3])
add_parts!(i, :L2, 0, src2=[], tgt2=[])

# Test Hamiltonian function

a = IsingModel()
add_parts!(a, :V1, 6)
add_parts!(a, :V2, 0)
add_parts!(a, :E, 0, p=[], q=[])
add_parts!(a, :L1, 7, src1=[1,2,3,6,5,4,5], tgt1=[2,3,6,5,4,1,2])
add_parts!(a, :L2, 0, src2=[], tgt2=[])

b = IsingModel()
add_parts!(b, :V1, 4)
add_parts!(b, :V2, 2)
add_parts!(b, :E, 2, p=[2,4], q=[1,2])
add_parts!(b, :L1, 4, src1=[1,2,4,3], tgt1=[2,4,3,1])
add_parts!(b, :L2, 1, src2=[1], tgt2=[2])

c = IsingModel()
add_parts!(c, :V1, 2)
add_parts!(c, :V2, 4)
add_parts!(c, :E, 2, p=[1,2], q=[1,3])
add_parts!(c, :L1, 1, src1=[1], tgt1=[2])
add_parts!(c, :L2, 4, src2=[1,2,4,3], tgt2=[2,4,3,1])

ising_states = [a,b,c]

h_vec = [calculate_hamiltonian(1.0, 0.1, state) for state in ising_states]

# Manual creation of one rewrite rule
R = @acset IsingModel begin
  V1=4
  V2=1
  E=4
  p=[1,2,3,4]
  q=[5,5,5,5]
  # no V2s, so no Es or L2s
end
I = @acset IsingModel begin
  V1=4
  V2=0
  # no edges at all in the interface
end
R = @acset IsingModel begin
  V1=5
  V2=0
  L1=4
  src1=[1,2,3,4]
  tgt1=[5,5,5,5]
  # no V2s, so no Es or L2s
end

l = ACSetTransformation(I, L, V1=[1,2,3,4]);
r = ACSetTransformation(I, R, V1=[1,2,3,4]);




