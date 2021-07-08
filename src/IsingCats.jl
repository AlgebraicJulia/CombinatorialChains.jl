module IsingCats

using Catlab, Catlab.Theories, Catlab.CategoricalAlgebra, Catlab.CategoricalAlgebra.CSets
using Catlab.CategoricalAlgebra.DPO, Catlab.Graphs, Catlab.Present, Catlab.Graphics
import Catlab.Graphics: to_graphviz


export IsingModel, SchemaIsingModel, calculate_hamiltonian, ising_state_accept

"Calculates Hamiltonian of Ising State"
function calculate_hamiltonian(ising_state::CSet, J::Number=1, μ::Number=0.1)
  return J * length(ising_state.tables.E) - μ * (length(ising_state.tables.V1) - length(ising_state.tables.V2))
end

"""
Calculates the probability of flipping state based on ΔH.  Takes two Ising States as arguments and returns a tuple 
consisting of a boolean and the value of ΔH.
"""
function ising_state_accept(state1::CSet, state2::CSet, β::Float64 = 0.42)
  ΔE = calculate_hamiltonian(state2) - calculate_hamiltonian(state1)
  𝜰 = exp(-β * ΔE)
  return ((rand() < 𝜰), ΔE)
end


# Schema for the (two state) Ising model 
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

# Create abstract and concrete types for IsingModel
const AbstractIsingModel = AbstractACSetType(SchemaIsingModel)
const IsingModel = ACSetType(SchemaIsingModel, index=[:src1,:tgt1,:src2,:tgt2]) 

function to_graphviz(j::AbstractIsingModel;
  prog::AbstractString="neato", graph_attrs::AbstractDict=Dict(),
  node_attrs::AbstractDict=Dict(), edge_attrs::AbstractDict=Dict(),
  node_labels::Bool=false, edge_labels::Bool=false,
  show_reflexive::Bool=false)
  pg = PropertyGraph{Any}(; prog = prog,
    graph=graph_attrs,
    #node = merge!(default_node_attrs(node_labels), node_attrs),
    node = merge!(Dict(:style => "filled"), node_attrs),
    edge = merge!(Dict(:arrowsize => "0.5"), edge_attrs),
  )

  shiftv2(v) = v + nparts(j, :V1)

  # blue vertices in spin 1
  for v₁ in parts(j, :V1)
    vₚ = add_vertex!(pg)
    set_vprop!(pg, vₚ, :color, "#6C9AC3")
  end

  # orange vertices in spin 2
  for v₂ in parts(j, :V2)
    vₚ = add_vertex!(pg)
    set_vprop!(pg, vₚ, :color, "#E28F41")
  end

  # grey edges in L1
  for e in parts(j, :L1)
    eₚ = add_edge!(pg, j[e, :src1], j[e,:tgt1])
    set_eprop!(pg, eₚ, :color,"#B1B3B6")
  end

  # grey edges in L2
  for e in parts(j, :L2)
    s = j[e, :src2] |> shiftv2
    t = j[e,:tgt2] |> shiftv2
    eₚ = add_edge!(pg, s, t)
    set_eprop!(pg, eₚ, :color,"#B1B3B6")
  end

  # bright orange edges in E
  for e in parts(j, :E)
    # adding spin violations, have to shift the targets (q) becayse V is V1 + V2
    s = j[e, :p]
    t = j[e,:q] |> shiftv2
    eₚ = add_edge!(pg, s, t)
    set_eprop!(pg, eₚ, :color,"#FA4616")
  end

  # call the underlying implementation for property graphs
  return to_graphviz(pg)
end

include("ConjuctionQueryHomomorphism.jl")
end
