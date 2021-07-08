module IsingCats

using Catlab, Catlab.Theories, Catlab.CategoricalAlgebra, Catlab.CategoricalAlgebra.CSets
using Catlab.CategoricalAlgebra.DPO, Catlab.Graphs, Catlab.Present, Catlab.Graphics
import Catlab.Graphics: to_graphviz

include("ConjuctionQueryHomomorphism.jl")

export IsingModel, SchemaIsingModel, calculate_hamiltonian

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

function calculate_hamiltonian(J::Number, μ::Number, ising_model::CSet)
  return J * length(ising_model.tables.E) - μ * (length(ising_model.tables.V1) - length(ising_model.tables.V2))
end

function ΔE(span)
end

function generate_state(N::int)
end

function accept_rewrite(rewrite_span, T)

    ΔE = ΔE(rewrite_span)

    if ΔE < 0:
        return true

    elif exp(-ΔE/(T))>rand(Float64, 2):
        return true
    else:
        return false
 end

function rewrite_ising(j::AbstractIsingModel, rewrite_rules, T):

  dpo_search = true

  while dpo_search:
    # Compute table of matches
    span = rewrite_rules[rand(1:length(reqrite_rules))]
    qₗ = homomorphism_query(codom(left(span)))
    matches = query(b, qₗ)
    α = make_homomorphism(matches[rand(1:length(rewrite_rules))], codom(left(span)), j)

    if valid_dpo(left(span), α) and accept_rewrite(span, T):
      # Rewrite
      rewrite_match(α, right(span), left(span))
      dpo_search = false

    end
  end


end

function find_equilibrium()
end

end
