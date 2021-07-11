module IsingCats

using Random
using Catlab, Catlab.Theories, Catlab.CategoricalAlgebra, Catlab.CategoricalAlgebra.CSets
using Catlab.CategoricalAlgebra.DPO, Catlab.Graphs, Catlab.Present, Catlab.Graphics
using Catlab.WiringDiagrams
import Catlab.Graphics: to_graphviz

include("ConjuctionQueryHomomorphism.jl")
using ..ConjunctiveQueryHomomorphism

export IsingModel, SchemaIsingModel,SymmetricIsingModel, SymmetricSchemaIsingModel,
calculate_hamiltonian, ising_state_accept, ΔE, p_transition, rulel, ruler, rule,
rewrite_ising, accept_rewrite, symmetrise, generate_state, run_ising


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

@present SymmetricSchemaIsingModel <: SchemaIsingModel begin
  invₑ::Hom(E,E)
  inv₁::Hom(L1,L1)
  inv₂::Hom(L2,L2)

  compose(invₑ,invₑ) == id(E)
  compose(invₑ,p) == q
  compose(invₑ,q) == p

  compose(inv₁,inv₁) == id(L1)
  compose(inv₁,src1) == tgt1
  compose(inv₁,tgt1) == src1

  compose(inv₂,inv₂) == id(L2)
  compose(inv₂,src2) == tgt2
  compose(inv₂,tgt2) == src2
end

""" Abstract type for symmetric Ising model.
"""
const AbstractSymmetricIsingModel = AbstractACSetType(SymmetricSchemaIsingModel)

""" A symmetric Ising model.
"""
const SymmetricIsingModel = CSetType(SymmetricSchemaIsingModel, index=[:src1,:tgt1,:src2,:tgt2])


function symmetrise(j::IsingModel)
  k = IsingModel()
  for v₁ in parts(j, :V1)
    add_parts!(k, :V1, 1)
  end
  for v₂ in parts(j, :V2)
    add_parts!(k, :V2, 1)
  end
  for e in parts(j, :E)
    add_parts!(k, :E, 1, p=j[e, :p], q=j[e,:q])
  end
  for e₁ in parts(j, :L1)
    add_parts!(k, :L1, 1, src1=j[e₁, :src1], tgt1=j[e₁,:tgt1])
    add_parts!(k, :L1, 1, src1=j[e₁, :tgt1], tgt1=j[e₁,:src1])
  end
  for e₂ in parts(j, :L2)
    add_parts!(k, :L2, 1, src2=j[e₂, :src2], tgt2=j[e₂,:tgt2])
    add_parts!(k, :L2, 1, src2=j[e₂, :tgt2], tgt2=j[e₂,:src2])
  end
  return k
end

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

function to_graphviz(j::AbstractSymmetricIsingModel;
  prog::AbstractString="neato", graph_attrs::AbstractDict=Dict(),
  node_attrs::AbstractDict=Dict(), edge_attrs::AbstractDict=Dict(),
  node_labels::Bool=false, edge_labels::Bool=false,
  show_reflexive::Bool=false)
  pg = SymmetricPropertyGraph{Any}(; prog = prog,
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

"""Calculates Hamiltonian of Ising State"""
function calculate_hamiltonian(ising_state::CSet, J::Number=1, μ::Number=0.1)
  return J * nparts(ising_state, :E) - μ * (nparts(ising_state, :V1) - nparts(ising_state, :V2))
end

# TODO: remove this method because it is redundant
function calculate_hamiltonian(J::Number, μ::Number, ising_state::CSet)
  return calculate_hamiltonian(ising_state, J, μ)
end

ΔE(state₁::AbstractIsingModel, state₂::AbstractIsingModel) = begin
  calculate_hamiltonian(state₂) - calculate_hamiltonian(state₁)
end

"""    p_transition(state₁, state₂, β)

Compute the probability of changing from state₁ to state₂ based on the change in energy
"""
function p_transition(state₁::ACSet, state₂::ACSet, β::Float64=0.42)
  𝜰 = exp(-β * ΔE(state₁, state₂))
end

"""    ising_state_accept(state₁, state₂, β)

Calculates the probability of flipping state based on ΔE.  Takes two Ising States as arguments and returns a tuple
consisting of a boolean and the value of ΔE. Where E is the Energy Hamiltonian (sometimes called H)
"""
function ising_state_accept(state1::CSet, state2::CSet, β::Float64 = 0.42)
  𝜰 = p_transition(state1, state2, β)
  return ((rand() < 𝜰), ΔE(state1, state2))
end

function generate_state(n::Int, m::Int)
  grid = bitrand(n,m)
  #grid = Bool[1 1 0 1; 0 1 1 1; 1 0 0 0; 0 1 0 1]
  #grid = Bool[1 1 1; 0 1 0; 1 0 1]
  V1 = sum(grid)
  V2 = (n*m) - V1
  L1 = 0
  L2 = 0
  E = 0
  src1 = []
  tgt1 = []
  src2 = []
  tgt2 = []
  p = []
  q = []
  blue = findall(x -> x == true, grid)
  orange = findall(x -> x == false, grid)
  for i=2:n
    for j=1:m
      if grid[i-1,j] == grid[i,j] == true
        L1 += 1
        push!(src1, findfirst(x -> x == CartesianIndex(i,j), blue))
        push!(tgt1, findfirst(x -> x == CartesianIndex(i-1,j), blue))
      elseif grid[i-1,j] == grid[i,j] == false
        L2 += 1
        push!(src2, findfirst(x -> x == CartesianIndex(i,j), orange))
        push!(tgt2, findfirst(x -> x == CartesianIndex(i-1,j), orange))
      elseif grid[i-1,j] != grid[i,j]
        E += 1
        if grid[i,j] == true
          push!(q, findfirst(x -> x == CartesianIndex(i-1,j), orange))
          push!(p, findfirst(x -> x == CartesianIndex(i,j), blue))
        elseif grid[i,j] == false
          push!(q, findfirst(x -> x == CartesianIndex(i,j), orange))
          push!(p, findfirst(x -> x == CartesianIndex(i-1,j), blue))
        end
     end
    end
  end

  for i=1:n
    for j=2:m
      if grid[i,j-1] == grid[i,j] == true
        L1 += 1
        push!(src1, findfirst(x -> x == CartesianIndex(i,j), blue))
        push!(tgt1, findfirst(x -> x == CartesianIndex(i,j-1), blue))
      elseif grid[i,j-1] == grid[i,j] == false
        L2 += 1
        push!(src2, findfirst(x -> x == CartesianIndex(i,j), orange))
        push!(tgt2, findfirst(x -> x == CartesianIndex(i,j-1), orange))
      elseif grid[i,j-1] != grid[i,j]
        E += 1
        if grid[i,j] == true
          push!(q, findfirst(x -> x == CartesianIndex(i,j-1), orange))
          push!(p, findfirst(x -> x == CartesianIndex(i,j), blue))
        elseif grid[i,j] == false
          push!(q, findfirst(x -> x == CartesianIndex(i,j), orange))
          push!(p, findfirst(x -> x == CartesianIndex(i,j-1), blue))
        end
      end
    end
  end
  J = IsingModel()
  add_parts!(J, :V1, V1)
  add_parts!(J, :V2, V2)
  add_parts!(J, :L1, L1)
  add_parts!(J, :L2, L2)
  add_parts!(J, :E,  E)
  J[:, :src1] = src1
  J[:, :tgt1] = tgt1
  J[:, :src2] = src2
  J[:, :tgt2] = tgt2
  J[:, :p] = p
  J[:, :q] = q
  return J
end
function accept_rewrite(rewrite_span, T)

    Δ = ΔE(codom(rewrite_span[1]), codom(rewrite_span[2]))
    if Δ < 0
        return true
    elseif exp(-Δ/(T))>rand(Float64)
        return true
    else
        return false
    end
 end

include("rewrite_rules.jl")

function rewrite_ising(j::IsingCats.AbstractIsingModel, T)
  # choose a random rule
  rules = Random.shuffle([(i,b) for i in 1:4, b in 1:2])
  for (k,b) in rules
    l,r = rule(k)
    if b==2
      r,l = l,r
    end

    qₗ = homomorphism_query(codom(l))
    matches = query(j, qₗ)

    αs = map(ρ -> make_homomorphism(matches[ρ], codom(l), j),
        1:length(matches))


    # quick hack for "monic on V1"
    αs = filter(αs) do α
      length(unique(collect(components(α).V1))) == length(collect(components(α).V1))
    end

    αs = filter(α->valid_dpo(l, α), αs)
    length(αs) > 0 || println("None of the matches for rule ($k,$b) are valid DPO")
    if length(αs) > 0
        # compute table of matches
        # pick random match
        ρ = rand(1:length(αs))
        α = αs[ρ]
        if accept_rewrite((l,r), T)
          # Rewrite
          return rewrite_match(l, r, α)
        end
        return j
    end
  end
  error("Could not find a valid match in all the rules")
end

function run_ising(j::IsingCats.AbstractIsingModel, T, n::Int, f)
  vals = Any[]
  for i in 1:n
    j = rewrite_ising(j, T)
    push!(vals, f(j))
  end
  return j, vals 
end
end
