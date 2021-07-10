#include("../src/IsingCats.jl")
using IsingCats
using IsingCats.ConjunctiveQueryHomomorphism
using Test

using Catlab, Catlab.Theories, Catlab.CategoricalAlgebra, Catlab.CategoricalAlgebra.CSets
using Catlab.WiringDiagrams
import Catlab.Graphics: to_graphviz


# Create an instance of that schema
@testset "Constructors" begin
i = IsingModel()
add_parts!(i, :V1, 3)
add_parts!(i, :V2, 1)
add_parts!(i, :E, 2, p=[1,3], q=[1,1])
add_parts!(i, :L1, 2, src1=[1,2], tgt1=[2,3])
add_parts!(i, :L2, 0, src2=[], tgt2=[])

@test nparts(i, :V1) == 3
@test nparts(i, :V2) == 1

i₂ = symmetrise(i)

@test nparts(i₂, :V1) == 3
@test nparts(i₂, :V2) == 1
@test nparts(i₂, :E) == 2
@test nparts(i₂, :L1) == 4
@test nparts(i₂, :L2) == 0
end

#Test Hamiltonian function

@testset "Hamiltonians" begin
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

@test [-0.6000000000000001, 1.8, 2.2] == [calculate_hamiltonian(state) for state in ising_states]
end


@testset "Rules" begin
# Manual creation of one rewrite rule
L = @acset IsingModel begin
  V1=4
  V2=1
  E=4
  p=[1,2,3,4]
  q=[1,1,1,1]
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

@test is_natural(l)
@test is_natural(r)

end


@testset "Conjunctive Queries" begin
  a = @acset IsingModel begin
    V1 = 2
    E = 2
    V2 = 2
    p = [1,2]
    q = [1,2]
  end

  b = @acset IsingModel begin
    V1 = 4
    E = 2
    V2 =2
    L1 = 4
    L2 = 1
    src1 = [1,2,4,3]
    tgt1 = [2,4,3,1]
    p = [4,2]
    q = [2,1]
    src2 = 1
    tgt2 = 2
  end

  qₐ = homomorphism_query(a)
  @test nparts(qₐ, :Box) == 6
  @test nparts(qₐ, :Port) == 10
  @test nparts(qₐ, :OuterPort) == 6
  @test nparts(qₐ, :Junction) == 6


  matches = query(b, qₐ)
  @test length(matches) == 4
  for i in 1:length(matches)
    @test is_natural(make_homomorphism(matches[i], a, b))
  end

end

@testset "Ising state acceptance" begin
  state_a = @acset IsingModel begin
    V1 = 2
    E = 3
    V2 =3
    L1 = 1
    L2 = 0
    src1 = [1]
    tgt1 = [2]
    p = [2,2,2]
    q = [1,2,3]
    src2 = []
    tgt2 = []
  end

  state_b = @acset IsingModel begin
    V1 = 1
    E = 1
    V2 = 4
    L1 = 0
    L2 = 3
    src1 = []
    tgt1 = []
    p = [1]
    q = [2]
    src2 = [2,2,2]
    tgt2 = [1,3,4]
  end

  @test calculate_hamiltonian(state_a) == 3.1
  @test calculate_hamiltonian(state_b) == 1.3
  @test ΔE(state_a, state_b) == -1.8
  @test p_transition(state_a, state_b) ≈ 2.1297401990391056
  @test p_transition(state_a, state_b, 0.5) > 2.1297401990391056
  @test p_transition(state_a, state_b, 0.5) < p_transition(state_a, state_b, 0.6)
  @test p_transition(state_a, state_b, 0.6) < p_transition(state_a, state_b, 0.7)
  @test p_transition(state_a, state_b, 0.5) < p_transition(state_a, state_b, 0.7)
end


a = @acset IsingModel begin
  V1 = 2
  E = 2
  V2 = 2
  p = [1,2]
  q = [1,2]
end

b = @acset IsingModel begin
  V1 = 4
  E = 2
  V2 =2
  L1 = 4
  L2 = 1
  src1 = [1,2,4,3]
  tgt1 = [2,4,3,1]
  p = [4,2]
  q = [2,1]
  src2 = 1
  tgt2 = 2
end

β = @acset SymmetricIsingModel begin
  V1 = 2
  E = 2
  V2 = 2
  p = [1,2]
  q = [1,2]
end

to_graphviz(a)
to_graphviz(b)
to_graphviz(β)
to_graphviz(homomorphism_query(a), junction_labels=:variable, box_labels=:name)
  # to_graphviz(homomorphism_query(homomorphism_query(a)), junction_labels=:variable, box_labels=:name)


@testset "Rule Generation" begin
  rulel(0)
  rulel(0) |> codom |> to_graphviz
  rulel(1) |> codom |> to_graphviz
  rulel(2) |> codom |> to_graphviz
  rulel(3) |> codom |> to_graphviz
  rulel(4) |> codom |> to_graphviz

  ruler(0) |> codom |> to_graphviz
  ruler(1) |> codom |> to_graphviz
  ruler(2) |> codom |> to_graphviz
  ruler(3) |> codom |> to_graphviz
  ruler(4) |> codom |> to_graphviz
end

# @testset "Sampling" begin
  J₀ = @acset IsingModel begin
    V1 = 9
    L1 = 12
    src1 = [1,2,5,5,5,5,7,8,3,6,1,4]
    tgt1 = [2,3,2,4,6,8,8,9,6,9,4,7]
  end

  J₀ˢ = symmetrise(J₀)

  m = homomorphisms(codom(rule(4)[1]), J₀ˢ, monic=true)[1]
  J₁ = rewrite_match(rule(4)..., m)
  to_graphviz(J₁)


  function rewrite_ising(j::IsingCats.AbstractIsingModel, T, maxrules=25, maxtries=100)
    # choose a random rule
    for k in 1:maxrules
      l,r = rule(rand(0:4))
      if rand(Bool)
        r,l = l,r
      end

      qₗ = homomorphism_query(codom(l))
      matches = query(j, qₗ)
      @show length(matches)

      αs = map(ρ -> make_homomorphism(matches[ρ], codom(l), j),
         1:length(matches))


      # quick hack for "monic on V1"
      αs = filter(αs) do α
        length(unique(collect(components(α).V1))) == length(collect(components(α).V1))
      end

      αs = filter(α->valid_dpo(l, α), αs)
      @show length(αs)
      if length(αs) > 0
        for i in 1:maxtries
          # compute table of matches
          # pick random match
          ρ = rand(1:length(αs))
          α = αs[ρ]
          @show valid_dpo(l, α)

          if valid_dpo(l, α) && accept_rewrite((l,r), T)
            # Rewrite
            j′ = rewrite_match(l, r, α)
            return j′
          end
        end
      end
    end

    error("Could not find a valid match in $maxtries attempts")
  end
#
  J₀ = @acset IsingModel begin
    V1 = 9
    L1 = 12
    src1 = [1,2,5,5,5,5,7,8,3,6,1,4]
    tgt1 = [2,3,2,4,6,8,8,9,6,9,4,7]
  end
  J₀ˢ = symmetrise(J₀)
  J = rewrite_ising(J₀ˢ, 2, 40, 100)
  to_graphviz(J₀)
  to_graphviz(J)
  @test nparts(J, :V2) == 1
# end


function run_ising(j::IsingCats.AbstractIsingModel, T, n::Int, f)
  vals = Any[]
  for i in 1:n
    j = rewrite_ising(j, T, 40, 100)
    push!(vals, f(j))
  end
  return j
end

@testset "sampler" begin
  J₀ = @acset IsingModel begin 
    V1 = 9
    L1 = 12
    src1 = [1,2,5,5,5,5,7,8,3,6,1,4]
    tgt1 = [2,3,2,4,6,8,8,9,6,9,4,7]
  end 
  J₀ = symmetrise(J₀)
  J = run_ising(J₀, 2, 6, calculate_hamiltonian)
  to_graphviz(J)
  @test nparts(J, :V2) <= 1
  J = run_ising(J₀, 2, 5, calculate_hamiltonian)
  @test nparts(J, :V2) <= 1
  to_graphviz(J)
end
