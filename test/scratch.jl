using IsingCats
using IsingCats.ConjunctiveQueryHomomorphism
using Test

using Catlab, Catlab.Theories, Catlab.CategoricalAlgebra, Catlab.CategoricalAlgebra.CSets
using Catlab.WiringDiagrams
import Catlab.Graphics: to_graphviz

function random_ising_state(height::Int, width::Int)
    lattice = bitrand(height, width) # 1 -> 'V2', 0 -> 'V1'
    i = IsingModel()
    add_parts!(i, :V1, 3)
    add_parts!(i, :V2, 1)
    add_parts!(i, :E, 2, p=[1,3], q=[1,1])
    add_parts!(i, :L1, 2, src1=[1,2], tgt1=[2,3])
    add_parts!(i, :L2, 0, src2=[], tgt2=[])
  end

  