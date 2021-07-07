"""
Create a undirected wiring diagram representing a query
which when applied to an C-set returns data representing every
homomorphism into that C-set.

Usage:
```julia
g = @acset Graph begin
    V = 3
    E = 2
    src = [1,2]
    tgt = [2,3]
end

g_query = homomorphism_query(g)

h = @acset Graph begin
    V = 4
    E = 4
    src = [1,1,2,3]
    tgt = [2,3,4,4]
end

q = query(h, g_query)

α = make_homomorphism(q[1],g,h)
```
"""
module ConjunctiveQueryHomomorphism

export homomorphism_query, make_homomorphism

using Catlab, Catlab.Programs, Catlab.CSetDataStructures, Catlab.CategoricalAlgebra.CSets
using Catlab.Programs.RelationalPrograms: UntypedNamedRelationDiagram

using Catlab.Theories: Schema, FreeSchema, SchemaType,
  CatDesc, CatDescType, ob, hom, dom, codom, codom_num

function symbolify(x::Symbol, i::Int)
  Symbol(string(x) * "_" * string(i))
end

function homomorphism_query(acs::ACSet{CD}) where {CD}
  reldia = UntypedNamedRelationDiagram{Symbol, Symbol}()
  junction_map = Dict{Tuple{Symbol,Int},Int}()

  for x in ob(CD)
    for i in parts(acs,x)
      name = symbolify(x,i)
      junction_map[(x,i)] = junction = add_part!(reldia, :Junction, variable=name)
      add_part!(reldia, :OuterPort, outer_junction=junction, outer_port_name=name)
    end
  end

  for x in ob(CD)
    outgoing = filter(f -> dom(CD,f) == x, hom(CD))
    for i in parts(acs,x)
      box = add_part!(reldia, :Box, name=x)
      add_part!(reldia, :Port, box=box, junction=junction_map[(x,i)], port_name=:_id)
      for f in outgoing
        y = codom(CD,f)
        add_part!(reldia, :Port, box=box, junction=junction_map[(y,acs[i,f])], port_name=f)
      end
    end
  end

  reldia
end

function make_homomorphism(row::NamedTuple, acs_src::ACSet{CD}, acs_tgt::ACSet{CD}) where {CD}
  components = NamedTuple{ob(CD)}([row[symbolify(x,i)] for i in parts(acs_src, x)] for x in ob(CD))
  ACSetTransformation(components, acs_src, acs_tgt)
end

end
