""" Graphviz support for Catlab's graph types.
"""
module GraphvizIsingGraphs
export parse_graphviz, to_graphviz

using StaticArrays: StaticVector, SVector

import Graphviz

using Catlab.CategoricalAlgebra.CSets
using Catlab.Present, Catlab.Graphics, Catlab.Graphs, Catlab.Theories

function to_graphviz(g::IsingGraph;
    prog::AbstractString="dot", graph_attrs::AbstractDict=Dict(),
    node_attrs::AbstractDict=Dict(), edge_attrs::AbstractDict=Dict(),
    node_labels::Bool=false, edge_labels::Bool=false,
    show_reflexive::Bool=false)
  pg = PropertyGraph{Any}(; prog = prog,
    graph = graph_attrs,
    node = merge!(default_node_attrs(node_labels), node_attrs),
    edge = merge!(Dict(:arrowsize => "0.5"), edge_attrs),
  )
  for v1 in parts(g, :V1)
    add_vertex!(pg, label=node_labels ? string(v) : "")
  end
  for v2 in parts(g, :V2)

  #Assume all edges are reflexive
  for e in parts(g, :E)
    if show_reflexive
      e′ = add_edge!(pg, p(g,e), q(g,e))
      if is_reflexive; set_eprop!(pg, e′, :style, "dashed") end
      if edge_labels; set_eprop!(pg, e′, :label, string(e)) end
    end
  end
  for e in parts(g, :L1)
    if show_reflexive
      e′ = add_edge!(pg, src1(g,e), tgt1(g,e))
      if is_reflexive; set_eprop!(pg, e′, :style, "dashed") end
      if edge_labels; set_eprop!(pg, e′, :label, string(e)) end
    end
  end
  for e in parts(g, :L2)
    if show_reflexive
      e′ = add_edge!(pg, src2(g,e), tgt2(g,e))
      if is_reflexive; set_eprop!(pg, e′, :style, "dashed") end
      if edge_labels; set_eprop!(pg, e′, :label, string(e)) end
    end
  end

  to_graphviz(pg)
end
