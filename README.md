# IsingCats

Ising model with C-Sets and DPO's.

C-Sets, or copresheaves are morphisms from a category C into FinSet, the category of finite sets and functions between these finite sets.

The category C can be thought of as some indexing category, with its objects and morphisms specifying some pattern within FinSet.

For example, directed graphs can be specified by an index category with two objects, E and V, and two morphisms, the source and target morphisms, both going from E to V.

By generalising this index category, the shape of which is often referred to as a 'schema', we control what objects and morphisms are selected by the functor F within FinSet. That is, F is isomorphic to a state over the schema.

In modelling systems, we would like to update states to reflect system updates, these updates being governed by the drive to minimise a Hamiltonian or some other objective function, or just by random environmental factors (e.g. gene mutations). These updates can be performed on the C-Set via a double pushout. The double pushout takes a rewrite rule, encoded as a span in FinSet, and replaces part of the state with an updated version.

For large schema, we can specify which section of the schema we would like to perform this rewrite on, which gives some computational gain.

The Ising Model
-----------------
-----------------

States of the Ising model are isomorphic to C-Sets over the following schema:

<p align="center">
<img src="https://github.com/aj-searle/IsingCats/main/src/_static/ising_schema.png"
title="Ising model" width="150"/>
</p>

And we restrict our attention to rewrite rules in which a single vertex is changed from vertex set 1 to vertex set 2, and edges updated accordingly. In more complex schema, certain update rules might be C-Sets from subschema, in which case one can restrict the algorithm to that subset, as mentioned before. (Here, we are not interested in update rules of subschema, because they are not sufficient for generating an equilibrium distribution).

Visualisation of Ising Model states
------------------------------------

An Ising State can be created and then visualised as follows

```julia
IsingState = @acset IsingModel begin
  V1 = 2
  E = 2
  V2 = 2
  p = [1,2]
  q = [1,2]
end

to_graphviz(IsingState)
```
