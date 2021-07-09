function rulel(i::Integer)
  0 <= i <= 4 || error("number of same color neighbors must be in [0,4]")
  I = IsingModel()
  add_parts!(I, :V1, i)
  add_parts!(I, :V2, 4-i)

  L = IsingModel()
  add_parts!(L, :V1, 1+i)
  add_parts!(L, :V2, 4-i)
  add_parts!(L, :E,  4-i)
  add_parts!(L, :L1,  i)

  L[:, :p] = 1
  L[:, :q] = parts(L,:V2)

  L[:, :src1] = 1
  L[:, :tgt1] = 2:nparts(L,:V1)

  h1 = i > 0 ? (2:nparts(L,:V1)) : Int[]
  h2 = 1:4-i

  @assert length(h1) == nparts(I, :V1)
  @assert length(h2) == nparts(I, :V2)

  l = ACSetTransformation(I, L, V1=h1, V2=h2);
  return l
end

function ruler(i::Integer)
  0 <= i <= 4 || error("number of same color neighbors must be in [0,4]")
  I = IsingModel()
  add_parts!(I, :V1, i)
  add_parts!(I, :V2, 4-i)

  R = IsingModel()
  add_parts!(R, :V1, i)
  add_parts!(R, :V2, 4-i+1)
  add_parts!(R, :E,  i)
  add_parts!(R, :L2, 4-i)

  R[:, :p] = parts(R,:V1)
  R[:, :q] = nparts(R,:V2)

  R[:, :src2] = nparts(R,:V2)
  R[:, :tgt2] = 1:nparts(R,:V2)-1

  h1 = 1:nparts(R,:V1)
  h2 = 1:4-i

  # @show h1
  # @show h2

  # @show nparts(I, :V1)
  # @show nparts(I, :V2)

  @assert length(h1) == nparts(I, :V1)
  @assert length(h2) == nparts(I, :V2)


  r = ACSetTransformation(I, R, V1=h1, V2=h2);
  return r
end

rule(i::Int) = rulel(i), ruler(i)
