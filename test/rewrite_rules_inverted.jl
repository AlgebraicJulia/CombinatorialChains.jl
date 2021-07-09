# Literally the same code, but with V1 and V2 etc. swapped to accomodate for the inverted center spin. Can this be compactified as well? (Not a CS person here, sorry, but happy to learn! :D )

L_array = []
I_array = []
R_array = []
l_array = []
r_array = []

for i = 0:4 # i counts number of vertices of opposite color
  L = @acset IsingModel begin
    V2=i
    V1=5-i
    E=i
    L1=4-i
    scr1 = []
    tgt1 = []
    p = []
    q = []
    for j in range(1, length=i)
      append!(p, [j])
      append!(q, [5-i])    #5-i takes indexing of central vertex into account
    for j in range(i+1,length=4-i) #Question: For i=4, this gives scr1 = tgt1 = Any[] instead of an empty array... is this a bad thing in Julia?
      append!(scr1, [j])
      append!(tgt1, [5-i]) #5-i takes indexing of central vertex into account
  end
  I = @acset IsingModel begin
    V2=i
    V1=4-i
    # no edges at all in the interface
  end
  R = @acset IsingModel begin
    V2=i+1
    V1=4-i
    E=4-i
    L2=i
    p=[]
    q=[]
    src2=[]
    tgt2=[]
    for j in range(i+1,length=4-i) 
      append!(p, [j])
      append!(q, [i+1])     #i+1 takes indexing of central vertex into account
    for j in range(1, length=i)
      append!(scr2, [j])
      append!(tgt2, [i+1])  #i+1 takes indexing of central vertex into account
  end

  # It turns out that in our setting, l and r are the same morphisms since the boundary does not change.
  h1 = []
  h2 = []
  for j in range(1,length=4-i)
    append!(h1,j)
  for j in range(1,length=i)
    append!(h2,j)
  l = ACSetTransformation(I, L, V1=h1, V2=h2);
  r = ACSetTransformation(I, R, V1=h1, V2=h2);

  append!(L_array, L)
  append!(I_array, I)
  append!(R_array, R)
  append!(l_array, l)
  append!(r_array, r)
end
