L1 = @acset IsingModel begin
  V = 2
  E = 1
  src = [1]
  tgt = [2]
  dec = ["e1"]
  label = ["a","b"]
end

R1 = @acset IsingModel begin
  V = 2
  E = 1
  src = [1]
  tgt = [2]
  dec = ["e1"]
  label = ["a","b"]
end

I1 =

L = ACSetTransformation(aI2, aarr, V=[1,2]);
R = ACSetTransformation(aI2, abiarr, V=[1,2]);
m = ACSetTransformation(aarr, aspan, V=[2,1], E=[1]);  # sends 'a'->'b' and 'b'->'a'
