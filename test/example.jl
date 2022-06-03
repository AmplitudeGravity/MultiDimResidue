#push!(LOAD_PATH, "/Users/gangchen");
using MultiResidue
using SymEngine
#using LinearAlgebra
@varss z 3
@vars a b c
@vars x y z
h1=a*z1*z2 + z1*z3;
h2=z1^3*z2-2b*z2^4;
h3=z1^2 - z2*z1 + z3^2;
ideal=[h1,h2,h3];
vars=[z1, z2, z3];
 
MultiResidue.homoEqn([h1,h2,h3],[z1,z2,z3],5)

MultiResidue.inhomoEqn([h1,h2,h3],[z1,z2,z3],5)

MultiResidue.eqnAnsatz([h1,h2,h3],[z1,z2,z3],5)

res2=multiResidue((2z2+1)/(z1^2-z3-1),[h1,h2,h3],[z1,z2,z3])

