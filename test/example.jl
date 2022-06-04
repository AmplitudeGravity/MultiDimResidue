#push!(LOAD_PATH, "/Users/gangchen");
using MultiResidue
using SymEngine
#using LinearAlgebra
@varss z 4
@vars a b c
@vars x y z
@funs p q
h1 = a*z1*z2 + z3*z1 + z2^2;
h2 = z1*z3 - 2b*z2^2;
h3 = z1^2 - z2*z1 + z3^2;
@time res1=multiResidue((2z2+1)/(z1^2-z3-1),[h1,h2,h3],[z1,z2,z3]);
res1


h1 = z1*z2 + z3*z1 + z1^2;
h2 = z1*z3 - 2b*z2^2;
h3 = z1^2 - z2*z1 + z1^2;
@time res2=multiResidue((2z2+1)/(z1^2-z3-1),[h1,h2,h3],[z1,z2,z3]);
res2


h1 = z1*z2 + z3*z4 + z2^2;
h2 = z1^3*z3 - 2z2^4;
h3 = z1^2 - z2*z1 + z3^2;
h4 = z2^2 - z2*z3 + z4^2;
@time res3=multiResidue((2z2+1)/(z1^2-z3-1),[h1,h2,h3,h4],[z1,z2,z3,z4]);
res3



#ideal=[h1,h2,h3];
#vars=[z1, z2, z3];
#MultiResidue.rankSym([h1 h2 h3; h3 h4 h2; h1+h2 h3*h3 h4*h1])
#MultiResidue.homoEqn([h1,h2,h3,h4],[z1,z2,z3,z4],3)
#subs([h1 h2; h3 h4],Dict(z1=>2,z2=>33))
#inEqn=MultiResidue.inhomoEqn([h1,h2,h3,h4],[z1,z2,z3,z4],3)


