push!(LOAD_PATH, "/Users/gangchen");
using MultiResidue
using SymEngine
@varss z 3
h1=z1*z2 + z1*z3;
h2=z1^3*z2-2z2^4;
h3=(z1^2 - z2*z1 + z3^2);


multiResidue((z2+1)/(z1^2-z3-1),[h1,h2,h3],[z1,z2,z3])

multiResidue((z2+1)/(z1^2-z3-1),[z1, z2,(z1 + z2)*(z1 + z2 + z3)],[z1,z2,z3])