# MultiDimResidue
This is the Julia package for Multidimensional Residues around zero-dimension intersection point. 
The main function is called "multiResidue" which is used to calculate $$\text{Res}_{h_1=\cdots=h_r=0} ~ {f dz_1\wedge dz_2\wedge \cdots \wedge dz_r\over h_1 h_2\cdots h_r},$$
where $h_i$ is homogeneous polynomial function and the solution of $h_1=\cdots=h_r=0$ is an isolated point $p$. The $f$ function is regular at $p$.

To use the package, you need install [SymEngine](https://github.com/symengine/SymEngine.jl) 

For general in-homogeneous polynomial function, it is easy to transform to the homogeneous cases by a trick using the global residue theorem. For more details, please
see the origine papers [arxiv 1609.07621](https://arxiv.org/pdf/1609.07621.pdf) and [arxiv 1709.08503](https://arxiv.org/pdf/1709.08503.pdf).

# Examples for the julia code
The package use the CAS SymEngine
```
using MultiResidue
using SymEngine
@vars x y z
```
## Example 1
A example is consider the residue around the intersection of three divisors 
```
h1=x;
h2=y*(x + 2y);
h3=x^2 + y*x + 3z*z;
```
and then calculate the counter intersection around the origin for the function f
```
f=(2x + 3y + 4z)/(z - 2)
multiResidue(f,[h1,h2,h3],[x,y,z])
```
The out put is
```
-1/8
```

## Example 2
Another example is 
```
h1=y - 2x + z;
h2=x^2 *y - x*z^2 + z^3;
h3=y^4 + x*y*z^2 - y^2*z^2 + z^4;
multiResidue((2x + 3y + 4z)/(z - 2),[h1,h2,h3],[x,y,z])
```
The out put is
```
93/1664
```
## Example 3
The third example, using macro @varss to generate multi variables
```
@varss z 3
h1=z1*z2 + z1*z3;
h2=z1^3*z2-2z2^4;
h3=(z1^2 - z2*z1 + z3^2);
multiResidue((z2+1)/(z1^2-z3-1),[h1,h2,h3],[z1,z2,z3])
```
The output is 
```
-5/3
```

## Example 4
The fourth example, for the equations with symbolic coefficients
```
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

res=multiResidue((2z2+1)/(z1^2-z3-1),[h1,h2,h3],[z1,z2,z3])
```
The output is 
```
-16*(-6 - 6*(4 + 8*b/a^2)/(4 - 4*a^(-2)))/(96*b/a^3 + (-3/2)*(128*b - 64*b/a^2)/a + (-1/2)*(-192*b/a^3 - 3*(64*b/a^2 + 64*a^2)/a)*(4 + 8*b/a^2)/(4 - 4*a^(-2))) 
+ 32*(6*(4 + 8*b/a^2)/(a*(4 - 4*a^(-2))) + 6*a^(-1))/(96*b/a^3 + (-3/2)*(128*b - 64*b/a^2)/a 
+ (-1/2)*(-192*b/a^3 - 3*(64*b/a^2 + 64*a^2)/a)*(4 + 8*b/a^2)/(4 - 4*a^(-2))) 
+ 384*b/(a^3*(96*b/a^3 + (-3/2)*(128*b - 64*b/a^2)/a + (-1/2)*(-192*b/a^3 - 3*(64*b/a^2 + 64*a^2)/a)*(4 + 8*b/a^2)/(4 - 4*a^(-2))))
+ 576*b/(a^2*(96*b/a^3 + (-3/2)*(128*b - 64*b/a^2)/a + (-1/2)*(-192*b/a^3 - 3*(64*b/a^2 + 64*a^2)/a)*(4 + 8*b/a^2)/(4 - 4*a^(-2)))) 
- 384*(4 + 8*b/a^2)/((4 - 4*a^(-2))*(96*b/a^3 + (-3/2)*(128*b - 64*b/a^2)/a + (-1/2)*(-192*b/a^3 - 3*(64*b/a^2 + 64*a^2)/a)*(4 + 8*b/a^2)/(4 - 4*a^(-2))))
- 576*(4 + 8*b/a^2)/(a*(4 - 4*a^(-2))*(96*b/a^3 + (-3/2)*(128*b - 64*b/a^2)/a + (-1/2)*(-192*b/a^3 - 3*(64*b/a^2 + 64*a^2)/a)*(4 + 8*b/a^2)/(4 - 4*a^(-2))))
```
Nowadays, we still do not have any efficient simplifcation code in Julia CAS. The SymEngine is updated so slow in recent years. However, by mathematica, the result can be simplified to 
```
-((a^3 + 2 a^4 + 2 a^5 + 4 b + 6 a b + a^2 (2 + 4 b))/(
 a^6 + 6 a^2 b + 2 b (-1 + 2 b))).
```



## Example of Non-zero dimensional residue
```
multiResidue((z2+1)/(z1^2-z3-1),[z1, z2,(z1 + z2)*(z1 + z2 + z3)],[z1,z2,z3])
```
The output is not fully determined and you will get 
```
false
```
The general math meaning or definition of such residue is still an open question. 

## Example for using the solve for linear equations
symbolic solution of linear equations
```
eqns=[x-y-1, 2(x + y)]
vars=[x,y]
solve(eqns,vars)
```

# Citation 
If you use **multiResidue.jl**, please cite the two papers [arxiv 1609.07621](https://arxiv.org/pdf/1609.07621.pdf) and [arxiv 1709.08503](https://arxiv.org/pdf/1709.08503.pdf) as following

```
@article{Chen:2016fgi,
    author = "Chen, Gang and Cheung, Yeuk-Kwan E. and Wang, Tianheng and Xu, Feng",
    title = "{A differential operator for integrating one-loop scattering equations}",
    eprint = "1609.07621",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    doi = "10.1007/JHEP01(2017)028",
    journal = "JHEP",
    volume = "01",
    pages = "028",
    year = "2017"
}
```

```
@article{Chen:2017bug,
    author = "Chen, Gang and Wang, Tianheng",
    title = "{BCJ Numerators from Differential Operator of Multidimensional Residue}",
    eprint = "1709.08503",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    doi = "10.1140/epjc/s10052-019-7604-8",
    journal = "Eur. Phys. J. C",
    volume = "80",
    number = "1",
    pages = "37",
    year = "2020"
}
```
