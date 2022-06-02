module MultiResidue
using SymEngine
using LinearAlgebra
using InvertedIndices
export MultiResidue, gd, 𝒟,degree, FrobeniusSolve, solve
function degree(f::Basic)
      if SymEngine.get_symengine_class(f)==:Add
             ls=get_args(f); 
             degvec = Array{Int64}(undef, length(ls))
             for i=1:length(ls)
                  degvec[i]=degree(ls[i])
             end
             maximum(degvec)
      elseif SymEngine.get_symengine_class(f)==:Mul;
            ls=get_args(f); 
             degvec = Array{Int64}(undef, length(ls))
             for i=1:length(ls)
                  degvec[i]=degree(ls[i])
             end
             sum(degvec)
      elseif SymEngine.get_symengine_class(f)==:Pow
            ls=get_args(f); 
            ls[2]*degree(ls[1])
      elseif SymEngine.get_symengine_class(f)==:Symbol
            1
      else 
            0
      end
end
function int(x) 
       floor(Int, x)
end
function FrobeniusSolve(vec::Vector{Int64},vs::Int64)
      res=[];
      js=vs ./vec.|>int;
      for i in Iterators.product((0:k for k in js)...)
             dot(vec,i)==vs ? push!(res,collect(i)) : nothing
      end
      res
end
function gd(f::Basic,var::Basic,order::Int64)
      if order==0
            f
      else
            xs=fill(var,order);
            diff(f,xs...)
      end
end
function gd(f::Basic,vars::Vector{Basic},orders::Vector{Int64})
      g=f
      for i=1:length(vars)
            g=gd(g,vars[i],orders[i])
      end
      g
end
function  𝒟(f::Basic,vars::Vector{Basic},deg::Int64)
      n=length(vars);
      par=FrobeniusSolve(fill(1,n),deg);
      a=SymFunction("a");
      res=0;
      for i=1:length(par)
           res+=a(par[i]...)gd(f,vars,par[i])
      end
      res
end  
function  𝒟(f::Basic,deg::Int64)
      vars=f|>get_args;
      n=length(vars);
      par=FrobeniusSolve(fill(1,n),deg);
      a=SymFunction("a");
      res=0;
      for i=1:length(par)
           res+=a(par[i]...)gd(f,vars,par[i])
      end
      res
end  
function monoGen(vars::Vector{Basic},orders::Vector{Int64})
      [vars[i]^orders[i] for i=1:length(vars)]|>prod
end
function localDual(f::Basic,vars::Vector{Basic},order::Int64)
   n=length(vars);
   deg=degree(f);
   diffdeg=order-deg;
   comb=FrobeniusSolve(fill(1,n),diffdeg);
   v=SymFunction("v");
   [monoGen(vars,comb[i]) for i=1:length(comb)].*f
end

function homoEqn(ideal::Vector{Basic},vars::Vector{Basic},order::Int64)
      f = [localDual(ideal[i], vars, order) for i=1:length(ideal)]|>Iterators.flatten|>collect|>unique;
      [𝒟(f[i],vars,order) for i=1:length(f)]
end

function inhomoEqn(ideal::Vector{Basic},vars::Vector{Basic},order::Int64)
      mat=[diff(ideal[i],vars[j]) for i in 1:length(ideal), j in 1:length(vars)];
      𝒟(det(mat),vars,order)
end

function eqnAnsatz(ideal::Vector{Basic},vars::Vector{Basic},order::Int64)
      intersectionNumber=[degree(ideal[i]) for i=1:length(ideal)]|>prod;
      homo = homoEqn(ideal, vars, order);
      inhomo = inhomoEqn(ideal, vars, order);
      push!(homo,inhomo-intersectionNumber);
      homo
end


function solve(ideal::Vector{Basic},vars::Vector{Basic})
      eqns=ideal;
      aVar=vars;
      aMat=[coeff(eqns[i],aVar[j])|>Int for i=1:length(eqns), j=1:length(aVar)];
      aMatrank=rank(aMat);
      aM=aMat;
      notlist=[];
      while length(aM[:,1])> aMatrank
            aMrank=rank(aM);
            for i=1:length(aM[:,1])
                  Mt=aM[Not(i), :];
                  if Mt|>rank == aMrank
                        push!(notlist,i+length(notlist))
                        aM=Mt;
                        break
                  end
            end
      end
      dimaM=size(aM)[1];
      aMSym=[Basic(aM[i,j]) for i=1:dimaM, j=1:dimaM];
      inhomTerm=eqns[Not(notlist)];
      for i=1:length(aVar)
            inhomTerm=[subs(inhomTerm[j],aVar[i],0)|>Basic for j=1:length(inhomTerm)]
      end
      inhomTerm=.-inhomTerm;
      invM=inv(aMSym);
      sol=[aVar[i]=>dot(invM[i,:],inhomTerm) for i=1:dimaM];
      sol=Dict(sol...)
end

function MultiResidue(num::Basic,homoideal::Vector{Basic},vars::Vector{Basic})
      dOrder=[degree(homoideal[i]) for i=1:length(homoideal)]|>sum;
      dOrder=dOrder-length(vars);
      coeqn=eqnAnsatz(homoideal, vars, dOrder);
      par=FrobeniusSolve(fill(1,length(vars)),dOrder);
      a=SymFunction("a");
      aVar=[a(par[i]...) for i=1:length(par)]; 
      varszero=Dict([vars[i]=>0 for i=1:length(vars)]...);
      sola=solve(coeqn,aVar);
      res=𝒟(num,vars,dOrder);
      res=subs(res,varszero);
      res=subs(res,sola)
end

end
