(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["Olsson`"]


(*Olsson7 is Olsson.wl*)


(*Olsson v0.2*)
(* only roc2 is added *)
(*ROC modules are removed*)
(*serrecog.m is attached *)
(*KdFQ,FTildeQ,seriesNvariable for serrecog are added*)
(*ToNumerics and ToNumericAC are added that rewrite the prefactors with if-else statement*)
(*ToNumericsAC is fixed, different from olsson5*)
(*serrecog2var is specialized to 2 variables*)


(* use psim3 as function of index. NOTE := instead of = 
psim3[index_]:=Flatten[Table[Pochhammer[a_,index[[i]] (a1_/;IntegerQ[a1]&&a1<0)]-> (-1)^(-index[[i]] a1)/Pochhammer[1-a,- a1  index[[i]]],{i,1,Length[index]}]];
sim1[index_]:=Flatten[Table[{Gamma[y_.(x_. index[[i]]+a_)]->   Pochhammer[y a,y x index[[i]]]Gamma[y a],Pochhammer[ x_.(z_. index[[i]]+a:_:0),n___]->   Pochhammer[  x a,n+x z  index[[i]]]/Pochhammer[x a,x z  index[[i]]]},{i,1,Length[index]}]];

ADDED
 *)
 
 (*work on serrecog.m
 (-1)^m, (-1)^n type factors can be multiplied to the x and y
 
 serrecog[{m,n},((-1)^m x^m (-y)^-b2 y^-n Gamma[a-b2] Gamma[c2] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+b2-c2,n])/(m! n! Gamma[a] Gamma[-b2+c2] Pochhammer[1-a+b2,-m+n] Pochhammer[c1,m])]
 WORKS
 
 BUT 
 serrecog[{n,m},((-1)^m x^m (-y)^-b2 y^-n Gamma[a-b2] Gamma[c2] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+b2-c2,n])/(m! n! Gamma[a] Gamma[-b2+c2] Pochhammer[1-a+b2,-m+n] Pochhammer[c1,m])]
 
 DOES NOT!!!*)


Print["Olsson.wl v1.0\n","Authors : B. Ananthanarayan, Souvik Bera, S. Friot & Tanay Pathak"];



(* ::Input::Initialization:: *)
Olsson::usage="Finds analytic continuations of general hypergeometric series 
Olsson[n_Integer, summationindices_list, expression, options]
n is the sum index wrt which analytic continuations are done ";
hypgeoinf::usage="gives the expression for pFp-1 at infinity
hypgeoinf[p_Integer,q_Integer,z_]";
hypgeoatinf::usage="gives the expression for pFp-1 at infinity
hypgeoatinf[a_List,b_List, x_]";
hypgeonum::usage="gives the definition for pFp-1
hypgeonum[p_Integer,q_Integer,z_,m_]";
hypgeopoch::usage="definition of pFp-1";
f21one::usage="gives the AC at 1";
f21et1::usage="uses et1";
f21et2::usage="uses et2";
f21et3::usage="uses et3";
callroc::usage="gives the roc of hypergeometric series
callroc[indices_List, expression]
";
characteristiclist::usage="gives the characteristic list of series. It takes input as
 characteristiclist[indices_List, expression]
Example : \n
characteristiclist[{m,p},\!\(\*FractionBox[\(Pochhammer[a, m + p] Pochhammer[b, m] Pochhammer[b1, p] \*SuperscriptBox[\(x\), \(m\)]\\\ \*SuperscriptBox[\(y\), \(p\)]\), \(Pochhammer[c, m + p] \(m!\) \(p!\)\)]\)]
o/p : {{{{m+p,m,p},{m+p}},{x,y}}}
";
Topoch::usage="converts a general expression containing 
gamma fn, prefactors, poch \[Rule] poch only ";
pochordering::usage="expresses the poch in a certain ordering such that the 
first summation index has always positive coefficient";
charlistrecog::usage="given a characteristic list, the output is the characteristic list with certain ordering
the first summation index has always positive coefficient";
serrecogn::usage="given an expression, it recognises the series";
serrecog::usage="given an expression, it recognises the series";
serrecog2var::usage="given an expression, it recognises the series";
KdFQ::usage="KDF series";
FTildeQ::usage="FTildeQ series";
PochhammerSimplify::usage="simplifies the Pochhammers
PochhammerSimplify[index_List,exp_,actions_List]
with actions : {Dimidiation,GammatoPochhammer,PositivePochhammer,PochhammertoGamma}";
PochSimplify2::usage="simplifies the Pochhammers";
pocsim::usage="simplifies the Pochhammers";
Pochdim::usage="apply Pochhammer dimidiation formulae
Examples : Pochhammer[a,2m]/.Pochdim \n o/p : \!\(\*SuperscriptBox[\(2\), \(2\\\ m\)]\) Pochhammer[\!\(\*FractionBox[\(1\), \(2\)]\)+\!\(\*FractionBox[\(a\), \(2\)]\),m] Pochhammer[\!\(\*FractionBox[\(a\), \(2\)]\),m]
Simplify//@(Pochhammer[a,3n+9m])/.Pochdim \n o/p : \!\(\*SuperscriptBox[\(3\), \(3\\\ \((3\\\ m + n)\)\)]\) Pochhammer[\!\(\*FractionBox[\(1\), \(3\)]\)+\!\(\*FractionBox[\(a\), \(3\)]\),3 m+n] Pochhammer[\!\(\*FractionBox[\(2\), \(3\)]\)+\!\(\*FractionBox[\(a\), \(3\)]\),3 m+n] Pochhammer[\!\(\*FractionBox[\(a\), \(3\)]\),3 m+n]";
gammatoPoch::usage="converts gamma functions to Pochhammer
Examples : Gamma[a+n]/.gammatoPoch[{n}][[1]] \n o/p : Gamma[a] Pochhammer[a,n]
Pochhammer[a+n,m]/.gammatoPoch[{n}][[2]] \n o/p : \!\(\*FractionBox[\(Pochhammer[a, m + n]\), \(Pochhammer[a, n]\)]\)
gammatoPoch contains both type of relations. It takes the list of summation indices as argument\n gammatoPoch[{m,n,p}]";
unsim::usage="converts Pochhammer to gamma functions by using Pochhammer[a,m+n]-> Pochhammer[a,n] Pochhammer[a+n,m]
Example : Pochhammer[a,2m-n]/.unsim[m]\n o/p : Pochhammer[a,-n] Pochhammer[a-n,2 m]";
positivepoch::usage="use the formula Pochhammer[a,-n]-> \!\(\*FractionBox[SuperscriptBox[\((\(-1\))\), \(n\)], \(Pochhammer[1 - a, n]\)]\)
Examples :  Pochhammer[a,-m]/.positivepoch[{m}]\n o/p : \!\(\*FractionBox[SuperscriptBox[\((\(-1\))\), \(m\)], \(Pochhammer[1 - a, m]\)]\)
Pochhammer[a,-2n-m]/.positivepoch[{n,m}]\n o/p : \!\(\*FractionBox[SuperscriptBox[\((\(-1\))\), \(m + 2\\\ n\)], \(Pochhammer[1 - a, m + 2\\\ n]\)]\)";
ToNumerics::usage="Rewrites the prefactor of a series with if-else form
ToNumerics[var_List,series_]";
ToNumericsAC::usage="Rewrites the prefactor of a AC with if-else form
ToNumerics[var_List,series_]";
condition::usage="finds the condition for a prefactor
condition[exp]";


(* ::Input::Initialization:: *)
Begin["`Private`"]


(*Get["series_list.m"];*)Get["ROC2.wl"];


(*Needs["serieslist`"];*)Needs["ROC2`"];


condition[exp_]:=Module[{var,num,denom,num1,denom1,complex,ii},
var=Variables[exp];
(*Print[var];*)
complex=ComplexExpand[exp,var]//Simplify;
{num,denom}={Numerator[complex],Denominator[complex]};
(*Print[{num,denom}];*)
num1=num* Conjugate[denom]//Simplify;
(*Print[num1];*)
denom1=denom*Conjugate[denom]//Simplify;
(*Print[{num1//Expand}];
Print[Coefficient[Expand[num1],ii]];*)
Return[Coefficient[Expand[num1/.I-> ii/.-I-> -ii],ii]]]


ToNumerics[var_List,series_]:=Module[{prefac1,prefac2,prefactorlist,newprefactors,serieslist,ifelse,prefactor,e},

prefactor=Expand[Simplify[series]]/.Table[var[[i]]-> 0,{i,1,Length[var]}];

serieslist={prefactor,Expand[Simplify[series/prefactor]]};
(*Print["serieslist",serieslist];*)
prefac1=serieslist[[1]];
(*Print["prefac1",prefac1];*)
prefac2=prefac1/(prefac1/.Gamma[z_]-> 1);
(*Print["{prefac2,serieslist[[2]]}",{prefac2,serieslist[[2]]}];*)
prefactorlist=If[Head[#]=!=List,{#,1},#]&/@(If[Head[#]===Times,List@@#,{#}]&@(prefac1/.Gamma[z_]-> 1)/.Power[x_,a_]-> List[x,a]);
(*Print["new",If[Head[#]=!=List,{#,1},#]&/@(If[Head[#]===Times,List@@#,{#}]&@(prefac1/.Gamma[z_]-> 1)/.Power[x_.,a_.]-> List[x,a])];*)
If[prefactorlist===1,newprefactors={1};Goto[end]];
If[Head[prefactorlist[[1]]]=!=List,prefactorlist={prefactorlist}];
(*Print["prefactorlist",prefactorlist];*)
ifelse=Table[{(Plus@@(Cases[If[Head[#]=== Plus,List@@#,{#}]&@condition[#[[1]]]/.Im[z_]-> -e,e x_.]/.e-> 1))>0,Power[#[[1]],#[[2]]],Power[1/#[[1]],-#[[2]]]}&@prefactorlist[[i]],{i,1,Length[prefactorlist]}];
(*Print["ifelse", ifelse];*)
newprefactors=Table[If@@ifelse[[i]],{i,1,Length[ifelse]}];
(*Print["newprefactors", (newprefactors)];*)
Label[end];
Return[Times@@newprefactors*prefac2*serieslist[[2]]]]


ToNumericsAC[var_List,series_]:=Module[{},Return[Plus@@(ToNumerics[var,#]&/@If[Head[#]===Plus,List@@#,{#}]&@Expand[Simplify[series]])]]


hypgeonum[p0_Integer,q0_Integer,z0_,m0_]:=Module[{p,q,list1,list2,term},If[p0=!= q0+1,Return["Not applicable",Module]];list1 = Table[a[i]=Subscript[Global`a, i],{i,1,p0}];list2 = Table[b[i]=Subscript[Global`b, i],{i,1,q0}];
term=( Times@@Table[a[i]=Pochhammer[list1[[i]],m0],{i,1,p0}] z0^m0)/(Times@@Table[b[i]=Pochhammer[list2[[i]],m0],{i,1,q0}] m0!);
Return[term]];


hypgeoinf[p0_Integer,q0_Integer,z0_]:=Module[{p,q,list1,list2,term},
If[p0=!= q0+1,Return["Not applicable",Module]];list1 = Table[a[i]=Subscript[Global`a, i],{i,1,p0}];list2 = Table[b[i]=Subscript[Global`b, i],{i,1,q0}];
term = (((-z0)^-list1[[1]] Times@@Gamma[list2])/(Times@@Gamma[list1]))*((Times[Gamma[list1[[1]]],Times@@(Table[a[i]=Gamma[list1[[i]]-list1[[1]]],{i,2,p0}])])/(Times@@(Table[a[i]=Gamma[list2[[i]]-list1[[1]]],{i,1,q0}])))*HypergeometricPFQ[Join[{list1[[1]]},Table[a[i]=1+list1[[1]]-list2[[i]],{i,1,q0}]],Table[a[i]=1+list1[[1]]-list1[[i]],{i,2,p0}],1/z0];
(*Print[term];*)
Return[Plus@@Table[a[i]=term/.{list1[[1]]-> list1[[i]],list1[[i]]-> list1[[1]]},{i,1,p0}]]];


pocsim[a0___,(n0:_Integer :1) m0_]:=Module[{a},
a={};
If[n0>0,For[i=0,i<= n0-1,i++,{a=Append[a,Pochhammer[a0/n0+i/n0,m0]]}];];
Return[If[n0>0,n0^(n0 m0) Times@@a,Pochhammer[a0,n0 m0]]]];


hypgeoatinf[a0_List,b0_List, x_] :=Module[{},Return[hypgeoinf[Length[a0], Length[b0], x]/. 
   Table[a1[i] = Subscript[Global`a, i] -> a0[[i]], {i, 1, Length[a0]}] /. 
  Table[b1[i] = Subscript[Global`b, i] -> b0[[i]], {i, 1, Length[b0]}]]];


hypgeopoch[a0_List,b0_List, x_,m_] :=Module[{},Return[hypgeonum[Length[a0], Length[b0], x,m]/. 
   Table[a1[i] = Subscript[Global`a, i] -> a0[[i]], {i, 1, Length[a0]}] /. 
  Table[b1[i] = Subscript[Global`b, i] -> b0[[i]], {i, 1, Length[b0]}]]];


(* ::Input::Initialization:: *)
f21one[{a_,b_},{c_},x_]:= (Gamma[c]Gamma[c-b-a])/(Gamma[c-b]Gamma[c-a]) HypergeometricPFQ[{a,b},{a+b-c+1},1-x]+ (Gamma[c]Gamma[a+b-c])/(Gamma[a]Gamma[b]) (1-x)^(c-a-b) HypergeometricPFQ[{c-a,c-b},{c-b-a+1},1-x];


(* ::Input::Initialization:: *)
f21et1[{a_,b_},{c_},x_]:=(1-x)^-a HypergeometricPFQ[{a,c-b},{c},x/(x-1)];


(* ::Input::Initialization:: *)
f21et2[{a_,b_},{c_},x_]:=(1-x)^-b HypergeometricPFQ[{c-a,b},{c},x/(x-1)];


(* ::Input::Initialization:: *)
f21et3[{a_,b_},{c_},x_]:=(1-x)^(c-a-b) HypergeometricPFQ[{c-a,c-b},{c},x];


(* ::Input::Initialization:: *)
(*sim1 = {Gamma[y_.(x_.m+a_)]\[RuleDelayed]   Pochhammer[y a,y x m]Gamma[y a],Gamma[y_.(x_.n+a_)]\[RuleDelayed]   Pochhammer[y a,y x n]Gamma[y a],Pochhammer[ x_.(z_.n+a_),m___]\[RuleDelayed]  Pochhammer[  x a,m+x z n]/Pochhammer[x a,x z n],Pochhammer[ x_.(z_.m+a_),n___]\[RuleDelayed]  Pochhammer[  x a,n+x z m]/Pochhammer[x a,x z m]}*)


(* ::Input::Initialization:: *)
(*psim3 =  {Pochhammer[a_,m (a1_/;IntegerQ[a1]&&a1<0)]\[RuleDelayed](-1)^(-m a1)/Pochhammer[1-a,- a1  m],Pochhammer[a_,(a1_/;IntegerQ[a1]&&a1<0)n]\[RuleDelayed](-1)^(- n a1 )/Pochhammer[1-a,- a1  n]}*)


gammatoPoch[index_]:=Flatten[Table[{Gamma[y_.(x_. index[[i]]+a_)]->   Pochhammer[y a,y x index[[i]]]Gamma[y a],Pochhammer[ x_.(z_. index[[i]]+a:_:0),n___]->   Pochhammer[  x a,n+x z  index[[i]]]/Pochhammer[x a,x z  index[[i]]]},{i,1,Length[index]}]];


unsim[m0_] :=  {Pochhammer[a_,x_. m0+y_. ]-> Pochhammer[y +a,x m0] Pochhammer[a,y ]};


positivepoch[index_]:= Flatten[Prepend[{positivepoch2[index]},Flatten[Table[Pochhammer[a_,index[[i]] (a1_/;IntegerQ[a1]&&a1<0)]->   (-1)^(-index[[i]] a1)/Pochhammer[1-a,- a1  index[[i]]],{i,1,Length[index]}]]]];


positivepoch2[index_]:=Pochhammer[a_,Plus@@(index*(\!\(\*
TagBox[
StyleBox[
RowBox[{"PatternTest", "[", 
RowBox[{
RowBox[{"Pattern", "[", 
RowBox[{"#", ",", 
RowBox[{"Blank", "[", "]"}]}], "]"}], ",", 
RowBox[{"Function", "[", 
RowBox[{"And", "[", 
RowBox[{
RowBox[{"IntegerQ", "[", "#", "]"}], ",", 
RowBox[{"Less", "[", 
RowBox[{"#", ",", "0"}], "]"}]}], "]"}], "]"}]}], "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\)&/@(ToExpression[#]&/@("a"<>ToString[#]&/@Range[Length[index]]))))]->  Power[-1,-Plus@@(index*(ToExpression[#1]&/@("a"<>ToString[#]&/@Range[Length[index]])))]/Pochhammer[1-a,-Plus@@(index*(ToExpression[#1]&/@("a"<>ToString[#]&/@Range[Length[index]])))]


(* ::Input::Initialization:: *)
Options[Olsson]={ PET1-> False,PET2-> False,PET3-> False,one-> False,inf -> False,sim-> False,roc-> False,(*hint-> False,*) sum-> False};
Olsson[args___] := Olssoncore[Olsson,args];
Olssoncore[Olsson,p_,index_List,poch0_,OptionsPattern[Olsson]]:=(*Olsson[p_,index_List,poch0_,actions_]:=*)Module[{poch,poch1,poch2,list,list1,list11,x,usim1,pocdec,pocdec2,result,rocresult,sim1,psim3,psim4,m0},
m0=index[[p]];


(*some useful commands*)


pocdec2=Pochhammer[a_,Plus@@(index*Table[\!\(\*
TagBox[
StyleBox[
RowBox[{"PatternTest", "[", 
RowBox[{
RowBox[{"Pattern", "[", 
RowBox[{"a1", ",", 
RowBox[{"Blank", "[", "]"}]}], "]"}], ",", 
RowBox[{"Function", "[", 
RowBox[{"And", "[", 
RowBox[{
RowBox[{"IntegerQ", "[", "#", "]"}], ",", 
RowBox[{"Greater", "[", 
RowBox[{"#", ",", "1"}], "]"}]}], "]"}], "]"}]}], "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\),{i,1,Length[index]}])]:>pocsim[a,a1 Plus@@index];
usim1=unsim[m0];
pocdec={Pochhammer[a___,m___]:>pocsim[a,m]};
sim1=gammatoPoch[index];
psim3=positivepoch[index];
psim4=positivepoch2[index];




poch1=Simplify//@((Expand//@(poch0/.HypergeometricPFQ[a_,b_,z_]->hypgeopoch[a,b,z,index[[p]]] ))//.usim1//.psim3//.pocdec//.pocdec2);
(*Print[poch1];*)

poch=({List@@Numerator[#],List@@Denominator[#]})&@(poch1/(poch1/.m0-> 0));
(*Print["poch",poch];*)
list=(Cases[#,Pochhammer[___,m0]]&/@poch)/.Pochhammer[z_,m0]-> z;
(*Print["list",list];*)
x=DeleteCases[((poch/.Pochhammer[z_,m_]-> 1)/(poch/.Pochhammer[z_,m_]-> 1/.m0-> 0))/.m0-> 1,1,2];
(*Print[x];*)
(*If[Length[index]===1,{list,x}=Flatten[Select[If[Head[#]===Times,List@@#,{#}]&@poch0/.Hypergeometric2F1[a_,b_,c_,z_]\[Rule] HypergeometricPFQ[{a,b},{c},z],Head[#]===HypergeometricPFQ&]/.HypergeometricPFQ[{a_,b_},{c_},z_]\[Rule]{ {{a,b},{c}},{{z},{}}},1]];*)
(*Print[{list,x}];
*)
(*Print[Simplify//@(poch0/.m0-> 0)];
Print[hypgeoatinf[list[[1]],list[[2]],(Times@@x[[1]])/(Times@@x[[2]])]];*)



result=poch0;

(*Print[list];*)
(*ETs*)
Which[OptionValue[PET1]===True,
(*If[Length[list[[1]]]=!=2||list[[2]]=!=1,Print[Length[list[[1]]],"F" ,Length[list[[2]]]," appeared"];Abort[];];*)result=Distribute[Simplify//@(poch1/.m0-> 0)*Simplify//@f21et1[list[[1]],list[[2]],(Times@@x[[1]])/(Times@@x[[2]])]];,
OptionValue[PET2]===True,
(*If[Length[list[[1]]]=!=2||list[[2]]=!=1,Print[Length[list[[1]]],"F" ,Length[list[[2]]]," appeared"];Abort[];];*)result=Distribute[Simplify//@(poch1/.m0-> 0)*Simplify//@f21et2[list[[1]],list[[2]],(Times@@x[[1]])/(Times@@x[[2]])]];,
OptionValue[PET3]===True,
(*If[Length[list[[1]]]=!=2||list[[2]]=!=1,Print[Length[list[[1]]],"F" ,Length[list[[2]]]," appeared"];Abort[];];*)result=Distribute[Simplify//@(poch1/.m0-> 0)*Simplify//@f21et3[list[[1]],list[[2]],(Times@@x[[1]])/(Times@@x[[2]])]];];





Which[
(*AC at 1*)
OptionValue[one]===True,
(*If[Length[list[[1]]]=!=2||list[[2]]=!=1,Print[Length[list[[1]]],"F" ,Length[list[[2]]]," appeared"];Abort[];];*)result=Distribute[Simplify//@(poch1/.m0-> 0)*Simplify//@f21one[list[[1]],list[[2]],(Times@@x[[1]])/(Times@@x[[2]])]];,

(*AC at inf*)
OptionValue[inf]===True,result=Distribute[Simplify//@(poch1/.m0-> 0)*Simplify//@hypgeoatinf[list[[1]],list[[2]],(Times@@x[[1]])/(Times@@x[[2]])]];,

(*True*)
OptionValue[sum]===True,
	result=Simplify//@(poch1/.m0-> 0)*HypergeometricPFQ[list[[1]],list[[2]],(Times@@x[[1]])/(Times@@x[[2]])];


];




(*~~~~~~~~~~~~~~~~~~~~~~~~simplification~~~~~~~~~~~~~~~~~~~~~~~*)
If[OptionValue[sim]===True,result=If[Head[#]=!=Plus,{#},List@@#]&@Expand[result];result=((Expand//@(result))/.HypergeometricPFQ[a1_List,b1_List,z_]:> hypgeopoch[a1,b1,z,m0])(*//.sim1//.psim3//.psim4*);
result=PochhammerSimplify[index,#,{"GammatoPochhammer","PositivePochhammer","Dimidiation"}]&/@result;result=Refine[#1,Element[index,Integers]&&#>0&/@index]&/@result;result=Plus@@result;];

(*Print["after sim",result];*)

(*~~~~~~~~~~~~~~~~~~~~~~~ roc2  ~~~~~~~~~~~~~~~~~~~~~~~~~~~*)

If[OptionValue[roc]===True,

rocresult=callroc[index,result];

];




Return[If[OptionValue[roc]===True,{And@@Flatten[rocresult],Simplify//@result},Simplify//@result]]];
(*Return[Simplify//@result]];*)


(*no need of pocdec1*)
Pochdim1[index_]:=Pochhammer[a_,Plus@@(index*Table[\!\(\*
TagBox[
StyleBox[
RowBox[{"PatternTest", "[", 
RowBox[{
RowBox[{"Pattern", "[", 
RowBox[{"a1", ",", 
RowBox[{"Blank", "[", "]"}]}], "]"}], ",", 
RowBox[{"Function", "[", 
RowBox[{"And", "[", 
RowBox[{
RowBox[{"IntegerQ", "[", "#", "]"}], ",", 
RowBox[{"Greater", "[", 
RowBox[{"#", ",", "1"}], "]"}]}], "]"}], "]"}]}], "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\),{i,1,Length[index]}])]:>  pocsim[a,a1 Plus@@index]


(*use Simplify//@exp and use Pochdim*)
Pochdim:=Pochhammer[a___,z___]:> pocsim[a,z];


Options[PochSimplify2]={ gatoPoch-> False,dimidiation -> False,positivepoch-> False,Pochtoga-> False};
PochSimplify2[args___] := PochSimplifycore[PochSimplify,args];
PochSimplifycore[PochSimplify2,p_,index_,exp_,OptionsPattern[PochSimplify2]]:=Module[{pocdec,sim1,psim3,psim4,pocdec2,usim1,exp1,optionlength},
Print["gatoPoch-> False,dimidiation -> False,positivepoch-> False,Pochtoga-> False"];
pocdec2=Pochhammer[a_,Plus@@(index*Table[\!\(\*
TagBox[
StyleBox[
RowBox[{"PatternTest", "[", 
RowBox[{
RowBox[{"Pattern", "[", 
RowBox[{"a1", ",", 
RowBox[{"Blank", "[", "]"}]}], "]"}], ",", 
RowBox[{"Function", "[", 
RowBox[{"And", "[", 
RowBox[{
RowBox[{"IntegerQ", "[", "#", "]"}], ",", 
RowBox[{"Greater", "[", 
RowBox[{"#", ",", "1"}], "]"}]}], "]"}], "]"}]}], "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\),{i,1,Length[index]}])]:>  pocsim[a,a1 Plus@@index];
pocdec={Pochhammer[a___,m___]:>  pocsim[a,m]};
sim1=gammatoPoch[index];
psim3=positivepoch[index];
psim4=positivepoch2[index];
usim1=unsim[index[[p]]];
exp1=Expand//@exp;
optionlength=Length[DeleteCases[{OptionValue[Pochtoga],OptionValue[dimidiation],OptionValue[positivepoch],OptionValue[gatoPoch]},False]];
Print[optionlength];
Do[
Which[OptionValue[Pochtoga],exp1=Expand//@(exp1//.usim1);Continue[];,
OptionValue[dimidiation], exp1=Expand//@(exp1//.pocdec//.pocdec2)Continue[];,
OptionValue[positivepoch],exp1=Expand//@(exp1//.psim3//.psim4)Continue[];,
OptionValue[gatoPoch],exp1=Expand//@(exp1//.sim1);Continue[];];,optionlength];
Return[Simplify//@exp1]]


PochhammerSimplify[(*p_,*)index_List,exp_,actions_List]:=Module[{pocdec,sim1,psim3,psim4,pocdec2,usim1,exp1},
(*Print["actions  :   ","{Dimidiation,GammatoPochhammer,PositivePochhammer,PochhammertoGamma}"];*)
pocdec2=Pochhammer[a_,Plus@@(index*Table[\!\(\*
TagBox[
StyleBox[
RowBox[{"PatternTest", "[", 
RowBox[{
RowBox[{"Pattern", "[", 
RowBox[{"a1", ",", 
RowBox[{"Blank", "[", "]"}]}], "]"}], ",", 
RowBox[{"Function", "[", 
RowBox[{"And", "[", 
RowBox[{
RowBox[{"IntegerQ", "[", "#", "]"}], ",", 
RowBox[{"Greater", "[", 
RowBox[{"#", ",", "1"}], "]"}]}], "]"}], "]"}]}], "]"}],
ShowSpecialCharacters->False,
ShowStringCharacters->True,
NumberMarks->True],
FullForm]\),{i,1,Length[index]}])]:>   pocsim[a,a1 Plus@@index];
pocdec={Pochhammer[a___,m___]:>    pocsim[a,m]};
sim1=gammatoPoch[index];
psim3=positivepoch[index];
psim4=positivepoch2[index];
(*usim1=gammatoPoch[{index[[p]]}];
Print[usim1];*)
(*usim1=unsim[index[[p]]];*)
(*Print[{pocdec2,pocdec,sim1,psim3,psim4,usim1}];*)
exp1=Expand//@(exp/.Pochhammer[a_,m_]-> Gamma[a+m]/Gamma[a]);
Table[Which[ToString[actions[[i]]]==="Dimidiation", exp1=Expand//@(Expand//@(Simplify//@(exp1)//.pocdec)//.pocdec2);(*Print[exp1];*),
ToString[actions[[i]]]==="PositivePochhammer",exp1=Expand//@(exp1//.psim3//.psim4);(*Print[exp1];*),
(*ToString[actions[[i]]]==="PochhammertoGamma",exp1=Expand//@(exp1//.usim1);Print[exp1];,*)
ToString[actions[[i]]]==="GammatoPochhammer",exp1=Expand//@(exp1//.sim1);(*Print[exp1];*)],{i,1,Length[actions]}];
(*Table[If[Length[actions[[i]]]===1,Which[ToString[actions[[i]]]==="Dimidiation", exp1=Expand//@(exp1//.pocdec//.pocdec2);,
ToString[actions[[i]]]==="PositivePochhammer",exp1=Expand//@(exp1//.psim3//.psim4);],usim1=unsim[index[[actions[[i]][[1]]]]];
Which[ToString[actions[[i]][[2]]]==="PochhammertoGamma",exp1=Expand//@(exp1//.usim1);,
ToString[actions[[i]][[2]]]==="GammatoPochhammer",exp1=Expand//@(exp1//.sim1);]],{i,1,Length[actions]}];*)
Return[Simplify//@exp1]]


(*     For Series Recognition     *)


characteristiclist[index_List,ser_]:=Module[{poch2,list1,list11},
poch2=(#/(#/.Table[index[[i]]-> 0,{i,1,Length[index]}]))&/@List@@If[Head[ser]=!=Plus,Flatten[{ser}],Expand[ser]];



(*Print[poch2];*)
list11=Table[{{Flatten[List@@If[Head[Numerator[#]]=!=Times,{Numerator[#]},Numerator[#]]/.Power[x_,n_]->Table[x,{i,1,2}]/.Pochhammer[z_,m_]-> m],Flatten[List@@If[Head[Denominator[#]]=!=Times,{Denominator[#]},Denominator[#]]/.Power[x_,n_]->Table[x,{i,1,2}]/.Pochhammer[z_,m_]-> m]}&@(poch2[[i]]/(poch2[[i]]/.Pochhammer[z_,m_]-> 1)),Table[((#/.Pochhammer[z_,m_]-> 1/.Table[If[k===j,index[[k]]-> 1,index[[k]]-> 0],{k,1,Length[index]}]))&@poch2[[i]],{j,1,Length[index]}]},{i,1,Length[poch2]}];
list1=Simplify//@DeleteDuplicates[list11];
(*Print[list1];*)
(*Print[Simplify//@list11];*)
Return[Simplify//@list11];]


Topoch[exp_]:=Module[
(* converts a general expression containing 
gamma fn, prefactors, poch \[Rule] poch only*)


{num,denom,pochnum,pochdenom},
num=Numerator[exp];denom=Denominator[exp];
(*Print[{num,denom}];*)

pochnum=Select[If[Head[#]===Times,List@@#,{#}]&@num,Head[#]===Pochhammer&];
pochdenom=Select[If[Head[#]===Times,List@@#,{#}]&@denom,Head[#]===Pochhammer&];

Return[{pochnum,pochdenom}]
]


pochordering[list_List]:=Module[{num,num1,denom,denom1,index},
(* expresses the poch in a certain ordering such that the 
first summation index has always positive coefficient*)


{num,denom}=list;
(*Print[Sort[Variables[Flatten[list][[All,2]]]]];*)
index=First@Sort[Variables[Flatten[list][[All,2]]]];

num1=If[Negative[Coefficient[Last@#,index]],(-1)^-#[[2]]/Pochhammer[1-#[[1]],-  #[[2]]],#]&/@num;
denom1=If[Negative[Coefficient[Last@#,index]],(-1)^-#[[2]]/Pochhammer[1-#[[1]],-  #[[2]]],#]&/@denom;
Return[Times@@num1/Times@@denom1]
]


charlistrecog[{num_,denom_}]:=Module[{a,n,x,var,topoch,num1,denom1,charlist},


(* given a characteristic list, the output is the characteristic list with certain ordering
the first summation index has always positive coefficient *)



var=Sort[Variables[{num,denom}]];
topoch={Table[Pochhammer[Subscript[a, i],num[[i]]],{i,1,Length[num]}] ,Table[Pochhammer[Subscript[b, i],denom[[i]]],{i,1,Length[denom]}]};
{num1,denom1}=(If[Head[#]===Pochhammer,List@#1,List@@#1]&/@{Numerator[#],Denominator[#]}//.Pochhammer-> List)&@pochordering[topoch];
charlist=Sort/@(If[Head[#]===List,#[[2]],#]&/@#1&/@{num1,denom1});
Return[charlist]]


serrecogn[exp_]:=Module[{prefactor,charlist,topoch,indices,num,denom,table,pos,var},
(* given an expression, it recognises the series *)

(*Get["C:\\ProgramData\\Mathematica\\Applications\\series_list.m"];*)
topoch=Topoch[exp];
(*Print[topoch];*)


indices=Sort[Variables[Flatten[topoch][[All,2]]]];
(*Print[indices];*)
prefactor=exp/.Pochhammer[z_,a_]-> 1/.Table[indices[[i]]-> 0,{i,1,Length[indices]}];
var=Table[((exp/.Pochhammer[z_,a_]-> 1)/(exp/.Pochhammer[z_,a_]-> 1/.indices[[i]]-> 0))/.indices[[i]]-> 1,{i,1,Length[indices]}];



{num,denom}=(If[Head[#]===Pochhammer,List@#1,List@@#1]&/@{Numerator[#],Denominator[#]}//.Pochhammer-> List)&@pochordering[topoch];
(*Print[{num,denom}];*)
(*
charlist=({List@@Numerator[#],List@@Denominator[#]}//.Pochhammer[z_,a_]\[Rule] a)&@pochordering[topoch];;*)
(*charlist=#[[All,2]]&/@{num,denom};*)
charlist=Sort/@(If[Head[#]===List,#[[2]],#]&/@#1&/@{num,denom});

Print[charlist];
(*charlist=Sort/@charlist;*)

(*Print["charlist    ",charlist];*)

(*table2={{{indices\[LeftDoubleBracket]1\[RightDoubleBracket]-indices\[LeftDoubleBracket]2\[RightDoubleBracket],indices\[LeftDoubleBracket]1\[RightDoubleBracket]+indices\[LeftDoubleBracket]2\[RightDoubleBracket],indices\[LeftDoubleBracket]2\[RightDoubleBracket]},{indices\[LeftDoubleBracket]1\[RightDoubleBracket]}} (*H1*),
{{(-1)^indices\[LeftDoubleBracket]1\[RightDoubleBracket],2 indices\[LeftDoubleBracket]1\[RightDoubleBracket]+indices\[LeftDoubleBracket]2\[RightDoubleBracket]},{(-1)^indices\[LeftDoubleBracket]2\[RightDoubleBracket],indices\[LeftDoubleBracket]1\[RightDoubleBracket]-indices\[LeftDoubleBracket]2\[RightDoubleBracket],indices\[LeftDoubleBracket]2\[RightDoubleBracket]}}(*H5*)};
table2a={H1,H5};*)

(*table2=table2var[indices];
Print[table2];*)

Which[Length[indices]===2,
table=table2var[indices];,
Length[indices]===3,table=table3var[indices];;
];


If[MemberQ[table[[All,1]],charlist],pos=Position[table[[All,1]],charlist];,Return[Print["Unknown series! Ordered characteristic list is : "];charlist]];


Print["Final result :    ",{"series name","prefactor", "expression for series"}];
Return[(*Which[Length[indices]===2,*)(*table2[[First@Flatten@pos,2]]*) {table[[First@Flatten@pos,2]],prefactor  , Times@@(If[Head[#]===List,#/.List-> Pochhammer,#]&/@num)/Times@@(If[Head[#]===List,#/.List-> Pochhammer,#]&/@denom) Times@@Table[var[[i]]^indices[[i]],{i,1,Length[indices]}]/Times@@Table[indices[[i]]!,{i,1,Length[indices]}]}]]


serrecog[ind_,exp_]:=Module[{indices,prefactor,charlist,charlist2,permutationscharlist,topoch,num,num1,denom1,denom,table,pos,var,name},
(* given an expression, it recognises the series *)




(*If[And[Length[indices]===2,MemberQ[exp1,Pochhammer[x_.,indices[[2]]-indices[[1]]]]],exp=exp1/.{indices[[1]]-> indices[[2]],indices[[2]]-> indices[[1]]},exp=exp1;];
*)
indices=Sort[ind];
(*indices=Sort[Variables[Flatten[topoch][[All,2]]]];*)
(*Print[indices];*)
prefactor=exp/.Table[indices[[i]]-> 0,{i,1,Length[indices]}];
var=Table[((exp/.Pochhammer[z_,a_]-> 1)/(exp/.Pochhammer[z_,a_]-> 1/.indices[[i]]-> 0))/.indices[[i]]-> 1,{i,1,Length[indices]}];
(*Print[var];*)
(*Print[prefactor];*)
topoch=Topoch[exp/prefactor];
(*Print["topoch   ",topoch];*)

{num,denom}=(If[Head[#]===Pochhammer,List@#1,List@@#1]&/@{Numerator[#],Denominator[#]}//.Pochhammer-> List)&@pochordering[topoch];

(*If[And[Length[indices]===2,Or[MemberQ[num,indices[[2]]-indices[[1]],2],MemberQ[denom,indices[[2]]-indices[[1]],2]]],{num,denom}={num,denom}/.{indices[[1]]\[Rule] indices[[2]],indices[[2]]\[Rule] indices[[1]]}];
*)
(*Print["num-denom   ",{num,denom}];*)
num1=Select[num,ContainsAny[indices,Variables[#[[2]]]]&];
denom1=Select[denom,ContainsAny[indices,Variables[#[[2]]]]&];
(*Print[{num1,denom1}];*)
(*
charlist=({List@@Numerator[#],List@@Denominator[#]}//.Pochhammer[z_,a_]\[Rule] a)&@pochordering[topoch];;*)
(*charlist=#[[All,2]]&/@{num,denom};*)
charlist=DeleteCases[Sort/@(If[Head[#]===List,#[[2]],#]&/@#1&/@{num1,denom1}),Power[(-1),a_.],2];


(*here is the new part*)
(*permutationscharlist=charlist/.((#/.List-> Rule)&/@#&/@(Transpose[#]&/@({indices,#}&/@Permutations[indices])));
Print[permutationscharlist];
*)
(*charlist=Sort/@charlist;*)

(*Print["charlist    ",charlist];*)


(* CHANGE HERE IF REQUIRED *)
Which[Length[indices]===2,
table=table2var[indices];(*Print[table[[All,1]]];*),
Length[indices]===3,
table=table3var[indices];,
Length[indices]> 3,
table=tableNvar[indices];
];
(*charlist2=Flatten[Intersection[table[[All,1]],(*permutationscharlist*)charlist],1];*)
(*Print["charlist2  :", charlist2];*)
(*Print[MemberQ[table[[All,1]],(*permutationscharlist*)charlist]];*)
If[MemberQ[table[[All,1]],(*permutationscharlist*)charlist],pos=Position[table[[All,1]],(*charlist2*)charlist];(*Print[pos];*)Goto[end];,pos={};];
name={};
If[Length[indices]===2,name=KdFQ[indices,charlist];If[name==={},name=(*If[#=!={},#[[1]],#]*)#&@(*Flatten[*)FTildeQ[indices,#]&@(*permutationscharlist*)charlist];(*];*)];
If[Length[indices]>=3,name=KdFQ[indices,charlist];];
(*Print[{pos,name}];*)
If[name==={}&&pos==={},Return[Print["Unknown series! \n ",charlist];{prefactor  , exp/prefactor}]];


Label[end];
(*Print["Final result :    ",{"series name","prefactor", "expression for series"}];*)
Return[ {If[pos==={},name,table[[First@Flatten@pos,2]]],prefactor  , exp/prefactor(*Times@@(If[Head[#]===List,#/.List-> Pochhammer,#]&/@num)/Times@@(If[Head[#]===List,#/.List-> Pochhammer,#]&/@denom) Times@@Table[var[[i]]^indices[[i]],{i,1,Length[indices]}]/Times@@Table[indices[[i]]!,{i,1,Length[indices]}]*)}]]


serrecog2var[ind_,exp0_]/;Length[ind]===2:=Module[{exp,var0,indices,list,num,denom,var,a,a1,a2,b,b1,b2,c,c1,c2,d,e},
indices=Sort[ind];
list=serrecog[indices,exp0];
If[Length[list]===2,Return[Times@@list]];
var0=Table[(#/.Pochhammer[z_,a_]-> 1)/.Table[If[i===j,indices[[i]]-> 1,indices[[i]]-> 0],{i,1,Length[indices]}],{j,1,Length[indices]}]&@Last[list];
exp=Times@@Power[var0,indices]*pochordering[Topoch[Last[list]]];
(*Print[exp];*)
(*{num,denom}=SortBy[#,Last]&/@((If[Head[#]===Pochhammer,List@#1,List@@#1]&/@{Numerator[#],Denominator[#]}//.Pochhammer-> List)&@pochordering[Topoch[Last[list]]]);*)
{num,denom}=SortBy[#,Last]&/@(Topoch[exp]/.Pochhammer[z_,a_]-> List[z,a]);
var=Table[(#/.Pochhammer[z_,a_]-> 1)/.Table[If[i===j,indices[[i]]-> 1,indices[[i]]-> 0],{i,1,Length[indices]}],{j,1,Length[indices]}]&@exp;

(*Print[{num,denom,var}];*)

Which[(*F1*)
First[list]==="F1",{b1,b2,a}=num[[All,1]];{c}=denom[[All,1]];Return[list[[2]]*Global`F1@@{a,b1,b2,c,var[[1]],var[[2]]}];,

(*F2*)
First[list]==="F2",{b1,b2,a}=num[[All,1]];{c1,c2}=denom[[All,1]];Return[list[[2]]*Global`F2@@{a,b1,b2,c1,c2,var[[1]],var[[2]]}];,

(*F3*)
First[list]==="F3",{a1,b1,a2,b2}=num[[All,1]];{c}=denom[[All,1]];Return[list[[2]]*Global`F3@@{a1,a2,b1,b2,c,var[[1]],var[[2]]}];,

(*F4*)
First[list]==="F4",{a,b}=num[[All,1]];{c1,c2}=denom[[All,1]];Return[list[[2]]*Global`F4@@{a,b,c1,c2,var[[1]],var[[2]]}];,

(*G1*)
First[list]==="G1",{b2,a}=num[[All,1]];{b1}=denom[[All,1]];Return[list[[2]]*Global`G1@@{a,1-b1,b2,-var[[2]],-var[[1]]}];,

(*G2*)
First[list]==="G2",{a,b2,a1}=num[[All,1]];{b1}=denom[[All,1]];Return[list[[2]]*Global`G2@@{a,a1,1-b1,b2,-var[[1]],-var[[2]]}];,

(*G3*)
First[list]==="G3",{a1}=num[[All,1]];{a}=denom[[All,1]];Return[list[[2]]*Global`G3@@{1-a,a1,-var[[1]],var[[2]]}];,

(*H1*)
First[list]==="H1",{a,c,b}=num[[All,1]];{d}=denom[[All,1]];Return[list[[2]]*Global`H1@@{a,b,c,d,var[[1]],var[[2]]}];,

(*H2*)
First[list]==="H2",{b,a,c,d}=num[[All,1]];{e}=denom[[All,1]];Return[list[[2]]*Global`H2@@{a,b,c,d,e,var[[1]],var[[2]]}];,

(*H3*)
First[list]==="H3",{b,a}=num[[All,1]];{c}=denom[[All,1]];Return[list[[2]]*Global`H3@@{a,b,c,var[[1]],var[[2]]}];,

(*H4*)
First[list]==="H4",{b,a}=num[[All,1]];{c,d}=denom[[All,1]];Return[list[[2]]*Global`H4@@{a,b,c,d,var[[1]],var[[2]]}];,

(*H5*)
First[list]==="H5",{a}=num[[All,1]];{b,c}=denom[[All,1]];Return[list[[2]]*Global`H5@@{a,1-b,c,-var[[1]],-var[[2]]}];,

(*H6*)
First[list]==="H6",{a,c}=num[[All,1]];{b}=denom[[All,1]];Return[list[[2]]*Global`H6@@{a,1-b,c,-var[[1]],-var[[2]]}];,

(*H7*)
First[list]==="H7",{a,b,c}=num[[All,1]];{d}=denom[[All,1]];Return[list[[2]]*Global`H7@@{a,b,c,d,var[[1]],var[[2]]}];,


(*KDF*)
StringStartsQ[First[list],"KdF"],
Return[list[[2]]*(*ToExpression[StringDelete[list[[1]], {" ","{","}",";",","}]]*)Global`KdF@@{Join[{Select[#,Last[#]===(Plus@@indices)&][[All,1]]},Table[Select[#,Last[#]===indices[[i]]&][[All,1]],{i,1,Length[indices]}]]&/@{num,denom},var}];,

(*FTilde*)
StringStartsQ[First[list],"FTilde"],
Return[list[[2]]*(*ToExpression[StringDelete[list[[1]], {" ","{","}",";",","}]]*)Global`FTilde@@{Join[{Select[#,Last[#]===indices[[1]]-indices[[2]]&][[All,1]]},Table[Select[#,Last[#]===indices[[i]]&][[All,1]],{i,1,Length[indices]}]]&/@{num,denom},var}];
];

]


variables[n0_]:="n"<>ToString[#]&/@Range[n0]


tableNvar[list_List]:=Module[{var,serieslist},
var=variables[Length[list]];
var=list;
serieslist={
{Sort/@{Append[list,Plus@@list],list},"FA("<>ToString[Length[list]]<>")"},
{Sort/@{Join[list,list],{Plus@@list}},"FB("<>ToString[Length[list]]<>")"},
{Sort/@{Join[{Plus@@list},{Plus@@list}],list},"FC("<>ToString[Length[list]]<>")"},
{Sort/@{Append[list,Plus@@list],{Plus@@list}},"FD("<>ToString[Length[list]]<>")"}

};

Return[serieslist]]


(*KdFQ[var_List,charlist_List]:=Module[{m,n,num,denom,numindex,denomindex},
numindex={};denomindex={};
{m,n}=var;
If[And@@(ContainsOnly[#,{m,n,m+n,1}]&/@charlist),{
{num,denom}=Tally[#]&/@charlist;
numindex={If[#==={},0,#[[1,2]]]&@Select[num,(#[[1]]==m+n)&],If[#==={},0,#[[1,2]]]&@Select[num,(#[[1]]==m)&],If[#==={},0,#[[1,2]]]&@Select[num,(#[[1]]==n)&]};
denomindex={If[#==={},0,#[[1,2]]]&@Select[denom,(#[[1]]==m+n)&],If[#==={},0,#[[1,2]]]&@Select[denom,(#[[1]]==m)&],If[#==={},0,#[[1,2]]]&@Select[denom,(#[[1]]==n)&]};}];
If[numindex==={}&&denomindex==={},Return[{}],Return["KdF"<>ToString[numindex]<>";"<>ToString[denomindex]]];
]*)


KdFQ[var_List,charlist_List]:=Module[{m,n,num,denom,numindex,denomindex},
numindex={};denomindex={};
(*{m,n}=var;*)
If[And@@(ContainsOnly[#,Join[var,{Plus@@var},{1}]]&/@charlist),{
{num,denom}=Tally[#]&/@charlist;
(*Print[{num,denom}];*)
numindex=Table[If[i===0,If[#==={},0,#[[1,2]]]&@Select[num,(#[[1]]==Plus@@var)&],If[#==={},0,#[[1,2]]]&@Select[num,(#[[1]]==var[[i]])&]],{i,0,Length[var]}];
denomindex=Table[If[i===0,If[#==={},0,#[[1,2]]]&@Select[denom,(#[[1]]==Plus@@var)&],If[#==={},0,#[[1,2]]]&@Select[denom,(#[[1]]==var[[i]])&]],{i,0,Length[var]}];}];
If[numindex==={}&&denomindex==={},Return[{}],Return["KdF"<>ToString[numindex]<>";"<>ToString[denomindex]]];
]


FTildeQ[var_List,charlist_List]:=Module[{m,n,num,denom,numindex,denomindex},
numindex={};denomindex={};
{m,n}=var;
If[And@@(ContainsOnly[#,{m,n,m-n,1}]&/@charlist),{
{num,denom}=Tally[#]&/@charlist;
numindex={If[#==={},0,#[[1,2]]]&@Select[num,(#[[1]]==m-n)&],If[#==={},0,#[[1,2]]]&@Select[num,(#[[1]]==m)&],If[#==={},0,#[[1,2]]]&@Select[num,(#[[1]]==n)&]};
denomindex={If[#==={},0,#[[1,2]]]&@Select[denom,(#[[1]]==m-n)&],If[#==={},0,#[[1,2]]]&@Select[denom,(#[[1]]==m)&],If[#==={},0,#[[1,2]]]&@Select[denom,(#[[1]]==n)&]};}];
If[numindex==={}&&denomindex==={},Return[{}],Return["FTilde"<>ToString[numindex]<>";"<>ToString[denomindex]]];
]


callroc[index_List, ser1_]:=Module[{poch3,list1,list11,rocresult},
poch3=(#/(#/.Table[index[[i]]-> 0,{i,1,Length[index]}]))&/@List@@If[Head[ser1]=!=Plus,Flatten[{ser1}],Expand[ser1]];



(*Print[poch2];*)
list11=Table[{{Flatten[List@@If[Head[Numerator[#]]=!=Times,{Numerator[#]},Numerator[#]]/.Power[x_,n_]->Table[x,{i,1,2}]/.Pochhammer[z_,m_]-> m],Flatten[List@@If[Head[Denominator[#]]=!=Times,{Denominator[#]},Denominator[#]]/.Power[x_,n_]->Table[x,{i,1,2}]/.Pochhammer[z_,m_]-> m]}&@(poch3[[i]]/(poch3[[i]]/.Pochhammer[z_,m_]-> 1)),Table[((#/.Pochhammer[z_,m_]-> 1/.Table[If[k===j,index[[k]]-> 1,index[[k]]-> 0],{k,1,Length[index]}]))&@poch3[[i]],{j,1,Length[index]}]},{i,1,Length[poch3]}];list1=Simplify//@DeleteDuplicates[list11];
Print[list1];

(*Which*)If[Length[index]===2,
rocresult=Table[Last[ROC2[index[[1]],index[[2]],list1[[i]][[1]][[1]],list1[[i]][[1]][[2]],list1[[i]][[2]][[1]],list1[[i]][[2]][[2]]]],{i,1,Length[list1]}]
(*,
Length[index]===3,
rocresult=Table[Last[ROC3[index[[1]],index[[2]],index[[3]],list1[[i]][[1]][[1]],list1[[i]][[1]][[2]],list1[[i]][[2]][[1]],list1[[i]][[2]][[2]],list1[[i]][[2]][[3]]]],{i,1,Length[list1]}];*)
];

Return[Flatten[rocresult]];



];


(* ::Chapter:: *)
(*Series recognition lists*)


(* ::Section::Closed:: *)
(*For 2 variables*)


table2var[{m_,n_}]:=Module[{list},

list={

{{{m,n,m+n},{m+n}},"F1"},

{{{m,n,m+n},{m,n}},"F2"},

{{{m,m,n,n},{m+n}},"F3"},

{{{m+n,m+n},{m,n}},"F4"},

{{{m-n,m+n},{m-n}},"G1"},

{{{m,m-n,n},{m-n}},"G2"},

{{{2 m-n},{m-2 n}},"G3"},

{{{m-n,n,m+n},{m}}, "H1"},

{{{m,m-n,n,n},{m}},"H2"},

{{{n,2 m+n},{m+n}},"H3"},

{{{n,2 m+n},{m,n}},"H4"},

{{{2 m+n},{m-n,n}}, "H5"},

{{{2 m-n,n},{m-n}}, "H6"},

{{{2 m-n,n,n},{m}}, "H7"},

{{{2 m,n,m+n},{m,m+n}},"new"}

(*{{{m,m,n,m+n},{m,m+n}},"KDF121;110"},

{{{m,m,m-n,n},{m,m-n}},"FTilde121;110"},

{{{m,m-n,m-n,n},{m-n,m-n}},"FTilde211;200"},

{{{m,m-n,n,n},{m-n,n}},"FTilde112;101"},

{{{m,n,n,m+n},{n,m+n}},"KDF112;101"},

{{{n,n,2 m+n},{n,m+n}},"BSeries NS"},

{{{m,m+n,m+n},{m,m+n}},"KDF210;110"},

{{{m-n,m-n,n,n},{m-n}},"FTilde202;100"},

{{{m,n,m+n,m+n},{m+n,m+n}},"KDF211;200"}*)

};


Return[list]]


(* ::Section::Closed:: *)
(*For 3 variables*)


table3var[{m_,n_,p_}]:=Module[{list},

list={

(*{{{2 m,n,m+n+p},{m,m+n}},"new3"},*)

{{{m,n,p,m+n+p},{m,n,p}},"FA(3)"},

{{{m,m,n,n,p,p},{m+n+p}},"FB(3)"},

{{{m+n+p,m+n+p},{m,n,p}},"FC(3)"},

{{{m,n,p,m+n+p},{m+n+p}},"FD(3)"}


};


Return[list]]


End[]


EndPackage[]
