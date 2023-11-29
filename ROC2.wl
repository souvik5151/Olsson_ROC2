(* ::Package:: *)

BeginPackage["ROC2`"]


Print["ROC2.wl v1.0\n","Authors : B. Ananthanarayan, Souvik Bera, S. Friot & Tanay Pathak"];
(*roc2v18 is ROC2.wl*)


ROC2::usage="ROC2[sum_index_1,sum_index_2,{Pochammers in the numerator},{Pochammers in the denominator},x,y] gives the roc of the double series. Do not include the factorials in the denominator. If there is no denominator put 1";


Begin["`Private`"]


ROC2[m0_,n0_,x0_,y0_,x1_,y1_]:=Module[{x=x0,y=y0,m=m0,n=n0,a,b,t,temp,temp1,rs,rs2,R,S,con,cases1,
cases2,pow1,pow2,pow,z1,z2,f1,num,denom,bound,bound0,bound1x,bound2x,bound3x,bound4x,bound5x,bound1y,
bound2y,bound3y,bound4y,bound5y,bound6y,bound6x,bound5,region0,region1,region2,region3,region4,region5,
region6,lrule,p1,p2,q1,q2,r,s1,s2,t1,t2,u,v,x2,y2,x3,y3,i,j,k,degree,maxdegree,mindegree,(*minfun,*)maxfun,pos,list,region11,region22,h3,u1,v1,gcd},
(*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~lrule takes the limit on Gamma function~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*)

lrule={Gamma[a_.+m x_.+n y_.]->(m x+y n)^a,Gamma[a_.+m x_.]->(m x)^a,Gamma[a_.+n y_.]->(y n)^a};
(*~~~~~~~~~~~~~~~~region initialization~~~~~~~~~~*)
region4[x2_,y2_]={};
region5[x2_,y2_]={};
region2[x2_,y2_]={};
region0[x2_,y2_]={};(*temp[m_,n_]={};*)
(* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~r and s calculation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*)

(*~~~~~~~~~~~~~~~~~~~~~~~seperating num and denom~~~~~~~~~~~~~~~~~~~~*)
num=Times@@(Table[a1[i]=Gamma[b+x[[i]]],{i,1,Length[x]}]);
denom=Times@@(Table[b1[i]=Gamma[b+y[[i]]],{i,1,Length[y]}]);
f1[m_,n_]=(num/(m! n! denom))/.{m0->m,n0->n};(*Print[f1[m0,n0]];*)(*Print[Limit[Expand//@(Simplify[{f1[m+1,n]/f1[m,n],f1[m,n+1]/f1[m,n]}])//.lrule/.m\[Rule]t m/.n\[Rule]t n,t\[Rule]\[Infinity]]//Factor];*)temp1=Factor[Limit[Expand//@(Simplify[{f1[m+1,n]/f1[m,n],f1[m,n+1]/f1[m,n]}])//.lrule/.m->t m/.n->t n,t->\[Infinity]]];(*Print[temp1];*)(*Print["0"];*)temp[m_,n_]=temp1;
(*temp[m_,n_]=Factor[Limit[Expand//@(Simplify[{f1[m+1,n]/f1[m,n],f1[m,n+1]/f1[m,n]}])//.lrule/.m\[Rule]t m/.n\[Rule]t n,t\[Rule]\[Infinity]]];*)(*Print[temp[m0,n0]];*)(*rs[m_,n_]=1/(Factor[Limit[Expand//@(Simplify[{f1[m+1,n]/f1[m,n],f1[m,n+1]/f1[m,n]}])//.lrule/.m\[Rule]t m/.n\[Rule]t n,t\[Rule]\[Infinity]]]);*)rs[m_,n_]=1/temp[m,n];
(*~~~~~~~~~~~~~~~R and S~~~~~~~~~~~~~~~~~~~~~~~~~~~*)
R=Abs[rs[m,n][[1]]/.{m->1,n->0}];
S=Abs[rs[m,n][[2]]/.{m->0,n->1}];
(*bound5[x2_,y2_]={};*)
(*Print[{R,S}];*)
(*Print[rs[m0,n0]];*)
(*~~~~~~~~~~~~~~~~~~~~~~conditions on m and n~~~~~~~~~~~~~~~~*)

(*list declaration*)
bound2x[y2_]={};
bound1x[y2_]={};
bound3x[y2_]={};
bound4x[y2_]={};
bound5x[y2_]={};
bound2y[x2_]={};
bound1y[x2_]={};
bound3y[x2_]={};
bound4y[x2_]={};
bound5y[x2_]={};
region3[x2_,y2_]={};
bound0[x2_,y2_]={};
region11[x2_,y2_]={};
region22[x2_,y2_]={};

(*IF1~~~~~~~if rs[m,n] is number~~~~~~~~~~~~~~~~~~~~~*)If[Or[NumberQ[rs[m,n][[1]]]&&NumberQ[rs[m,n][[2]]],R===\[Infinity]&&S===\[Infinity]]===True,{region0[x2_,y2_]={(Abs[x2]<R)&&(Abs[y2]<S)},bound0[x2_,y2_]={}},
(*~~~~ELSE1~~~~,else general treatment*)
(*finding the roots of numerator and denominator of r and s*)
{(*Print["not number"];*)
p1=Flatten[Solve[Simplify[Numerator[rs[m,n][[1]]/.m->n t1]]==0,t1]];
p2=Flatten[Solve[Simplify[Denominator[rs[m,n][[1]]/.m->n t1]]==0,t1]];q1=Flatten[Solve[Simplify[Numerator[rs[m,n][[2]]/.m->n t1]]==0,t1]];
q2=Flatten[Solve[Simplify[Denominator[rs[m,n][[2]]/.m->n t1]]==0,t1]];r=DeleteDuplicates[Join[p1,p2,q1,q2]];
s1=DeleteDuplicates[Append[Table[a[i]=r[[i]][[2]],{i,1,Length[r]}],0]];
s2=Sort[Select[s1,#>=0&],#1<#2&];
con=Append[Table[a[i]=t1>s2[[i]]&&t1<s2[[i+1]],{i,1,Length[s2]-1}],t1>s2[[Length[s2]]]];
(*Print[con];*)
rs2[m_,n_]={};
pow={};
For[i=1,i<= 2,i++,{{gcd=GCD@@DeleteCases[Cases[Simplify[Abs[#],m>0&&n>0],(a:_Integer:1)x_^p_.:>p,1],1];If[gcd=!=0,pow1=gcd,pow1=1];(*Print[pow1];*)pow=AppendTo[pow,pow1];rs2[m_,n_]=AppendTo[rs2[m,n],PowerExpand[(Simplify[Numerator[#]/Denominator[#],m>0&&n>0])^(1/pow1)]]}&@rs[m,n][[i]]}];(*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~r and s calculation for different conditions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*)
(*Print[rs2[m0,n0]];*)
cases1=Factor//@(DeleteDuplicates[Table[t2[i]=Refine[Abs[Simplify[(rs2[m,n]/.m-> t1 n),m>0&&n>0]]/.m-> t1 n,con[[i]]&&t1>= 0],{i,1,Length[con]}]]);
(*Print[cases1];*)
(*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~boundaries~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*)

bound[u1_,v1_]=Simplify[Flatten[Table[Simplify[GroebnerBasis[{u^(1/pow[[1]])==cases1[[i]][[1]],v^(1/pow[[2]])==cases1[[i]][[2]]},{u^(1/pow[[1]]),v^(1/pow[[2]])},{t1}]],{i,1,Length[cases1]}]/.u->u1^pow[[1]]/.v->v1^pow[[2]]],u1>0&&v1>0];(*Print[pow];*)(*Print[bound[x1,y1]];*)(*minfun={Min@@((List@@@Flatten[#][[All,2]])[[All,1]]),DeleteDuplicates[(List@@@Flatten[#][[All,2]])[[All,2]]][[1]]}&;*)
(*maxfun={Max@@((List@@@Flatten[#][[All,2]])[[All,1]]),DeleteDuplicates[(List@@@Flatten[#][[All,2]])[[All,2]]][[1]]}&;*)
(*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*)
For[i=1,i<= Length[bound[u1,v1]],i++,{degree = {Exponent[#,u1],Exponent[#,v1]}&@Expand[bound[u1,v1][[i]]];,(*Print[degree];*)
mindegree = Min[degree];If[Position[degree,mindegree]==={{1}},{(*Print["x"];*)region3[x2_,y2_]=AppendTo[region3[x2,y2],(region[u1,v1,bound[u1,v1][[i]],R^(1/pow[[1]]),S^(1/pow[[2]]),1][[2]])/.u1-> Abs[x2]^(1/pow[[1]])/.v1-> Abs[y2]^(1/pow[[2]])];bound0[x2_,y2_]=AppendTo[bound0[x2,y2],(First[region[u1,v1,bound[u1,v1][[i]],R^(1/pow[[1]]),S^(1/pow[[2]]),1]])/.u1-> Abs[x2]^(1/pow[[1]])/.v1-> Abs[y2]^(1/pow[[2]])];bound3x[y2_]=AppendTo[bound3y[x2],Last[region[u1,v1,bound[u1,v1][[i]],R^(1/pow[[1]]),S^(1/pow[[2]]),1]]/.v1-> Abs[y2]^(1/pow[[2]])];},
{(*Print["y"];*)region3[x2_,y2_]=AppendTo[region3[x2,y2],(region[u1,v1,bound[u1,v1][[i]],R^(1/pow[[1]]),S^(1/pow[[2]]),2][[2]])/.u1-> Abs[x2]^(1/pow[[1]])/.v1-> Abs[y2]^(1/pow[[2]])];bound0[x2_,y2_]=AppendTo[bound0[x2,y2],(First[region[u1,v1,bound[u1,v1][[i]],R^(1/pow[[1]]),S^(1/pow[[2]]),2]])/.u1-> Abs[x2]^(1/pow[[1]])/.v1-> Abs[y2]^(1/pow[[2]])];bound3y[x2_]=AppendTo[bound3y[x2],Last[region[u1,v1,bound[u1,v1][[i]],R^(1/pow[[1]]),S^(1/pow[[2]]),2]]/.u1-> Abs[x2]^(1/pow[[1]])];}];
}];



(*Print[region3[x1,y1]];
Print["3x",bound3x[y1]];
Print["3y",bound3y[x1]];*)
If[Flatten[bound3x[y2]]=!=List[],{list=DeleteCases[bound3x[y2],{}];
i=1;While[True,If[list===List[],Break[],{(*Print[i];,*)pos ={#1}&/@Flatten[Position[Simplify[LogicalExpand/@#[[All,2]]],Simplify[LogicalExpand[Simplify[And@@#[[All,2]]]]]]&@list];,bound4x[y3_]=({Min@@#[[All,1]],LogicalExpand[Simplify[And@@#[[All,2]]]]}&@list)/.y2-> y3;,(*Print[bound4x[y1]];*)bound5x[y2_]=AppendTo[bound5x[y2],bound4x[y2]];,list=Delete[list,pos];}];i++];}];
(*Print[bound3y[x1]];*)
If[Flatten[bound3y[x2]]=!=List[],{list=DeleteCases[bound3y[x2],{}];
i=1;While[True,If[list===List[],Break[],{(*Print[i];,*)pos ={#1}&/@Flatten[Position[Simplify[LogicalExpand/@#[[All,2]]],Simplify[LogicalExpand[Simplify[And@@#[[All,2]]]]]]&@list];,bound4y[x3_]=({Min@@#[[All,1]],LogicalExpand[Simplify[And@@#[[All,2]]]]}&@list)/.x2-> x3;,(*Print[bound4y[x1]];,*)bound5y[x2_]=AppendTo[bound5y[x2],bound4y[x2]];,list=Delete[list,pos];}];i++];}]}];
(*Print[region3[x1,y1]];*)
(*Print[bound5y[x1]];*)
(*Print[region0[x1,y1]];*)
(*Print["all ok"];*)
(*Print["bound5",bound5y[x1]];*)(*Print[region3[x1,y1]];*)
(*If[region0[x2,y2]=!=List[],region0[x2_,y2_]=Simplify[And@@region0[x2,y2]]];*)
If[region3[x2,y2]=!=List[],region0[x2_,y2_]=FullSimplify[And@@Flatten[region3[x2,y2]]]];(*Print[region0[x1,y1]];;*)
If[bound5y[x2]=!=List[],{region1[x2_,y2_]=Refine[Piecewise[bound5y[x2],Undefined],0<Abs[x2]<R&&0<Abs[y2]<S];(*Print["region1",region1[x1,y1]];*)
region11[x2_,y2_]=0<Abs[x2]<R&&0<Abs[y2]<S&&Normal[Abs[y2]^(1/pow[[2]])<region1[x2,y2]];,
region0[x2_,y2_]={And@@DeleteCases[Flatten[List@@@LogicalExpand/@#],0<Abs[x_]]}&@Flatten[{region11[x2,y2]}]}];
If[bound5x[y2]=!=List[],{region2[x2_,y2_]=Refine[Piecewise[bound5x[y2],Undefined],0<Abs[x2]<R&&0<Abs[y2]<S];
region22[x2_,y2_]=0<Abs[x2]<R&&0<Abs[y2]<S&&Normal[Abs[x2]^(1/pow[[1]])<region2[x2,y2]];,region0[x2_,y2_]={And@@DeleteCases[Flatten[List@@@LogicalExpand/@#],0<Abs[x_]]}&@Flatten[{region22[x2,y2]}]}];
If[region0[x2,y2]===List[],region0[x2_,y2_]=Append[region0[x2,y2],Abs[x2]<R&&Abs[y2]<S]];
Return[{"{R,S}, Cartesian Curve, ROC -> ",{R,S},(*rs[m0,n0],con,cases1,bound3[x1,y1],*)(*,rs2[m0,n0],bound[x1,y1],pow,con,cases1,bound[x1,y1],region5[x1,y1],*)(*cases1,pow1,pow2,bound[x1,y1]*)(*,bound[x1,y1],bound2[x1,y1],bound3[x1,y1]*)Simplify[Flatten[bound0[x1,y1]]](*region0[x1,y1],region4[x1,y1],*)(*region4[x1,y1],region5[x1,y1]*),Flatten[{region0[x1,y1]}]}]]



(*~~~~~~~~~~~~~~~~~~~~~~~~~~~ minfun~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*)


minfun[x0_,bound0_]:=Module[{x4,x5,bound1y,bound4y,bound5y,pos,list},bound1y[x5_]=bound0/.x0-> x5;(*Print[bound1y[x5]];*)bound4y[x5_]={};bound5y[x5_]={};
If[Flatten[bound1y[x5]]=!=List[],{list=List@@@Flatten[bound1y[x5]][[All,2]];
i=1;While[True,If[list===List[],Break[],{If[MemberQ[Head/@list[[All,2]],Or]===True,
{pos=#1&/@Position[Simplify[LogicalExpand/@#[[All,2]]],Simplify[LogicalExpand[Simplify[And@@#[[All,2]]]]]]&@list;
bound4y[x4_]=({Min@@#[[First/@pos,1]],Simplify[And@@#[[All,2]]]}&@list)/.x5-> x4;bound5y[x4_]=AppendTo[bound5y[x4],bound4y[x4]];
list={{list[[First[Flatten[Cases[pos,{_,_}]]],1]],Delete[list[[All,2]],pos][[1]]}};(*Print["or"];*)},
{(*Print[i];*)pos ={#1}&/@Flatten[Position[Simplify[LogicalExpand/@#[[All,2]]],Simplify[LogicalExpand[Simplify[And@@#[[All,2]]]]]]&@list];,
bound4y[x4_]=({Min@@#[[Flatten[pos],1]],Simplify[And@@#[[All,2]]]}&@list)/.x5-> x4;,(*Print[bound4y[x1]];,*)
bound5y[x4_]=AppendTo[bound5y[x4],bound4y[x4]];,list=Delete[list,pos];}]}];i++];}];(*Print[bound5y[x0]];*)Return[bound5y[x0]]]


(* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  region  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*)



region[x0_,y0_,bound0_,R0_,S0_,n0_]:=Module[{x,y,bound,R,S,n,i,j,k,bound1y,bound6y,bound1x,bound6x,bound1,region3,region2,region1,bound3y,bound2y,bound3x,bound2x(*,minfun*),orbound,andbound,region0,xmax,ymax,result},
R=R0;S=S0;n=n0;
Off[MaxValue::infeas];
bound[x_,y_]=bound0/.x0-> x/.y0-> y;
(*Print[bound[x,y]];*)
result[x_,y_]={};
bound1[x_,y_]={};
andbound[x_]={};
orbound[x_]={};
bound3y[x_]={};bound3x[y_]={};
region3[x_,y_]={};region2[x_,y_]={};region1[x_,y_]={};(*region0[x_,y_]={};*)
(*minfun={Min@@((List@@@Flatten[#][[All,2]])[[All,1]]),DeleteDuplicates[(List@@@Flatten[#][[All,2]])[[All,2]]][[1]]}&;*)
Which[n===2,
{bound1y[x_]=Simplify[Solve[{bound[x,y]==0,R>x> 0&&y> 0},y,Reals]];(*Print["bound1y",bound1y[x0]];*)If[bound1y[x]=!=List[],{bound6y[x_]=(List@@@Flatten[Normal[bound1y[x]]])[[All,2]];(*Print["here"];*)(*Print["bound6y ",bound6y[x0]];*)For[k=1,k<= Length[bound6y[x]],k++,{If[Reduce[{y==bound6y[x][[k]],0<x<R&&0<y<S}]=!= False,bound1[x_,y_]=AppendTo[bound1[x,y],y==bound6y[x][[k]]]]}];,(*bound00[x_,y_]={Refine[{y==bound6y[x],0<x<R&&0<y<S}]};*)If[MemberQ[Limit[bound6y[x],{x-> 0}],0]===True,{For[j=1,j<= Length[bound6y[x]],j++,{(*Print[bound6y[x]];*)If[Limit[bound6y[x][[j]],{x-> 0}]===0,
{ymax=Simplify[MaxValue[{bound6y[x][[j]],bound6y[x][[j]]<= S},x]];
xmax=Max[Simplify[MaxValue[{#,0<=  #<= R,0<=  y<= S},y,Reals]&/@DeleteCases[Flatten[Simplify[Normal[Solve[{y==bound6y[x][[j]],x>=  0&&y>=   0},x]]]/.Rule-> List],x]]];
result[x_,y_]=Or[(y<bound6y[x][[j]]&&0<= x<= xmax),y>bound6y[x][[j]]&&0<=y<=ymax&&Flatten[Solve[ymax==bound6y[x][[j]],x]/.Rule-> Less][[1]]]//Simplify;
region3[x_,y_]=AppendTo[region3[x,y],result[x,y]];(*((y<bound6y[x][[j]])&&(x<R&&y<S))||((y<MaxValue[{bound6y[x][[j]],And[bound6y[x][[j]]<S,((List@@@Flatten[bound1y[x]][[All,2]])[[All,2]])[[1]]]},x])&&(x<Reduce[{bound6y[x][[j]]==MaxValue[{bound6y[x][[j]],And[bound6y[x][[j]],((List@@@Flatten[bound1y[x]][[All,2]])[[All,2]])[[1]]]},x],x\[GreaterEqual] 0},x,Reals][[2]]))*)(*((y<bound6y[x][[j]])&&(x<R&&y<S))||((y<MaxValue[{bound6y[x][[j]],((List@@@Flatten[bound1y[x]][[All,2]])[[All,2]])[[1]]},x])&&(x<Reduce[{bound6y[x][[j]]==MaxValue[{bound6y[x][[j]],((List@@@Flatten[bound1y[x]][[All,2]])[[All,2]])[[1]]},x],x\[GreaterEqual] 0},x,Reals][[2]]))*)(*Print["region3 ",region3[x0,y0]];*)},{If[And[Flatten[Solve[{bound6y[x][[j]]==S,0< x<= R},x]]==={},Limit[bound6y[x][[j]],x-> 0]<= S,Solve[{bound6y[x][[j]]==y&&x==R,0<y<= S},{y},Reals]==={}],{orbound[x_]=AppendTo[orbound[x],bound6y[x][[j]]];(*Print["orbound",orbound[x0]];*)},{andbound[x_]=AppendTo[andbound[x],bound6y[x][[j]]];(*Print["andbound",andbound[x0]];*)}]}];
}];(*Print["here"];*)If[andbound[x]=!={},region1[x_,y_]=AppendTo[region1[x,y],(And@@region3[x,y])&&(And@@Table[b[i]=y<andbound[x][[i]],{i,1,Length[andbound[x]]}])],region1[x_,y_]=AppendTo[region1[x,y],And@@region3[x,y]]];(*Print[region1[x0,y0]];*)
(*Print[1];*)If[orbound[x]=!={},region1[x_,y_]={And@@region1[x,y]||(And@@Table[b[i]=y<orbound[x][[i]],{i,1,Length[orbound[x]]}])};(*region1[x_,y_]=And@@region1[x,y];*)];
(*Print["reg1",region1[x0,y0]];*)},{bound2y[x_]=Refine[minfun[x,bound1y[x]]];,bound3y[x_]=AppendTo[bound3y[x],bound2y[x]]}]}]},
n===1,
{bound1x[y_]=Simplify[Solve[{bound[x,y]==0,S>y> 0&&x> 0},x,Reals]];(*Print[bound1x[y0]];*)If[bound1x[y]=!=List[],{bound6x[y_]=(List@@@Flatten[Normal[bound1x[y]]])[[All,2]];(*Print["bou6x ",bound6x[y0]];*)For[k=1,k<= Length[bound6x[y]],k++,{If[Reduce[{x==bound6x[y][[k]],0<x<R&&0<y<S}]=!= False,bound1[x_,y_]=AppendTo[bound1[x,y],x==bound6x[y][[k]]]]}];(*Print[bound1[x0,y0]];*)(*bound00[x_,y_]={Refine[{x==bound6x[y],0<x<R&&0<y<S}]};*)(*here !*)(*Print[bound00[x0,y0]];*)If[MemberQ[Limit[bound6x[y],{y-> 0}],0]===True,{For[j=1,j<= Length[bound6x[y]],j++,{If[Limit[bound6x[y][[j]],{y-> 0}]===0,{xmax=Simplify[MaxValue[{bound6x[y][[j]],bound6x[y][[j]]<= R},y]];
ymax=Max[Simplify[MaxValue[{#,0<=  #<= S,0<=  x<= R},x,Reals]&/@DeleteCases[Flatten[Simplify[Normal[Solve[{x==bound6x[y][[j]],y>=  0&&x>=   0},y]]]/.Rule-> List],y]]];
result[x_,y_]=Or[(x<bound6x[y][[j]]&&0<= y<= ymax),x>bound6x[y][[j]]&&0<=x<=xmax&&Flatten[Solve[xmax==bound6x[y][[j]],y]/.Rule-> Less][[1]]]//Simplify;
region3[x_,y_]=AppendTo[region3[x,y],result[x,y]];(*region3[x_,y_]=AppendTo[region3[x,y],((x<bound6x[y][[j]])&&(x<R&&y<S))||((x<MaxValue[{bound6x[y][[j]],((List@@@Flatten[bound1x[y]][[All,2]])[[All,2]])[[1]]},y])&&(y<Reduce[{bound6x[y][[j]]==MaxValue[{bound6x[y][[j]],((List@@@Flatten[bound1x[y]][[All,2]])[[All,2]])[[1]]},y],y\[GreaterEqual] 0},y,Reals][[2]]))];*)(*Print[region3[x0,y0]];*)},{If[And[Flatten[Solve[{bound6x[y][[j]]==R,0< y<= S},y]]==={},Limit[bound6x[y][[j]],y-> 0]<= R,Solve[{bound6x[y][[j]]==x&&y==S,0<x<= R},{x},Reals]==={}],{orbound[y_]=AppendTo[orbound[y],bound6x[y][[j]]];(*Print["orbound",orbound[y0]];*)},{andbound[y_]=AppendTo[andbound[y],bound6x[y][[j]]];(*Print["andbound",andbound[y0]];*)}]}];
}];If[andbound[y]=!={},region1[x_,y_]=AppendTo[region1[x,y],(And@@region3[x,y])&&(And@@Table[b[i]=x<andbound[y][[i]],{i,1,Length[andbound[y]]}])],region1[x_,y_]=AppendTo[region1[x,y],And@@region3[x,y]]];(*Print[region1[x0,y0]];*)
(*Print[1];*)If[orbound[y]=!={},region1[x_,y_]=And@@region1[x,y]||(And@@Table[b[i]=x<orbound[y][[i]],{i,1,Length[orbound[y]]}]);(*region1[x_,y_]=And@@region1[x,y]*)];
(*Print[region1[x0,y0]];*)},{bound2x[y_]=Refine[minfun[y,bound1x[y]]];,bound3x[y_]=AppendTo[bound3x[y],bound2x[y]]}]}]}
];
(*Print["bound2x ", bound2x[y0]];*)(*Print["reg1"];Print[region1[x0,y0]];Print[And[Abs[x0]<R&&Abs[y0]<S,And@@Simplify[region1[x0,y0]]]];*)
If[region1[x,y]=!=List[],{region0[x_,y_]:=And[x<R&&y<S,(And@@(Simplify[region1[x,y]]))];(*region0[x_,y_]=region0/.x\[Rule] x0/.y\[Rule] y0;*)(*Print[region0[x0,y0](*/.x\[Rule] x0/.y\[Rule] y0*)];*)},region0[x_,y_]=x<R&&y<S];
Return[{Flatten[bound1[x0,y0]],Flatten[{region0[x0,y0]}],If[n===1,Flatten[bound3x[y0]],Flatten[bound3y[x0]]]}]] 


End[]

EndPackage[]

