(* ::Package:: *)

(* ::Title:: *)
(*This script was used to perform calculations for the paper:*)
(*Seyed Mahmoud Arzideh, Andr\[EAcute]s C\[OAcute]rdoba, Jeffrey G. Ethier, Jay D. Schieber, David C. Venerus; Equibiaxial elongation of entangled polyisobutylene melts: Experiments and theoretical predictions. J. Rheol. 1 May 2024; 68 (3): 341\[Dash]353. https://doi.org/10.1122/8.0000809.*)


(* ::Subtitle:: *)
(*The paper was submitted for publication in the Journal of Rheology on 21-Dec-2023.*)
(*email: andcorduri@gmail.com (Andr\[EAcute]s C\[OAcute]rdoba)*)


(* ::Section:: *)
(*Definition of Functions*)


(* ::Input::Initialization:: *)
(*This function solves the Rolie-Double-Poly (RDP) model of Boudara et. al. [J.Rheol. 63(1), 71\[Dash]91] and calculates the entropy generation rate using the method of Schieber and C\[OAcute]rdoba [Phys. Fluids 33(8)]*)
Clear[RDPmodel]
RDPmodel[\[Epsilon]d0_?NumericQ,m_?NumericQ,MWDist_?MatrixQ,Me_?NumericQ,\[Tau]e_?NumericQ,GN0_?NumericQ,\[Lambda]max_?NumericQ,rang_?ArrayQ,ICin_?ArrayQ]:=Module[{\[Kappa],c,d\[Gamma],cxx,cyx,cxy,cyy,czz,RHS,ec1,ec2,ec3,IC,crelax,sol,cond,nm,Z,C1=1.69,C2=4.17,C3=-1.55,\[Tau]d,\[Tau]s,Rrep,RCR,Rret,\[Beta]th=1,\[Beta]CCR=1,\[Delta]=-1/2,RCCR,res1,res2,res3,sampl,dfdc,entropyg,R=8.3145,A,res1b,res1c,resRrep,resRCR,resRret,resRCCR,rescrelax,t},
Z=#[[1]]/Me&/@MWDist;
nm=Length[MWDist];
\[Tau]d=3#^3 \[Tau]e (1-(2C1)/#^(1/2)+C2/#+C3/#^(3/2))&/@Z;
\[Tau]s=#^2 \[Tau]e&/@Z;
c=Table[{{cxx[i,j][t],0,0},{0,cyy[i,j][t],0},{0,0,czz[i,j][t]}},{i,1,nm},{j,nm}];
A=Table[Sum[MWDist[[j,2]]*c[[i,j]],{j,1,nm}],{i,1,nm}];
\[Kappa]=\[Epsilon]d0*{{1,0,0},{0,m,0},{0,0,-(1+m)}};
Rrep=Table[(1/\[Tau]d[[i]])(c[[i,j]]-IdentityMatrix[3]),{i,1,nm},{j,1,nm}];
RCR=Table[(\[Beta]th/\[Tau]d[[j]])(c[[i,j]]-IdentityMatrix[3]),{i,1,nm},{j,1,nm}];
Rret=Table[((2(1-(Tr[A[[i]]]/3)^(-1/2)))/\[Tau]s[[i]])((1-\[Lambda]max^-2)/(1-(Tr[A[[i]]]/3)\[Lambda]max^-2))c[[i,j]],{i,1,nm},
{j,1nm}];
RCCR=Table[2 \[Beta]CCR ((1-(Tr[A[[j]]]/3)^(-1/2))/\[Tau]s[[j]])((1-\[Lambda]max^-2)/(1-(Tr[A[[j]]]/3)\[Lambda]max^-2))(Tr[A[[i]]]/3)^\[Delta] (c[[i,j]]-IdentityMatrix[3]),{i,1,nm},{j,1,nm}];
crelax=-(Rrep+RCR+Rret+RCCR);
RHS=Table[\[Kappa].(c[[i,j]])+(c[[i,j]]).Transpose[\[Kappa]],{i,1,nm},{j,1,nm}]+crelax;
ec1=Flatten[Table[Drop[Sort[Flatten[MapThread[#1==#2&,{D[c[[i,j]],t],RHS[[i,j]]},2]]],6],{i,1,nm},{j,1,nm}]];
IC=Flatten[Table[{(cxx[i,j][t]/.t->0)==ICin[[i,j,1]],(cyy[i,j][t]/.t->0)==ICin[[i,j,2]],(czz[i,j][t]/.t->0)==ICin[[i,j,3]]},{i,1,nm},{j,1,nm}]];
sol=NDSolve[Join[ec1,IC],Flatten[Diagonal[#]&/@Flatten[c,1]],{t,0,rang[[2]]}][[1]];
sampl=Exp[Table[i,{i,N[Log[rang[[1]]]],N[Log[rang[[2]]]],(N[Log[rang[[2]]]]-N[Log[rang[[1]]]])/100}]];
res1=A/.sol/.t->#&/@sampl;
res1b=Table[Inverse[#[[i]]],{i,1,nm}]&/@res1;
res1c=Table[Tr[#[[i]]],{i,1,nm}]&/@res1;
resRrep=Table[(1/\[Tau]d[[i]])((c[[i,j]]/.sol)-IdentityMatrix[3]),{i,1,nm},{j,1,nm}]/.t->#&/@sampl;
resRCR=Table[(\[Beta]th/\[Tau]d[[j]])((c[[i,j]]/.sol)-IdentityMatrix[3]),{i,1,nm},{j,1,nm}]/.t->#&/@sampl;
resRret=MapThread[Table[((2(1-(#2[[i]]/3)^(-1/2)))/\[Tau]s[[i]])((1-\[Lambda]max^-2)/(1-(#2[[i]]/3)\[Lambda]max^-2))(c[[i,j]]/.sol),{i,1,nm},
{j,1nm}]/.t->#1&,{sampl,res1c}];
resRCCR=MapThread[Table[2 \[Beta]CCR ((1-(#2[[j]]/3)^(-1/2))/\[Tau]s[[j]])((1-\[Lambda]max^-2)/(1-(#2[[j]]/3)\[Lambda]max^-2))(#2[[i]]/3)^\[Delta] ((c[[i,j]]/.sol)-IdentityMatrix[3]),{i,1,nm},{j,1,nm}]/.t->#1&,{sampl,res1c}];
rescrelax=MapThread[-(#1+#2+#3+#4)&,{resRrep,resRCR,resRret,resRCCR}];
res2=Sum[MWDist[[i,2]]((1-\[Lambda]max^-2)/(1-(Tr[#[[i]]]/3)\[Lambda]max^-2))#[[i]],{i,1,nm}]&/@res1 ;
res3=MapThread[{#1,If[\[Epsilon]d0>0,GN0 (#2[[1,1]]-#2[[3,3]])/\[Epsilon]d0,0]}&,{sampl,res2}];
entropyg=MapThread[{#1,-(1/2)Sum[Sum[((Z[[i]]R MWDist[[i,2]])/MWDist[[i,1]])((((1-\[Lambda]max^-2)/(1-(#3[[i]]/3)\[Lambda]max^-2))KroneckerDelta[k,l])-#2[[i,k,l]])*If[i==j,MWDist[[i,2]],MWDist[[j,2]]]*#4[[j,i,k,l]],{i,1,nm},{j,1,nm}],{k,1,3},{l,1,3}] }&,{sampl,res1b,res1c,rescrelax}];
{res3,entropyg,Table[{(cxx[i,j][t]/.sol/.t->Last[sampl]),(cyy[i,j][t]/.sol/.t->Last[sampl]),(czz[i,j][t]/.sol/.t->Last[sampl])},{i,1,nm},{j,nm}]}]


(* ::Input::Initialization:: *)
(*This function solves the Rolie-Double-Poly (RDP) model of Boudara et. al. [J.Rheol. 63(1), 71\[Dash]91] and calculates the entropy generation rate using the method of Schieber and C\[OAcute]rdoba [Phys. Fluids 33(8)]*)
Clear[RDPmodelC2]
RDPmodelC2[trc_?NumericQ,cxx_?NumericQ,MWDist_?MatrixQ,Me_?NumericQ,\[Tau]e_?NumericQ,GN0_?NumericQ,\[Lambda]max_?NumericQ]:=Module[{\[Kappa],c,d\[Gamma],cyx,cxy,cyy,czz,RHS,ec1,ec2,ec3,IC,crelax,sol,cond,nm,Z,C1=1.69,C2=4.17,C3=-1.55,\[Tau]d,\[Tau]s,Rrep,RCR,Rret,\[Beta]th=1,\[Beta]CCR=1,\[Delta]=-1/2,RCCR,res1,res2,res3,sampl,dfdc,entropyg,R=8.3145,A,res1b,res1c,resRrep,resRCR,resRret,resRCCR,rescrelax,t},
Z=#[[1]]/Me&/@MWDist;
nm=Length[MWDist];
\[Tau]d=3#^3 \[Tau]e (1-(2C1)/#^(1/2)+C2/#+C3/#^(3/2))&/@Z;
\[Tau]s=#^2 \[Tau]e&/@Z;
sol=Table[NSolve[{trc==cxx+cyy+czz,cyy==czz},{cyy,czz}],{i,1,nm},{j,1,nm}];
c=Table[{{cxx,0,0},{0,cyy,0},{0,0,czz}}/.sol[[i,j,1]],{i,1,nm},{j,nm}];
A=Table[MWDist[[i,2]]*c[[i,i]]+Sum[If[i!=j,MWDist[[j,2]]*c[[i,j]],0],{j,1,nm}],{i,1,nm}];
res1b=Table[Inverse[A[[i]]],{i,1,nm}];
res1c=Table[Tr[A[[i]]],{i,1,nm}];
resRrep=Table[(1/\[Tau]d[[i]])(c[[i,j]]-IdentityMatrix[3]),{i,1,nm},{j,1,nm}];
resRCR=Table[(\[Beta]th/\[Tau]d[[j]])(c[[i,j]]-IdentityMatrix[3]),{i,1,nm},{j,1,nm}];
resRret=Table[((2(1-(res1c[[i]]/3)^(-1/2)))/\[Tau]s[[i]])((1-\[Lambda]max^-2)/(1-(res1c[[i]]/3)\[Lambda]max^-2))c[[i,j]],{i,1,nm},{j,1nm}];
resRCCR=Table[2 \[Beta]CCR ((1-(res1c[[j]]/3)^(-1/2))/\[Tau]s[[j]])((1-\[Lambda]max^-2)/(1-(res1c[[j]]/3)\[Lambda]max^-2))(res1c[[i]]/3)^\[Delta] (c[[i,j]]-IdentityMatrix[3]),{i,1,nm},{j,1,nm}];
rescrelax=MapThread[-(#1+#2+#3+#4)&,{resRrep,resRCR,resRret,resRCCR}];
res2=Sum[MWDist[[i,2]]((1-\[Lambda]max^-2)/(1-(Tr[A[[i]]]/3)\[Lambda]max^-2))A[[i]],{i,1,nm}] ;
entropyg=-(1/2)Sum[Sum[((Z[[i]]R MWDist[[i,2]])/MWDist[[i,1]])((((1-\[Lambda]max^-2)/(1-(res1c[[i]]/3)\[Lambda]max^-2))KroneckerDelta[k,l])-res1b[[i,k,l]])*If[i==j,MWDist[[i,2]],MWDist[[j,2]]]*rescrelax[[j,i,k,l]],{i,1,nm},{j,1,nm}],{k,1,3},{l,1,3}]]


(* ::Input::Initialization:: *)
Clear[PowerTicks]
PowerTicks[label_][min_,max_]:=Block[{min10,max10},min10=Floor[Log10[min]];
max10=Ceiling[Log10[max]];
Join[Table[{10^i,If[label,Superscript["10",i],Spacer[{0,0}]]},{i,min10,max10}],Flatten[Table[{k 10^i,Spacer[{0,0}],{0.005,0.`},{Thickness[0.001`]}},{i,min10,max10},{k,9}],1]]]


(* ::Input::Initialization:: *)
Clear[ErrBarP]
(*The function "ErrBarP" is for plotting the error bars in a log-log plot*)
ErrBarP[data_/;MatrixQ[data,NumericQ],err_/;VectorQ[err,NumericQ],del_?NumericQ,color_]:=
Module[{data2,data3,err2},
data2=Select[MapThread[{#1[[1]],#1[[2]],#2}&,{data,err}],#[[2]]-#[[3]]>0&&#[[2]]+#[[3]]>0&];
data3={#[[1]],#[[2]]}&/@data2;err2=#[[3]]&/@data2;
{{color,MapThread[Line[{Log[{#1[[1]],(#1[[2]]-#2)}],Log[{#1[[1]],(#1[[2]]+#2)}]}]&,{data3,err2}]},{color,MapThread[Line[{{Log[#1[[1]]]-del,Log[(#1[[2]]-#2)]},{Log[#1[[1]]]+del,Log[(#1[[2]]-#2)]}}]&,{data3,err2}]},{color,MapThread[Line[{{Log[#1[[1]]]-del,Log[(#1[[2]]+ #2)]},{Log[#1[[1]]]+del,Log[(#1[[2]]+#2)]}}]&,{data3,err2}]}}]


(* ::Section:: *)
(*Sample calculations for the MW = 108 kDa and PDI = 3.3 system*)


(* ::Input::Initialization:: *)
(*Specifying the Molecular weight distribution*)
Mn=32.5;Mw=108;
Mmp=Mn^(5/2)/Mw^(3/2);soldist={\[Mu]->Log[Mn^(3/2)/Sqrt[Mw]],\[Sigma]->Sqrt[2] Sqrt[Log[Sqrt[Mw]/Sqrt[Mn]]],C[1]->0,C[2]->0};
cutoffdist=FindRoot[10^-4 (PDF[LogNormalDistribution[\[Mu],\[Sigma]],Mmp]/.soldist)==PDF[LogNormalDistribution[\[Mu],\[Sigma]],xc]/.soldist,{xc,Mw}];
MCv=3000;
MkPIB=274.2;
NCsAll={10,14,20,28,34,40};
NKsAll={20,30,40,60,80};
NCsAllLimits=Append[Prepend[Join[Drop[NKsAll*MkPIB/MCv,1],NCsAll],0],xc*1000/MCv/.cutoffdist];
MwRDP=Table[{Mean[{NCsAllLimits[[i]]*MCv/1000,NCsAllLimits[[i+1]]*MCv/1000}],NIntegrate[PDF[LogNormalDistribution[\[Mu],\[Sigma]],x]/.soldist,{x,NCsAllLimits[[i]]*MCv/1000,NCsAllLimits[[i+1]]*MCv/1000}]},{i,1,Length[NCsAllLimits]-1}];


(* ::Input::Initialization:: *)
(*Discretized form of the molecular weight distribution used in the calculations below*)
Prepend[MwRDP,{"M,kDa","fraction"}]//TableForm


(* ::Input::Initialization:: *)
Quiet[Needs["PlotLegends`"]];
distD=EmpiricalDistribution[(#[[2]]&/@MwRDP)->(#[[1]]&/@MwRDP)];
Legend1={{{Graphics[{Black, Thick,Line[{{0,0},{3,0}}]}],Text[Style["Log-Normal distribution",FontSize->15,FontFamily->"Arial",Black]]},{Graphics[{Red,Dashed, Thick,Line[{{0,0},{3,0}}]}],Text[Style["Discretized version used in RDP",FontSize->15,FontFamily->"Arial",Black]]}},LegendShadow->False,LegendSize->{0.95,0.2},LegendPosition->{-0.1,-0.5},LegendTextSpace->6,LegendBorder->False};
ShowLegend[Plot[{CDF[LogNormalDistribution[\[Mu],\[Sigma]],x]/.soldist/.{C[1]->0,C[2]->0}/.{Mn->113.950,Mw->130.414},CDF[distD,x]},{x,0,130},Frame->True,ImageSize->600,BaseStyle->{FontSize->21,FontFamily->"Arial"},FrameLabel->{"M, kDa","Cumulative probability"},PlotStyle->{Black,{Red,Dashed}},AspectRatio->0.8,Exclusions->None],Legend1]


(* ::Input::Initialization:: *)
(*Specifying the Parameters*)
mIn=1;(*This parameter defined the type of elongational flow, 1 is for biaxial, simple elongation is -1/2 and planar is 0. See Meissner et. al., J. Nonnewton. Fluid. Mech. ,11 (1982) 221-237*)
MeIn=6.690; (*Entanglement molecular weight in kDa, taken from Fetters et. al.*)
\[Tau]eIn=0.08;(*Rouse relaxation time of one entanglement
segment in seconds, Subscript[\[Tau], e]. I chose this to be Subscript[\[Tau], e]=0.502*Subscript[tau, c] based on Macromolecules,vol.54,pp.8033\[Dash]8042,2021. Where Subscript[\[Tau], c] is the cluster-shuffling characteristic time in CFSM*)
GN0In=253.035;(*The plateau modulus in kPa.*)
\[Lambda]maxIn=2(*Maximum stretch ratio*);


(* ::Subsection:: *)
(*Inception of biaxial elongation*)


(* ::Input::Initialization:: *)
sol\[Epsilon]dp3=RDPmodel[0.3,mIn,MwRDP,MeIn,\[Tau]eIn,GN0In,\[Lambda]maxIn,{10^-2,2*10^1},Table[{1,1,1},{i,1,Length[MwRDP]},{j,1,Length[MwRDP]}]];


(* ::Input::Initialization:: *)
sol\[Epsilon]dp1=RDPmodel[0.1,mIn,MwRDP,MeIn,\[Tau]eIn,GN0In,\[Lambda]maxIn,{10^-2,5*10^1},Table[{1,1,1},{i,1,Length[MwRDP]},{j,1,Length[MwRDP]}]];


(* ::Input::Initialization:: *)
sol\[Epsilon]dp03=RDPmodel[0.03,mIn,MwRDP,MeIn,\[Tau]eIn,GN0In,\[Lambda]maxIn,{10^-2,2*10^2},Table[{1,1,1},{i,1,Length[MwRDP]},{j,1,Length[MwRDP]}]];


(* ::Input::Initialization:: *)
sol\[Epsilon]dp01=RDPmodel[0.01,mIn,MwRDP,MeIn,\[Tau]eIn,GN0In,\[Lambda]maxIn,{10^-2,5*10^2},Table[{1,1,1},{i,1,Length[MwRDP]},{j,1,Length[MwRDP]}]];


(* ::Input::Initialization:: *)
Quiet[Needs["PlotLegends`"]];
Legend1={{{Graphics[{Blue, Thick,Line[{{0,0},{6,0}}]}],Text[Style["RDP,\!\(\*OverscriptBox[\(\(\\\ \)\*SubscriptBox[\(\[Epsilon]\), \(B\)]\), \(.\)]\)=0.01 \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",FontSize->16,FontFamily->"Arial",Blue]]},{Graphics[{Darker[Green], Thick,Line[{{0,0},{6,0}}]}],Text[Style["RDP,\!\(\*OverscriptBox[\(\(\\\ \)\*SubscriptBox[\(\[Epsilon]\), \(B\)]\), \(.\)]\)=0.03 \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",FontSize->16,FontFamily->"Arial",Darker[Green]]]},{Graphics[{Darker[Cyan], Thick,Line[{{0,0},{6,0}}]}],Text[Style["RDP,\!\(\*OverscriptBox[\(\(\\\ \)\*SubscriptBox[\(\[Epsilon]\), \(B\)]\), \(.\)]\)=0.1 \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",FontSize->16,FontFamily->"Arial",Darker[Cyan]]]},{Graphics[{Magenta, Thick,Line[{{0,0},{6,0}}]}],Text[Style["RDP,\!\(\*OverscriptBox[\(\(\\\ \)\*SubscriptBox[\(\[Epsilon]\), \(B\)]\), \(.\)]\)=0.3 \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",FontSize->16,FontFamily->"Arial",Magenta]]}},LegendShadow->False,LegendSize->{0.63,0.3},LegendPosition->{-0.76,0.35},LegendTextSpace->4,LegendBorder->False};
ShowLegend[ListLogLogPlot[{sol\[Epsilon]dp01[[1]],sol\[Epsilon]dp03[[1]],sol\[Epsilon]dp1[[1]],sol\[Epsilon]dp3[[1]]},PlotRange->{{0.5*10^-2,4*10^2},{1*10^1,3*10^4}},Frame->True,Axes->False,BaseStyle->{FontSize->22,FontFamily->"Arial"},AspectRatio->0.7,FrameLabel->{"t, s","\!\(\*SubscriptBox[\(\[Eta]\), \(B\)]\), kPa\[CenterDot]s"},ImageSize->700,Joined->True,PlotStyle->{{Blue},{Darker[Green]},{Darker[Cyan]},Magenta},FrameTicks->{{PowerTicks[True],PowerTicks[False]},{PowerTicks[True],PowerTicks[False]}}],Legend1]


(* ::Input::Initialization:: *)
Quiet[Needs["PlotLegends`"]];
Legend1={{{Graphics[{Blue, Thick,Line[{{0,0},{6,0}}]}],Text[Style["RDP,\!\(\*OverscriptBox[\(\(\\\ \)\*SubscriptBox[\(\[Epsilon]\), \(B\)]\), \(.\)]\)=0.01 \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",FontSize->16,FontFamily->"Arial",Blue]]},{Graphics[{Darker[Green], Thick,Line[{{0,0},{6,0}}]}],Text[Style["RDP,\!\(\*OverscriptBox[\(\(\\\ \)\*SubscriptBox[\(\[Epsilon]\), \(B\)]\), \(.\)]\)=0.03 \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",FontSize->16,FontFamily->"Arial",Darker[Green]]]},{Graphics[{Darker[Cyan], Thick,Line[{{0,0},{6,0}}]}],Text[Style["RDP,\!\(\*OverscriptBox[\(\(\\\ \)\*SubscriptBox[\(\[Epsilon]\), \(B\)]\), \(.\)]\)=0.1 \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",FontSize->16,FontFamily->"Arial",Darker[Cyan]]]},{Graphics[{Magenta, Thick,Line[{{0,0},{6,0}}]}],Text[Style["RDP,\!\(\*OverscriptBox[\(\(\\\ \)\*SubscriptBox[\(\[Epsilon]\), \(B\)]\), \(.\)]\)=0.3 \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",FontSize->16,FontFamily->"Arial",Magenta]]}},LegendShadow->False,LegendSize->{0.63,0.3},LegendPosition->{-0.74,0.35},LegendTextSpace->4,LegendBorder->False};
InsetP1=ListLogLinearPlot[{sol\[Epsilon]dp01[[2]],sol\[Epsilon]dp03[[2]]},PlotRange->{{0.5*10^-2,4*10^2},All},Frame->True,Axes->False,BaseStyle->{FontSize->10,FontFamily->"Arial"},AspectRatio->0.7,FrameLabel->{"t, seconds","Entropy generation rate, J/(K s)"},ImageSize->230,Joined->True,PlotStyle->{{Blue},{Darker[Green]},{Darker[Cyan]},Magenta},FrameTicks->{{Automatic,Automatic},{PowerTicks[True],PowerTicks[False]}}];
ShowLegend[ListLogLinearPlot[{sol\[Epsilon]dp01[[2]],sol\[Epsilon]dp03[[2]],sol\[Epsilon]dp1[[2]],sol\[Epsilon]dp3[[2]]},PlotRange->{{0.5*10^-2,4*10^2},All},Frame->True,Axes->False,BaseStyle->{FontSize->22,FontFamily->"Arial"},AspectRatio->0.7,FrameLabel->{"t, s","Entropy generation rate, J/(K s)"},ImageSize->700,Joined->True,PlotStyle->{{Blue},{Darker[Green]},{Darker[Cyan]},Magenta},Epilog->Inset[InsetP1,{Log[1*10^-1],6}],FrameTicks->{{Automatic,Automatic},{PowerTicks[True],PowerTicks[False]}}],Legend1]


(* ::Subsection:: *)
(*Cessation after biaxial elongation*)


(* ::Input::Initialization:: *)
sol\[Epsilon]dp3c=RDPmodel[0,mIn,MwRDP,MeIn,\[Tau]eIn,GN0In,\[Lambda]maxIn,{10^-2,2*10^2},Last[sol\[Epsilon]dp3]];


(* ::Input::Initialization:: *)
sol\[Epsilon]dp1c=RDPmodel[0,mIn,MwRDP,MeIn,\[Tau]eIn,GN0In,\[Lambda]maxIn,{10^-2,2*10^2},Last[sol\[Epsilon]dp1]];


(* ::Input::Initialization:: *)
sol\[Epsilon]dp03c=RDPmodel[0,mIn,MwRDP,MeIn,\[Tau]eIn,GN0In,\[Lambda]maxIn,{10^-2,2*10^2},Last[sol\[Epsilon]dp03]];


(* ::Input::Initialization:: *)
sol\[Epsilon]dp01c=RDPmodel[0,mIn,MwRDP,MeIn,\[Tau]eIn,GN0In,\[Lambda]maxIn,{10^-2,5*10^2},Last[sol\[Epsilon]dp01]];


(* ::Input::Initialization:: *)
Quiet[Needs["PlotLegends`"]];
Legend1={{{Graphics[{Blue, Thick,Line[{{0,0},{6,0}}]}],Text[Style["RDP,\!\(\*OverscriptBox[\(\(\\\ \)\*SubscriptBox[\(\[Epsilon]\), \(B\)]\), \(.\)]\)=0.01 \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",FontSize->16,FontFamily->"Arial",Blue]]},{Graphics[{Darker[Green], Thick,Line[{{0,0},{6,0}}]}],Text[Style["RDP,\!\(\*OverscriptBox[\(\(\\\ \)\*SubscriptBox[\(\[Epsilon]\), \(B\)]\), \(.\)]\)=0.03 \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",FontSize->16,FontFamily->"Arial",Darker[Green]]]},{Graphics[{Darker[Cyan], Thick,Line[{{0,0},{6,0}}]}],Text[Style["RDP,\!\(\*OverscriptBox[\(\(\\\ \)\*SubscriptBox[\(\[Epsilon]\), \(B\)]\), \(.\)]\)=0.1 \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",FontSize->16,FontFamily->"Arial",Darker[Cyan]]]},{Graphics[{Magenta, Thick,Line[{{0,0},{6,0}}]}],Text[Style["RDP,\!\(\*OverscriptBox[\(\(\\\ \)\*SubscriptBox[\(\[Epsilon]\), \(B\)]\), \(.\)]\)=0.3 \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",FontSize->16,FontFamily->"Arial",Magenta]]}},LegendShadow->False,LegendSize->{0.63,0.3},LegendPosition->{-0.3,0.35},LegendTextSpace->4,LegendBorder->False};
InsetP1=ListLogLinearPlot[{sol\[Epsilon]dp01c[[2]],sol\[Epsilon]dp03c[[2]]},PlotRange->{{0.5*10^-2,4*10^2},All},Frame->True,Axes->False,BaseStyle->{FontSize->10,FontFamily->"Arial"},AspectRatio->0.7,FrameLabel->{"t, seconds","Entropy generation rate, J/(K s)"},ImageSize->230,Joined->True,PlotStyle->{{Blue},{Darker[Green]},{Darker[Cyan]},Magenta},FrameTicks->{{Automatic,Automatic},{PowerTicks[True],PowerTicks[False]}}];
ShowLegend[ListLogLinearPlot[{sol\[Epsilon]dp01c[[2]],sol\[Epsilon]dp03c[[2]],sol\[Epsilon]dp1c[[2]],sol\[Epsilon]dp3c[[2]]},PlotRange->{{0.5*10^-2,4*10^2},All},Frame->True,Axes->False,BaseStyle->{FontSize->22,FontFamily->"Arial"},AspectRatio->0.7,FrameLabel->{"t, s","Entropy generation rate, J/(K s)"},ImageSize->700,Joined->True,PlotStyle->{{Blue},{Darker[Green]},{Darker[Cyan]},Magenta},Epilog->Inset[InsetP1,{Log[3*10^1],3}],FrameTicks->{{Automatic,Automatic},{PowerTicks[True],PowerTicks[False]}}],Legend1]


(* ::Subsection:: *)
(*Sudden reversal after biaxial elongation*)


(* ::Input::Initialization:: *)
sol\[Epsilon]dp3r=RDPmodel[-0.3,mIn,MwRDP,MeIn,\[Tau]eIn,GN0In,\[Lambda]maxIn,{10^-2,2*10^2},Last[sol\[Epsilon]dp3]];


(* ::Input::Initialization:: *)
sol\[Epsilon]dp1r=RDPmodel[-0.1,mIn,MwRDP,MeIn,\[Tau]eIn,GN0In,\[Lambda]maxIn,{10^-2,2*10^2},Last[sol\[Epsilon]dp1]];


(* ::Input::Initialization:: *)
sol\[Epsilon]dp03r=RDPmodel[-0.03,mIn,MwRDP,MeIn,\[Tau]eIn,GN0In,\[Lambda]maxIn,{10^-2,2*10^2},Last[sol\[Epsilon]dp03]];


(* ::Input::Initialization:: *)
sol\[Epsilon]dp01r=RDPmodel[-0.01,mIn,MwRDP,MeIn,\[Tau]eIn,GN0In,\[Lambda]maxIn,{10^-2,5*10^2},Last[sol\[Epsilon]dp01]];


(* ::Input::Initialization:: *)
Quiet[Needs["PlotLegends`"]];
Legend1={{{Graphics[{Blue, Thick,Line[{{0,0},{6,0}}]}],Text[Style["RDP,\!\(\*OverscriptBox[\(\(\\\ \)\*SubscriptBox[\(\[Epsilon]\), \(B\)]\), \(.\)]\)=-0.01 \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",FontSize->16,FontFamily->"Arial",Blue]]},{Graphics[{Darker[Green], Thick,Line[{{0,0},{6,0}}]}],Text[Style["RDP,\!\(\*OverscriptBox[\(\(\\\ \)\*SubscriptBox[\(\[Epsilon]\), \(B\)]\), \(.\)]\)=-0.03 \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",FontSize->16,FontFamily->"Arial",Darker[Green]]]},{Graphics[{Darker[Cyan], Thick,Line[{{0,0},{6,0}}]}],Text[Style["RDP,\!\(\*OverscriptBox[\(\(\\\ \)\*SubscriptBox[\(\[Epsilon]\), \(B\)]\), \(.\)]\)=-0.1 \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",FontSize->16,FontFamily->"Arial",Darker[Cyan]]]},{Graphics[{Magenta, Thick,Line[{{0,0},{6,0}}]}],Text[Style["RDP,\!\(\*OverscriptBox[\(\(\\\ \)\*SubscriptBox[\(\[Epsilon]\), \(B\)]\), \(.\)]\)=-0.3 \!\(\*SuperscriptBox[\(s\), \(-1\)]\)",FontSize->16,FontFamily->"Arial",Magenta]]}},LegendShadow->False,LegendSize->{0.63,0.3},LegendPosition->{-0.7,0.35},LegendTextSpace->4,LegendBorder->False};
InsetP1=ListLogLinearPlot[{sol\[Epsilon]dp01r[[2]],sol\[Epsilon]dp03r[[2]]},PlotRange->{{0.5*10^-2,4*10^2},All},Frame->True,Axes->False,BaseStyle->{FontSize->10,FontFamily->"Arial"},AspectRatio->0.7,FrameLabel->{"t, seconds","Entropy generation rate, J/(K s)"},ImageSize->230,Joined->True,PlotStyle->{{Blue},{Darker[Green]},{Darker[Cyan]},Magenta},FrameTicks->{{Automatic,Automatic},{PowerTicks[True],PowerTicks[False]}}];
ShowLegend[ListLogLinearPlot[{sol\[Epsilon]dp01r[[2]],sol\[Epsilon]dp03r[[2]],sol\[Epsilon]dp1r[[2]],sol\[Epsilon]dp3r[[2]]},PlotRange->{{0.5*10^-2,4*10^2},All},Frame->True,Axes->False,BaseStyle->{FontSize->22,FontFamily->"Arial"},AspectRatio->0.7,FrameLabel->{"t, s","Entropy generation rate, J/(K s)"},ImageSize->700,Joined->True,PlotStyle->{{Blue},{Darker[Green]},{Darker[Cyan]},Magenta},Epilog->Inset[InsetP1,{Log[1*10^-1],25}],FrameTicks->{{Automatic,Automatic},{PowerTicks[True],PowerTicks[False]}}],Legend1]


(* ::Subsection:: *)
(*Check for generally allowed values of the conformation tensor (As a function of the smallest eigenvalue of the conformation tensor)*)


(* ::Input::Initialization:: *)
soltr1=Table[{cxx,RDPmodelC2[1,cxx,MwRDP,MeIn,\[Tau]eIn,GN0In,\[Lambda]maxIn]},{cxx,0.01,0.3,0.005}];
soltrp8=Table[{cxx,RDPmodelC2[0.8,cxx,MwRDP,MeIn,\[Tau]eIn,GN0In,\[Lambda]maxIn]},{cxx,0.01,0.3,0.005}];
soltrp6=Table[{cxx,RDPmodelC2[0.6,cxx,MwRDP,MeIn,\[Tau]eIn,GN0In,\[Lambda]maxIn]},{cxx,0.01,0.3,0.005}];
soltrp5=Table[{cxx,RDPmodelC2[0.5,cxx,MwRDP,MeIn,\[Tau]eIn,GN0In,\[Lambda]maxIn]},{cxx,0.01,0.3,0.005}];


(* ::Input::Initialization:: *)
Quiet[Needs["PlotLegends`"]];
Clear[Legend0]
Legend0[size_,pos_]:={{{Graphics[{Black,Thick,Line[{{0,0},{3,0}}]}],Text[Style["trc=1",FontSize->18,FontFamily->"Arial",Black]]},{Graphics[{Blue,Dashed,Thick,Line[{{0,0},{3,0}}]}],Text[Style["trc=0.8",FontSize->18,FontFamily->"Arial",Blue]]},{Graphics[{Red,DotDashed,Thick,Line[{{0,0},{3,0}}]}],Text[Style["trc=0.6",FontSize->18,FontFamily->"Arial",Red]]},{Graphics[{Magenta,Dotted,Thick,Line[{{0,0},{3,0}}]}],Text[Style["trc=0.5",FontSize->18,FontFamily->"Arial",Magenta]]}},LegendShadow->False,LegendSize->size,LegendPosition->pos,LegendTextSpace->3,LegendBorder->False};
PA=ShowLegend[ListPlot[{soltr1,soltrp8,soltrp6,soltrp5},PlotRange->All,Frame->True,Axes->True,BaseStyle->22,Joined->True,AspectRatio->0.7,PlotStyle->{{Thick,Black},{Thick,Blue,Dashed},{Thick,Red,DotDashed},{Thick,Magenta,Dotted},{Thick,Darker[Green]},{Thick,Darker[Yellow]},{Thick,Brown},{Thick,Pink}},ImageSize->700,FrameLabel->{"\!\(\*SubscriptBox[\(\[Lambda]\), \(min\)]\)","Entropy generation rate, J/(K s)"},AspectRatio->0.7],Legend0[{0.4,0.33},{0.2,0.3}]]
