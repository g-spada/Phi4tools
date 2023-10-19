(* ::Package:: *)

(* ::Title:: *)
(*Phi4tools*)


BeginPackage["GSberveglieri`Phi4tools`"];


(* ::Subtitle:: *)
(*Public symbols*)


NickelIndex::usage = "NickelIndex[n, v3, v4, d] gives the Nickel Index for the d-th diagram of the n-point function with v3 cubic vertices and v4 quatic vertices.
NickelIndex[n, v3, v4] gives the Nickel Indices for the diagrams of the n-point function with v3 cubic vertices and v4 quatic vertices."


SymmetryFactorDiagram::usage = "SymmetryFactorDiagram[n, v3, v4, d] gives the tensorial symmetry factor for the d-th diagram of the n-point function with v3 cubic vertices and v4 quartic vertices.
SymmetryFactorDiagram[n, v3, v4] gives the list of tensorial symmetry factors for the diagrams of the n-point function with v3 cubic vertices and v4 quartic vertices."
Options[SymmetryFactorDiagram]={"Tensor"->"O(N)"}


VisualizeDiagram::usage = "VisualizeDiagram[n, v3, v4, d] draws d-th Feynman diagram of the n-point function with v3 cubic vertices and v4 quartic vertices
VisualizeDiagram[n, v3, v4] draws all Feynman diagrams of the n-point function with v3 cubic vertices and v4 quartic vertices"
Options[VisualizeDiagram]={"Substitutions"->"Nothing"}


DrawGraph::usage = "DrawGraph[str] draws the graph corresponding to the Nickel index written as input."


(*FindGraph::usage = "FindGraph[graph] finds the values {n, {v3,v4},{d}} corresponding to graph"*)


(*FindIndex::usage = "FindIndex[graph] finds the values {n, {v3,v4},{d}} corresponding to the Nickel index"*)


IntegrandDiagram::usage = "IntegrandDiagram[n, v3, v4, d] prints the integrand for the d-th Feynman diagram of the n-point function with v3 cubic vertices and v4 quartic vertices
IntegrandDiagram[n, v3, v4] prints the list of integrands for all the Feynman diagrams of the n-point function with v3 cubic vertices and v4 quartic vertices"
Options[IntegrandDiagram]={"Substitutions"->"Nothing","ExternalMomentum"->False}


ValueDiagram::usage = "ValueDiagram[n, v3, v4, d] gives the three-dimensional integrated value of the d-th Feynman diagram of the n-point function with v3 cubic vertices and v4 quartic vertices
ValueDiagram[n, v3, v4, d] the list of three-dimensional integrated values of the Feynman diagram for the n-point function with v3 cubic vertices and v4 quartic vertices"
Options[ValueDiagram]={"Derivative"->False}


InformationDiagram::usage = "InformationDiagram[n, v3, v4, d] gives details about the d-th diagram of the n-point function with v3 cubic vertices and v4 quartic vertices at zero external momentum
InformationDiagram[index] gives the details of the diagram having a specific Nickel index.
InformationDiagram[graph] gives the details of the diagram matching the provided graph."
Options[InformationDiagram]={"Substitutions"->"Nothing", "Tensor"->None, "ShowIntegrand"->False}


WriteExplicit::usage = "WriteExplicit[integrand] writes integrand in the form ready to be integrated in three dimensions.";
Options[WriteExplicit] = {"Simplification"->"Refine"}


DeriveAndWriteExplicit::usage = "DeriveAndWriteExplicit[integrand] derives integrand with respect to the external momentum squared p^2 and writes the result at p=0 in the form ready to be integrated in three dimensions."


CountLoops::usage = "CountLoops[integrand] gives the number of effective loop of a integrand, i.e. the number of internal momenta on which the integrand has to be integrated on.";


MomVars::usage = "MomVars[expr] gives the list of Momentum variables appearing in expr."


(* ::Subsubsection:: *)
(*Diagram elements*)


Propagator::usage = "Propagator[q] represents the scalar propagator of momentum q."
\[ScriptCapitalG]::usage = "Propagator[q] represents the scalar propagator of momentum q."


BubbleSubdiagram::usage = "BubbleSubdiagram[q] represents the one-loop bubble subdiagram with total external momentum q."
\[ScriptCapitalB]::usage = "BubbleSubdiagram[q] represents the one-loop bubble subdiagram with total external momentum q."


SunsetSubdiagram::usage = "SunsetSubdiagram[q] represents the regularized two-loop sunset subdiagram with total external momentum q."
\[ScriptCapitalS]::usage = "SunsetSubdiagram[q] represents the regularized two-loop sunset subdiagram with total external momentum q."


TriangleSubdiagram::usage = "TriangleSubdiagram[q1,q2,q3] represents the one-loop triangle subdiagram."
\[ScriptCapitalT]::usage = "TriangleSubdiagram[q1,q2,q3] represents the one-loop triangle subdiagram."


SquareSubdiagram::usage = "SquareSubdiagram[q1,q2,q3,q4] represents the one-loop square subdiagram."
\[ScriptCapitalQ]::usage = "SquareSubdiagram[q1,q2,q3,q4] represents the one-loop square subdiagram."


TadSunsetSubdiagram::usage = "TadSunsetSubdiagram[] represents a tadpole-like subdiagram containing the regularized sunset."
\[ScriptT]\[ScriptCapitalS]::usage = "TadSunsetSubdiagram[] represents a tadpole-like subdiagram containing the regularized sunset."


TadTriangleBubblesSubdiagram::usage = "TadTriangleBubblesSubdiagram[] represents a tadpole-like subdiagram containing a triangle and two bubbles."
\[ScriptT]\[ScriptCapitalT]\[ScriptCapitalB]::usage = "TadTriangleBubblesSubdiagram[] represents a tadpole-like subdiagram containing a triangle and two bubbles."


(* ::Subsubsection:: *)
(*Momenta representation*)


Momentum::usage = "Momentum[i] represents the i-th internal momentum variable.
Momentum[i,sub] represents the sub component of the i-th internal momentum variable."
\[ScriptQ]::usage = "Momentum[i] represents the i-th internal momentum variable.
Momentum[i,sub] represents the sub component of the i-th internal momentum variable."


ExternalMomentum::usage = "ExternalMomentum[i] represents the i-th external momentum variable.
ExternalMomentum[i,sub] represents the sub component of the i-th external momentum variable."
\[ScriptP]::usage = "ExternalMomentum[i] represents the i-th external momentum variable.
ExternalMomentum[i,sub] represents the sub component of the i-th external momentum variable."


(* ::Subsubsection:: *)
(*O(N) and cubic factors*)


NComponents::usage = "NComponents[] represents the number N of field components."
\[CapitalNu]::usage = "NComponents[] represents the number N of field components."


XCubicRatio::usage = "XCubicRatio[] represents the ratio between the cubic and O(N)-symmetric coupling constants."
\[CapitalChi]::usage = "XCubicRatio[] represents the ratio between the cubic and O(N)-symmetric coupling constants."


(* ::Subtitle:: *)
(*Private*)


Begin["`Private`"];


(* ::Chapter:: *)
(*Import from .txt*)


(* ::Subsection::Closed:: *)
(*Define internal momenta and their Momentum representation*)


(*qaux = {q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12};*)


qauxify[i_Integer /; i > 0] :=
    ToExpression[StringJoin[ToString[Context[qauxify]],"q", ToString[i]]]


qaux = Table[qauxify[i],{i,1,20}]


(*momentumfy[qq_Symbol] :=
    Module[{str = ToString[qq]},
        If[StringTake[str, 1] == "q",
            Momentum[ToExpression[StringTake[str, {2, -1}]]]
        ]
    ]*)


(*momentumfy[qqs_Subscript] :=
    Module[{str = ToString[qqs[[1]]]},
        If[StringTake[str, 1] == "q",
            Momentum[ToExpression[StringTake[str, {2, -1}]], qqs[[2]]
                ]
        ]
    ]*)


(* ::Text:: *)
(*Define external momentum*)


pextify[i_Integer /; i > 0] :=
    ToExpression[StringJoin[ToString[Context[pextify]],"p", ToString[i]]]


pextlist = Table[pextify[i],{i,10}];


(* ::Subsection:: *)
(*Functions*)


subUndirUNDO[obj_]:=obj/. {a_,b_}:> UndirectedEdge[a,b]


(*changeEO={e1->o1,e2->o2,e3->o3,e4->o4};*)


splitsimport[imp_]:=Map[StringSplit[ToString@#,"\t"]&,imp]


formatsimport[table_]:={ToExpression[table[[1]],StandardForm],table[[2]],ToExpression[table[[3]],StandardForm],ToExpression[table[[4]],StandardForm]}


(*importtxtEdges[text_]:=Module[{tmp1},
tmp1=Import[text,"List"];
formatsimport/@splitsimport@tmp1
]*)


(*importEdges1PI[text_]:=Module[{tmp},
tmp=importtxtEdges[text];
Transpose[{tmp[[All,-1]],(subUndirUNDO/@tmp[[All,-2]])/.changeEO}]
]*)


importLables1PI[text_]:=Module[{tmp},
tmp=importtxtEdges[text];
tmp[[All,2]]
]


(* ::Text:: *)
(*Alternative import method that explicitly converts external points to internal private variables*)


extVertexReplaceNames = {
	"e1"->ToString[Context[extVertexReplaceNames]]<>"o1",
	"e2"->ToString[Context[extVertexReplaceNames]]<>"o2",
	"e3"->ToString[Context[extVertexReplaceNames]]<>"o3",
	"e4"->ToString[Context[extVertexReplaceNames]]<>"o4"}

importtxtEdges[textfile_]:=Module[{str,strtab,nid,nlab,edges,symf},
	str=Import[textfile,"Text"];
	strtab=Partition[StringSplit[str,{"\n","\t"}],4];
	nid=ToExpression[strtab[[All,1]]];
	nlab=strtab[[All,2]];
	edges=ToExpression[StringReplace[strtab[[All,3]],extVertexReplaceNames]];
	symf=ToExpression[strtab[[All,4]]];
	Transpose[{nid,nlab,edges,symf}]
]

importEdges1PI[filename_]:=Module[{tmp},
	If[
		FileExistsQ[filename],
		tmp=importtxtEdges[filename]; Transpose[{tmp[[All,-1]],(subUndirUNDO/@tmp[[All,-2]])}],
		Message[VisualizeDiagram::nvert]
	]
]


(* ::Subsection:: *)
(*Import edges*)


(*Table[fourpt1PI[0,oo4]=importEdges1PI[FileNameJoin[{PacletManager`PacletResource["GSberveglieri/Phi4tools","EDGES"],"4pt0"<>ToString[oo4]<>".txt"}]];,{oo4,2,8}];*)


(*Table[twopt1PI[0,oo4]=importEdges1PI[FileNameJoin[{PacletManager`PacletResource["GSberveglieri/Phi4tools","EDGES"],"2pt0"<>ToString[oo4]<>".txt"}]];,{oo4,2,8}];*)


(*Table[zeropt1PI[0,oo4]=importEdges1PI[FileNameJoin[{PacletManager`PacletResource["GSberveglieri/Phi4tools","EDGES"],"0pt0"<>ToString[oo4]<>".txt"}]];,{oo4,2,8}];*)


(* ::Subsection:: *)
(*Import label*)


(*oo3=0;
Table[nickelIndices[nn,oo3,oo4]=importLables1PI[FileNameJoin[{PacletManager`PacletResource["GSberveglieri/Phi4tools","EDGES"],ToString[nn]<>"pt"<>ToString[oo3]<>ToString[oo4]<>".txt"}]];(*Print[nn,oo3,oo4];*),{nn,0,4,2},{oo4,2,8}];*)


(* ::Subsection:: *)
(*Import symmetry factors*)


(*importsymno[n_,o_]:=Module[{facraw},
facraw=Import[FileNameJoin[{PacletManager`PacletResource["GSberveglieri/Phi4tools","FACTORS"],"cubfactorG"<>ToString@n<>"o"<>ToString@o<>".txt"}],"Table"];
If[n==0||n==2,ToExpression[facraw[[All,-1]]],
ToExpression/@(facraw[[All,-4;;-1]])]
]*)


(* ::Text:: *)
(*Alternative import method that explicitly converts {n,x} to internal private variables*)


nxFactorsReplaceNames = {"n"->ToString[Context[nxFactorsReplaceNames]]<>"n","x"->ToString[Context[nxFactorsReplaceNames]]<>"x"};
importtxtCubicFactors[textfile_]:=Module[{str,strnwl,strtab,cubf},
	str=Import[textfile,"Text"];
	strnwl=StringSplit[str,{"\n"}];
    strtab=StringSplit[strnwl,{"\t"}];
	cubf=ToExpression[StringReplace[nxFactorsReplaceNames]/@strtab[[All,3;;-1]]];
     If[Dimensions[cubf][[2]]==1, Return[Flatten[cubf]],Return[cubf]]
]

importsymno[n_,o_]:=Module[{filename},
	filename = FileNameJoin[{PacletManager`PacletResource["GSberveglieri/Phi4tools","FACTORS"],"cubfactorG"<>ToString@n<>"o"<>ToString@o<>".txt"}];
	If[
		FileExistsQ[filename],
		importtxtCubicFactors[filename],
		Message[SymmetryFactorDiagram::nvert]
	]
]


(*Table[symfactorGno[0,ii]=importsymno[0,ii],{ii,4,8}];*)


(*Table[symfactorGno[2,ii]=importsymno[2,ii],{ii,2,8}];*)


(*Table[symfactorGno[4,ii]=importsymno[4,ii],{ii,1,8}];*)


ImportedSymmetryFactorQ[n_, o4_] := SameQ[ Head[symFactorData[n,o4]], List]


symfactorGno[n_,o4_]:=Module[{},
	If[
		Not[ImportedSymmetryFactorQ[n,o4]],
		symFactorData[n,o4]=importsymno[n,o4]
	];
	Return[symFactorData[n,o4]]
]


(* ::Subsection:: *)
(*Import integrands*)


(*
importtxtIntegrands[text_]:=Module[{tmp1},
	tmp1=Import[text,"List"];
	formatsimport/@splitsimport@tmp1
]
*)


(* ::Text:: *)
(*Alternative import method that explicitly converts momentum variables to internal private variables*)


momentaReplaceNames = Table[StringJoin["q", ToString[i]] -> StringJoin["GSberveglieri`Phi4tools`Private`q", ToString[i]],{i,1,9}]
symbolsReplaceNames = Table[symb -> "GSberveglieri`Phi4tools`Private`"<>symb,{symb,{"prop","bubble","sunset","triangle","square","tadSunset","tadTrianBub"}}]
importtxtIntegrands[textfile_]:=Module[{str,strtab,nid,nlab,loop,integr},
	str=Import[textfile,"Text"];
	strtab=Partition[StringSplit[str,{"\n","\t"}],4];
	nid=ToExpression[strtab[[All,1]]];
	nlab=strtab[[All,2]];
    loop=ToExpression[strtab[[All,3]]];
	integr=ToExpression[StringReplace[strtab[[All,4]],Join[momentaReplaceNames,symbolsReplaceNames]]];
	Transpose[{nid,nlab,loop,integr}]
]


IntegrandFileName[nn_,oo4_]:=FileNameJoin[{PacletManager`PacletResource["GSberveglieri/Phi4tools","INTEGRANDS"],"int"<>ToString[nn]<>"pt0"<>ToString[oo4]<>".txt"}]


IntegrandFileExistsQ[nn_,oo4_]:=FileExistsQ[IntegrandFileName[nn,oo4]]


(*nn=2;*)


(*Table[integrand1PIno[nn,oo4]=importtxtIntegrands[
	FileNameJoin[{PacletManager`PacletResource["GSberveglieri/Phi4tools","INTEGRANDS"],"int"<>ToString[nn]<>"pt0"<>ToString[oo4]<>".txt"}]],{oo4,2,8}];*)


(*nn=4;*)


(*Table[integrand1PIno[nn,oo4]=importtxtIntegrands[
	FileNameJoin[{PacletManager`PacletResource["GSberveglieri/Phi4tools","INTEGRANDS"],"int"<>ToString[nn]<>"pt0"<>ToString[oo4]<>".txt"}]];,{oo4,2,8}];*)


(*nn=0;*)


(*Table[integrand1PIno[nn,oo4]=importtxtIntegrands[
	FileNameJoin[{PacletManager`PacletResource["GSberveglieri/Phi4tools","INTEGRANDS"],"int"<>ToString[nn]<>"pt0"<>ToString[oo4]<>".txt"}]],{oo4,3,8}];*)


(*Clear@nn;*)


(* ::Subsection:: *)
(*Import values*)


importvalno[n_,o_]:=Module[{filename,valraw},
	filename = FileNameJoin[{PacletManager`PacletResource["GSberveglieri/Phi4tools","VALUES"],"G"<>n<>"o"<>ToString@o<>".txt"}];
	If[
		FileExistsQ[filename],
		valraw=Import[filename,"Table"]; Transpose[{valraw[[All,1]],Table[Around[ToExpression@valraw[[i,3]],valraw[[i,4]]],{i,Length@valraw}]}],
		Message[ValueDiagram::nvert]
	]
]


(*Table[resultGIno["0",ii]=importvalno["0",ii],{ii,4,8}];*)


(*Table[resultGIno["2",ii]=importvalno["2",ii],{ii,2,8}];*)


(*Table[resultGIno["2D",ii]=importvalno["2d",ii],{ii,2,7}];*)


(*Table[resultGIno["4",ii]=importvalno["4",ii],{ii,2,8}];*)


{resultGIno["2",0],resultGIno["2",1]}={{{1,-1}},{{1,0}}};


{resultGIno["2D",0],resultGIno["2D",1]}={{{1,-1}},{{1,0}}};


{resultGIno["4",0],resultGIno["4",1]}={{{1,0}},{{1,1/3}}};


ImportedValueQ[str_, o4_] := SameQ[ Head[valueData[str,o4]], List]


resultGIno[str_,o4_]:=Module[{strR = StringReplace[str,"2D"->"2d"]},
	If[
		Not[ImportedValueQ[str,o4]],
		valueData[str,o4]=importvalno[strR,o4]
	];
	Return[valueData[str,o4]]
]


(* ::Chapter:: *)
(*Nickel labels*)


(* ::Section::Closed:: *)
(*Basic functions*)


(* ::Subsubsection::Closed:: *)
(*Basic functions*)


(* ::Text:: *)
(*I kept the labels a more handy forms at the last moment will convert the in strings.*)
(*Legend: the label is in the form of a list of a string, the "1/2" have to be substituted with "e" and the "b" with "|"*)


(* ::Text:: *)
(*Some useful functions*)


nearvertex[list_,vx_]:=Cases[list,vx\[UndirectedEdge]_|_\[UndirectedEdge]vx]/.vx\[UndirectedEdge]a_:>a/.a_\[UndirectedEdge]vx:>a


selecttoassign[listv_]:=DeleteCases[listv,Alternatives@@Join[Select[listv,NumericQ],{e}]]


assignrule[list1_,list2_]:=Table[list1[[i]]-> list2[[i]],{i,Length@list1}]


delteb[list_]:=DeleteCases[list,b]


(* ::Text:: *)
(*Identify the change of vertex, i.e. add the bars in the right spot*)


addBarstolist[list_,v_]:=Module[{i=1,tmp=list,ll},
ll=Length@list;
While[i<(ll-1+v),If[Length@tmp=!=i,If[(tmp[[i]]/.a_\[UndirectedEdge]b_:> a)==(tmp[[i+1]]/.a_\[UndirectedEdge]b_:> a),i++(*;Print[tmp]*),tmp=Insert[tmp,b,i+1];i=i+2(*;Print[tmp]*)],tmp=Insert[tmp,b,i+1];i=i+2(*;Print[tmp]*)]];
tmp=Join[tmp,ConstantArray[b,-Length@tmp+ll+v-1]];
tmp]


(* ::Subsubsection::Closed:: *)
(*Brute force*)


(* ::Text:: *)
(*labelNickelBruteforce: brute force function to write the Nickel label of a given diagram.*)
(*All the possible numerations of the vertices are tried and the first in numerical order (and correct one) of the possible labels is chosen.*)
(*((Since the code is pretty simple it should be quite solid, but it is slow and becomes too time consuming for the orders ~ 7))*)


labelNickelBruteForce[graph_]:=Module[{tmp,vert,assign0,rules,assigned0,assigned,labels},
tmp=graph/.{o1-> e,o2-> e,o3-> e,o4-> e,0-> f};
vert=DeleteCases[VertexList@tmp,e];
assign0=Permutations@Range[0,Length@vert-1];
rules=Table[assignrule[vert,assign0[[i]]],{i,Length@assign0}];
assigned0=Table[tmp/.rules[[i]],{i,Length@assign0}];
assigned=Select[assigned0,MemberQ[#,0\[UndirectedEdge]1|1\[UndirectedEdge]0]&];
(*Print[tmp,vert,assign0,rules,assigned];*)
labels=Table[addBarstolist[Sort[(Sort/@assigned[[i]])/.e-> 1/2],Length@vert]/.a_\[UndirectedEdge]b_:> b,{i,Length@assigned}];
labels[[Position[delteb/@labels,SortBy[delteb/@labels,First][[1]]][[1,1]]]]
]


(* ::Input:: *)
(*(*AbsoluteTiming[labelbruteforcev2/@(twopt1PI[0,6][[All,2]])]\[LeftDoubleBracket]1\[RightDoubleBracket]*)*)


(* ::Input:: *)
(*(*Table[AbsoluteTiming[labelbruteforcev2/@(fourpt1PI[0,i][[All,2]])]\[LeftDoubleBracket]1\[RightDoubleBracket],{i,5,6}]*)*)


(* ::Section::Closed:: *)
(*Labelling diagrams*)


(* ::Subsubsection::Closed:: *)
(*General functions*)


(* ::Text:: *)
(*addBars: adds the bars of the Nickel labels.*)
(*It starts from the lists of the vertices at the left and at the right of the edges of the graphs and from the  number of vertices*)


addBars[list_,list2_,v_]:=Module[{skips,positioninsert,llbinside},
skips=Flatten@{0,list[[2;;-1]]-list[[1;;-2]]};
positioninsert=Flatten[{Position[skips,1],Join[Position[skips,2],Position[skips,2]],Join[Position[skips,3],Position[skips,3],Position[skips,3]]},1];
llbinside=Length@positioninsert;
Join[Insert[list2,b,positioninsert],ConstantArray[b,v-llbinside]]
]


(* ::Subsubsection::Closed:: *)
(*labellingG0Tree*)


(* ::Text:: *)
(*labellingG0Tree: a tree graph of all possible assignations of the vertices 0, 1, and the adjacent ones is scanned and all the possible labels are in output. *)
(*Works for 0-point functions*)


labellingG0Tree[list_]:=Module[{list1,vertices,llv,comvertexedouble,rule0,list2,vertex0,configuration0,ruleA,j,priority0ext,jj,list3,vertex1all,vertex1numb,vertex1,lastassigned,configuration1,ruleB,priority1ext,names,jjj,assignement,remainedv,lrv,lastassignedC,ruleC,y,assigned,label,sortass0,sortass},
list1=Sort[list/.{o1-> e,o2-> e,o3-> e,o4-> e,0-> f}];
(*Print[list1];*)
vertices=VertexList@list1;
llv=Length@vertices;
rule0=Table[vertices[[i]]-> 0,{i,llv}];
j=1;
While[j<= Length@rule0,
list2=Sort[Sort/@(list1/.rule0[[j]])];
vertex0=selecttoassign[nearvertex[list2,0]];
(*Print[list2,vertex0];*)
configuration0=ReverseSortBy[Tally@vertex0,#[[2]]&];
ruleA[j]=Which[configuration0[[All,2]]=={3,1},
(*Print["OK at 0 2 fix"];*){{configuration0[[1,1]]-> 1,configuration0[[2,1]]-> 2}},

configuration0[[All,2]]=={2,2},
(*Print["OK at 0 2 poss"];*){{configuration0[[1,1]]-> 1,configuration0[[2,1]]-> 2},{configuration0[[2,1]]-> 1,configuration0[[1,1]]-> 2}},

configuration0[[All,2]]=={2,1,1},
(*Print["OK at 0 3 poss"];*){{configuration0[[1,1]]-> 1,configuration0[[2,1]]-> 2,configuration0[[3,1]]-> 3},{configuration0[[1,1]]-> 1,configuration0[[3,1]]-> 2,configuration0[[2,1]]-> 3}},

configuration0[[All,2]]=={2,1},
(*Print["OK at 0 2 fix"];*){{configuration0[[1,1]]-> 1,configuration0[[2,1]]-> 2}},

configuration0[[All,2]]=={1,1,1},
(*Print["OK at 0 4 poss"];*)Table[assignrule[Permutations[vertex0][[i]],Range[1,3]],{i,6}],

configuration0[[All,2]]=={1,1,1,1},
(*Print["OK at 0 4 poss"];*)Table[assignrule[Permutations[vertex0][[i]],Range[1,4]],{i,24}]
];
(*Print[{rule0\[LeftDoubleBracket]j\[RightDoubleBracket],ruleA[j]}];*)j++];
j=1;
While[j<= Length@rule0,jj=1;While[jj<= Length@ruleA[j],
list3=Sort[Sort/@(list1/.(rule0[[j]])/.(ruleA[j][[jj]]))];
(*Print[list3];*)
vertex1all=DeleteCases[nearvertex[list3,1],0];
vertex1numb=Sort@Select[vertex1all/.(e-> 1/2),NumericQ];
vertex1=selecttoassign[vertex1all];
(*Print[vertex1];*)

lastassigned[j]=((ruleA[j][[jj,-1]])/.(a_-> b_):> b);
(*Print[lastassigned[j]];*)

configuration1=ReverseSortBy[Tally@vertex1,#[[2]]&];
ruleB[j,jj]=Which[configuration1[[All,2]]=={2,1},
(*Print["OK at 1 2 fix"];*){{configuration1[[1,1]]-> lastassigned[j]+1,configuration1[[2,1]]-> lastassigned[j]+2}},

configuration1[[All,2]]=={1,1,1},
(*Print["OK at 1 3 poss"];*)Table[assignrule[Permutations[vertex1][[i]],lastassigned[j]+Range[1,3]],{i,6}],

configuration1[[All,2]]=={1,1},
(*Print["OK at 1 2 poss"];*)Table[assignrule[Permutations[vertex1][[i]], lastassigned[j]+Range[1,2]],{i,2}],

Length@configuration1== 1,
(*Print["OK at 1 1 fix"];*){{configuration1[[1,1]]-> lastassigned[j]+1}},

Length@configuration1== 0,{}
];
(*Print[{j,jj,rule0\[LeftDoubleBracket]j\[RightDoubleBracket],ruleA[j]\[LeftDoubleBracket]jj\[RightDoubleBracket],ruleB[j,jj]}];*)jj++];j++];
names={};
j=1;
While[j<= Length@rule0,jj=1;While[jj<= Length@ruleA[j],jjj=1;
If[ruleB[j,jj]=={},
assignement=Flatten@{rule0[[j]],ruleA[j][[jj]]};
(*Print[assignement];*)
remainedv=selecttoassign[VertexList[(list1/.assignement)]];
lrv=Length@remainedv;
(*Print[lastassigned[j]];*)
ruleC=Table[assignrule[Permutations[remainedv][[i]],lastassigned[j]+Range[1,lrv]],{i,lrv!}];
y=1;While[y<= Length@ruleC,
assigned=list1/.assignement/.ruleC[[y]];
sortass0=Sort[(Sort/@assigned)/.e-> 1/2];
label=addBars[sortass0/.a_\[UndirectedEdge]b_:> a,sortass0/.a_\[UndirectedEdge]b_:> b,Length@VertexList@list1];
AppendTo[names,label];y++
];];

While[jjj<= Length@ruleB[j,jj],
assignement=Flatten@{rule0[[j]],ruleA[j][[jj]],ruleB[j,jj][[jjj]]};
(*Print[assignement];*)
remainedv=selecttoassign[VertexList[(list1/.assignement)]];
lrv=Length@remainedv;
lastassignedC[j,jj,jjj]=((ruleB[j,jj][[jjj,-1]])/.(a_-> b_):> b);
(*Print[lastassignedC[j,jj,jjj]];*)
ruleC=Table[assignrule[Permutations[remainedv][[i]],lastassignedC[j,jj,jjj]+Range[1,lrv]],{i,lrv!}];
y=1;While[y<= Length@ruleC,
assigned=list1/.assignement/.ruleC[[y]];
sortass=Sort[(Sort/@assigned)/.e-> 1/2];
label=addBars[sortass/.a_\[UndirectedEdge]b_:> a,sortass/.a_\[UndirectedEdge]b_:> b,Length@VertexList@list1];
AppendTo[names,label];y++
];jjj++];jj++];j++];
names
]


(* ::Text:: *)
(*The first in numerical order (and correct one) of the possible labels of labellingTree01 is chosen.*)


labelNickelG0[list_]:=Module[{tmp},
tmp=labellingG0Tree@list;
tmp[[Position[delteb/@tmp,SortBy[delteb/@tmp,First][[1]]][[1,1]]]]]


(* ::Subsubsection::Closed:: *)
(*labellingTree G2 & G4*)


(* ::Text:: *)
(*labellingTree: a tree graph of all possible assignations of the vertices 0, 1, and the adjacent ones is scanned and all the possible labels are in output. *)
(*Works for 2- and 4-point functions*)


labellingTree[list_]:=Module[{list1,vertexext,comvertexe,llcve,comvertexedouble,rule0,list2,vertex0,configuration0,ruleA,j,priority0ext,jj,list3,vertex1all,vertex1numb,vertex1,lastassigned,configuration1,ruleB,priority1ext,names,jjj,assignement,remainedv,lrv,lastassignedC,ruleC,y,assigned,label,sortass0,sortass},
list1=Sort[list/.{o1-> e,o2-> e,o3-> e,o4-> e,0-> f}];
(*Print[list1];*)
vertexext=Cases[list1,e\[UndirectedEdge]_]/.e\[UndirectedEdge]a_:>a;
(*Print[vertexext];*)
comvertexe=Commonest[vertexext];llcve=Length@comvertexe;
rule0=If[llcve==1,{comvertexe[[1]]-> 0},
comvertexedouble=Flatten@Position[Table[DuplicateFreeQ@nearvertex[list1,comvertexe[[i]]],{i,llcve}],False];
If[Length@comvertexedouble==1,{comvertexe[[comvertexedouble[[1]]]]-> 0},
If[comvertexedouble=!={},Table[comvertexe[[comvertexedouble[[i]]]]-> 0,{i,Length@comvertexedouble}],
Table[comvertexe[[i]]-> 0,{i,Length@comvertexe}]]]];
(*Print[comvertexe,Length@rule0,rule0];*)

(*2nd level of scan*)
j=1;
While[j<= Length@rule0,
list2=Sort[Sort/@(list1/.rule0[[j]])];
vertex0=selecttoassign[nearvertex[list2,0]];
(*Print[list2,vertex0];*)
configuration0=ReverseSortBy[Tally@vertex0,#[[2]]&];
ruleA[j]=Which[configuration0[[All,2]]=={2,1},
(*Print["OK at 0 2 fix"];*){{configuration0[[1,1]]-> 1,configuration0[[2,1]]-> 2}},

configuration0[[All,2]]=={1,1,1},
(*Print["OK at 0 3 poss"];*)priority0ext=Select[vertex0,MemberQ[vertexext,#]&];Which[Length@priority0ext==1,Table[Join[{priority0ext[[1]]-> 1},assignrule[Permutations[DeleteCases[vertex0,priority0ext[[1]]]][[i]],Range[2,3]]],{i,2}],
Length@priority0ext==2,Join[Table[Join[{priority0ext[[1]]-> 1},assignrule[Permutations[{priority0ext[[2]],DeleteCases[vertex0,Alternatives@@priority0ext][[1]]}][[i]],Range[2,3]]],{i,2}],
Table[Join[{priority0ext[[2]]-> 1},assignrule[Permutations[{priority0ext[[1]],DeleteCases[vertex0,Alternatives@@priority0ext][[1]]}][[i]],Range[2,3]]],{i,2}]],Length@priority0ext==0||Length@priority0ext==3,Table[assignrule[Permutations[vertex0][[i]],Range[1,3]],{i,6}]],

configuration0[[All,2]]=={1,1},
(*Print["OK at 0 2 poss"];*)priority0ext=Select[vertex0,MemberQ[vertexext,#]&];If[Length@priority0ext==1,{{priority0ext[[1]]-> 1,DeleteCases[vertex0,priority0ext[[1]]][[1]]-> 2}},Table[assignrule[Permutations[vertex0][[i]],Range[1,2]],{i,2}]],

Length@configuration0== 1,
(*Print["OK at 0 1 fix"];*){{configuration0[[1,1]]-> 1}}
];
(*Print[{rule0\[LeftDoubleBracket]j\[RightDoubleBracket],ruleA[j]}];*)j++];

(*3rd level of scan*)
j=1;
While[j<= Length@rule0,jj=1;While[jj<= Length@ruleA[j],
list3=Sort[Sort/@(list1/.(rule0[[j]])/.(ruleA[j][[jj]]))];
(*Print[list3];*)
vertex1all=DeleteCases[nearvertex[list3,1],0];
vertex1numb=Sort@Select[vertex1all/.(e-> 1/2),NumericQ];
vertex1=selecttoassign[vertex1all];
(*Print[vertex1];*)

lastassigned[j]=((ruleA[j][[jj,-1]])/.(a_-> b_):> b);
(*Print[lastassigned[j]];*)

configuration1=ReverseSortBy[Tally@vertex1,#[[2]]&];
ruleB[j,jj]=Which[configuration1[[All,2]]=={2,1},
(*Print["OK at 1 2 fix"];*){{configuration1[[1,1]]-> lastassigned[j]+1,configuration1[[2,1]]-> lastassigned[j]+2}},

configuration1[[All,2]]=={1,1,1},
(*Print["OK at 1 3 poss"];*)Table[assignrule[Permutations[vertex1][[i]],lastassigned[j]+Range[1,3]],{i,6}],

configuration1[[All,2]]=={1,1},
(*Print["OK at 1 2 poss"];*)Table[assignrule[Permutations[vertex1][[i]], lastassigned[j]+Range[1,2]],{i,2}],

Length@configuration1== 1,
(*Print["OK at 1 1 fix"];*){{configuration1[[1,1]]-> lastassigned[j]+1}},

Length@configuration1== 0,{}
];
(*Print[{j,jj,rule0\[LeftDoubleBracket]j\[RightDoubleBracket],ruleA[j]\[LeftDoubleBracket]jj\[RightDoubleBracket],ruleB[j,jj]}];*)jj++];j++];

(*4th level of scan*)
names={};
j=1;
While[j<= Length@rule0,jj=1;While[jj<= Length@ruleA[j],jjj=1;
If[ruleB[j,jj]=={},
assignement=Flatten@{rule0[[j]],ruleA[j][[jj]]};
(*Print[assignement];*)
remainedv=selecttoassign[VertexList[(list1/.assignement)]];
lrv=Length@remainedv;
(*Print[lastassigned[j]];*)
ruleC=Table[assignrule[Permutations[remainedv][[i]],lastassigned[j]+Range[1,lrv]],{i,lrv!}];
y=1;While[y<= Length@ruleC,
assigned=list1/.assignement/.ruleC[[y]];
sortass0=Sort[(Sort/@assigned)/.e-> 1/2];
label=addBars[sortass0/.a_\[UndirectedEdge]b_:> a,sortass0/.a_\[UndirectedEdge]b_:> b,Length@VertexList@list1-1];
AppendTo[names,label];y++
];];

While[jjj<= Length@ruleB[j,jj],
assignement=Flatten@{rule0[[j]],ruleA[j][[jj]],ruleB[j,jj][[jjj]]};
(*Print[assignement];*)
remainedv=selecttoassign[VertexList[(list1/.assignement)]];
lrv=Length@remainedv;
lastassignedC[j,jj,jjj]=((ruleB[j,jj][[jjj,-1]])/.(a_-> b_):> b);
(*Print[lastassignedC[j,jj,jjj]];*)
ruleC=Table[assignrule[Permutations[remainedv][[i]],lastassignedC[j,jj,jjj]+Range[1,lrv]],{i,lrv!}];
y=1;While[y<= Length@ruleC,
assigned=list1/.assignement/.ruleC[[y]];
sortass=Sort[(Sort/@assigned)/.e-> 1/2];
label=addBars[sortass/.a_\[UndirectedEdge]b_:> a,sortass/.a_\[UndirectedEdge]b_:> b,Length@VertexList@list1-1];
AppendTo[names,label];y++
];jjj++];jj++];j++];
names
]


(* ::Text:: *)
(*labelNickel: the first in numerical order (and correct one) of the possible labels of labellingTree is chosen.*)


labelNickel[list_]:=Module[{tmp},
tmp=labellingTree@list;
tmp[[Position[delteb/@tmp,SortBy[delteb/@tmp,First][[1]]][[1,1]]]]]


(* ::Section::Closed:: *)
(*Ordering, labelling and finding*)


(* ::Subsection::Closed:: *)
(*Ordering*)


(* ::Text:: *)
(*orderingn1PI: gives the map between the list of diagrams that we have and the ordered sequence with the Nickel labels*)


orderingn1PI[n_,o_]:=Module[{tmp},
If[n==0,
tmp=zeropt1PI[0,o][[All,2]];OrderingBy[delteb/@labelNickelG0/@tmp,First],
tmp=Which[n==2,
twopt1PI[0,o][[All,2]],
n==4,
fourpt1PI[0,o][[All,2]]];
OrderingBy[delteb/@labelNickel/@tmp,First]]]


(* ::Text:: *)
(*Checks*)


(* ::Input:: *)
(*(*orderingn1PI[4,6]==order\[CapitalGamma]4o[6]*)*)


(* ::Input:: *)
(*(*orderingn1PI[4,7]==order\[CapitalGamma]4o[7]*)*)


(* ::Input:: *)
(*(*orderingn1PI[2,7]==order\[CapitalGamma]2o[7]*)*)


(* ::Input:: *)
(*(*orderingn1PI[0,6]*)*)


stringNickel[list_]:=StringJoin[ToString/@(list/.1/2-> "e"/.b-> "|")]


(* ::Subsection::Closed:: *)
(*Labelling*)


(* ::Text:: *)
(*labelsNickelGo: gives the label for all the diagrams at a given order for \[CapitalGamma]^(2) & \[CapitalGamma]^(4)*)


labelsNickelGo[pt_,o_]:=Module[{tmplist,labelsraw,order,labelsordered},
Which[pt==2,
tmplist=twopt1PI[0,o];,
pt==4,
tmplist=fourpt1PI[0,o];,
pt=!=2&&pt=!=4, "It works just for \!\(\*SuperscriptBox[\(\[CapitalGamma]\), \((2)\)]\) & \!\(\*SuperscriptBox[\(\[CapitalGamma]\), \((4)\)]\)"
];
labelsraw=labelNickel/@(tmplist[[All,2]]);
(*Print[labelsraw];*)
order=OrderingBy[delteb/@labelsraw,First];
labelsordered=labelsraw[[order]];
(*Print[order,labelsordered];*)
stringNickel/@labelsordered
]


(* ::Text:: *)
(*labelsNickelallGo: gives the label for all the diagrams at a given order for  \[CapitalGamma]^(0), \[CapitalGamma]^(2) & \[CapitalGamma]^(4)*)


labelsNickelallGo[pt_,o_]:=Module[{tmplist,labelsraw,order,labelsordered},
Which[pt==0,
tmplist=zeropt1PI[0,o];,
pt==2,
tmplist=twopt1PI[0,o];,
pt==4,
tmplist=fourpt1PI[0,o];,
pt=!=0&&pt=!=2&&pt=!=4, "It works just for \!\(\*SuperscriptBox[\(\[CapitalGamma]\), \((0)\)]\), \!\(\*SuperscriptBox[\(\[CapitalGamma]\), \((2)\)]\) & \!\(\*SuperscriptBox[\(\[CapitalGamma]\), \((4)\)]\)"
];
labelsraw=If[pt==0,labelNickelG0/@(tmplist[[All,2]]),labelNickel/@(tmplist[[All,2]])];
(*Print[labelsraw];*)
order=OrderingBy[delteb/@labelsraw,First];
labelsordered=labelsraw[[order]];
(*Print[order,labelsordered];*)
stringNickel/@labelsordered
]


(* ::Text:: *)
(*Examples*)


(* ::Input:: *)
(*(*labelsNickelGo[2,5]*)*)


(* ::Input:: *)
(*(*labelsNickelGo[4,5]*)*)


(* ::Input:: *)
(*(*AbsoluteTiming@labelsNickelGo[4,7]*)*)


(* ::Input:: *)
(*(*labelsNickelallGo[0,7]*)*)


(* ::Text:: *)
(*labelsNickelallGo: gives the label for all the diagrams at a given order for  \[CapitalGamma]^(0), \[CapitalGamma]^(2) & \[CapitalGamma]^(4) for both cubic and quartic vertices*)


labelsNickelallGo3o4[pt_,o3_,o4_]:=Module[{tmplist,labelsraw,order,labelsordered},
Which[pt==0,
tmplist=zeropt1PI[o3,o4];,
pt==2,
tmplist=twopt1PI[o3,o4];,
pt==4,
tmplist=fourpt1PI[o3,o4];,
pt=!=0&&pt=!=2&&pt=!=4, "It works just for \!\(\*SuperscriptBox[\(\[CapitalGamma]\), \((0)\)]\), \!\(\*SuperscriptBox[\(\[CapitalGamma]\), \((2)\)]\) & \!\(\*SuperscriptBox[\(\[CapitalGamma]\), \((4)\)]\)"
];
labelsraw=If[pt==0,labelNickelG0/@(tmplist[[All,2]]),labelNickel/@(tmplist[[All,2]])];
(*Print[labelsraw];*)
order=OrderingBy[delteb/@labelsraw,First];
labelsordered=labelsraw[[order]];
(*Print[order,labelsordered];*)
stringNickel/@labelsordered
]


(* ::Input:: *)
(*(*labelsNickelallGo3o4[4,2,4]*)*)


(* ::Subsection::Closed:: *)
(*Finding*)


(* ::Text:: *)
(*fromlabelCountnptorder: given a Nickel label, it gives out n = {0, 2, 4} for respectively the \[CapitalGamma]^(0), \[CapitalGamma]^(2) & \[CapitalGamma]^(4) and the order*)


fromlabelCountnptorder[string_]:={StringCount[string,"e"],StringCount[string,"|"]}


(* ::Text:: *)
(*fromlabelCountnptorderB: given a Nickel label, it gives out n = {0, 2, 4} for respectively the \[CapitalGamma]^(0), \[CapitalGamma]^(2) & \[CapitalGamma]^(4) and the number of cubic and quartic vertices (respectively o3 and o4) in the form {n, {o3, o4}}*)


fromlabelCountnptorderB[string_]:=Module[{nn,vv,ll,o3},
nn=StringCount[string,"e"];
vv=StringCount[string,"|"];
ll=StringLength[string];
o3=nn+6vv-2ll;
(*Print[{nn,vv,ll}];*)
{nn,{o3,vv-o3}}
]


(* ::Text:: *)
(*fromlabelFindposition: given a Nickel label, it gives the n = {0, 2, 4} for respectively the \[CapitalGamma]^(0), \[CapitalGamma]^(2) & \[CapitalGamma]^(4), the order and its position in the ordered list*)


fromlabelFindposition[string_]:=Flatten[{fromlabelCountnptorder[string],Position[labelsNickelGo@@fromlabelCountnptorder[string],string]},1]


(* ::Text:: *)
(*Example*)


(* ::Input:: *)
(*(*fromlabelFindposition["e123|e24|34|44||"]*)*)


(* ::Text:: *)
(*fromlabelFindpositionB: generalization of fromlabelFindposition to cubic and quartic vertices and to n = 0.*)


fromlabelFindpositionB[string_]:=Flatten[{fromlabelCountnptorderB[string],Position[labelsNickelallGo3o4@@Flatten[fromlabelCountnptorderB[string]],string]},1]


(* ::Input:: *)
(*(*fromlabelFindpositionB@"ee12|e33|44|45|5|e|"*)*)


(* ::Chapter:: *)
(*Replacements (and visualizations) (from Simple graphs to Complex graphs)*)


(* ::Section:: *)
(*The 4 simple substitutions: the bubble (b), the sunset (s), the triangle (t), the square (c) plus the triangle with bubble (\[Psi]) and the kite (\[Kappa]) (NEW)*)


(* ::Subsection::Closed:: *)
(*Substitutions*)


(* ::Subsubsection::Closed:: *)
(*Bubbles & Sunsets*)


(* ::Text:: *)
(*findbns: Gather in sub-list the same edges (connecting the same vertices) in our theories we can just have cases with 2 or 3 edges connecting the same vertices (that will become bubble (b) or sunset (s))*)


(* ::Input:: *)
(*(*findbns[objgraph_]:=GatherBy[Range[Length@objgraph],objgraph\[LeftDoubleBracket]#\[RightDoubleBracket]&]/.{x_,y_}\[RuleDelayed]{ objgraph\[LeftDoubleBracket]x\[RightDoubleBracket],objgraph\[LeftDoubleBracket]y\[RightDoubleBracket]}/.{x_,y_,z_}\[RuleDelayed]{objgraph\[LeftDoubleBracket]x\[RightDoubleBracket],objgraph\[LeftDoubleBracket]y\[RightDoubleBracket],objgraph\[LeftDoubleBracket]z\[RightDoubleBracket]}/.{x_}\[RuleDelayed]objgraph\[LeftDoubleBracket]x\[RightDoubleBracket]/;NumericQ[x]*)*)


findbns[objgraph_]:=GatherBy[Range[Length@objgraph],objgraph[[#]]&]/.{x_,y_}:>{ objgraph[[x]],objgraph[[y]]}/.{x_,y_,z_}/;(NumericQ[x]&&NumericQ[y]&&NumericQ[z]):>{objgraph[[x]],objgraph[[y]],objgraph[[z]]}/.{x_}:>objgraph[[x]]/;NumericQ[x]


(* ::Text:: *)
(*Substitution of bubbles (b) and sunsets (s). Each special vertex is indicated with the letter(s) that identify the type and in the subscript the name of the normal vertices the it connects to.*)


subsbns[grp_]:=findbns[grp]/.{a_\[UndirectedEdge]c_,a_\[UndirectedEdge]c_}:> {a\[UndirectedEdge]\!\(\*SubscriptBox[\(b\), \({a, c}\)]\),\!\(\*SubscriptBox[\(b\), \({a, c}\)]\)\[UndirectedEdge]c}/.{a_\[UndirectedEdge]c_,a_\[UndirectedEdge]c_,a_\[UndirectedEdge]c_}:> {a\[UndirectedEdge]\!\(\*SubscriptBox[\(s\), \({a, c}\)]\),\!\(\*SubscriptBox[\(s\), \({a, c}\)]\)\[UndirectedEdge]c}


subsjusts[grp_]:=findbns[grp]/.{a_\[UndirectedEdge]c_,a_\[UndirectedEdge]c_,a_\[UndirectedEdge]c_}:> {a\[UndirectedEdge]\!\(\*SubscriptBox[\(s\), \({a, c}\)]\),\!\(\*SubscriptBox[\(s\), \({a, c}\)]\)\[UndirectedEdge]c}


(* ::Subsubsection::Closed:: *)
(*Assess better substitutions (old sequence)*)


(* ::Text:: *)
(*maximcycles: starting from of the possible cycles found, it creates the groups (list of list) of not adjacent one (without any edge in common)*)
(*extractlongest: simply exacts the longest list from a list of list (in our case the biggest group of non adjacent cycles)*)


maximcycles[graph_]:=Module[{origlist,list,ll,i,props,j,seen,group},
origlist=graph;
list=Map[Sort,graph,2];
ll=Length@graph;seen={};group=Table[{},{p,ll}];i=1;
While[i<= ll,props=list[[i]];AppendTo[group[[i]],origlist[[i]]];j=1;
While[j<= ll&&!MemberQ[seen,i],If[Table[MemberQ[list[[j]],props[[k]]],{k,Length@props}]/.List->Or,j++,props=Flatten@AppendTo[props,list[[j]]];seen=AppendTo[seen,j];
AppendTo[group[[i]],origlist[[j]]];j++]];i++];group
]


extractlongest[listlist_]:=Module[{leng,max,pos},
leng=Length/@listlist;max=Max[leng];
pos=Position[leng,max];
If[Extract[listlist,pos]=={},{},Extract[listlist,pos][[1]]]]


(* ::Subsubsection::Closed:: *)
(*Triangles*)


(* ::Text:: *)
(*The special vertices (b) and (s) are removed since they can be part of a triangles or squares*)


keeprealvetex[grp_]:=Delete[grp,Position[grp,Alternatives@@Select[grp,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),Subscript[s, _]|Subscript[b, _],2]&]]]


(* ::Text:: *)
(*In this new graph the cycles (closed loop) of length 3 (edges) are found*)


findtriansafe[grp_]:=If[grp=={},{},FindCycle[grp,{3},All]]


(* ::Text:: *)
(*The best combinations of substitution is selected (see maximcyles)*)


findbesttrian[grp_]:=extractlongest@maximcycles@findtriansafe@keeprealvetex@grp


(* ::Text:: *)
(*Replacement of the special vertex t*)


trianforinsert[trianglist_]:=trianglist/.{a_\[UndirectedEdge]b_,b_\[UndirectedEdge]c_,c_\[UndirectedEdge]a_}:> {\!\(\*SubscriptBox[\(t\), \({a, b, c}\)]\)\[UndirectedEdge]a,\!\(\*SubscriptBox[\(t\), \({a, b, c}\)]\)\[UndirectedEdge]b,\!\(\*SubscriptBox[\(t\), \({a, b, c}\)]\)\[UndirectedEdge]c}


(* ::Text:: *)
(*Substitution of the triangles: old edges are removed ad addition of new vertices*)


subsTrian[grp_]:=Module[{bestTriang},
bestTriang=findbesttrian[grp];
Join[Delete[grp,Position[Sort/@grp,Alternatives@@Sort/@Flatten@bestTriang]],trianforinsert@bestTriang]
]


(* ::Subsubsection::Closed:: *)
(*Triangles and squares*)


(* ::Text:: *)
(*See subsection before, everything is generalized to cycles of length 4 that are substituted with the special vertex c (square). NB now maximcyles takes in considerations triangles and squares and they have the same weight, i.e. the function wants to maximize the number of substitutions (and the gain in number of loop), even if the complexity of the square is higher that the one of the triangles, but it seems a better procedure to give always priority to the triangles.*)


keeprealvetex[grp_]:=Delete[grp,Position[grp,Alternatives@@Select[grp,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),Subscript[s, _]|Subscript[b, _],2]&]]]


find1loop34safe[grp_]:=If[grp=={},{},FindCycle[grp,{3,4},All]]


findbest34group[grp_]:=extractlongest@maximcycles@find1loop34safe@keeprealvetex@grp


(* ::Text:: *)
(*The subscript and the order of the edges in the square is important. Even if in the graphical representation (c) is just a vertex, it accounts also for the order in which it connect the normal vertices around it (i.e. in the graphical representation different graphs may appear the same)*)


triansquarforinsert[trianglist_]:=trianglist/.{a_\[UndirectedEdge]b_,b_\[UndirectedEdge]c_,c_\[UndirectedEdge]a_}:> {\!\(\*SubscriptBox[\(t\), \({a, b, c}\)]\)\[UndirectedEdge]a,\!\(\*SubscriptBox[\(t\), \({a, b, c}\)]\)\[UndirectedEdge]b,\!\(\*SubscriptBox[\(t\), \({a, b, c}\)]\)\[UndirectedEdge]c}/.{a_\[UndirectedEdge]b_,b_\[UndirectedEdge]d_,d_\[UndirectedEdge]e_,e_\[UndirectedEdge]a_}:> {\!\(\*SubscriptBox[\(c\), \({a, b, d, e}\)]\)\[UndirectedEdge]a,\!\(\*SubscriptBox[\(c\), \({a, b, d, e}\)]\)\[UndirectedEdge]b,\!\(\*SubscriptBox[\(c\), \({a, b, d, e}\)]\)\[UndirectedEdge]d,\!\(\*SubscriptBox[\(c\), \({a, b, d, e}\)]\)\[UndirectedEdge]e}


subsTrianSquares[grp_]:=Module[{bestTriangSqu},
bestTriangSqu=findbest34group[grp];
Join[Delete[grp,Position[Sort/@grp,Alternatives@@Sort/@Flatten@bestTriangSqu]],triansquarforinsert@bestTriangSqu]
]


(* ::Subsubsection:: *)
(*\[Psi] and \[Kappa]*)


(* ::Subsubsection::Closed:: *)
(*\[Psi]*)


(* ::Text:: *)
(*find\[Psi]1noprint: given the list of edges, it finds the groups of propagators for the complex vertices \[Psi]a & \[Psi]b*)


find\[Psi]1noprint[grp_]:=Module[{cycletri,cyclebub,ltri,lbub,cycleadj,pos\[Psi]0,pos\[Psi],struct,lstruc,vert,vertbub,vertpiv,structab,struct\[Psi]a,struct\[Psi]b,ij,ji},
cycletri=Map[Sort,If[grp=={},{},FindCycle[grp,{3},All]],2];
cyclebub=Map[Sort,If[grp=={},{},FindCycle[grp,{2},All]],2];
{ltri,lbub}=Length/@{cycletri,cyclebub};
cycleadj=Table[ContainsAll[cycletri[[i]],cyclebub[[j]]],{i,ltri},{j,lbub}];
pos\[Psi]0=Position[cycleadj,True];
pos\[Psi]=Flatten[DeleteCases[GatherBy[pos\[Psi]0,First],x_/;Length@x>1],1];
If[pos\[Psi]==={},{{},{}},
lstruc=Length@pos\[Psi];
struct=Table[{cycletri[[pos\[Psi][[ii,1]]]],cyclebub[[pos\[Psi][[ii,2]]]]},{ii,lstruc}];
(*Print[struct];*)
vert=VertexList/@(struct[[All,1]]);
vertbub=VertexList/@(struct[[All,2]]);
vertpiv=Table[Complement[vert[[iii]],vertbub[[iii]]],{iii,lstruc}];
(*Print[{vert,vertbub,vertpiv}];*)
(*looking for \[Psi]a class*)
structab=Table[{ContainsAll[Join[vert[[jj]],{o1,o2,o3,o4}], AdjacencyList[grp,vertpiv[[jj]]]],Table[ContainsAll[Join[vert[[jj]],{o1,o2,o3,o4}], AdjacencyList[grp,vertbub[[jj,k]]]],{k,2}]},{jj,lstruc}];
(*Print[structab];*)
struct\[Psi]a={};struct\[Psi]b={};ij=1;ji=1;
(*Print[structab\[LeftDoubleBracket]1,1\[RightDoubleBracket]||And@@(structab\[LeftDoubleBracket]1,2\[RightDoubleBracket])];*)
(*Print[!(structab\[LeftDoubleBracket]1,1\[RightDoubleBracket])&&Xor@@(structab\[LeftDoubleBracket]1,2\[RightDoubleBracket])];*)
While[ij<=lstruc, If[structab[[ij,1]]||And@@(structab[[ij,2]]),AppendTo[struct\[Psi]a,struct[[ij]]]];ij++;];
While[ji<=lstruc, If[!(structab[[ji,1]])&&Xor@@(structab[[ji,2]]),AppendTo[struct\[Psi]b,struct[[ji]]]];ji++;];
{struct\[Psi]a,struct\[Psi]b}
]
]


(* ::Text:: *)
(*replace\[Psi]ok: replaces the groups of propagators with the complex vertices \[Psi]a & \[Psi]b*)


replace\[Psi]ok[grp_]:=Module[{struca,strucb,grpouta,grpoutab,lla,llb,propstructa,propstructb,sub\[Psi]a,sub\[Psi]b,tmp,j,test,grptmp},
struca=find\[Psi]1noprint[grp][[1]];
strucb=find\[Psi]1noprint[grp][[2]];
If[struca=={},grpouta=grp;sub\[Psi]a={};,
lla=Length@struca;
propstructa=Table[Join[struca[[i,1]],{struca[[i,2,1]]}],{i,lla}];
(*Print[{struca,propstructa}];*)
j=1;test=0;grptmp=grp;sub\[Psi]a={};
While[j<= lla,tmp=cango[Sort/@grptmp,Sort/@propstructa[[j]]];If[tmp[[1]]==0,grptmp=tmp[[2]];AppendTo[sub\[Psi]a,struca[[j]]];];j++;];
grpouta=grptmp;];
If[strucb=={},grpoutab=grpouta;sub\[Psi]b={};,
llb=Length@strucb;
propstructb=Table[Join[strucb[[i,1]],{strucb[[i,2,1]]}],{i,llb}];
(*Print[{strucb,propstructb}];*)
j=1;test=0;grptmp=grpouta;sub\[Psi]b={};
While[j<= llb,tmp=cango[Sort/@grptmp,Sort/@propstructb[[j]]];If[tmp[[1]]==0,grptmp=tmp[[2]];AppendTo[sub\[Psi]b,strucb[[j]]];];j++;];
grpoutab=grptmp;];
{grpoutab,{sub\[Psi]a,sub\[Psi]b}}
]


(* ::Text:: *)
(*\[Psi]forinsert: last stage, it inserts the the complex vertices \[Psi]a & \[Psi]b*)


\[Psi]forinsert[triangblist_]:=Module[{struct,suba,lstruca,verta,vertbuba,vertpiva,verticesa,subb,lstrucb,vertb,vertbubb,vertpivb,verticesb},
struct=triangblist;
(*Print[struct];*)
If[struct[[1]]=={},suba={},
lstruca=Length@struct[[1]];
verta=VertexList/@(struct[[1,All,1]]);
vertbuba=VertexList/@(struct[[1,All,2]]);
vertpiva=Table[Complement[verta[[iii]],vertbuba[[iii]]],{iii,lstruca}];
verticesa=Table[Flatten@{vertpiva[[i]],vertbuba[[i]]},{i,lstruca}];
(*Print[verticesa];*)
suba=Table[ {Subscript[\[Psi]a, verticesa[[i]]]\[UndirectedEdge]verticesa[[i,1]],Subscript[\[Psi]a, verticesa[[i]]]\[UndirectedEdge]verticesa[[i,2]],Subscript[\[Psi]a, verticesa[[i]]]\[UndirectedEdge]verticesa[[i,3]]},{i,lstruca}];
];
If[struct[[2]]=={},subb={},
lstrucb=Length@struct[[2]];
vertb=VertexList/@(struct[[2,All,1]]);
vertbubb=VertexList/@(struct[[2,All,2]]);
vertpivb=Table[Complement[vertb[[iii]],vertbubb[[iii]]],{iii,lstrucb}];
verticesb=Table[Flatten@{vertpivb[[i]],vertbubb[[i]]},{i,lstrucb}];
(*Print[verticesb];*)
subb=Table[ {Subscript[\[Psi]b, verticesb[[i]]]\[UndirectedEdge]verticesb[[i,1]],Subscript[\[Psi]b, verticesb[[i]]]\[UndirectedEdge]verticesb[[i,2]],Subscript[\[Psi]b, verticesb[[i]]]\[UndirectedEdge]verticesb[[i,3]]},{i,lstrucb}];
];
Join[suba,subb]
]


(* ::Text:: *)
(*subs\[Psi]ab: given the list of edges, it does the substitutions for the complex vertices \[Psi]a & \[Psi]b*)


subs\[Psi]ab[grp_]:=Module[{tmp,sunsDone,suns,sunsDonedelet},
sunsDone=subsjusts@findbns@grp;
suns=Cases[sunsDone,(x_/;MemberQ[x,Subscript[s, _],2])];
sunsDonedelet=keeprealvetexs@Flatten@subsjusts@findbns@grp;
tmp=replace\[Psi]ok[sunsDonedelet];
(Join@@{tmp[[1]],suns,\[Psi]forinsert[tmp[[2]]]})/.a_\[UndirectedEdge]o3:>o3\[UndirectedEdge]a/.a_\[UndirectedEdge]o4:>o4\[UndirectedEdge]a
]


(* ::Input:: *)
(**)


(* ::Text:: *)
(*All the same functions of the rest of the section but just for \[Psi]a (that is analytic differently than \[Psi]b that is numeric)*)


replace\[Psi]a[grp_]:=Module[{struca,strucb,grpouta,grpoutab,lla,llb,propstructa,propstructb,sub\[Psi]a,sub\[Psi]b,tmp,j,test,grptmp},
struca=find\[Psi]1noprint[grp][[1]];
strucb=find\[Psi]1noprint[grp][[2]];
If[struca=={},grpouta=grp;sub\[Psi]a={};,
lla=Length@struca;
propstructa=Table[Join[struca[[i,1]],{struca[[i,2,1]]}],{i,lla}];
(*Print[{struca,propstructa}];*)
j=1;test=0;grptmp=grp;sub\[Psi]a={};
While[j<= lla,tmp=cango[Sort/@grptmp,Sort/@propstructa[[j]]];If[tmp[[1]]==0,grptmp=tmp[[2]];AppendTo[sub\[Psi]a,struca[[j]]];];j++;];
grpouta=grptmp;];
{grpouta,{sub\[Psi]a}}
]


\[Psi]aforinsert[triangblist_]:=Module[{struct,suba,lstruca,verta,vertbuba,vertpiva,verticesa,subb,lstrucb,vertb,vertbubb,vertpivb,verticesb},
struct=triangblist;
(*Print[struct];*)
If[struct[[1]]=={},suba={},
lstruca=Length@struct[[1]];
verta=VertexList/@(struct[[1,All,1]]);
vertbuba=VertexList/@(struct[[1,All,2]]);
vertpiva=Table[Complement[verta[[iii]],vertbuba[[iii]]],{iii,lstruca}];
verticesa=Table[Flatten@{vertpiva[[i]],vertbuba[[i]]},{i,lstruca}];
(*Print[verticesa];*)
suba=Table[ {Subscript[\[Psi]a, verticesa[[i]]]\[UndirectedEdge]verticesa[[i,1]],Subscript[\[Psi]a, verticesa[[i]]]\[UndirectedEdge]verticesa[[i,2]],Subscript[\[Psi]a, verticesa[[i]]]\[UndirectedEdge]verticesa[[i,3]]},{i,lstruca}];
];
suba
]


subs\[Psi]a[grp_]:=Module[{tmp,sunsDone,suns,sunsDonedelet},
sunsDone=subsjusts@findbns@grp;
suns=Cases[sunsDone,(x_/;MemberQ[x,Subscript[s, _],2])];
sunsDonedelet=keeprealvetexs@Flatten@subsjusts@findbns@grp;
tmp=replace\[Psi]a[sunsDonedelet];
(Join@@{tmp[[1]],suns,\[Psi]aforinsert[tmp[[2]]]})/.a_\[UndirectedEdge]o3:>o3\[UndirectedEdge]a/.a_\[UndirectedEdge]o4:>o4\[UndirectedEdge]a
]


(* ::Subsubsection::Closed:: *)
(*\[Kappa]*)


(* ::Text:: *)
(*find\[Kappa]1noprint: given the list of edges, it finds the groups of propagators for the complex vertex \[Kappa]*)


find\[Kappa]1noprint[grp_]:=Module[{cycletri,ltri,sortedgrp,cycleadj,pos\[Kappa],struct,lstruc,vert1,vert2,vertpiv,vertall,structmom,struct\[Kappa],ij},
cycletri=Map[Sort,If[grp=={},{},FindCycle[grp,{3},All]],2];
ltri=Length@cycletri;pos\[Kappa]={};
(*Print[cycletri];*)
sortedgrp=Sort/@grp;
cycleadj=Flatten[Table[Table[If[(Length@Intersection[cycletri[[i]],cycletri[[i+j]]])==1&&And@@Table[Count[sortedgrp,(*Intersection*)(Sort/@Union[cycletri[[i]],cycletri[[i+j]]])[[k]]]==1,{k,Length@Union[cycletri[[i]],cycletri[[i+j]]]}],AppendTo[pos\[Kappa],{i,i+j}]],{j,ltri-i}],{i,ltri}],1];
If[pos\[Kappa]==={},{},
lstruc=Length@pos\[Kappa];
(*Print[pos\[Kappa]];*)
struct=Table[{cycletri[[pos\[Kappa][[ii,1]]]],cycletri[[pos\[Kappa][[ii,2]]]]},{ii,lstruc}];
(*Print[struct];*)
vert1=VertexList/@(struct[[All,1]]);
vert2=VertexList/@(struct[[All,2]]);
vertpiv=Table[Intersection[vert1[[iii]],vert2[[iii]]],{iii,lstruc}];
vertall=Table[Union[vert1[[iii]],vert2[[iii]]],{iii,lstruc}];
(*Print[{vert1,vert2,vertpiv}];*)
(*looking for \[Psi]a class*)
structmom=Table[Table[ContainsAll[Join[vertall[[jj]],{o1,o2,o3,o4}], AdjacencyList[grp,vertpiv[[jj,k]]]],{k,2}],{jj,lstruc}];
(*Print[structmom];*)
struct\[Kappa]={};ij=1;
(*Print[And@@(structmom\[LeftDoubleBracket]1\[RightDoubleBracket])];*)
While[ij<=lstruc, If[And@@(structmom[[ij]]),AppendTo[struct\[Kappa],struct[[ij]]]];ij++;];
struct\[Kappa]
]
]


(* ::Text:: *)
(*replace\[Kappa]: replaces the groups of propagators with the complex vertices \[Kappa]*)


replace\[Kappa][grp_]:=Module[{struc,grpout,ll,propstruct,sub\[Kappa],tmp,j,test,grptmp},
struc=find\[Kappa]1noprint[grp];
If[struc=={},grpout=grp;sub\[Kappa]={};,
ll=Length@struc;
propstruct=Table[Union[struc[[i,1]],struc[[i,2]]],{i,ll}];
(*Print[{struc,propstruct}];*)
j=1;test=0;grptmp=grp;sub\[Kappa]={};
While[j<= ll,tmp=cango[Sort/@grptmp,Sort/@propstruct[[j]]];If[tmp[[1]]==0,grptmp=tmp[[2]];AppendTo[sub\[Kappa],struc[[j]]];];
(*Print[grptmp];*)j++;];
grpout=grptmp;];
{grpout,sub\[Kappa]}
]


(* ::Text:: *)
(*\[Kappa]forinsert: last stage, it inserts the the complex vertices \[Kappa]*)


\[Kappa]forinsert[triangblist_]:=Module[{struct,sub,lstruc,vert,vertbuba,vertpiv,vertnopiv,vertices},
struct=triangblist;
(*Print[struct];*)
If[struct=={},sub={},
lstruc=Length@struct;
vert=Table[VertexList/@(struct[[j]]),{j,lstruc}];
vertpiv=Table[Intersection@@vert[[j]],{j,lstruc}];
vertnopiv=Table[Complement[Union@@vert[[j]],vertpiv[[j]]],{j,lstruc}];
vertices=Table[{vertnopiv[[j,1]],vertpiv[[j,1]],vertnopiv[[j,2]],vertpiv[[j,2]]},{j,lstruc}];
(*Print[vert,vertpiv,vertnopiv];*)
sub=Table[ {Subscript[\[Kappa], vertices[[i]]]\[UndirectedEdge]vertices[[i,1]],Subscript[\[Kappa], vertices[[i]]]\[UndirectedEdge]vertices[[i,2]],Subscript[\[Kappa], vertices[[i]]]\[UndirectedEdge]vertices[[i,3]],Subscript[\[Kappa], vertices[[i]]]\[UndirectedEdge]vertices[[i,4]]},{i,lstruc}];
];
sub
]


(* ::Text:: *)
(*subs\[Kappa]: given the list of edges, it does the substitutions for the complex vertices \[Kappa]*)


subs\[Kappa][grp_]:=Module[{tmp,sunsDone,suns,sunsDonedelet},
sunsDone=subsjusts@findbns@grp;
suns=Cases[sunsDone,(x_/;MemberQ[x,Subscript[s, _],2])];
sunsDonedelet=keeprealvetexs@Flatten@subsjusts@findbns@grp;
tmp=replace\[Kappa][sunsDonedelet];
(Join@@{tmp[[1]],suns,\[Kappa]forinsert[tmp[[2]]]})/.a_\[UndirectedEdge]o3:>o3\[UndirectedEdge]a/.a_\[UndirectedEdge]o4:>o4\[UndirectedEdge]a
]


(* ::Subsubsection::Closed:: *)
(*All together*)


(* ::Text:: *)
(*replace\[Psi]ab\[Kappa]: replaces the groups of propagators with the complex vertices \[Psi]a, \[Psi]b and \[Kappa]*)


replace\[Psi]ab\[Kappa][grp_]:=Module[{struca,strucb,struc,grpouta,grpoutab,grpout,lla,llb,ll,propstructa,propstructb,propstruct,sub\[Psi]a,sub\[Psi]b,sub\[Kappa],tmp,j,test,grptmp},
struca=find\[Psi]1noprint[grp][[1]];
strucb=find\[Psi]1noprint[grp][[2]];
struc=find\[Kappa]1noprint[grp];
If[struca=={},grpouta=grp;sub\[Psi]a={};,
lla=Length@struca;
propstructa=Table[Join[struca[[i,1]],{struca[[i,2,1]]}],{i,lla}];
(*Print[{struca,propstructa}];*)
j=1;test=0;grptmp=grp;sub\[Psi]a={};
While[j<= lla,tmp=cango[Sort/@grptmp,Sort/@propstructa[[j]]];If[tmp[[1]]==0,grptmp=tmp[[2]];AppendTo[sub\[Psi]a,struca[[j]]];];j++;];
grpouta=grptmp;];
If[strucb=={},grpoutab=grpouta;sub\[Psi]b={};,
llb=Length@strucb;
propstructb=Table[Join[strucb[[i,1]],{strucb[[i,2,1]]}],{i,llb}];
(*Print[{strucb,propstructb}];*)
j=1;test=0;grptmp=grpouta;sub\[Psi]b={};
While[j<= llb,tmp=cango[Sort/@grptmp,Sort/@propstructb[[j]]];If[tmp[[1]]==0,grptmp=tmp[[2]];AppendTo[sub\[Psi]b,strucb[[j]]];];j++;];
grpoutab=grptmp;];
If[struc=={},grpout=grpoutab;sub\[Kappa]={};,
ll=Length@struc;
propstruct=Table[Union[struc[[i,1]],struc[[i,2]]],{i,ll}];
(*Print[{struc,propstruct}];*)
j=1;test=0;grptmp=grpoutab;sub\[Kappa]={};
While[j<= ll,tmp=cango[Sort/@grptmp,Sort/@propstruct[[j]]];If[tmp[[1]]==0,grptmp=tmp[[2]];AppendTo[sub\[Kappa],struc[[j]]];];
(*Print[grptmp];*)j++;];
grpout=grptmp;];
{grpout,{sub\[Psi]a,sub\[Psi]b,sub\[Kappa]}}
]


(* ::Text:: *)
(*subs\[Psi]ab\[Kappa]: given the list of edges, it does the substitutions for the complex vertices \[Psi]a, \[Psi]b and \[Kappa]*)


subs\[Psi]ab\[Kappa][grp_]:=Module[{tmp,sunsDone,suns,sunsDonedelet},
sunsDone=subsjusts@findbns@grp;
suns=Cases[sunsDone,(x_/;MemberQ[x,Subscript[s, _],2])];
sunsDonedelet=keeprealvetexs@Flatten@subsjusts@findbns@grp;
tmp=replace\[Psi]ab\[Kappa][sunsDonedelet];
(Join@@{tmp[[1]],suns,\[Psi]forinsert[{tmp[[2,1]],tmp[[2,2]]}],\[Kappa]forinsert[tmp[[2,3]]]})/.a_\[UndirectedEdge]o3:>o3\[UndirectedEdge]a/.a_\[UndirectedEdge]o4:>o4\[UndirectedEdge]a
]


(* ::Subsubsection::Closed:: *)
(*Bubbles, triangles, and squares (old and new sequence) (after removing and substituting \[Psi] and/or \[Kappa])*)


(* ::Text:: *)
(*See subsection before, everything is generalized to cycles of length 4 that are substituted with the special vertex c (square). NB now maximcyles takes in considerations triangles and squares and they have the same weight, i.e. the function wants to maximize the number of substitutions (and the gain in number of loop), even if the complexity of the square is higher that the one of the triangles, but it seems a better procedure to give always priority to the triangles.*)


keeprealvetexs[grp_]:=DeleteCases[grp,Subscript[s, _]\[UndirectedEdge]a_|a_\[UndirectedEdge]Subscript[s, _]]


keeprealvetexs\[Psi]a[grp_]:=DeleteCases[grp,Subscript[s, _]\[UndirectedEdge]a_|a_\[UndirectedEdge]Subscript[s, _]|Subscript[\[Psi]a, _]\[UndirectedEdge]a_|a_\[UndirectedEdge]Subscript[\[Psi]a, _]]


keeprealvetexs\[Psi][grp_]:=DeleteCases[grp,Subscript[s, _]\[UndirectedEdge]a_|a_\[UndirectedEdge]Subscript[s, _]|Subscript[\[Psi]a, _]\[UndirectedEdge]a_|a_\[UndirectedEdge]Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]\[UndirectedEdge]a_|a_\[UndirectedEdge]Subscript[\[Psi]b, _]]


keeprealvetexs\[Psi]\[Kappa][grp_]:=DeleteCases[grp,Subscript[s, _]\[UndirectedEdge]a_|a_\[UndirectedEdge]Subscript[s, _]|Subscript[\[Psi]a, _]\[UndirectedEdge]a_|a_\[UndirectedEdge]Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]\[UndirectedEdge]a_|a_\[UndirectedEdge]Subscript[\[Psi]b, _]|Subscript[\[Kappa], _]\[UndirectedEdge]a_|a_\[UndirectedEdge]Subscript[\[Kappa], _]]


find1loop234safe[grp_]:=If[grp=={},{},FindCycle[grp,4,All]]


(* ::Input:: *)
(*(*findbest234group[grp_]:=extractlongest@maximcycles@find1loop234safe@keeprealvetexs@grp*)*)


cango[list_,check_]:=Module[{listtmp,i,test,j,cc},
i=1;listtmp=list;test=0;cc=Length@check;
While[i<=cc&&test==0,If[MemberQ[listtmp,check[[i]]],listtmp=Delete[listtmp,Position[listtmp,check[[i]]][[1]]];,test=1;];i++;];
{test,listtmp}]


maximcyclestransitionsnoinside[subgraph_,fullgraph_]:=Module[{origlist,list,listprop,ll,i,listprop1,listpropold,listprop1tmp,j,seen,group},
origlist=subgraph;
list=Map[Sort,subgraph,2];
listprop=Sort/@fullgraph;
ll=Length@origlist;seen={};group=Table[{},{p,ll}];i=1;
While[i<= ll(*&&!MemberQ[seen,i]*),listprop1=cango[listprop,list[[i]]];AppendTo[group[[i]],origlist[[i]]];j=1;
While[j<= ll,If[i==j,Null;,listpropold=listprop1;listprop1tmp=cango[listprop1[[2]],list[[j]]];
If[listprop1tmp[[1]]==1,listprop1=listpropold,listprop1=listprop1tmp;AppendTo[group[[i]],origlist[[j]]]]];j++];i++];
group
]


extractlongestweight[listlist_]:=Module[{leng,max,pos,bestones,weighttot},
leng=Length/@listlist;max=Max[leng];
pos=Position[leng,max];
(*Print[{max,pos}];*)
bestones=Extract[listlist,pos];
If[bestones=={},{},(*Print[bestones];*)(*weighttot=Total/@(Table[Length/@bestones\[LeftDoubleBracket]i\[RightDoubleBracket],{i,Length@bestones}]/.{2\[Rule] -4,3\[Rule] -2,4\[Rule] -1});
Print[weighttot];(bestones\[LeftDoubleBracket]Flatten@Position[weighttot,Min[weighttot]]\[RightDoubleBracket])\[LeftDoubleBracket]1\[RightDoubleBracket];*)
bestones[[1]]]]


(* ::Text:: *)
(*The subscript and the order of the edges in the square is important. Even if in the graphical representation (c) is just a vertex, it accounts also for the order in which it connect the normal vertices around it (i.e. in teh graphical representation different graphs may appear the same)*)


bubtriansquarforinsert[bubtriangsqlist_]:=bubtriangsqlist/.{a_\[UndirectedEdge]c_,c_\[UndirectedEdge]a_}:> {a\[UndirectedEdge]\!\(\*SubscriptBox[\(b\), \({a, c}\)]\),\!\(\*SubscriptBox[\(b\), \({a, c}\)]\)\[UndirectedEdge]c}/.{a_\[UndirectedEdge]b_,b_\[UndirectedEdge]c_,c_\[UndirectedEdge]a_}:> {\!\(\*SubscriptBox[\(t\), \({a, b, c}\)]\)\[UndirectedEdge]a,\!\(\*SubscriptBox[\(t\), \({a, b, c}\)]\)\[UndirectedEdge]b,\!\(\*SubscriptBox[\(t\), \({a, b, c}\)]\)\[UndirectedEdge]c}/.{a_\[UndirectedEdge]b_,b_\[UndirectedEdge]d_,d_\[UndirectedEdge]e_,e_\[UndirectedEdge]a_}:> {\!\(\*SubscriptBox[\(c\), \({a, b, d, e}\)]\)\[UndirectedEdge]a,\!\(\*SubscriptBox[\(c\), \({a, b, d, e}\)]\)\[UndirectedEdge]b,\!\(\*SubscriptBox[\(c\), \({a, b, d, e}\)]\)\[UndirectedEdge]d,\!\(\*SubscriptBox[\(c\), \({a, b, d, e}\)]\)\[UndirectedEdge]e}


subsTrianSquares[grp_]:=Module[{bestTriangSqu},
bestTriangSqu=findbest34group[grp];
Join[Delete[grp,Position[Sort/@grp,Alternatives@@Sort/@Flatten@bestTriangSqu]],triansquarforinsert@bestTriangSqu]
]


(* ::Text:: *)
(*sunsSub: given the list of edges, it does the substitutions of the sunsets*)


sunsSub[grp_]:=Module[{sunsDone,suns,sunsDonedelet},
sunsDone=subsjusts@findbns@grp;
suns=Cases[sunsDone,(x_/;MemberQ[x,Subscript[s, _],2])];
sunsDonedelet=keeprealvetexs@Flatten@subsjusts@findbns@grp;
Join[sunsDonedelet,suns]/.a_\[UndirectedEdge]o1:>o1\[UndirectedEdge]a/.a_\[UndirectedEdge]o2:>o2\[UndirectedEdge]a/.a_\[UndirectedEdge]o3:>o3\[UndirectedEdge]a/.a_\[UndirectedEdge]o4:>o4\[UndirectedEdge]a
]


subsBubTrianSquares[grp_]:=Module[{sunsDone,suns,sunsDonedelet,bestBubTriangSqu,toinsert,deleted},
sunsDone=subsjusts@findbns@grp;
suns=Cases[sunsDone,(x_/;MemberQ[x,Subscript[s, _],2])];
sunsDonedelet=keeprealvetexs@Flatten@subsjusts@findbns@grp;
bestBubTriangSqu=extractlongestweight@maximcyclestransitionsnoinside[find1loop234safe@sunsDonedelet,sunsDonedelet];
toinsert=bubtriansquarforinsert@bestBubTriangSqu;
deleted=cango[Sort/@sunsDonedelet,Sort/@Join@@bestBubTriangSqu];
(*Print[sunsDone];Print[bestBubTriangSqu];Print[toinsert];Print[deleted];*)
If[deleted[[1]]==1,Message[subsBubTrianSquares::substsubdiagr];,
Join[deleted[[2]],suns,toinsert]]/.a_\[UndirectedEdge]o1:>o1\[UndirectedEdge]a/.a_\[UndirectedEdge]o2:>o2\[UndirectedEdge]a/.a_\[UndirectedEdge]o3:>o3\[UndirectedEdge]a/.a_\[UndirectedEdge]o4:>o4\[UndirectedEdge]a
]


subs\[Psi]aBubTrianSquares[grp_]:=Module[{suns\[Psi]Done,suns\[Psi],suns\[Psi]Donedelet,bestBubTriangSqu,toinsert,deleted},
suns\[Psi]Done=subs\[Psi]a@grp;
suns\[Psi]=Cases[suns\[Psi]Done,(x_/;MemberQ[x,Subscript[s, _],2])|(x_/;MemberQ[x,Subscript[\[Psi]a, _],2])];
suns\[Psi]Donedelet=keeprealvetexs\[Psi]a@Flatten@subs\[Psi]a@grp;
bestBubTriangSqu=extractlongestweight@maximcyclestransitionsnoinside[find1loop234safe@suns\[Psi]Donedelet,suns\[Psi]Donedelet];
toinsert=bubtriansquarforinsert@bestBubTriangSqu;
deleted=cango[Sort/@suns\[Psi]Donedelet,Sort/@Join@@bestBubTriangSqu];
(*Print[sunsDone];Print[bestBubTriangSqu];Print[toinsert];Print[deleted];*)
If[deleted[[1]]==1,Message[subs\[Psi]aBubTrianSquares::substsubdiagr];,
Join[deleted[[2]],suns\[Psi],toinsert]]/.a_\[UndirectedEdge]o1:>o1\[UndirectedEdge]a/.a_\[UndirectedEdge]o2:>o2\[UndirectedEdge]a/.a_\[UndirectedEdge]o3:>o3\[UndirectedEdge]a/.a_\[UndirectedEdge]o4:>o4\[UndirectedEdge]a
]


subs\[Psi]BubTrianSquares[grp_]:=Module[{suns\[Psi]Done,suns\[Psi],suns\[Psi]Donedelet,bestBubTriangSqu,toinsert,deleted},
suns\[Psi]Done=subs\[Psi]ab@grp;
suns\[Psi]=Cases[suns\[Psi]Done,(x_/;MemberQ[x,Subscript[s, _],2])|(x_/;MemberQ[x,Subscript[\[Psi]a, _],2])|(x_/;MemberQ[x,Subscript[\[Psi]b, _],2])];
suns\[Psi]Donedelet=keeprealvetexs\[Psi]@Flatten@subs\[Psi]ab@grp;
bestBubTriangSqu=extractlongestweight@maximcyclestransitionsnoinside[find1loop234safe@suns\[Psi]Donedelet,suns\[Psi]Donedelet];
toinsert=bubtriansquarforinsert@bestBubTriangSqu;
deleted=cango[Sort/@suns\[Psi]Donedelet,Sort/@Join@@bestBubTriangSqu];
(*Print[sunsDone];Print[bestBubTriangSqu];Print[toinsert];Print[deleted];*)
If[deleted[[1]]==1,Message[subs\[Psi]BubTrianSquares::substsubdiagr];,
Join[deleted[[2]],suns\[Psi],toinsert]]/.a_\[UndirectedEdge]o1:>o1\[UndirectedEdge]a/.a_\[UndirectedEdge]o2:>o2\[UndirectedEdge]a/.a_\[UndirectedEdge]o3:>o3\[UndirectedEdge]a/.a_\[UndirectedEdge]o4:>o4\[UndirectedEdge]a
]


subs\[Psi]\[Kappa]BubTrianSquares[grp_]:=Module[{suns\[Psi]\[Kappa]Done,suns\[Psi]\[Kappa],suns\[Psi]\[Kappa]Donedelet,bestBubTriangSqu,toinsert,deleted},
suns\[Psi]\[Kappa]Done=subs\[Psi]ab\[Kappa]@grp;
suns\[Psi]\[Kappa]=Cases[suns\[Psi]\[Kappa]Done,(x_/;MemberQ[x,Subscript[s, _],2])|(x_/;MemberQ[x,Subscript[\[Psi]a, _],2])|(x_/;MemberQ[x,Subscript[\[Psi]b, _],2])|(x_/;MemberQ[x,Subscript[\[Kappa], _],2])];
suns\[Psi]\[Kappa]Donedelet=keeprealvetexs\[Psi]\[Kappa]@Flatten@subs\[Psi]ab\[Kappa]@grp;
bestBubTriangSqu=extractlongestweight@maximcyclestransitionsnoinside[find1loop234safe@suns\[Psi]\[Kappa]Donedelet,suns\[Psi]\[Kappa]Donedelet];
toinsert=bubtriansquarforinsert@bestBubTriangSqu;
deleted=cango[Sort/@suns\[Psi]\[Kappa]Donedelet,Sort/@Join@@bestBubTriangSqu];
(*Print[sunsDone];Print[bestBubTriangSqu];Print[toinsert];Print[deleted];*)
If[deleted[[1]]==1,Message[subs\[Psi]\[Kappa]BubTrianSquares::substsubdiagr];,
Join[deleted[[2]],suns\[Psi]\[Kappa],toinsert]]/.a_\[UndirectedEdge]o1:>o1\[UndirectedEdge]a/.a_\[UndirectedEdge]o2:>o2\[UndirectedEdge]a/.a_\[UndirectedEdge]o3:>o3\[UndirectedEdge]a/.a_\[UndirectedEdge]o4:>o4\[UndirectedEdge]a
]


subsBubTrianSquares::substsubdiagr = "Error."


subs\[Psi]aBubTrianSquares::substsubdiagr = "Error."


subs\[Psi]BubTrianSquares::substsubdiagr = "Error."


subs\[Psi]\[Kappa]BubTrianSquares::substsubdiagr = "Error."


(* ::Subsection:: *)
(*Application of substitution for 0, 2, 4 point functions (up to level 3)*)


(* ::Text:: *)
(*Taking the original lists and applying the basic substitutions:*)
(**)
(*s just the sunsets,*)
(*bst for bubbles, sunsets and triangles  with the old algorithm,*)
(*bstc for bubbles sunsets, triangles and squares with the old algorithm,*)
(*bstc2 for bubbles sunsets, triangles and squares *)
(*(NB: this new list are not "flat" but have list corresponding to the special vertices.)*)
(*bstc2a for bubbles sunsets, triangles and squares and \[Psi]a*)
(*bstc3a for bubbles sunsets, triangles and squares, \[Psi]a and \[Psi]b*)
(*bstc3 for bubbles sunsets, triangles and squares, \[Psi]a, \[Psi]b and \[Kappa]*)


(* ::Subsubsection:: *)
(*Sunsets*)


zeropt1PIs[n3_,n4_]:=Module[{tmp=zeropt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsjusts@tmp[[i,2]])},{i,tt}]]


twopt1PIs[n3_,n4_]:=Module[{tmp=twopt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsjusts@tmp[[i,2]])},{i,tt}]]


fourpt1PIs[n3_,n4_]:=Module[{tmp=fourpt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsjusts@tmp[[i,2]])},{i,tt}]]


(* ::Subsubsection::Closed:: *)
(*Sunsets, Bubbles, and Triangles*)


zeropt1PIbst[n3_,n4_]:=Module[{tmp=zeropt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsTrian@subsbns@tmp[[i,2]])},{i,tt}]]


twopt1PIbst[n3_,n4_]:=Module[{tmp=twopt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsTrian@subsbns@tmp[[i,2]])},{i,tt}]]


fourpt1PIbst[n3_,n4_]:=Module[{tmp=fourpt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsTrian@subsbns@tmp[[i,2]])},{i,tt}]]


(* ::Subsubsection::Closed:: *)
(*Sunsets, Bubbles, Triangles, and Squares*)


zeropt1PIbstc[n3_,n4_]:=Module[{tmp=zeropt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsTrianSquares@subsbns@tmp[[i,2]])},{i,tt}]]


twopt1PIbstc[n3_,n4_]:=Module[{tmp=twopt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsTrianSquares@subsbns@tmp[[i,2]])},{i,tt}]]


fourpt1PIbstc[n3_,n4_]:=Module[{tmp=fourpt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsTrianSquares@subsbns@tmp[[i,2]])},{i,tt}]]


(* ::Subsubsection::Closed:: *)
(*Sunsets, Bubbles, Triangles, and Squares NEW SEQUENCE*)


zeropt1PIbstc2[n3_,n4_]:=Module[{tmp=zeropt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsBubTrianSquares@tmp[[i,2]])},{i,tt}]]


twopt1PIbstc2[n3_,n4_]:=Module[{tmp=twopt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsBubTrianSquares@tmp[[i,2]])},{i,tt}]]


fourpt1PIbstc2[n3_,n4_]:=Module[{tmp=fourpt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsBubTrianSquares@tmp[[i,2]])},{i,tt}]]


(* ::Subsubsection::Closed:: *)
(*\[Psi] and \[Psi] \[Kappa]*)


zeropt1PI\[Psi][n3_,n4_]:=Module[{tmp=zeropt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subs\[Psi]ab@tmp[[i,2]])},{i,tt}]]


twopt1PI\[Psi][n3_,n4_]:=Module[{tmp=twopt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subs\[Psi]ab@tmp[[i,2]])},{i,tt}]]


fourpt1PI\[Psi][n3_,n4_]:=Module[{tmp=fourpt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subs\[Psi]ab@tmp[[i,2]])},{i,tt}]]


(* ::Input:: *)
(**)


zeropt1PI\[Psi]\[Kappa][n3_,n4_]:=Module[{tmp=zeropt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subs\[Psi]ab\[Kappa]@tmp[[i,2]])},{i,tt}]]


twopt1PI\[Psi]\[Kappa][n3_,n4_]:=Module[{tmp=twopt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subs\[Psi]ab\[Kappa]@tmp[[i,2]])},{i,tt}]]


fourpt1PI\[Psi]\[Kappa][n3_,n4_]:=Module[{tmp=fourpt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subs\[Psi]ab\[Kappa]@tmp[[i,2]])},{i,tt}]]


(* ::Subsubsection::Closed:: *)
(*Sunsets, Bubbles, Triangles, and Squares NEW SEQUENCE + \[Psi] & \[Kappa]*)


zeropt1PIbstc2a[n3_,n4_]:=Module[{tmp=zeropt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subs\[Psi]aBubTrianSquares@tmp[[i,2]])},{i,tt}]]


twopt1PIbstc2a[n3_,n4_]:=Module[{tmp=twopt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subs\[Psi]aBubTrianSquares@tmp[[i,2]])},{i,tt}]]


fourpt1PIbstc2a[n3_,n4_]:=Module[{tmp=fourpt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subs\[Psi]aBubTrianSquares@tmp[[i,2]])},{i,tt}]]


(* ::Input:: *)
(**)


zeropt1PIbstc3a[n3_,n4_]:=Module[{tmp=zeropt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subs\[Psi]BubTrianSquares@tmp[[i,2]])},{i,tt}]]


twopt1PIbstc3a[n3_,n4_]:=Module[{tmp=twopt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subs\[Psi]BubTrianSquares@tmp[[i,2]])},{i,tt}]]


fourpt1PIbstc3a[n3_,n4_]:=Module[{tmp=fourpt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subs\[Psi]BubTrianSquares@tmp[[i,2]])},{i,tt}]]


(* ::Input:: *)
(**)


zeropt1PIbstc3[n3_,n4_]:=Module[{tmp=zeropt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subs\[Psi]\[Kappa]BubTrianSquares@tmp[[i,2]])},{i,tt}]]


twopt1PIbstc3[n3_,n4_]:=Module[{tmp=twopt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subs\[Psi]\[Kappa]BubTrianSquares@tmp[[i,2]])},{i,tt}]]


fourpt1PIbstc3[n3_,n4_]:=Module[{tmp=fourpt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subs\[Psi]\[Kappa]BubTrianSquares@tmp[[i,2]])},{i,tt}]]


(* ::Section::Closed:: *)
(*The substitutions involving special vertices*)


(* ::Text:: *)
(*The substitutions at this order flatten the list to be able to identify cycles*)
(*In this subsection we identify cycles (or combinations of cycles) that contain also special vertices, these will be substituted with new special vertices whose value we know (special tadpoles) or that we can compute numerically.*)


(* ::Subsection::Closed:: *)
(*Complex tadpoles*)


(* ::Subsubsection::Closed:: *)
(*Tadpole with sunset (\[Tau]) (Simple loop with just the sunset)*)


find4safe[grp_]:=If[grp=={},{},FindCycle[grp,{4},All]]


(* ::Text:: *)
(*Identifying the tadpole with the sunset (\[Tau])*)


specialverticesNos=Flatten@{{Subscript[b, _],Subscript[t, _],Subscript[c, _],Subscript[\[Psi]a, _]},{Subscript[\[Sigma], _],Subscript[\[Sigma]2, _],Subscript[\[Theta], _],Subscript[\[Theta]2, _],Subscript[\[Gamma], _],Subscript[\[Gamma]\[Kappa], _],Subscript[\[Gamma]c, _],Subscript[\[Psi]b, _],Subscript[\[Kappa], _]}};


findsunloop[grp_]:=Select[If[grp=={},{},FindCycle[grp,{4},All]],(Count[#,\!\(\*SubscriptBox[\(s\), \({a_, b_}\)]\)\[UndirectedEdge]a_]+Count[#,\!\(\*SubscriptBox[\(s\), \({a_, b_}\)]\)\[UndirectedEdge]b_]==1&&!MemberQ[VertexList@Flatten@#,Alternatives@@specialverticesNos])&]


(* ::Text:: *)
(*Replacement: identifying the pivot to be substituted with the special vertex*)
(*NB the "tadpole type" special vertices may have real propagators as edge (solid lines)*)
(*NB2 Be careful in the zeropoint graphs*)


replaceSunTad\[Tau][grp_]:=Module[{tmpgrp=Flatten@grp,listsunslop,vetexsl,pivots,rulesrename},
listsunslop=findsunloop@tmpgrp;
If[listsunslop=={},tmpgrp,
vetexsl=VertexList/@listsunslop;
pivots=Table[DeleteCases[vetexsl[[i]],Alternatives@@Flatten[Cases[vetexsl[[i]],\!\(\*SubscriptBox[\(s\), \({a_, b_}\)]\)]/.\!\(\*SubscriptBox[\(s\), \({a_, b_}\)]\):> {\!\(\*SubscriptBox[\(s\), \({a, b}\)]\),a,b}]],{i,Length@vetexsl}];
rulesrename=Table[pivots[[j,1]]-> Subscript[\[Tau],pivots[[j,1]]],{j,Length@pivots}];
(*Print[listsunslop]; Print[vetexsl]; Print[pivots]; Print[rulesrename];*)
Delete[tmpgrp,Position[Sort/@tmpgrp,Alternatives@@Sort/@Flatten@listsunslop]]/.rulesrename
]
]


replaceSunTad\[Tau]case0pt[grp_]:=Module[{tmp=replaceSunTad\[Tau][grp]},
If[tmp==={},If[(Length@findsunloop@Flatten@grp)==2,{Subscript[\[Tau], \[Phi]]\[UndirectedEdge]Subscript[\[Tau], \[Psi]]},"newcase"],tmp]]


(* ::Subsubsection::Closed:: *)
(*Tadpole with 2 bubbles and a triangle (\[Beta])*)


find6safe[grp_]:=If[grp=={},{},FindCycle[grp,{6},All]]


(* ::Text:: *)
(*Identifying the tadpole with 2 bubbles and a triangle (\[Beta])*)


specialverticesNobt=Flatten@{{Subscript[s, _],Subscript[c, _],Subscript[\[Psi]a, _]},{Subscript[\[Sigma], _],Subscript[\[Sigma]2, _],Subscript[\[Theta], _],Subscript[\[Theta]2, _],Subscript[\[Gamma], _],Subscript[\[Gamma]\[Kappa], _],Subscript[\[Gamma]c, _],Subscript[\[Psi]b, _],Subscript[\[Kappa], _]}};


findbubbletrloop[grp_]:=Select[If[grp=={},{},FindCycle[grp,{6},All]],(Count[#,\!\(\*SubscriptBox[\(b\), \({a_, b_}\)]\)\[UndirectedEdge]a_]+Count[#,\!\(\*SubscriptBox[\(b\), \({a_, b_}\)]\)\[UndirectedEdge]b_]==2&&Count[#,\!\(\*SubscriptBox[\(t\), \({a_, b_, c_}\)]\)\[UndirectedEdge]a_]+Count[#,\!\(\*SubscriptBox[\(t\), \({a_, b_, c_}\)]\)\[UndirectedEdge]b_]+Count[#,\!\(\*SubscriptBox[\(t\), \({a_, b_, c_}\)]\)\[UndirectedEdge]c_]==1&&!MemberQ[VertexList@Flatten@#,Alternatives@@specialverticesNobt])&]


replacebbt\[Beta][grp_]:=Module[{tmpgrp=Flatten@grp,listbubltrlop,vetexsl,pivots,notpivots,rulesrename,almostdelet,delposib},
listbubltrlop=findbubbletrloop@tmpgrp;
If[listbubltrlop=={},tmpgrp,
vetexsl=VertexList/@listbubltrlop;
pivots=Table[Select[Flatten[Cases[vetexsl[[i]],\!\(\*SubscriptBox[\(t\), \({a_, b_, c_}\)]\)]/.\!\(\*SubscriptBox[\(t\), \({a_, b_, c_}\)]\):> {a,b,c}],!MemberQ[vetexsl[[i]],#]&],{i,Length@vetexsl}];
notpivots=Table[Select[Flatten[Cases[vetexsl[[i]],\!\(\*SubscriptBox[\(t\), \({a_, b_, c_}\)]\)]/.\!\(\*SubscriptBox[\(t\), \({a_, b_, c_}\)]\):> {a,b,c}],MemberQ[vetexsl[[i]],#]&],{i,Length@vetexsl}];
rulesrename=Table[pivots[[j,1]]-> Subscript[\[Beta],pivots[[j,1]]],{j,Length@pivots}];
almostdelet=Delete[tmpgrp,Position[Sort/@tmpgrp,Alternatives@@Sort/@Flatten@listbubltrlop]]/.rulesrename;
delposib=Table[\!\(\*SubscriptBox[\(t\), \({Subscript[\[Beta], pivots[\([k, 1]\)]], notpivots[\([k, 1]\)], notpivots[\([k, 2]\)]}\)]\)\[UndirectedEdge]Subscript[\[Beta],pivots[[k,1]]]|\!\(\*SubscriptBox[\(t\), \({notpivots[\([k, 1]\)], Subscript[\[Beta], pivots[\([k, 1]\)]], notpivots[\([k, 2]\)]}\)]\)\[UndirectedEdge]Subscript[\[Beta],pivots[[k,1]]]|\!\(\*SubscriptBox[\(t\), \({notpivots[\([k, 1]\)], notpivots[\([k, 2]\)], Subscript[\[Beta], pivots[\([k, 1]\)]]}\)]\)\[UndirectedEdge]Subscript[\[Beta],pivots[[k,1]]],{k,Length@pivots}];
(*Print[listbubltrlop]; Print[vetexsl]; Print[pivots]; Print[notpivots]; Print[rulesrename];Print[delposib];*)
DeleteCases[almostdelet,Alternatives@@delposib]
]
]


(* ::Text:: *)
(*Composition of Complex Tadpoles*)


(* ::Input:: *)
(*(*replacecomplextadflat[grp_]:=replacebubtrloopnflat@replacesunloopnflat@grp*)*)


replaceTad\[Tau]\[Beta][grp_]:=Module[{tmp=replaceSunTad\[Tau]case0pt@grp,usualtmp},
usualtmp=replacebbt\[Beta]@tmp;
If[usualtmp==={},
If[(Length@findbubbletrloop@tmp)==1&&MemberQ[VertexList@tmp,Subscript[\[Tau], _]],{Subscript[\[Beta], \[Phi]]\[UndirectedEdge]Cases[VertexList@tmp,Subscript[\[Tau], _]][[1]]},If[(Length@findbubbletrloop@Flatten@grp)==2,{Subscript[\[Beta], \[Phi]]\[UndirectedEdge]Subscript[\[Beta], \[Psi]]},"newcase"]],usualtmp]
]


(* ::Subsection::Closed:: *)
(*Numerical vertices*)


(* ::Subsubsection::Closed:: *)
(*Sunset in square (\[Sigma])*)


(* ::Text:: *)
(*Identifying the sunset in square (\[Sigma])*)


specialverticesNos=Flatten@{{Subscript[b, _],Subscript[t, _],Subscript[c, _],Subscript[\[Psi]a, _]},{Subscript[\[Sigma], _],Subscript[\[Sigma]2, _],Subscript[\[Theta], _],Subscript[\[Theta]2, _],Subscript[\[Gamma], _],Subscript[\[Gamma]\[Kappa], _],Subscript[\[Gamma]c, _],Subscript[\[Psi]b, _],Subscript[\[Kappa], _]}};


findsunsq[grp_]:=Select[If[grp=={},{},FindCycle[grp,{5},All]],(Count[#,\!\(\*SubscriptBox[\(s\), \({a_, b_}\)]\)\[UndirectedEdge]a_]+Count[#,\!\(\*SubscriptBox[\(s\), \({a_, b_}\)]\)\[UndirectedEdge]b_]==1&&!MemberQ[VertexList@Flatten@#,Alternatives@@specialverticesNos])&]


(* ::Text:: *)
(*Replacement: identifying the pivots to be connected with the special vertex*)


replaceSunSq\[Sigma][grp_]:=Module[{tmpgrp=Flatten@grp,listsunsqall,listsunsq,vetexsl,pivots,rulesrenamever,rulesrenameinver,listsunsqsave},
listsunsqall=findsunsq@tmpgrp;
If[listsunsqall=={},tmpgrp,
listsunsq=extractlongest@maximcycles@listsunsqall;
vetexsl=VertexList/@listsunsq;
pivots=Table[DeleteCases[vetexsl[[i]],Alternatives@@Flatten[Cases[vetexsl[[i]],\!\(\*SubscriptBox[\(s\), \({a_, b_}\)]\)]/.\!\(\*SubscriptBox[\(s\), \({a_, b_}\)]\):> {\!\(\*SubscriptBox[\(s\), \({a, b}\)]\),a,b}]],{i,Length@vetexsl}];
rulesrenamever=Table[pivots[[j,1]]\[UndirectedEdge]pivots[[j,2]]->{pivots[[j,1]]\[UndirectedEdge]\!\(\*SubscriptBox[\(\[Sigma]\), \({pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket]}\)]\),\!\(\*SubscriptBox[\(\[Sigma]\), \({pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket]}\)]\)\[UndirectedEdge]pivots[[j,2]]},{j,Length@pivots}];
rulesrenameinver=Table[pivots[[j,2]]\[UndirectedEdge]pivots[[j,1]]->{pivots[[j,2]]\[UndirectedEdge]\!\(\*SubscriptBox[\(\[Sigma]\), \({pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket]}\)]\),\!\(\*SubscriptBox[\(\[Sigma]\), \({pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket]}\)]\)\[UndirectedEdge]pivots[[j,1]]},{j,Length@pivots}];
listsunsqsave=Table[DeleteCases[listsunsq[[j]],pivots[[j,1]]\[UndirectedEdge]pivots[[j,2]]|pivots[[j,2]]\[UndirectedEdge]pivots[[j,1]]],{j,Length@pivots}];
(*Print[listsunsq];Print[pivots];Print[listsunsqsave];Print[rulesrenamever];Print[rulesrenameinver];*)
Flatten@(Delete[tmpgrp,Position[Sort/@tmpgrp,Alternatives@@Sort/@Flatten@listsunsqsave]]/.rulesrenamever/.rulesrenameinver)
]
]


(* ::Subsubsection::Closed:: *)
(*Sunset in square with bubble (\[Sigma]2) (or sunset in sunset)*)


(* ::Text:: *)
(*Identifying the sunset in square with bubble (\[Sigma]2)*)


specialverticesNobs=Flatten@{{Subscript[t, _],Subscript[c, _],Subscript[\[Psi]a, _]},{Subscript[\[Sigma], _],Subscript[\[Sigma]2, _],Subscript[\[Theta], _],Subscript[\[Theta]2, _],Subscript[\[Gamma], _],Subscript[\[Gamma]\[Kappa], _],Subscript[\[Gamma]c, _],Subscript[\[Psi]b, _],Subscript[\[Kappa], _]}};


findsunsqb[grp_]:=Select[If[grp=={},{},FindCycle[grp,{6},All]],(Count[#,\!\(\*SubscriptBox[\(s\), \({a_, b_}\)]\)\[UndirectedEdge]a_]+Count[#,\!\(\*SubscriptBox[\(s\), \({a_, b_}\)]\)\[UndirectedEdge]b_]==1&&Count[#,\!\(\*SubscriptBox[\(b\), \({a_, b_}\)]\)\[UndirectedEdge]a_]+Count[#,\!\(\*SubscriptBox[\(b\), \({a_, b_}\)]\)\[UndirectedEdge]b_]==1&&!MemberQ[VertexList@Flatten@#,Alternatives@@specialverticesNobs])&]


(* ::Text:: *)
(*Replacement: identifying the pivots to be connected with the special vertex*)


replaceSunSqnB\[Sigma]2[grp_]:=Module[{tmpgrp=Flatten@grp,listsunsqall,listsunsq,vetexsl,pivots,rulesrename1,rulesrename2,rulesrename3,rulesrename4,rulesrename5,rulesrename6,rulesrename7,rulesrename8,listsunsqbsave},
listsunsqall=findsunsqb@tmpgrp;
If[listsunsqall=={},tmpgrp,
listsunsq=extractlongest@maximcycles@listsunsqall;
vetexsl=VertexList/@listsunsq;
pivots=Table[DeleteCases[vetexsl[[i]],Alternatives@@Flatten[Join[Cases[vetexsl[[i]],\!\(\*SubscriptBox[\(s\), \({a_, b_}\)]\)]/.\!\(\*SubscriptBox[\(s\), \({a_, b_}\)]\):> {\!\(\*SubscriptBox[\(s\), \({a, b}\)]\),a,b},{Subscript[b, _]}]]],{i,Length@vetexsl}];
rulesrename1=Table[\!\(\*SubscriptBox[\(b\), \({pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket]}\)]\)\[UndirectedEdge]pivots[[j,1]]->\!\(\*SubscriptBox[\(\[Sigma]2\), \({pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket]}\)]\)\[UndirectedEdge]pivots[[j,1]],{j,Length@pivots}];
rulesrename2=Table[\!\(\*SubscriptBox[\(b\), \({pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket]}\)]\)\[UndirectedEdge]pivots[[j,2]]->\!\(\*SubscriptBox[\(\[Sigma]2\), \({pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket]}\)]\)\[UndirectedEdge]pivots[[j,2]],{j,Length@pivots}];
rulesrename3=Table[pivots[[j,1]]\[UndirectedEdge]\!\(\*SubscriptBox[\(b\), \({pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket]}\)]\)->pivots[[j,1]]\[UndirectedEdge]\!\(\*SubscriptBox[\(\[Sigma]2\), \({pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket]}\)]\),{j,Length@pivots}];
rulesrename4=Table[pivots[[j,2]]\[UndirectedEdge]\!\(\*SubscriptBox[\(b\), \({pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket]}\)]\)->pivots[[j,2]]\[UndirectedEdge]\!\(\*SubscriptBox[\(\[Sigma]2\), \({pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket]}\)]\),{j,Length@pivots}];
rulesrename5=Table[\!\(\*SubscriptBox[\(b\), \({pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket]}\)]\)\[UndirectedEdge]pivots[[j,1]]->\!\(\*SubscriptBox[\(\[Sigma]2\), \({pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket]}\)]\)\[UndirectedEdge]pivots[[j,1]],{j,Length@pivots}];
rulesrename6=Table[\!\(\*SubscriptBox[\(b\), \({pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket]}\)]\)\[UndirectedEdge]pivots[[j,2]]->\!\(\*SubscriptBox[\(\[Sigma]2\), \({pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket]}\)]\)\[UndirectedEdge]pivots[[j,2]],{j,Length@pivots}];
rulesrename7=Table[pivots[[j,1]]\[UndirectedEdge]\!\(\*SubscriptBox[\(b\), \({pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket]}\)]\)->pivots[[j,1]]\[UndirectedEdge]\!\(\*SubscriptBox[\(\[Sigma]2\), 
StyleBox[
RowBox[{"{", 
RowBox[{
RowBox[{"pivots", "\[LeftDoubleBracket]", 
RowBox[{"j", ",", "2"}], "\[RightDoubleBracket]"}], ",", 
RowBox[{"pivots", "\[LeftDoubleBracket]", 
RowBox[{"j", ",", "1"}], "\[RightDoubleBracket]"}]}], "}"}],
FontWeight->"Bold"]]\),{j,Length@pivots}];
rulesrename8=Table[pivots[[j,2]]\[UndirectedEdge]\!\(\*SubscriptBox[\(b\), \({pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket]}\)]\)->pivots[[j,2]]\[UndirectedEdge]\!\(\*SubscriptBox[\(\[Sigma]2\), \({pivots\[LeftDoubleBracket]j, 2\[RightDoubleBracket], pivots\[LeftDoubleBracket]j, 1\[RightDoubleBracket]}\)]\),{j,Length@pivots}];
listsunsqbsave=Table[DeleteCases[listsunsq[[j]],\!\(\*SubscriptBox[\(b\), \({a_, e_}\)]\)\[UndirectedEdge]a_|\!\(\*SubscriptBox[\(b\), \({a_, e_}\)]\)\[UndirectedEdge]e_|a_\[UndirectedEdge]\!\(\*SubscriptBox[\(b\), \({a_, e_}\)]\)|e_\[UndirectedEdge]\!\(\*SubscriptBox[\(b\), \({a_, e_}\)]\)],{j,Length@listsunsq}];
(*Print[listsunsq];Print[pivots];Print[listsunsqsave];Print[rulesrenamever];Print[rulesrenameinver];*)
Flatten@(Delete[tmpgrp,Position[Sort/@tmpgrp,Alternatives@@Sort/@Flatten@listsunsqbsave]]/.rulesrename1/.rulesrename2/.rulesrename3/.rulesrename4/.rulesrename5/.rulesrename6/.rulesrename7/.rulesrename8)
]
]


(* ::Subsubsection::Closed:: *)
(*Triangle with 2 bubbles (\[Theta])*)


(* ::Text:: *)
(*Identifying the triangle with 2 bubbles (\[Theta])*)


specialverticesNob=Flatten@{{Subscript[s, _],Subscript[t, _],Subscript[c, _],Subscript[\[Psi]a, _]},{Subscript[\[Sigma], _],Subscript[\[Sigma]2, _],Subscript[\[Theta], _],Subscript[\[Theta]2, _],Subscript[\[Gamma], _],Subscript[\[Gamma]\[Kappa], _],Subscript[\[Gamma]c, _],Subscript[\[Psi]b, _],Subscript[\[Kappa], _]}};


findbbptr[grp_]:=Select[If[grp=={},{},FindCycle[grp,{5},All]],(Count[#,\!\(\*SubscriptBox[\(b\), \({a_, b_}\)]\)\[UndirectedEdge]a_]+Count[#,\!\(\*SubscriptBox[\(b\), \({a_, b_}\)]\)\[UndirectedEdge]b_]==2&&!MemberQ[VertexList@Flatten@#,Alternatives@@specialverticesNob])&]


(* ::Text:: *)
(*Replacement: identifying the pivot edge and substitute it with two connecting the new special vertex*)


replacebbpT\[Theta][grp_]:=Module[{tmpgrp=Flatten@grp,listbbptrall,listbbptr,vetexsl,pivotsprp,pivotsprpinv,replacedpivotsprp,replacedpivotsprpinv,rulereplace,rulereplaceinv},
listbbptrall=findbbptr@tmpgrp;
If[listbbptrall=={},tmpgrp,
listbbptr=extractlongest@maximcycles@listbbptrall;
vetexsl=VertexList/@listbbptr;
pivotsprp=Table[DeleteCases[listbbptr[[i]],\!\(\*SubscriptBox[\(b\), \({a_, e_}\)]\)\[UndirectedEdge]a_|\!\(\*SubscriptBox[\(b\), \({a_, e_}\)]\)\[UndirectedEdge]e_|a_\[UndirectedEdge]\!\(\*SubscriptBox[\(b\), \({a_, e_}\)]\)|e_\[UndirectedEdge]\!\(\*SubscriptBox[\(b\), \({a_, e_}\)]\)],{i,Length@listbbptr}];
pivotsprpinv=Table[pivotsprp[[j]]/.c_\[UndirectedEdge]d_:>d\[UndirectedEdge]c,{j,Length@pivotsprp}];
replacedpivotsprp=Table[pivotsprp[[j]]/.c_\[UndirectedEdge]d_:>{c\[UndirectedEdge]\!\(\*SubscriptBox[\(\[Theta]\), \({c, d}\)]\),\!\(\*SubscriptBox[\(\[Theta]\), \({c, d}\)]\)\[UndirectedEdge]d},{j,Length@pivotsprp}];
replacedpivotsprpinv=Table[pivotsprpinv[[j]]/.c_\[UndirectedEdge]d_:>{c\[UndirectedEdge]\!\(\*SubscriptBox[\(\[Theta]\), \({c, d}\)]\),\!\(\*SubscriptBox[\(\[Theta]\), \({c, d}\)]\)\[UndirectedEdge]d},{j,Length@pivotsprp}];
rulereplace=Table[pivotsprp[[k,1]]-> replacedpivotsprp[[k,1]],{k,Length@pivotsprp}];
rulereplaceinv=Table[pivotsprpinv[[k,1]]-> replacedpivotsprpinv[[k,1]],{k,Length@pivotsprp}];
(*Print[listbbptr];Print[pivotsprp];Print[rulereplace];*)
Flatten@(Delete[tmpgrp,Position[Sort/@tmpgrp,Alternatives@@Sort/@Flatten@Table[DeleteCases[listbbptr[[i]],pivotsprp[[i,1]]],{i,Length@listbbptr}]]]/.rulereplace/.rulereplaceinv)
]
]


(* ::Subsubsection::Closed:: *)
(*Square with 3 bubbles (\[Theta]2)*)


(* ::Text:: *)
(*Identifying the square with 3 bubbles (\[Theta]2)*)


findbbbpsq[grp_]:=Select[If[grp=={},{},FindCycle[grp,{7},All]],(Count[#,\!\(\*SubscriptBox[\(b\), \({a_, b_}\)]\)\[UndirectedEdge]a_]+Count[#,\!\(\*SubscriptBox[\(b\), \({a_, b_}\)]\)\[UndirectedEdge]b_]==3&&!MemberQ[VertexList@Flatten@#,Alternatives@@specialverticesNob])&]


(* ::Text:: *)
(*Replacement: identifying the pivot edge and substitute it with two connecting the new special vertex*)


replacebbbpSq\[Theta]2[grp_]:=Module[{tmpgrp=Flatten@grp,listbbptrall,listbbptr,vetexsl,pivotsprp,pivotsprpinv,replacedpivotsprp,replacedpivotsprpinv,rulereplace,rulereplaceinv},
listbbptrall=findbbbpsq@tmpgrp;
If[listbbptrall=={},tmpgrp,
listbbptr=extractlongest@maximcycles@listbbptrall;
vetexsl=VertexList/@listbbptr;
pivotsprp=Table[DeleteCases[listbbptr[[i]],\!\(\*SubscriptBox[\(b\), \({a_, e_}\)]\)\[UndirectedEdge]a_|\!\(\*SubscriptBox[\(b\), \({a_, e_}\)]\)\[UndirectedEdge]e_|a_\[UndirectedEdge]\!\(\*SubscriptBox[\(b\), \({a_, e_}\)]\)|e_\[UndirectedEdge]\!\(\*SubscriptBox[\(b\), \({a_, e_}\)]\)],{i,Length@listbbptr}];
pivotsprpinv=Table[pivotsprp[[j]]/.c_\[UndirectedEdge]d_:>d\[UndirectedEdge]c,{j,Length@pivotsprp}];
replacedpivotsprp=Table[pivotsprp[[j]]/.c_\[UndirectedEdge]d_:>{c\[UndirectedEdge]\!\(\*SubscriptBox[\(\[Theta]2\), \({c, d}\)]\),\!\(\*SubscriptBox[\(\[Theta]2\), \({c, d}\)]\)\[UndirectedEdge]d},{j,Length@pivotsprp}];
replacedpivotsprpinv=Table[pivotsprpinv[[j]]/.c_\[UndirectedEdge]d_:>{c\[UndirectedEdge]\!\(\*SubscriptBox[\(\[Theta]2\), \({c, d}\)]\),\!\(\*SubscriptBox[\(\[Theta]2\), \({c, d}\)]\)\[UndirectedEdge]d},{j,Length@pivotsprp}];
rulereplace=Table[pivotsprp[[k,1]]-> replacedpivotsprp[[k,1]],{k,Length@pivotsprp}];
rulereplaceinv=Table[pivotsprpinv[[k,1]]-> replacedpivotsprpinv[[k,1]],{k,Length@pivotsprp}];
(*Print[listbbptr];Print[pivotsprp];Print[rulereplace];*)
Flatten@(Delete[tmpgrp,Position[Sort/@tmpgrp,Alternatives@@Sort/@Flatten@Table[DeleteCases[listbbptr[[i]],pivotsprp[[i,1]]],{i,Length@listbbptr}]]]/.rulereplace/.rulereplaceinv)
]
]


(* ::Subsection::Closed:: *)
(*Composition on the substitutions*)


replaceComplexV\[Theta]2\[Theta]\[Sigma]2\[Sigma][grp_]:=replacebbbpSq\[Theta]2@replacebbpT\[Theta]@replaceSunSqnB\[Sigma]2@replaceSunSq\[Sigma][grp]


(* ::Subsubsection::Closed:: *)
(*Sunsets, Bubbles, Triangles, and Squares NEW SEQUENCE + \[Psi] & \[Kappa] + \[Tau], \[Beta], \[Theta], \[Theta]2, \[Sigma], \[Sigma]2*)


(* ::Text:: *)
(*Taking the original lists and applying the basic substitutions:*)
(**)
(*bstc4a for bubbles sunsets, triangles and squares, \[Psi]a, \[Psi]b and \[Tau], \[Beta], \[Theta], \[Theta]2, \[Sigma], \[Sigma]2*)
(*bstc4 for bubbles sunsets, triangles and squares, \[Psi]a, \[Psi]b, \[Kappa] and \[Tau], \[Beta], \[Theta], \[Theta]2, \[Sigma], \[Sigma]2*)


zeropt1PIbstc4a[n3_,n4_]:=Module[{tmp=zeropt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(replaceComplexV\[Theta]2\[Theta]\[Sigma]2\[Sigma]@replaceTad\[Tau]\[Beta]@subs\[Psi]BubTrianSquares@tmp[[i,2]])},{i,tt}]]


twopt1PIbstc4a[n3_,n4_]:=Module[{tmp=twopt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(replaceComplexV\[Theta]2\[Theta]\[Sigma]2\[Sigma]@replaceTad\[Tau]\[Beta]@subs\[Psi]BubTrianSquares@tmp[[i,2]])},{i,tt}]]


fourpt1PIbstc4a[n3_,n4_]:=Module[{tmp=fourpt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(replaceComplexV\[Theta]2\[Theta]\[Sigma]2\[Sigma]@replaceTad\[Tau]\[Beta]@subs\[Psi]BubTrianSquares@tmp[[i,2]])},{i,tt}]]


(* ::Input:: *)
(**)


zeropt1PIbstc4[n3_,n4_]:=Module[{tmp=zeropt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(replaceComplexV\[Theta]2\[Theta]\[Sigma]2\[Sigma]@replaceTad\[Tau]\[Beta]@subs\[Psi]\[Kappa]BubTrianSquares@tmp[[i,2]])},{i,tt}]]


twopt1PIbstc4[n3_,n4_]:=Module[{tmp=twopt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(replaceComplexV\[Theta]2\[Theta]\[Sigma]2\[Sigma]@replaceTad\[Tau]\[Beta]@subs\[Psi]\[Kappa]BubTrianSquares@tmp[[i,2]])},{i,tt}]]


fourpt1PIbstc4[n3_,n4_]:=Module[{tmp=fourpt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(replaceComplexV\[Theta]2\[Theta]\[Sigma]2\[Sigma]@replaceTad\[Tau]\[Beta]@subs\[Psi]\[Kappa]BubTrianSquares@tmp[[i,2]])},{i,tt}]]


(* ::Section:: *)
(*Visualization*)


(* ::Subsection::Closed:: *)
(*General functions*)


renameExternalPoints[edgelist_] := ReplaceAll[edgelist, 
	{GSberveglieri`Phi4tools`Private`o1->"o1",GSberveglieri`Phi4tools`Private`o2->"o2",GSberveglieri`Phi4tools`Private`o3->"o3",GSberveglieri`Phi4tools`Private`o4->"o4"}
]


replaceComplexStrings[grp_]:=ReplaceAll[
	grp,
	{ b -> "b", s -> "s", t -> "t", c -> "c", \[Tau] -> "\[Tau]", \[Beta] -> "\[Beta]", \[Sigma] -> "\[Sigma]", \[Sigma]2 -> "\[Sigma]2", \[Theta] -> "\[Theta]", \[Theta]2 -> "\[Theta]2",
	 \[Gamma] -> "\[Gamma]", \[Gamma]\[Kappa] -> "\[Gamma]\[Kappa]", \[Gamma]c -> "\[Gamma]c", \[Psi]a -> "\[Psi]a", \[Psi]b -> "\[Psi]b", \[Kappa] -> "\[Kappa]"}
]


(* ::Text:: *)
(*Visualization of the graphs, legend and examples follow*)


(* ::Input:: *)
(*(*specialverticessimple*)*)


specialverticessimple = {
	{Subscript["b", _], Subscript["s", _], Subscript["t", _], Subscript["c", _], Subscript["\[Psi]a", _]},
	{Subscript["\[Sigma]", _],Subscript["\[Sigma]2", _],Subscript["\[Theta]", _],Subscript["\[Theta]2", _],Subscript["\[Gamma]", _],
	 Subscript["\[Gamma]\[Kappa]", _],Subscript["\[Gamma]c", _],Subscript["\[Psi]b", _],Subscript["\[Kappa]", _]},
	{Subscript["\[Tau]", _],Subscript["\[Beta]", _]}};


egdeofv[x_]:={x\[UndirectedEdge]_,_\[UndirectedEdge]x}


visualizeGraph[grp_]:=Module[
	{graph=replaceComplexStrings[renameExternalPoints[Flatten@grp]], notprop, externalprop},
	notprop=Part[graph,Flatten@Position[graph,Alternatives@@Flatten[egdeofv/@Flatten@Join[specialverticessimple[[1;;2]],{Subscript["\[Tau]", "\[Phi]"]|Subscript["\[Beta]", "\[Phi]"]}]]]];
	externalprop=Part[graph,Flatten@Position[graph,"o1"\[UndirectedEdge]_|"o2"\[UndirectedEdge]_|"o3"\[UndirectedEdge]_|"o4"\[UndirectedEdge]_|_\[UndirectedEdge]"o1"|_\[UndirectedEdge]"o2"|_\[UndirectedEdge]"o3"|_\[UndirectedEdge]"o4"]];
	Graph[
		graph,
		VertexStyle -> Table[(Alternatives@@specialverticessimple[[i]])->{Green,Red,Blue}[[i]],{i,3}],
		VertexShapeFunction -> {
			Subscript["b", _]-> "Capsule",
			Subscript["s", _]-> "Star",
			Subscript["t", _]-> "Triangle",
			Subscript["c", _]-> "Square",
			Subscript["\[Tau]",_]-> "ConcaveDiamond",
			Subscript["\[Beta]",_]-> "ConcaveTriangle",
			Subscript["\[Sigma]", _]-> "ConcaveSquare",
			Subscript["\[Sigma]2", _]-> "ConcavePentagon",
			Subscript["\[Theta]", _]-> "RoundedRectangle",
			Subscript["\[Theta]2", _]-> "RoundedHexagon",
			Subscript["\[Gamma]", _]-> "Diamond",
			Subscript["\[Gamma]\[Kappa]", _]-> "Pentagon",
			Subscript["\[Gamma]c", _]-> "DownTrapezoid",
			Subscript["\[Psi]a", _]-> "RoundedTriangle",
			Subscript["\[Psi]b", _]-> "RoundedTriangle",
			Subscript["\[Kappa]", _]-> "ConcaveDiamond" },
		VertexSize -> 0.15,
		EdgeStyle -> Join[Table[notprop[[i]]->Dashed,{i,Length@notprop}],Table[externalprop[[i]]->Directive[Dotted,Thick],{i,Length@externalprop}]]
	]
]


(* ::Text:: *)
(*One function to visualize the graph of p-point function with no-substitutions, just the one already used or all those written*)


(* ::Subsection::Closed:: *)
(*Applications*)


visualizer[p_,n_,m_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PI[n,m],p==2,twopt1PI[n,m],p==4,fourpt1PI[n,m]][[All,2]];
visualizeGraph/@Flatten/@graphsnomult
]


visualizerSub1[p_,n_,m_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PIbstc[n,m],p==2,twopt1PIbstc[n,m],p==4,fourpt1PIbstc[n,m]][[All,2]];
visualizeGraph/@replaceTad\[Tau]\[Beta]/@graphsnomult
]


visualizerNew[p_,n_,m_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PIbstc[n,m],p==2,twopt1PIbstc[n,m],p==4,fourpt1PIbstc[n,m]][[All,2]];
visualizeGraph/@replaceComplexV\[Theta]2\[Theta]\[Gamma]c\[Gamma]\[Kappa]\[Gamma]\[Sigma]2\[Sigma]/@replaceTad\[Tau]\[Beta]/@graphsnomult
]


visualizerS[p_,n_,m_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PIs[n,m],p==2,twopt1PIs[n,m],p==4,fourpt1PIs[n,m]][[All,2]];
visualizeGraph/@graphsnomult
]


visualizerSub2[p_,n_,m_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PIbstc2[n,m],p==2,twopt1PIbstc2[n,m],p==4,fourpt1PIbstc2[n,m]][[All,2]];
visualizeGraph/@replaceTad\[Tau]\[Beta]/@graphsnomult
]


visualizerSub2a[p_,n_,m_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PIbstc2a[n,m],p==2,twopt1PIbstc2a[n,m],p==4,fourpt1PIbstc2a[n,m]][[All,2]];
visualizeGraph/@replaceTad\[Tau]\[Beta]/@graphsnomult
]


visualizer\[Psi][p_,n_,m_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PI\[Psi][n,m],p==2,twopt1PI\[Psi][n,m],p==4,fourpt1PI\[Psi][n,m]][[All,2]];
visualizeGraph/@graphsnomult
]


visualizer\[Psi]\[Kappa][p_,n_,m_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PI\[Psi]\[Kappa][n,m],p==2,twopt1PI\[Psi]\[Kappa][n,m],p==4,fourpt1PI\[Psi]\[Kappa][n,m]][[All,2]];
visualizeGraph/@graphsnomult
]


visualizerSub3a[p_,n_,m_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PIbstc3a[n,m],p==2,twopt1PIbstc3a[n,m],p==4,fourpt1PIbstc3a[n,m]][[All,2]];
visualizeGraph/@replaceTad\[Tau]\[Beta]/@graphsnomult
]


visualizerSub3[p_,n_,m_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PIbstc3[n,m],p==2,twopt1PIbstc3[n,m],p==4,fourpt1PIbstc3[n,m]][[All,2]];
visualizeGraph/@replaceTad\[Tau]\[Beta]/@graphsnomult
]


visualizerSub4a[p_,n_,m_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PIbstc4a[n,m],p==2,twopt1PIbstc4a[n,m],p==4,fourpt1PIbstc4a[n,m]][[All,2]];
visualizeGraph/@graphsnomult
]


visualizerSub4[p_,n_,m_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PIbstc3[n,m],p==2,twopt1PIbstc3[n,m],p==4,fourpt1PIbstc3[n,m]][[All,2]];
visualizeGraph/@replaceComplexV\[Theta]2\[Theta]\[Sigma]2\[Sigma]/@replaceTad\[Tau]\[Beta]/@graphsnomult
]


(* ::Subsection:: *)
(*Different layout*)


visualizeGraphSpring[grp_]:=Module[
	{graph=replaceComplexStrings[renameExternalPoints[Flatten@grp]], notprop, externalprop},
	notprop=Part[graph,Flatten@Position[graph, Alternatives@@Flatten[egdeofv/@Flatten@Join[specialverticessimple[[1;;2]],{Subscript["\[Tau]", "\[Phi]"]|Subscript["\[Beta]", "\[Phi]"]}]]]];
	externalprop=Part[graph,Flatten@Position[graph,"o1"\[UndirectedEdge]_|"o2"\[UndirectedEdge]_|"o3"\[UndirectedEdge]_|"o4"\[UndirectedEdge]_|_\[UndirectedEdge]"o1"|_\[UndirectedEdge]"o2"|_\[UndirectedEdge]"o3"|_\[UndirectedEdge]"o4"]];
	Graph[
		graph,
		VertexStyle -> Table[(Alternatives@@specialverticessimple[[i]])->{Green,Red,Blue}[[i]],{i,3}],
		VertexShapeFunction -> {
			Subscript["b", _]-> "Capsule",
			Subscript["s", _]-> "Star",
			Subscript["t", _]-> "Triangle",
			Subscript["c", _]-> "Square",
			Subscript["\[Tau]",_]-> "ConcaveDiamond",
			Subscript["\[Beta]",_]-> "ConcaveTriangle",
			Subscript["\[Sigma]", _]-> "ConcaveSquare",
			Subscript["\[Sigma]2", _]-> "ConcavePentagon",
			Subscript["\[Theta]", _]-> "RoundedRectangle",
			Subscript["\[Theta]2", _]-> "RoundedHexagon",
			Subscript["\[Gamma]", _]-> "Diamond",
			Subscript["\[Gamma]\[Kappa]", _]-> "Pentagon",
			Subscript["\[Gamma]c", _]-> "DownTrapezoid",
			Subscript["\[Psi]a", _]-> "RoundedTriangle",
			Subscript["\[Psi]b", _]-> "RoundedTriangle",
			Subscript["\[Kappa]", _]-> "ConcaveDiamond" },
		VertexSize -> 0.15,
		EdgeStyle -> Join[Table[notprop[[i]]->Dashed,{i,Length@notprop}],Table[externalprop[[i]]->Directive[Dotted,Thick],{i,Length@externalprop}]],
		GraphLayout->"SpringEmbedding"
	]
]


visualizeGraphsimpleSpring[grp_]:=Module[{extedges},
extedges=Select[grp,StringContainsQ[#[[2]],"e"]&];
Graph[grp,VertexSize->0.15,EdgeStyle->Table[extedges[[i]]->Directive[Dotted,Thick],{i,Length@extedges}],GraphLayout->"SpringEmbedding"]
]


visualizersp[p_,n_,m_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PI[n,m],p==2,twopt1PI[n,m],p==4,fourpt1PI[n,m]][[All,2]];
visualizeGraphSpring/@Flatten/@graphsnomult
]


visualizerspS[p_,n_,m_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PIs[n,m],p==2,twopt1PIs[n,m],p==4,fourpt1PIs[n,m]][[All,2]];
visualizeGraphSpring/@graphsnomult
]


visualizerspSub2[p_,n_,m_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PIbstc2[n,m],p==2,twopt1PIbstc2[n,m],p==4,fourpt1PIbstc2[n,m]][[All,2]];
visualizeGraphSpring/@replaceTad\[Tau]\[Beta]/@graphsnomult
]


visualizerspSub2a[p_,n_,m_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PIbstc2a[n,m],p==2,twopt1PIbstc2a[n,m],p==4,fourpt1PIbstc2a[n,m]][[All,2]];
visualizeGraphSpring/@replaceTad\[Tau]\[Beta]/@graphsnomult
]


visualizerspSub2W[p_,n_,m_,w_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PIbstc2W[n,m,w],p==2,twopt1PIbstc2W[n,m,w],p==4,fourpt1PIbstc2W[n,m,w]][[All,2]];
visualizeGraphSpring/@replaceTad\[Tau]\[Beta]/@graphsnomult
]


(* ::Section::Closed:: *)
(*Counting loops and Order diagrams*)


loopnumber[grp_]:= EdgeCount[grp]-VertexCount[grp]+1


(* ::Input:: *)
(*(*loopnumberwithsun[grp_]:= EdgeCount[grp]-VertexCount[grp]+1+VertexCount[grp,Subscript[s, _]]*)*)


(* ::Section::Closed:: *)
(*Weighted substitution*)


(* ::Text:: *)
(*In substitution of the main analytic blocks are preferred whose with the highest weight.*)
(*The weighted blocks are {b, t, c} (in this order).*)
(*Using W with w = {1,1,1} (neutralweights) it is equivalent to the functions without W*)
(*Other commonly used set of weight are invertedweights = {2,3,4} and standardweights = {4,2,1}*)


(* ::Subsubsection::Closed:: *)
(*Functions*)


extractlongestweightWeights[listlist_,ww_]:=Module[{leng,max,pos,bestones,weighttot},
leng=Length/@listlist;max=Max[leng];
pos=Position[leng,max];
(*Print[{max,pos}];*)
bestones=Extract[listlist,pos];
If[bestones=={},{},(*Print[bestones];*)weighttot=Total/@(Table[Length/@bestones[[i]],{i,Length@bestones}]/.{2-> -ww[[1]],3->-ww[[2]],4->-ww[[3]]});
(*Print[weighttot];*)(bestones[[Flatten@Position[weighttot,Min[weighttot]]]])[[1]]
(*bestones\[LeftDoubleBracket]1\[RightDoubleBracket]*)]]


extractlongestweightWeightsChoice[listlist_,ww_,cc_]:=Module[{leng,max,pos,bestones,weighttot,finalpos},
leng=Length/@listlist;max=Max[leng];
pos=Position[leng,max];
(*Print[{max,pos}];*)
bestones=Extract[listlist,pos];
If[bestones=={},{},(*Print[bestones];*)weighttot=Total/@(Table[Length/@bestones[[i]],{i,Length@bestones}]/.{2-> -ww[[1]],3->-ww[[2]],4->-ww[[3]]});
(*Print[weighttot];*)finalpos=bestones[[Flatten@Position[weighttot,Min[weighttot]]]];If[cc>Length@finalpos,"Choice not available",finalpos[[cc]]]
(*bestones\[LeftDoubleBracket]1\[RightDoubleBracket]*)]]


neutralweights={1,1,1};
standardweights={4,2,1};
invertedweights={2,3,4};


subsBubTrianSquaresW[grp_,w_]:=Module[{sunsDone,suns,sunsDonedelet,bestBubTriangSqu,toinsert,deleted},
sunsDone=subsjusts@findbns@grp;
suns=Cases[sunsDone,(x_/;MemberQ[x,Subscript[s, _],2])];
sunsDonedelet=keeprealvetexs@Flatten@subsjusts@findbns@grp;
bestBubTriangSqu=extractlongestweightWeights[maximcyclestransitionsnoinside[find1loop234safe@sunsDonedelet,sunsDonedelet],w];
toinsert=bubtriansquarforinsert@bestBubTriangSqu;
deleted=cango[Sort/@sunsDonedelet,Sort/@Join@@bestBubTriangSqu];
(*Print[sunsDone];Print[bestBubTriangSqu];Print[toinsert];Print[deleted];*)
If[deleted[[1]]==1,Message[subsBubTrianSquaresW::substsubdiagr];,
Join[deleted[[2]],suns,toinsert]]/.a_\[UndirectedEdge]o1:>o1\[UndirectedEdge]a/.a_\[UndirectedEdge]o2:>o2\[UndirectedEdge]a/.a_\[UndirectedEdge]o3:>o3\[UndirectedEdge]a/.a_\[UndirectedEdge]o4:>o4\[UndirectedEdge]a
]


subs\[Psi]aBubTrianSquaresW[grp_,w_]:=Module[{suns\[Psi]Done,suns\[Psi],suns\[Psi]Donedelet,bestBubTriangSqu,toinsert,deleted},
suns\[Psi]Done=subs\[Psi]a@grp;
suns\[Psi]=Cases[suns\[Psi]Done,(x_/;MemberQ[x,Subscript[s, _],2])|(x_/;MemberQ[x,Subscript[\[Psi]a, _],2])];
suns\[Psi]Donedelet=keeprealvetexs\[Psi]a@Flatten@subs\[Psi]a@grp;
bestBubTriangSqu=extractlongestweightWeights[maximcyclestransitionsnoinside[find1loop234safe@suns\[Psi]Donedelet,suns\[Psi]Donedelet],w];
toinsert=bubtriansquarforinsert@bestBubTriangSqu;
deleted=cango[Sort/@suns\[Psi]Donedelet,Sort/@Join@@bestBubTriangSqu];
(*Print[sunsDone];Print[bestBubTriangSqu];Print[toinsert];Print[deleted];*)
If[deleted[[1]]==1,Message[subs\[Psi]aBubTrianSquaresW::substsubdiagr];,
Join[deleted[[2]],suns\[Psi],toinsert]]/.a_\[UndirectedEdge]o1:>o1\[UndirectedEdge]a/.a_\[UndirectedEdge]o2:>o2\[UndirectedEdge]a/.a_\[UndirectedEdge]o3:>o3\[UndirectedEdge]a/.a_\[UndirectedEdge]o4:>o4\[UndirectedEdge]a
]


subs\[Psi]BubTrianSquaresW[grp_,w_]:=Module[{suns\[Psi]Done,suns\[Psi],suns\[Psi]Donedelet,bestBubTriangSqu,toinsert,deleted},
suns\[Psi]Done=subs\[Psi]ab@grp;
suns\[Psi]=Cases[suns\[Psi]Done,(x_/;MemberQ[x,Subscript[s, _],2])|(x_/;MemberQ[x,Subscript[\[Psi]a, _],2])|(x_/;MemberQ[x,Subscript[\[Psi]b, _],2])];
suns\[Psi]Donedelet=keeprealvetexs\[Psi]@Flatten@subs\[Psi]ab@grp;
bestBubTriangSqu=extractlongestweightWeights[maximcyclestransitionsnoinside[find1loop234safe@suns\[Psi]Donedelet,suns\[Psi]Donedelet],w];
toinsert=bubtriansquarforinsert@bestBubTriangSqu;
deleted=cango[Sort/@suns\[Psi]Donedelet,Sort/@Join@@bestBubTriangSqu];
(*Print[sunsDone];Print[bestBubTriangSqu];Print[toinsert];Print[deleted];*)
If[deleted[[1]]==1,Message[subs\[Psi]BubTrianSquaresW::substsubdiagr];,
Join[deleted[[2]],suns\[Psi],toinsert]]/.a_\[UndirectedEdge]o1:>o1\[UndirectedEdge]a/.a_\[UndirectedEdge]o2:>o2\[UndirectedEdge]a/.a_\[UndirectedEdge]o3:>o3\[UndirectedEdge]a/.a_\[UndirectedEdge]o4:>o4\[UndirectedEdge]a
]


subs\[Psi]\[Kappa]BubTrianSquaresW[grp_,w_]:=Module[{suns\[Psi]\[Kappa]Done,suns\[Psi]\[Kappa],suns\[Psi]\[Kappa]Donedelet,bestBubTriangSqu,toinsert,deleted},
suns\[Psi]\[Kappa]Done=subs\[Psi]ab\[Kappa]@grp;
suns\[Psi]\[Kappa]=Cases[suns\[Psi]\[Kappa]Done,(x_/;MemberQ[x,Subscript[s, _],2])|(x_/;MemberQ[x,Subscript[\[Psi]a, _],2])|(x_/;MemberQ[x,Subscript[\[Psi]b, _],2])|(x_/;MemberQ[x,Subscript[\[Kappa], _],2])];
suns\[Psi]\[Kappa]Donedelet=keeprealvetexs\[Psi]\[Kappa]@Flatten@subs\[Psi]ab\[Kappa]@grp;
bestBubTriangSqu=extractlongestweightWeights[maximcyclestransitionsnoinside[find1loop234safe@suns\[Psi]\[Kappa]Donedelet,suns\[Psi]\[Kappa]Donedelet],w];
toinsert=bubtriansquarforinsert@bestBubTriangSqu;
deleted=cango[Sort/@suns\[Psi]\[Kappa]Donedelet,Sort/@Join@@bestBubTriangSqu];
(*Print[sunsDone];Print[bestBubTriangSqu];Print[toinsert];Print[deleted];*)
If[deleted[[1]]==1,Message[subs\[Psi]\[Kappa]BubTrianSquaresW::substsubdiagr];,
Join[deleted[[2]],suns\[Psi]\[Kappa],toinsert]]/.a_\[UndirectedEdge]o1:>o1\[UndirectedEdge]a/.a_\[UndirectedEdge]o2:>o2\[UndirectedEdge]a/.a_\[UndirectedEdge]o3:>o3\[UndirectedEdge]a/.a_\[UndirectedEdge]o4:>o4\[UndirectedEdge]a
]


subsBubTrianSquaresW::substsubdiagr = "Error."


subs\[Psi]aBubTrianSquaresW::substsubdiagr = "Error."


subs\[Psi]BubTrianSquaresW::substsubdiagr = "Error."


subs\[Psi]\[Kappa]BubTrianSquaresW::substsubdiagr = "Error."


(* ::Subsubsection::Closed:: *)
(*Application level 2*)


zeropt1PIbstc2W[n3_,n4_,w_]:=Module[{tmp=zeropt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsBubTrianSquaresW[tmp[[i,2]],w])},{i,tt}]]


twopt1PIbstc2W[n3_,n4_,w_]:=Module[{tmp=twopt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsBubTrianSquaresW[tmp[[i,2]],w])},{i,tt}]]


fourpt1PIbstc2W[n3_,n4_,w_]:=Module[{tmp=fourpt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsBubTrianSquaresW[tmp[[i,2]],w])},{i,tt}]]


visualizerSub2W[p_,n_,m_,w_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PIbstc2W[n,m,w],p==2,twopt1PIbstc2W[n,m,w],p==4,fourpt1PIbstc2W[n,m,w]][[All,2]];
visualizeGraph/@replaceTad\[Tau]\[Beta]/@graphsnomult
]


(* ::Input:: *)
(*(*twopt1PIbstc2[0,6]==twopt1PIbstc2W[0,6,neutralweights]*)*)


(* ::Input:: *)
(*(*visualizerSub2[2,0,4]*)
(*visualizerSub2W[2,0,4,standardweights]*)*)


(* ::Input:: *)
(*(*{visualizer[0,0,4],visualizerSub2[0,0,4]}*)*)


(* ::Input:: *)
(*(*zeropt1PIbstc2W[0,4,{2,2,4}]*)
(*zeropt1PIbstc2W[0,4,{2,2,1}]*)*)


(* ::Input:: *)
(*(*Transpose[{visualizerSub2W[2,0,4,invertedweights],visualizerSub2W[2,0,4,standardweights]}]*)*)


(* ::Subsubsection::Closed:: *)
(*Application level 4*)


zeropt1PIbstc4W[n3_,n4_,w_]:=Module[{tmp=zeropt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(replaceComplexV\[Theta]2\[Theta]\[Sigma]2\[Sigma]@replaceTad\[Tau]\[Beta]@subs\[Psi]\[Kappa]BubTrianSquaresW[tmp[[i,2]],w])},{i,tt}]]


twopt1PIbstc4W[n3_,n4_,w_]:=Module[{tmp=twopt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(replaceComplexV\[Theta]2\[Theta]\[Sigma]2\[Sigma]@replaceTad\[Tau]\[Beta]@subs\[Psi]\[Kappa]BubTrianSquaresW[tmp[[i,2]],w])},{i,tt}]]


fourpt1PIbstc4W[n3_,n4_,w_]:=Module[{tmp=fourpt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(replaceComplexV\[Theta]2\[Theta]\[Sigma]2\[Sigma]@replaceTad\[Tau]\[Beta]@subs\[Psi]\[Kappa]BubTrianSquaresW[tmp[[i,2]],w])},{i,tt}]]


visualizerSub4W[p_,n_,m_,w_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PIbstc4W[n,m,w],p==2,twopt1PIbstc4W[n,m,w],p==4,fourpt1PIbstc4W[n,m,w]][[All,2]];
visualizeGraph/@replaceComplexV\[Theta]2\[Theta]\[Sigma]2\[Sigma]/@replaceTad\[Tau]\[Beta]/@graphsnomult
]


(* ::Input:: *)
(*(*twopt1PIbstc4[0,6]==twopt1PIbstc4W[0,6,neutralweights]*)*)


(* ::Input:: *)
(*(*Transpose@{visualizer[2,0,5],visualizerSub4W[2,0,5,neutralweights],visualizerSub4W[2,0,5,standardweights],visualizerSub4W[2,0,5,invertedweights]}*)*)


(* ::Subsubsection::Closed:: *)
(*Weights and choices (possibility to choose other configurations if there are several minima)*)


extractlongestweightWeightsChoice[listlist_,ww_,cc_]:=Module[{leng,max,pos,bestones,weighttot,finalpos},
leng=Length/@listlist;max=Max[leng];
pos=Position[leng,max];
(*Print[{max,pos}];*)
bestones=Extract[listlist,pos];
If[bestones=={},{},(*Print[bestones];*)weighttot=Total/@(Table[Length/@bestones[[i]],{i,Length@bestones}]/.{2-> -ww[[1]],3->-ww[[2]],4->-ww[[3]]});
(*Print[weighttot];*)finalpos=bestones[[Flatten@Position[weighttot,Min[weighttot]]]];If[cc>Length@finalpos,"Choice not available",finalpos[[cc]]]
(*bestones\[LeftDoubleBracket]1\[RightDoubleBracket]*)]]


subsBubTrianSquaresWC[grp_,w_,cc_]:=Module[{sunsDone,suns,sunsDonedelet,bestBubTriangSqu,toinsert,deleted},
sunsDone=subsjusts@findbns@grp;
suns=Cases[sunsDone,(x_/;MemberQ[x,Subscript[s, _],2])];
sunsDonedelet=keeprealvetexs@Flatten@subsjusts@findbns@grp;
bestBubTriangSqu=extractlongestweightWeightsChoice[maximcyclestransitionsnoinside[find1loop234safe@sunsDonedelet,sunsDonedelet],w,cc];
toinsert=bubtriansquarforinsert@bestBubTriangSqu;
deleted=cango[Sort/@sunsDonedelet,Sort/@Join@@bestBubTriangSqu];
(*Print[sunsDone];Print[bestBubTriangSqu];Print[toinsert];Print[deleted];*)
If[deleted[[1]]==1,Message[subsBubTrianSquaresWC::substsubdiagr];,
Join[deleted[[2]],suns,toinsert]]/.a_\[UndirectedEdge]o1:>o1\[UndirectedEdge]a/.a_\[UndirectedEdge]o2:>o2\[UndirectedEdge]a/.a_\[UndirectedEdge]o3:>o3\[UndirectedEdge]a/.a_\[UndirectedEdge]o4:>o4\[UndirectedEdge]a
]


subs\[Psi]BubTrianSquaresWC[grp_,w_,cc_]:=Module[{suns\[Psi]Done,suns\[Psi],suns\[Psi]Donedelet,bestBubTriangSqu,toinsert,deleted},
suns\[Psi]Done=subs\[Psi]ab@grp;
suns\[Psi]=Cases[suns\[Psi]Done,(x_/;MemberQ[x,Subscript[s, _],2])|(x_/;MemberQ[x,Subscript[\[Psi]a, _],2])|(x_/;MemberQ[x,Subscript[\[Psi]b, _],2])];
suns\[Psi]Donedelet=keeprealvetexs\[Psi]@Flatten@subs\[Psi]ab@grp;
bestBubTriangSqu=extractlongestweightWeightsChoice[maximcyclestransitionsnoinside[find1loop234safe@suns\[Psi]Donedelet,suns\[Psi]Donedelet],w,cc];
toinsert=bubtriansquarforinsert@bestBubTriangSqu;
deleted=cango[Sort/@suns\[Psi]Donedelet,Sort/@Join@@bestBubTriangSqu];
(*Print[sunsDone];Print[bestBubTriangSqu];Print[toinsert];Print[deleted];*)
If[deleted[[1]]==1,Message[subs\[Psi]BubTrianSquaresWC::substsubdiagr];,
Join[deleted[[2]],suns\[Psi],toinsert]]/.a_\[UndirectedEdge]o1:>o1\[UndirectedEdge]a/.a_\[UndirectedEdge]o2:>o2\[UndirectedEdge]a/.a_\[UndirectedEdge]o3:>o3\[UndirectedEdge]a/.a_\[UndirectedEdge]o4:>o4\[UndirectedEdge]a
]


subs\[Psi]\[Kappa]BubTrianSquaresWC[grp_,w_,cc_]:=Module[{suns\[Psi]\[Kappa]Done,suns\[Psi]\[Kappa],suns\[Psi]\[Kappa]Donedelet,bestBubTriangSqu,toinsert,deleted},
suns\[Psi]\[Kappa]Done=subs\[Psi]ab\[Kappa]@grp;
suns\[Psi]\[Kappa]=Cases[suns\[Psi]\[Kappa]Done,(x_/;MemberQ[x,Subscript[s, _],2])|(x_/;MemberQ[x,Subscript[\[Psi]a, _],2])|(x_/;MemberQ[x,Subscript[\[Psi]b, _],2])|(x_/;MemberQ[x,Subscript[\[Kappa], _],2])];
suns\[Psi]\[Kappa]Donedelet=keeprealvetexs\[Psi]\[Kappa]@Flatten@subs\[Psi]ab\[Kappa]@grp;
bestBubTriangSqu=extractlongestweightWeightsChoice[maximcyclestransitionsnoinside[find1loop234safe@suns\[Psi]\[Kappa]Donedelet,suns\[Psi]\[Kappa]Donedelet],w,cc];
toinsert=bubtriansquarforinsert@bestBubTriangSqu;
deleted=cango[Sort/@suns\[Psi]\[Kappa]Donedelet,Sort/@Join@@bestBubTriangSqu];
(*Print[sunsDone];Print[bestBubTriangSqu];Print[toinsert];Print[deleted];*)
If[deleted[[1]]==1,Message[subs\[Psi]\[Kappa]BubTrianSquaresWC::substsubdiagr];,
Join[deleted[[2]],suns\[Psi]\[Kappa],toinsert]]/.a_\[UndirectedEdge]o1:>o1\[UndirectedEdge]a/.a_\[UndirectedEdge]o2:>o2\[UndirectedEdge]a/.a_\[UndirectedEdge]o3:>o3\[UndirectedEdge]a/.a_\[UndirectedEdge]o4:>o4\[UndirectedEdge]a
]


zeropt1PIbstc2WC[n3_,n4_,w_,cc_]:=Module[{tmp=zeropt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsBubTrianSquaresWC[tmp[[i,2]],w,cc])},{i,tt}]]


twopt1PIbstc2WC[n3_,n4_,w_,cc_]:=Module[{tmp=twopt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsBubTrianSquaresWC[tmp[[i,2]],w,cc])},{i,tt}]]


fourpt1PIbstc2WC[n3_,n4_,w_,cc_]:=Module[{tmp=fourpt1PI[n3,n4],tt},
tt=Length@tmp;
Table[{tmp[[i,1]],(subsBubTrianSquaresWC[tmp[[i,2]],w,cc])},{i,tt}]]


(* ::Input:: *)
(*(*twopt1PIbstc2W[0,6,neutralweights]==twopt1PIbstc2WC[0,6,neutralweights,1]*)*)


(* ::Input:: *)
(*(*twopt1PIbstc2W[0,4,standardweights]*)*)


(* ::Input:: *)
(*(*twopt1PIbstc2WC[0,4,standardweights,2]*)*)


visualizerSub2WC[p_,n_,m_,w_,cc_]:=Module[{graphsnomult},
graphsnomult=Which[p==0,zeropt1PIbstc2WC[n,m,w,cc],p==2,twopt1PIbstc2WC[n,m,w,cc],p==4,fourpt1PIbstc2WC[n,m,w,cc]][[All,2]];
visualizeGraph/@replaceTad\[Tau]\[Beta]/@graphsnomult
]


(* ::Input:: *)
(*(*visualizerSub2WC[4,0,7,{2,2,3},1]\[LeftDoubleBracket]order\[CapitalGamma]4o7\[RightDoubleBracket]\[LeftDoubleBracket]162\[RightDoubleBracket]*)*)


(* ::Input:: *)
(*(*visualizerSub2WC[4,0,7,{2,3,3},1]\[LeftDoubleBracket]order\[CapitalGamma]4o7\[RightDoubleBracket]\[LeftDoubleBracket]162\[RightDoubleBracket]*)*)


(* ::Input:: *)
(*(*visualizerSub2WC[4,0,7,{3,3,4},1]\[LeftDoubleBracket]order\[CapitalGamma]4o7\[RightDoubleBracket]\[LeftDoubleBracket]162\[RightDoubleBracket]*)*)


subsBubTrianSquaresWC::substsubdiagr = "Error."


subs\[Psi]BubTrianSquaresWC::substsubdiagr = "Error."


subs\[Psi]\[Kappa]BubTrianSquaresWC::substsubdiagr = "Error."


(* ::Chapter:: *)
(*Integrand (from Complex graphs to Integrands)*)


(* ::Section::Closed:: *)
(*Momentum assignation (from Complex graphs to Integrands with q-momenta)*)


(* ::Text:: *)
(*Different level indicate different substitutions, this is valid for the writing and visualizing the graph as well*)
(*Level 1: {b, s, t, c, \[Tau], \[Beta]} (with the old algorithm)*)
(*Level 2: {b, s, t, c, \[Tau], \[Beta]} (with the new algorithm)*)
(*Level 2W: {b, s, t, c, \[Tau], \[Beta]}*)
(*Level 2a: {b, s, t, c, \[Tau], \[Beta], \[Psi]a} *)
(*Level 3a: {b, s, t, c, \[Tau], \[Beta], \[Psi]a, \[Psi]b} *)
(*Level 3: {b, s, t, c, \[Tau], \[Beta], \[Psi]a, \[Psi]b, \[Kappa]}*)
(*Level 4a: {b, s, t, c, \[Tau], \[Beta], \[Psi]a, \[Psi]b, \[Theta], \[Theta]2, \[Sigma], \[Sigma]2} *)
(*Level 4: {b, s, t, c, \[Tau], \[Beta], \[Psi]a, \[Psi]b, \[Kappa], \[Theta], \[Theta]2, \[Sigma], \[Sigma]2} *)
(*Level 4W: {b, s, t, c, \[Tau], \[Beta], \[Psi]a, \[Psi]b, \[Kappa], \[Theta], \[Theta]2, \[Sigma], \[Sigma]2} *)


(* ::Subsubsection::Closed:: *)
(*Functions*)


countingpro2[list_]:=Table[Min[(Dimensions@(list[[i]]))/.{}/;(list[[i]]===0)->{0}/.{}-> 1,(Dimensions@(-list[[i]]))/.{}/;(list[[i]]===0)->{0}/.{}-> 1],{i,Length@list}]


shiftmomok[list_]:=Module[{lst=list,countlst},
countlst=Table[Min[(Dimensions@(lst[[i]]))/.{}/;(lst[[i]]===0)->{0}/.{}-> 1,(Dimensions@(-lst[[i]]))/.{}/;(lst[[i]]===0)->{0}/.{}-> 1],{i,Length@lst}];
RotateLeft[lst,If[MemberQ[countlst,0],0,(Position[countlst, Max[countlst]])[[-1]]]]
]


(* ::Subsection::Closed:: *)
(*\[CapitalGamma]^(4)*)


(* ::Subsubsection::Closed:: *)
(*Level given*)


GtP4ptslevelgivennocompltad[n3_,n4_,edges_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=edges;
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=Flatten@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-4]];
paux=Insert[paux,pext1,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,pext2,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
paux=Insert[paux,pext3,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o3]&][[1]]][[1,1]]];
paux=Insert[paux,-pext1-pext2-pext3,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o4]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext1|pext2|pext3];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}]
),{lj,1,dl}]]


GtP4ptslevelgiven[n3_,n4_,edges_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=edges;
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-4]];
paux=Insert[paux,pext1,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,pext2,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
paux=Insert[paux,pext3,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o3]&][[1]]][[1,1]]];
paux=Insert[paux,-pext1-pext2-pext3,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o4]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext1|pext2|pext3];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}]
),{lj,1,dl}]]


(* ::Subsubsection::Closed:: *)
(*Level 2, 2W, 2Wass and 2a*)


GtP4ptslevel2[n3_,n4_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=fourpt1PIbstc2[n3,n4];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-4]];
paux=Insert[paux,pext1,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,pext2,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
paux=Insert[paux,pext3,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o3]&][[1]]][[1,1]]];
paux=Insert[paux,-pext1-pext2-pext3,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o4]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext1|pext2|pext3];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}]
),{lj,1,dl}]]


GtP4ptslevel2W[n3_,n4_,w_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=fourpt1PIbstc2W[n3,n4,w];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-4]];
paux=Insert[paux,pext1,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,pext2,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
paux=Insert[paux,pext3,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o3]&][[1]]][[1,1]]];
paux=Insert[paux,-pext1-pext2-pext3,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o4]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext1|pext2|pext3];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}]
),{lj,1,dl}]]


GtP4ptslevel2Wass[n3_,n4_,w_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=fourpt1PIbstc2W[n3,n4,w];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-4]];
paux=Insert[paux,pext1,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,pext2,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
paux=Insert[paux,pext3,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o3]&][[1]]][[1,1]]];
paux=Insert[paux,-pext1-pext2-pext3,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o4]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas,paux[[5;;-1]]])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext1|pext2|pext3];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}]
),{lj,1,dl}]]


GtP4ptslevel2a[n3_,n4_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=fourpt1PIbstc2a[n3,n4];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-4]];
paux=Insert[paux,pext1,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,pext2,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
paux=Insert[paux,pext3,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o3]&][[1]]][[1,1]]];
paux=Insert[paux,-pext1-pext2-pext3,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o4]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext1|pext2|pext3];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}]
),{lj,1,dl}]]


(* ::Subsubsection::Closed:: *)
(*Level 3a, 3 and 4a*)


GtP4ptslevel3a[n3_,n4_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=fourpt1PIbstc3a[n3,n4];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-4]];
paux=Insert[paux,pext1,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,pext2,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
paux=Insert[paux,pext3,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o3]&][[1]]][[1,1]]];
paux=Insert[paux,-pext1-pext2-pext3,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o4]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext1|pext2|pext3];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}]
),{lj,1,dl}]]


GtP4ptslevel3[n3_,n4_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,\[Kappa]vert,vertex\[Kappa]pos,analytic\[Kappa],squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=fourpt1PIbstc3[n3,n4];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-4]];
paux=Insert[paux,pext1,Position[daux,Select[daux,MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,pext2,Position[daux,Select[daux,MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),o2]&][[1]]][[1,1]]];
paux=Insert[paux,pext3,Position[daux,Select[daux,MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),o3]&][[1]]][[1,1]]];
paux=Insert[paux,-pext1-pext2-pext3,Position[daux,Select[daux,MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),o4]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext1|pext2|pext3];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|Subscript[\[Kappa], _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

\[Kappa]vert=Cases[pos,Subscript[\[Kappa], _]];
vertex\[Kappa]pos=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),\[Kappa]vert[[i]]]&)],{i,Length@\[Kappa]vert}];
analytic\[Kappa]=(\[Kappa]fnc@@@Table[Part[paux3,vertex\[Kappa]pos[[j]]],{j,Length@vertex\[Kappa]pos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Kappa] analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}]
),{lj,1,dl}]]


GtP4ptslevel4a[n3_,n4_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,\[Kappa]vert,vertex\[Kappa]pos,analytic\[Kappa],squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=fourpt1PIbstc4a[n3,n4];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-4]];
paux=Insert[paux,pext1,Position[daux,Select[daux,MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,pext2,Position[daux,Select[daux,MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),o2]&][[1]]][[1,1]]];
paux=Insert[paux,pext3,Position[daux,Select[daux,MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),o3]&][[1]]][[1,1]]];
paux=Insert[paux,-pext1-pext2-pext3,Position[daux,Select[daux,MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),o4]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext1|pext2|pext3];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|Subscript[\[Kappa], _]|Subscript[\[Theta], _]|Subscript[\[Theta]2, _]|Subscript[\[Sigma], _]|Subscript[\[Sigma]2, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Theta], _]\[UndirectedEdge]_],(t2b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Theta], _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Theta]2, _]\[UndirectedEdge]_],(c3b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Theta]2, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Sigma], _]\[UndirectedEdge]_],(c1s/@(Part[paux3,Flatten@Position[daux,Subscript[\[Sigma], _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Sigma]2, _]\[UndirectedEdge]_],(c1s1b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Sigma]2, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother * analytic\[Psi]b * analytic\[Psi]a * analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}]
),{lj,1,dl}]]


(* ::Subsubsection::Closed:: *)
(*Level 4 and 4W*)


GtP4ptslevel4[n3_,n4_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,\[Kappa]vert,vertex\[Kappa]pos,analytic\[Kappa],squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=fourpt1PIbstc4[n3,n4];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-4]];
paux=Insert[paux,pext1,Position[daux,Select[daux,MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,pext2,Position[daux,Select[daux,MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),o2]&][[1]]][[1,1]]];
paux=Insert[paux,pext3,Position[daux,Select[daux,MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),o3]&][[1]]][[1,1]]];
paux=Insert[paux,-pext1-pext2-pext3,Position[daux,Select[daux,MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),o4]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext1|pext2|pext3];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|Subscript[\[Kappa], _]|Subscript[\[Theta], _]|Subscript[\[Theta]2, _]|Subscript[\[Sigma], _]|Subscript[\[Sigma]2, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

\[Kappa]vert=Cases[pos,Subscript[\[Kappa], _]];
vertex\[Kappa]pos=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),\[Kappa]vert[[i]]]&)],{i,Length@\[Kappa]vert}];
analytic\[Kappa]=(\[Kappa]fnc@@@Table[Part[paux3,vertex\[Kappa]pos[[j]]],{j,Length@vertex\[Kappa]pos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Theta], _]\[UndirectedEdge]_],(t2b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Theta], _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Theta]2, _]\[UndirectedEdge]_],(c3b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Theta]2, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Sigma], _]\[UndirectedEdge]_],(c1s/@(Part[paux3,Flatten@Position[daux,Subscript[\[Sigma], _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Sigma]2, _]\[UndirectedEdge]_],(c1s1b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Sigma]2, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Kappa] analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}]
),{lj,1,dl}]]


GtP4ptslevel4W[n3_,n4_,w_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,\[Kappa]vert,vertex\[Kappa]pos,analytic\[Kappa],squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=fourpt1PIbstc4W[n3,n4,w];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-4]];
paux=Insert[paux,pext1,Position[daux,Select[daux,MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,pext2,Position[daux,Select[daux,MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),o2]&][[1]]][[1,1]]];
paux=Insert[paux,pext3,Position[daux,Select[daux,MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),o3]&][[1]]][[1,1]]];
paux=Insert[paux,-pext1-pext2-pext3,Position[daux,Select[daux,MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),o4]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext1|pext2|pext3];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|Subscript[\[Kappa], _]|Subscript[\[Theta], _]|Subscript[\[Theta]2, _]|Subscript[\[Sigma], _]|Subscript[\[Sigma]2, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

\[Kappa]vert=Cases[pos,Subscript[\[Kappa], _]];
vertex\[Kappa]pos=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),\[Kappa]vert[[i]]]&)],{i,Length@\[Kappa]vert}];
analytic\[Kappa]=(\[Kappa]fnc@@@Table[Part[paux3,vertex\[Kappa]pos[[j]]],{j,Length@vertex\[Kappa]pos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Theta], _]\[UndirectedEdge]_],(t2b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Theta], _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Theta]2, _]\[UndirectedEdge]_],(c3b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Theta]2, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Sigma], _]\[UndirectedEdge]_],(c1s/@(Part[paux3,Flatten@Position[daux,Subscript[\[Sigma], _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Sigma]2, _]\[UndirectedEdge]_],(c1s1b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Sigma]2, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Kappa] analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}]
),{lj,1,dl}]]


(* ::Subsection::Closed:: *)
(*\[CapitalGamma]^(2)*)


(* ::Subsubsection::Closed:: *)
(*Level given*)


GtP2ptslevelgivennocompltad[n3_,n4_,edges_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=edges;
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=Flatten@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


GtP2ptslevelgiven[n3_,n4_,edges_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=edges;
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


(* ::Subsubsection::Closed:: *)
(*Level 2, 2W and 2a*)


GtP2ptslevel2[n3_,n4_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,qaux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=twopt1PIbstc2[n3,n4];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


GtP2ptslevel2W[n3_,n4_,w_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=twopt1PIbstc2W[n3,n4,w];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


GtP2ptslevel2Wass[n3_,n4_,w_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=twopt1PIbstc2W[n3,n4,w];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas,paux[[3;;-1]]])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


GtP2ptslevel2a[n3_,n4_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=twopt1PIbstc2a[n3,n4];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas,paux[[3;;-1]]])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


(* ::Subsubsection::Closed:: *)
(*Level 3a, 3 and 4a*)


GtP2ptslevel3a[n3_,n4_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=twopt1PIbstc3a[n3,n4];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


GtP2ptslevel3[n3_,n4_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,\[Kappa]vert,vertex\[Kappa]pos,analytic\[Kappa],squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=twopt1PIbstc3[n3,n4];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|Subscript[\[Kappa], _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

\[Kappa]vert=Cases[pos,Subscript[\[Kappa], _]];
vertex\[Kappa]pos=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),\[Kappa]vert[[i]]]&)],{i,Length@\[Kappa]vert}];
analytic\[Kappa]=(\[Kappa]fnc@@@Table[Part[paux3,vertex\[Kappa]pos[[j]]],{j,Length@vertex\[Kappa]pos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Kappa] analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


GtP2ptslevel4a[n3_,n4_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,\[Kappa]vert,vertex\[Kappa]pos,analytic\[Kappa],squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=twopt1PIbstc4a[n3,n4];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|Subscript[\[Kappa], _]|Subscript[\[Theta], _]|Subscript[\[Theta]2, _]|Subscript[\[Sigma], _]|Subscript[\[Sigma]2, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Theta], _]\[UndirectedEdge]_],(t2b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Theta], _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Theta]2, _]\[UndirectedEdge]_],(c3b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Theta]2, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Sigma], _]\[UndirectedEdge]_],(c1s/@(Part[paux3,Flatten@Position[daux,Subscript[\[Sigma], _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Sigma]2, _]\[UndirectedEdge]_],(c1s1b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Sigma]2, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother* analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


(* ::Subsubsection::Closed:: *)
(*Level 4 and 4W*)


GtP2ptslevel4[n3_,n4_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,\[Kappa]vert,vertex\[Kappa]pos,analytic\[Kappa],squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=twopt1PIbstc4[n3,n4];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|Subscript[\[Kappa], _]|Subscript[\[Theta], _]|Subscript[\[Theta]2, _]|Subscript[\[Sigma], _]|Subscript[\[Sigma]2, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

\[Kappa]vert=Cases[pos,Subscript[\[Kappa], _]];
vertex\[Kappa]pos=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),\[Kappa]vert[[i]]]&)],{i,Length@\[Kappa]vert}];
analytic\[Kappa]=(\[Kappa]fnc@@@Table[Part[paux3,vertex\[Kappa]pos[[j]]],{j,Length@vertex\[Kappa]pos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Theta], _]\[UndirectedEdge]_],(t2b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Theta], _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Theta]2, _]\[UndirectedEdge]_],(c3b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Theta]2, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Sigma], _]\[UndirectedEdge]_],(c1s/@(Part[paux3,Flatten@Position[daux,Subscript[\[Sigma], _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Sigma]2, _]\[UndirectedEdge]_],(c1s1b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Sigma]2, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Kappa] analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


GtP2ptslevel4W[n3_,n4_,w_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,\[Kappa]vert,vertex\[Kappa]pos,analytic\[Kappa],squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=twopt1PIbstc4W[n3,n4,w];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|Subscript[\[Kappa], _]|Subscript[\[Theta], _]|Subscript[\[Theta]2, _]|Subscript[\[Sigma], _]|Subscript[\[Sigma]2, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

\[Kappa]vert=Cases[pos,Subscript[\[Kappa], _]];
vertex\[Kappa]pos=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),\[Kappa]vert[[i]]]&)],{i,Length@\[Kappa]vert}];
analytic\[Kappa]=(\[Kappa]fnc@@@Table[Part[paux3,vertex\[Kappa]pos[[j]]],{j,Length@vertex\[Kappa]pos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Theta], _]\[UndirectedEdge]_],(t2b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Theta], _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Theta]2, _]\[UndirectedEdge]_],(c3b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Theta]2, _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Sigma], _]\[UndirectedEdge]_],(c1s/@(Part[paux3,Flatten@Position[daux,Subscript[\[Sigma], _]\[UndirectedEdge]_]]))/.List-> Times,1]*
If[MemberQ[daux,Subscript[\[Sigma]2, _]\[UndirectedEdge]_],(c1s1b/@(Part[paux3,Flatten@Position[daux,Subscript[\[Sigma]2, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Kappa] analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


(* ::Subsection::Closed:: *)
(*\[CapitalGamma]^(2)'*)


(* ::Subsubsection::Closed:: *)
(*General functions*)


(* ::Text:: *)
(*Chose the right path*)


weighting[listedg_]:=listedg/.a_/;MatchQ[a,o1\[UndirectedEdge]_|o2\[UndirectedEdge]_]:> 1/.a_/;MatchQ[a,Subscript[b, _]\[UndirectedEdge]_|_\[UndirectedEdge]Subscript[b, _]]:> 2/5/.a_/;MatchQ[a,Subscript[s, _]\[UndirectedEdge]_|_\[UndirectedEdge]Subscript[s, _]]:>3/.a_/;MatchQ[a,Subscript[t, _]\[UndirectedEdge]_|_\[UndirectedEdge]Subscript[t, _]]:> 5/3/.a_/;MatchQ[a,Subscript[c, _]\[UndirectedEdge]_|_\[UndirectedEdge]Subscript[c, _]]:> 40/.a_/;MatchQ[a,x_\[UndirectedEdge]y_]:> 1


buildingedges[listvert_]:=Table[UndirectedEdge[listvert[[i]],listvert[[i+1]]],{i,Length@listvert-1}]


reshuffle[listedg_,listpath_]:=Module[{llp= Length@listpath,tmpedg=listedg},
For[i=1,i<=llp,i++,
	If[MemberQ[tmpedg,listpath[[llp+1-i]]],
	tmpedg=Permute[tmpedg,Cycles[{Table[n,{n,3,Position[tmpedg,listpath[[llp+1-i]]][[1,1]]}]}]],
	tmpedg=Permute[tmpedg,Cycles[{Table[n,{n,3,Position[tmpedg,listpath[[llp+1-i]]/.a_\[UndirectedEdge]b_:>b\[UndirectedEdge]a][[1,1]]}]}]]];(*Print[tmpedg]*)
	];
tmpedg]


(*choicepext[listedges_]:=Module[{lstdg=listedges,weightlist,shortvertex,shortedges,llp,tmpedg},
weighting[listedg_]:=listedg/.a_/;MatchQ[a,o1\[UndirectedEdge]_|o2\[UndirectedEdge]_]:> 1/.a_/;MatchQ[a,Subscript[b, _]\[UndirectedEdge]_|_\[UndirectedEdge]Subscript[b, _]]:> 2/5/.a_/;MatchQ[a,Subscript[s, _]\[UndirectedEdge]_|_\[UndirectedEdge]Subscript[s, _]]:>3/.a_/;MatchQ[a,Subscript[t, _]\[UndirectedEdge]_|_\[UndirectedEdge]Subscript[t, _]]:> 5/3/.a_/;MatchQ[a,Subscript[c, _]\[UndirectedEdge]_|_\[UndirectedEdge]Subscript[c, _]]:> 40/.a_/;MatchQ[a,x_\[UndirectedEdge]y_]:> 1;
weightlist=weighting[lstdg];
shortvertex=FindShortestPath[Graph[lstdg,EdgeWeight->weightlist],x,0];
buildingedges[listvert_]:=Table[UndirectedEdge[listvert[[i]],listvert[[i+1]]],{i,Length@listvert-1}];
shortedges=buildingedges[shortvertex];
(*Print[weightlist];
Print[shortedges];*)
tmpedg=lstdg;
llp= Length@shortedges;
For[i=1,i<=llp,i++,
	If[MemberQ[tmpedg,shortedges[[llp+1-i]]],
	tmpedg=Permute[tmpedg,Cycles[{Table[n,{n,3,Position[tmpedg,shortedges[[llp+1-i]]][[1,1]]}]}]],
	tmpedg=Permute[tmpedg,Cycles[{Table[n,{n,3,Position[tmpedg,shortedges[[llp+1-i]]/.a_\[UndirectedEdge]b_:>b\[UndirectedEdge]a][[1,1]]}]}]]];(*Print[tmpedg]*)
	];
tmpedg]*)


(* ::Input:: *)
(*(*twopt1PIbstc[0,4][[2,2]]*)*)


(* ::Input:: *)
(*(*visualizer[2,0,4][[2]]*)*)


(* ::Input:: *)
(*(*choicepext[Flatten[twopt1PIbstc[0,4][[2,2]]]]*)*)


choicepextnop[listedges_]:=Module[{lstdg=listedges,weightlist,shortvertex,shortedges,llp,tmpedg},
weighting[listedg_]:=listedg/.a_/;MatchQ[a,o1\[UndirectedEdge]_|o2\[UndirectedEdge]_]:> 1/.a_/;MatchQ[a,Subscript[b, _]\[UndirectedEdge]_|_\[UndirectedEdge]Subscript[b, _]]:> 2/5/.a_/;MatchQ[a,Subscript[s, _]\[UndirectedEdge]_|_\[UndirectedEdge]Subscript[s, _]]:>3/.a_/;MatchQ[a,Subscript[t, _]\[UndirectedEdge]_|_\[UndirectedEdge]Subscript[t, _]]:> 5/3/.a_/;MatchQ[a,Subscript[c, _]\[UndirectedEdge]_|_\[UndirectedEdge]Subscript[c, _]]:> 40/.a_/;MatchQ[a,x_\[UndirectedEdge]y_]:> 1;
weightlist=weighting[lstdg];
shortvertex=FindShortestPath[Graph[lstdg,EdgeWeight->weightlist],x,0];
buildingedges[listvert_]:=Table[UndirectedEdge[listvert[[i]],listvert[[i+1]]],{i,Length@listvert-1}];
shortedges=buildingedges[shortvertex];
(*Print[weightlist];
Print[shortedges];*)
tmpedg=lstdg;
llp= Length@shortedges;
For[i=1,i<=llp,i++,
	If[MemberQ[tmpedg,shortedges[[llp+1-i]]],
	tmpedg=Permute[tmpedg,Cycles[{Table[n,{n,3,Position[tmpedg,shortedges[[llp+1-i]]][[1,1]]}]}]],
	tmpedg=Permute[tmpedg,Cycles[{Table[n,{n,3,Position[tmpedg,shortedges[[llp+1-i]]/.a_\[UndirectedEdge]b_:>b\[UndirectedEdge]a][[1,1]]}]}]]];(*Print[tmpedg]*)];
tmpedg]


choicepextnop[listedges_]:=Module[{lstdg=listedges,weightlist,shortvertex,shortedges,llp,tmpedg},
weighting[listedg_]:=listedg/.a_/;MatchQ[a,o1\[UndirectedEdge]_|o2\[UndirectedEdge]_]:> 1/.a_/;MatchQ[a,Subscript[b, _]\[UndirectedEdge]_|_\[UndirectedEdge]Subscript[b, _]]:> 2/5/.a_/;MatchQ[a,Subscript[s, _]\[UndirectedEdge]_|_\[UndirectedEdge]Subscript[s, _]]:>3/.a_/;MatchQ[a,Subscript[t, _]\[UndirectedEdge]_|_\[UndirectedEdge]Subscript[t, _]]:> 2/5/.a_/;MatchQ[a,Subscript[c, _]\[UndirectedEdge]_|_\[UndirectedEdge]Subscript[c, _]]:> 40/.a_/;MatchQ[a,x_\[UndirectedEdge]y_]:> 1;
weightlist=weighting[lstdg];
shortvertex=FindShortestPath[Graph[lstdg,EdgeWeight->weightlist],o2,o1];
buildingedges[listvert_]:=Table[UndirectedEdge[listvert[[i]],listvert[[i+1]]],{i,Length@listvert-1}];
shortedges=buildingedges[shortvertex];
(*Print[weightlist];
Print[shortedges];*)
tmpedg=lstdg;
llp= Length@shortedges;
For[i=1,i<=llp,i++,If[MemberQ[tmpedg,shortedges[[llp+1-i]]],tmpedg=Permute[tmpedg,Cycles[{Table[n,{n,3,Position[tmpedg,shortedges[[llp+1-i]]][[1,1]]}]}]],tmpedg=Permute[tmpedg,Cycles[{Table[n,{n,3,Position[tmpedg,shortedges[[llp+1-i]]/.a_\[UndirectedEdge]b_:>b\[UndirectedEdge]a][[1,1]]}]}]]];(*Print[tmpedg]*)];
tmpedg]


(* ::Text:: *)
(*Derive*)


loopinintegrand[x_] :=
    Module[{tmp, test, i, l},
        tmp = x;
        test = 1;
        l = 0;
        i = 1;
        While[
            i <= Length[tmp] && test == 1
            ,
            If[MemberQ[tmp, qaux[[i]], All],
                l++;
                i++
                ,
                test = 0
            ]
        ];
        l
    ]


loopinfunc[funcint_]:=Table[loopinintegrand[funcint[[i]]],{i,Length@funcint}]


loopinintegrandok[x_] :=
    Module[{tmp, test, i, l},
        tmp = x;
        test = 1;
        l = 0;
        i = 1;
        While[
            i <= Length @ qaux && test == 1
            ,
            If[MemberQ[tmp, qaux[[i]], All],
                l++;
                i++
                ,
                test = 0
            ]
        ];
        l
    ]


loopinfuncok[funcint_]:=Table[loopinintegrandok[funcint[[i]]],{i,Length@funcint}]


(* ::Input:: *)
(*(*functionder[intg_]:=intg/.tder1[a_,b_,c_]\[RuleDelayed]td1@@(vec3dsq/@{a,b,c})/.tder2[a_,b_,c_]\[RuleDelayed]td2@@(vec3dsq/@{a,b,c})/.tder11[a_,b_,c_]\[RuleDelayed]td11@@(vec3dsq/@{a,b,c})*)*)


(* ::Input:: *)
(*(*functionder2[intg_]:=intg/.tder1[a_,b_,c_]\[RuleDelayed]1/4trider1@@(Sqrt/@vec3dsq/@{a,b,c})/.tder2[a_,b_,c_]\[RuleDelayed]1/4trider2@@(Sqrt/@vec3dsq/@{a,b,c})/.tder11[a_,b_,c_]\[RuleDelayed]1/4trider11@@(Sqrt/@vec3dsq/@{a,b,c})*)*)


functionderMath[intg_]:=intg/.tder1[a_,b_,c_]:>1/2 trider1@@(Sqrt/@vec3dsq/@{a,b,c})/.tder2[a_,b_,c_]:>1/4 trider2@@(Sqrt/@vec3dsq/@{a,b,c})/.tder11[a_,b_,c_]:>1/4 trider11@@(Sqrt/@vec3dsq/@{a,b,c})


(* ::Input:: *)
(*(*writeGbubblder[int_]:=int/.suns2[p_,q_]\[RuleDelayed]((bubbl[p](Gt[p+q]-Gt[p]))(*-(((*If[p===q2,*)2\[Pi](*,((2\[Pi])^2)]*)2cren)/(vec3dsq[p](1+vec3dnorm[p])^2Sin[Subscript[p, 4]]))*))/.Gt[p_]\[RuleDelayed] 1/(vec3dsq[p]+1)/.bubbl[x_]\[RuleDelayed] ArcTan[ Sqrt[vec3dsq[x]]/2]/(4\[Pi] Sqrt[vec3dsq[x]])*)*)


(* ::Subsubsection:: *)
(*Level given*)


GtP2ptslevelgivendernocompltad[n3_,n4_,edges_]:=Module[
	{pos,predaux,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,
	propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,
	analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
	
	diaglist=edges;
	dl=Length@diaglist;
	(*Print[diaglist];*)
	Table[
		(
			deltas={};
			predaux=diaglist[[lj]][[2]];
			daux=If[MemberQ[predaux,o1\[UndirectedEdge]0],Flatten@predaux,choicepextnop@Flatten@predaux];
			pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
			If[
				(*Length@pos==1*)False,
				diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
				(
					paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
					paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
					paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
					For[idx=1,idx<=Length@pos,idx++,
						xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
						xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
						xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
						AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];
					];
					deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
					paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
					paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

					propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

					trinvert=Cases[pos,Subscript[t, _]];
					vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
					analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

					\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
					vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
					analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

					\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
					vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
					analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

					squarevert=Cases[pos,Subscript[c, _]];
					vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
					(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
					analyticsquar=(squarfnc@@@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]))/.List-> Times;

					analyticother=If[
						MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],
						(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,
						1
						]*If[
							MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],
							(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,
							1
						];
					analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];
					diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad
				)
			]
		),
		{lj,1,dl}]]


GtP2ptslevelgivender[n3_,n4_,edges_]:=Module[
	{pos,predaux,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,
	trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,
	vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
	
	diaglist=edges;
	dl=Length@diaglist;
	Table[
		(
			deltas={};
			predaux=diaglist[[lj]][[2]];
			daux=If[MemberQ[predaux,o1\[UndirectedEdge]0],replaceTad\[Tau]\[Beta]@predaux,choicepextnop@replaceTad\[Tau]\[Beta]@predaux];
			pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
			paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
			paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
			paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
			For[
				idx=1, idx<=Length@pos, idx++,
				xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
				xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
				xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
				AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];
			];
			deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
			paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
			paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];
			propagators = Table[
				If[
					MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],
					1,
					Gt[paux3[[k]]]
				],
				{k,Length@daux}
			]/.List-> Times;
			trinvert=Cases[pos,Subscript[t, _]];
			vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
			analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;
			\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
			vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
			analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;
			\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
			vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
			analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;
			squarevert=Cases[pos,Subscript[c, _]];
			vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
			(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
			analyticsquar=(squarfnc@@@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]))/.List-> Times;
			
			analyticother=If[
				MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],
				(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,
				1
			]*If[
				MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],
				(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,
			1];
			analyticcmplxtad = suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]];
			diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b*analytic\[Psi]a*analytictri*analyticsquar*propagators*analyticcmplxtad
		),
		{lj,1,dl}
	]
]


(* ::Subsubsection::Closed:: *)
(*Level 2, 2W*)


GtP2ptslevel2v1der[n3_,n4_]:=Module[{pos,predaux,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=twopt1PIbstc2[n3,n4];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
predaux=diaglist[[lj]][[2]];
daux=If[MemberQ[predaux,o1\[UndirectedEdge]0],replaceTad\[Tau]\[Beta]@predaux,choicepextnop@replaceTad\[Tau]\[Beta]@predaux];
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}])/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad]
),{lj,1,dl}]]


GtP2ptslevel2v1derW[n3_,n4_,w_]:=Module[{pos,predaux,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=twopt1PIbstc2W[n3,n4,w];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
predaux=diaglist[[lj]][[2]];
daux=If[MemberQ[predaux,o1\[UndirectedEdge]0],replaceTad\[Tau]\[Beta]@predaux,choicepextnop@replaceTad\[Tau]\[Beta]@predaux];
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad]
),{lj,1,dl}]]


GtP2ptslevel2v1derWass[n3_,n4_,w_]:=Module[{pos,predaux,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=twopt1PIbstc2W[n3,n4,w];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
predaux=diaglist[[lj]][[2]];
(* WHY THE IF?? BOTH BRANCHES ARE THE SAME
daux=If[MemberQ[predaux,o1\[UndirectedEdge]0],replaceTad\[Tau]\[Beta]@predaux,replaceTad\[Tau]\[Beta]@predaux];*)
daux=replaceTad\[Tau]\[Beta]@predaux;
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas,paux[[3;;-1]]])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad]
),{lj,1,dl}]]


(* ::Subsection::Closed:: *)
(*\[CapitalGamma]^(0)*)


(* ::Subsubsection::Closed:: *)
(*Level given*)


GtP0ptslevelgivennocompltad[n3_,n4_,edges_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=edges;
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=Flatten@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,o1|o2|o3|o4];
If[Length@pos==1(*||Length@pos\[Equal]2*),diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|Subscript[\[Tau], \[Phi]]|Subscript[\[Beta], \[Phi]]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


GtP0ptslevelgiven[n3_,n4_,edges_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=edges;
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,o1|o2|o3|o4];
If[Length@pos==1(*||Length@pos\[Equal]2*),diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|Subscript[\[Tau], \[Phi]]|Subscript[\[Beta], \[Phi]]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


(* ::Subsubsection::Closed:: *)
(*Level 2, 2W and 2a*)


GtP0ptslevel2[n3_,n4_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=zeropt1PIbstc2[n3,n4];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1(*||Length@pos\[Equal]2*),diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|Subscript[\[Tau], \[Phi]]|Subscript[\[Beta], \[Phi]]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


GtP0ptslevel2W[n3_,n4_,w_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=zeropt1PIbstc2W[n3,n4,w];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|Subscript[\[Tau], \[Phi]]|Subscript[\[Beta], \[Phi]]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


GtP0ptslevel2Wass[n3_,n4_,w_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=zeropt1PIbstc2W[n3,n4,w];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas,paux])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|Subscript[\[Tau], \[Phi]]|Subscript[\[Beta], \[Phi]]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


GtP0ptslevel2a[n3_,n4_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,propagators,trinvert,vertexpos,analytictri,\[Psi]avert,vertex\[Psi]apos,analytic\[Psi]a,\[Psi]bvert,vertex\[Psi]bpos,analytic\[Psi]b,squarevert,vertexposquad,analyticsquar,analyticother,analyticcmplxtad(*,analyticpiec*)},
diaglist=zeropt1PIbstc2a[n3,n4];
dl=Length@diaglist;
(*Print[diaglist];*)
Table[( 
deltas={};
daux=replaceTad\[Tau]\[Beta]@(diaglist[[lj]][[2]]);
pos=DeleteCases[VertexList@daux,0|o1|o2|o3|o4];
If[Length@pos==1,diaglist[[lj]][[1]]*suntad^Count[pos,Subscript[\[Tau],_]]*bubbl2tr^Count[pos,Subscript[\[Beta],_]],
paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux]];
For[idx=1,idx<=Length@pos,idx++,
xlist=Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),pos[[idx]]]&];
xpos=DeleteDuplicates@Flatten@Table[Position[daux,xlist[[k]]],{k,1,Length@xlist}];
xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length@xlist}];
AppendTo[deltas,Sum[xpossign[[k,2]]*paux[[xpossign[[k,1]]]],{k,1,Length@xpossign}]==0];];
deltasol=(Reduce[deltas])/.And->List/.Equal:>Rule;
paux2=DeleteCases[DeleteDuplicates[Flatten[(paux/.deltasol)/.Plus->List]/.Times[n_,qq_]:>qq/;NumericQ[n]],0];
paux3=paux/.deltasol/.Table[If[paux2=={},Null;,paux2[[k]]->qauxify[k]],{k,1,Length@paux2}];

propagators=Table[If[MemberQ[daux[[k]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|Subscript[\[Tau], \[Phi]]|Subscript[\[Beta], \[Phi]]|o1|o2|o3|o4],1,Gt[paux3[[k]]]],{k,Length@daux}]/.List-> Times;

trinvert=Cases[pos,Subscript[t, _]];
vertexpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),trinvert[[i]]]&)],{i,Length@trinvert}];
analytictri=(trinfnc@@@Table[Part[paux3,vertexpos[[j]]],{j,Length@vertexpos}])/.List-> Times;

\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]avert[[i]]]&)],{i,Length@\[Psi]avert}];
analytic\[Psi]a=(\[Psi]afnc@@@Table[Part[paux3,vertex\[Psi]apos[[j]]],{j,Length@vertex\[Psi]apos}])/.List-> Times;

\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten@Position[daux,_?(MemberQ[(#/.a_\[UndirectedEdge]b_:>{a,b}),\[Psi]bvert[[i]]]&)],{i,Length@\[Psi]bvert}];
analytic\[Psi]b=(\[Psi]bfnc@@@Table[Part[paux3,vertex\[Psi]bpos[[j]]],{j,Length@vertex\[Psi]bpos}])/.List-> Times;

squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten@Position[daux,_?(MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),squarevert[[i]]]&)],{i,Length@squarevert}];
(*analyticsquar=(squarfnc@@@Table[Part[paux3,vertexposquad\[LeftDoubleBracket]j\[RightDoubleBracket]],{j,Length@vertexposquad}])/.List\[Rule] Times;*)
analyticsquar=(squarfnc@@@shiftmomok/@(Table[Part[paux3,vertexposquad[[j]]],{j,Length@vertexposquad}]/.{pext->0,-pext->0}))/.List-> Times;

analyticother=If[MemberQ[daux,Subscript[b, _]\[UndirectedEdge]_],(bubbl/@(Part[paux3,Flatten@Position[daux,Subscript[b, _]\[UndirectedEdge]_]]))/.List-> Times,1]*If[MemberQ[daux,Subscript[s, _]\[UndirectedEdge]_],(suns/@(Part[paux3,Flatten@Position[daux,Subscript[s, _]\[UndirectedEdge]_]]))/.List-> Times,1];
analyticcmplxtad=suntad^Count[pos,Subscript[\[Tau],_]] bubbl2tr^Count[pos,Subscript[\[Beta],_]];

diaglist[[lj]][[1]]*analyticother*analytic\[Psi]b analytic\[Psi]a analytictri*analyticsquar*propagators*analyticcmplxtad/.{pext->0,-pext->0}]
),{lj,1,dl}]]


(* ::Subsubsection:: *)
(*Level 4 and 4W (no need)*)


(* ::Section::Closed:: *)
(*Re-assignation of momenta*)


(* ::Subsection:: *)
(*Generalization functions written for \[CapitalGamma]^(4)*)


(* ::Subsection::Closed:: *)
(*Function momenta*)


(* ::Text:: *)
(*Level 3a (or lower)*)


GtP4ptsdConfig[n3_,n4_,fourpt1PIod_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,trinvert,\[Psi]avert,vertex\[Psi]apos,\[Psi]bvert,vertex\[Psi]bpos,squarevert,vertexposquad,propagatorspos,deauxOK,vertexTripos,bubvert,vertexposbub,sunvert,vertexpossuns},
diaglist=fourpt1PIod;dl=Length[diaglist];Table[deltas={};daux=replaceTad\[Tau]\[Beta][diaglist[[lj]][[2]]];pos=DeleteCases[VertexList[daux],0|o1|o2|o3|o4];
If[Length[pos]==1,diaglist[[lj]][[1]] suntad^Count[pos,Subscript[\[Tau], _]] bubbl2tr^Count[pos,Subscript[\[Beta], _]],paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length[daux]-4]];paux=Insert[paux,pext1,Position[daux,Select[daux,MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},o1]&][[1]]][[1,1]]];paux=Insert[paux,pext2,Position[daux,Select[daux,MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},o2]&][[1]]][[1,1]]];paux=Insert[paux,pext3,Position[daux,Select[daux,MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},o3]&][[1]]][[1,1]]];paux=Insert[paux,-pext1-pext2-pext3,Position[daux,Select[daux,MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},o4]&][[1]]][[1,1]]];
For[idx=1,idx<=Length[pos],idx++,xlist=Select[daux,MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},pos[[idx]]]&];xpos=DeleteDuplicates[Flatten[Table[Position[daux,xlist[[k]]],{k,1,Length[xlist]}]]];xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length[xlist]}];AppendTo[deltas,\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(Length[xpossign]\)]\(xpossign[\([k, 2]\)]\ paux[\([xpossign[\([k, 1]\)]]\)]\)\)==0];];deltasol=Reduce[deltas]/. And->List/. Equal:>Rule;(*Print[deltasol];*)
paux2=DeleteCases[DeleteDuplicates[Flatten[paux/. deltasol/. Plus->List]/. Times[n_,qq_]:>qq/;NumericQ[n]],0|pext1|pext2|pext3];
paux3=paux/. deltasol/. Table[If[paux2=={},Null,paux2[[k]]->qauxify[k]],{k,1,Length[paux2]}];(*Print[paux3];*)
deauxOK=daux[[5;;-1]];
propagatorspos=Select[Range[Length@deauxOK],!MemberQ[deauxOK[[#]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4]&];
trinvert=Cases[pos,Subscript[t, _]];
vertexTripos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},trinvert[[i]]]&)]],{i,Length[trinvert]}];
\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Psi]avert[[i]]]&)]],{i,Length[\[Psi]avert]}];
\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Psi]bvert[[i]]]&)]],{i,Length[\[Psi]bvert]}];
squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},squarevert[[i]]]&)]],{i,Length[squarevert]}];
bubvert=Cases[pos,Subscript[b, _]];
vertexposbub=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},bubvert[[i]]]&)]],{i,Length[bubvert]}];
sunvert=Cases[pos,Subscript[s, _]];
vertexpossuns=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},sunvert[[i]]]&)]],{i,Length[sunvert]}];

{{propagatorspos,vertexposbub,vertexpossuns,vertex\[Psi]apos,vertex\[Psi]bpos,vertexTripos,vertexposquad},
paux3[[5;;-1]]/. {pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}}],{lj,1,dl}]]


(* ::Text:: *)
(*Level 4*)


GtP4ptsdConfiglevel4[n3_,n4_,fourpt1PIod_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,trinvert,\[Psi]avert,vertex\[Psi]apos,\[Psi]bvert,vertex\[Psi]bpos,squarevert,vertexposquad,propagatorspos,deauxOK,vertexTripos,bubvert,vertexposbub,sunvert,vertexpossuns,vertex\[Kappa]pos,vertexpostheta,vertexpostheta2,vertexpossigma,vertexpossigma2,\[Kappa]vert,thetavert,theta2vert,sigmavert,sigma2vert},
diaglist=fourpt1PIod;
dl=Length[diaglist];
Table[deltas={};daux=replaceTad\[Tau]\[Beta][diaglist[[lj]][[2]]];pos=DeleteCases[VertexList[daux],0|o1|o2|o3|o4];If[Length[pos]==1,diaglist[[lj]][[1]] suntad^Count[pos,Subscript[\[Tau], _]] bubbl2tr^Count[pos,Subscript[\[Beta], _]],paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length[daux]-4]];paux=Insert[paux,pext1,Position[daux,Select[daux,MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},o1]&][[1]]][[1,1]]];paux=Insert[paux,pext2,Position[daux,Select[daux,MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},o2]&][[1]]][[1,1]]];paux=Insert[paux,pext3,Position[daux,Select[daux,MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},o3]&][[1]]][[1,1]]];paux=Insert[paux,-pext1-pext2-pext3,Position[daux,Select[daux,MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},o4]&][[1]]][[1,1]]];
For[idx=1,idx<=Length[pos],idx++,xlist=Select[daux,MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},pos[[idx]]]&];xpos=DeleteDuplicates[Flatten[Table[Position[daux,xlist[[k]]],{k,1,Length[xlist]}]]];xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length[xlist]}];AppendTo[deltas,\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(Length[xpossign]\)]\(xpossign[\([k, 2]\)]\ paux[\([xpossign[\([k, 1]\)]]\)]\)\)==0];];deltasol=Reduce[deltas]/. And->List/. Equal:>Rule;(*Print[deltasol];*)
paux2=DeleteCases[DeleteDuplicates[Flatten[paux/. deltasol/. Plus->List]/. Times[n_,qq_]:>qq/;NumericQ[n]],0|pext1|pext2|pext3];
paux3=paux/. deltasol/. Table[If[paux2=={},Null,paux2[[k]]->qauxify[k]],{k,1,Length[paux2]}];(*Print[paux3];*)
deauxOK=daux[[5;;-1]];
propagatorspos=Select[Range[Length@deauxOK],!MemberQ[deauxOK[[#]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|Subscript[\[Kappa], _]|Subscript[\[Theta], _]|Subscript[\[Theta]2, _]|Subscript[\[Sigma], _]|Subscript[\[Sigma]2, _]|o1|o2|o3|o4]&];
trinvert=Cases[pos,Subscript[t, _]];
vertexTripos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},trinvert[[i]]]&)]],{i,Length[trinvert]}];
\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Psi]avert[[i]]]&)]],{i,Length[\[Psi]avert]}];
\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Psi]bvert[[i]]]&)]],{i,Length[\[Psi]bvert]}];
\[Kappa]vert=Cases[pos,Subscript[\[Kappa], _]];
vertex\[Kappa]pos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Kappa]vert[[i]]]&)]],{i,Length[\[Kappa]vert]}];
squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},squarevert[[i]]]&)]],{i,Length[squarevert]}];
bubvert=Cases[pos,Subscript[b, _]];
vertexposbub=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},bubvert[[i]]]&)]],{i,Length[bubvert]}];
sunvert=Cases[pos,Subscript[s, _]];
vertexpossuns=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},sunvert[[i]]]&)]],{i,Length[sunvert]}];
thetavert=Cases[pos,Subscript[\[Theta], _]];
vertexpostheta=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},thetavert[[i]]]&)]],{i,Length[thetavert]}];
theta2vert=Cases[pos,Subscript[\[Theta]2, _]];
vertexpostheta2=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},theta2vert[[i]]]&)]],{i,Length[theta2vert]}];
sigmavert=Cases[pos,Subscript[\[Sigma], _]];
vertexpossigma=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},sigmavert[[i]]]&)]],{i,Length[sigmavert]}];
sigma2vert=Cases[pos,Subscript[\[Sigma]2, _]];
vertexpossigma2=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},sigma2vert[[i]]]&)]],{i,Length[sigma2vert]}];

{{propagatorspos,vertexposbub,vertexpossuns,vertex\[Psi]apos,vertex\[Psi]bpos,vertexTripos,vertexposquad,vertex\[Kappa]pos,vertexpostheta,vertexpostheta2,vertexpossigma,vertexpossigma2},
paux3[[5;;-1]]/. {pext1->0,-pext1->0,pext2->0,-pext2->0,pext3->0,-pext3->0}}],{lj,1,dl}]]


(* ::Subsection::Closed:: *)
(*Basic functions*)


tovector2[x_]:={If[MemberQ[x,q1,All],Coefficient[x,q1],0],If[MemberQ[x,q2,All],Coefficient[x,q2],0]}


tovector3[x_]:={If[MemberQ[x,q1,All],Coefficient[x,q1],0],If[MemberQ[x,q2,All],Coefficient[x,q2],0],If[MemberQ[x,q3,All],Coefficient[x,q3],0]}


fun={gw,bbw,ssw,vv\[Psi]aw,vv\[Psi]bw,ttw,ssqw};


weightprov[conf_,grp_]:=Module[{tabtmp,j},tabtmp={};
For[j=1,j<= Length@conf,j++,AppendTo[tabtmp,Table[fun[[j]]@grp[[conf[[j,i]]]],{i,Length@conf[[j]]}]]];
tabtmp
]


bbw[{x_,y_}]:=bw[x]
ssw[{x_,y_}]:=sw[x]
vv\[Psi]aw[{x_,y_,z_}]:=v\[Psi]aw[SortBy[{x,y,z},Total@Abs@#&][[-1]]]
vv\[Psi]bw[{x_,y_,z_}]:=v\[Psi]bw[SortBy[{x,y,z},Total@Abs@#&][[-1]]]
vv\[Kappa]w[{x_,y_,z_,w_}]:=v\[Kappa]w[SortBy[{x,y,z,w},Total@Abs@#&][[-1]]]
(*ttw[{x_,y_,z_}]:=tw[SortBy[{x,y,z},Total@Abs@#&]\[LeftDoubleBracket]-1\[RightDoubleBracket]]*)
ttw[{x_,y_,z_}]:=tw[Sort@{x,y,z}]
ssqw[{x_,y_,w_,z_}]:=sqw[Sort@{x,y,w,z}]
vv\[Theta]w[{x_,y_}]:=v\[Theta]w[x]
vv\[Theta]2w[{x_,y_}]:=v\[Theta]2w[x]
vv\[Sigma]w[{x_,y_}]:=v\[Sigma]w[x]
vv\[Sigma]2w[{x_,y_}]:=v\[Sigma]2w[x]





funlev4={gw,bbw,ssw,vv\[Psi]aw,vv\[Psi]bw,ttw,ssqw,vv\[Kappa]w,vv\[Theta]w,vv\[Theta]2w,vv\[Sigma]w,vv\[Sigma]2w};


weightprovlev4[conf_,grp_]:=Module[{tabtmp,j},tabtmp={};
For[j=1,j<= Length@conf,j++,AppendTo[tabtmp,Table[funlev4[[j]]@grp[[conf[[j,i]]]],{i,Length@conf[[j]]}]]];
tabtmp
]


(* ::Subsection::Closed:: *)
(*Matrices (linear transformations on the internal momenta)*)


(* ::Subsubsection::Closed:: *)
(*Matrices 2d*)


matricesTransOK2d=Select[Flatten[Table[{Tuples[Range[-1,1],2][[i]],Tuples[Range[-1,1],2][[j]]},{i,1,9},{j,1,9}],1],Abs[Det@#]==1&];


Length@matricesTransOK2d


testVMtf2d[v_,m_]:=Module[{tmpp},
tmpp=v . m;
Abs[tmpp[[1]]]<2&&Abs[tmpp[[2]]]<2
]


testNo22d[vectlist_,x_]:=Table[testVMtf2d[vectlist[[i]],x],{i,Length@vectlist}]/.List-> And


matricesTransOKconf2d[vectlist_]:=Select[matricesTransOK2d,testNo22d[vectlist,#]&]


(* ::Subsubsection::Closed:: *)
(*Matrices 3d*)


matricesTransOK=Select[Flatten[Table[{Tuples[Range[-1,1],3][[i]],Tuples[Range[-1,1],3][[j]],Tuples[Range[-1,1],3][[y]]},{i,1,27},{j,1,27},{y,1,27}],2],Abs[Det@#]==1&];


Length@matricesTransOK


testVMtf[v_,m_]:=Module[{tmpp},
tmpp=v . m;
Abs[tmpp[[1]]]<2&&Abs[tmpp[[2]]]<2&&Abs[tmpp[[3]]]<2
]


testNo2[vectlist_,x_]:=Table[testVMtf[vectlist[[i]],x],{i,Length@vectlist}]/.List-> And


matricesTransOKconf[vectlist_]:=Select[matricesTransOK,testNo2[vectlist,#]&]


(* ::Subsection::Closed:: *)
(*Weights*)


(* ::Subsubsection::Closed:: *)
(*General for visualization*)


reverse2[x_]:=Reverse/@x


tableform2[x_]:=TableForm/@x


matrixform2[x_]:=MatrixForm/@x


(* ::Text:: *)
(**)


showbetter32d[x_]:=Table[x[[i]]{q1,q2},{i,3}]


showbetter42d[x_]:=Table[x[[i]]{q1,q2},{i,4}]


showbetter3[x_]:=Table[x[[i]]{q1,q2,q3},{i,3}]


showbetter4[x_]:=Table[x[[i]]{q1,q2,q3},{i,4}]


(* ::Subsection::Closed:: *)
(*w 2d*)


(* ::Subsubsection::Closed:: *)
(*Simpler pieces ({gw, bw, sw, v\[Psi]aw, v\[Psi]bw}) 2d*)


basicsequnce2d={{0,0},{1,0},{0,1},{1,1}};


{{gww[{0,0}],gww[{1,0}],gww[{0,1}],gww[{1,1}]},{bww[{0,0}],bww[{1,0}],bww[{0,1}],bww[{1,1}]},{sww[{0,0}],sww[{1,0}],sww[{0,1}],sww[{1,1}]},{v\[Psi]aww[{0,0}],v\[Psi]aww[{1,0}],v\[Psi]aww[{0,1}],v\[Psi]aww[{1,1}]},{v\[Psi]bww[{0,0}],v\[Psi]bww[{1,0}],v\[Psi]bww[{0,1}],v\[Psi]bww[{1,1}]}}={{0,1,13/10,12/5},{0,1,6/5,8/5},{0,1,6/5,9/5},{0,1,6/5,8/5},{0,1,6/5,10}};


(* ::Subsubsection::Closed:: *)
(*Simpler pieces ({gw, bw, sw, v\[Psi]aw, v\[Psi]bw, v\[Kappa]ww, v\[Theta]ww, v\[Theta]2ww, v\[Sigma]ww, v\[Sigma]2ww}) 2d*)


basicsequnce2d={{0,0},{1,0},{0,1},{1,1}};


{{v\[Kappa]ww[{0,0}],v\[Kappa]ww[{1,0}],v\[Kappa]ww[{0,1}],v\[Kappa]ww[{1,1}]},{v\[Theta]ww[{0,0}],v\[Theta]ww[{1,0}],v\[Theta]ww[{0,1}],v\[Theta]ww[{1,1}]},{v\[Theta]2ww[{0,0}],v\[Theta]2ww[{1,0}],v\[Theta]2ww[{0,1}],v\[Theta]2ww[{1,1}]},{v\[Sigma]ww[{0,0}],v\[Sigma]ww[{1,0}],v\[Sigma]ww[{0,1}],v\[Sigma]ww[{1,1}]},{v\[Sigma]2ww[{0,0}],v\[Sigma]2ww[{1,0}],v\[Sigma]2ww[{0,1}],v\[Sigma]2ww[{1,1}]}}={{0,1,6/5,10},{0,1,6/5,10},{0,1,6/5,10},{0,1,6/5,10},{0,1,6/5,10}};


(* ::Subsubsection::Closed:: *)
(*Triangles 2d*)


(* ::Input:: *)
(*Sort@DeleteDuplicates[Sort/@Abs/@Select[Flatten[Table[{Tuples[Range[-1,1],2][[i]],Tuples[Range[-1,1],2][[j]],Tuples[Range[-1,1],2][[i]]+Tuples[Range[-1,1],2][[j]]},{i,1,9},{j,1,9}],1],!MemberQ[Abs/@Flatten@#,2]&]]*)
(*Length@%*)


(* ::Input:: *)
(*Sort/@reverse2/@%%*)


(* ::Text:: *)
(*New order*)


orderokTrPar2d={{{0,0},{0,0},{0,0}},{{0,0},{1,0},{1,0}},{{0,0},{0,1},{0,1}},{{0,0},{1,1},{1,1}},{{0,1},{1,0},{1,1}}};


(* ::Input:: *)
(*(*Rationalize[{0,1,1.3,1.7,2.4,2.8,3.4,5}+0.1]*)*)


weightTr2d=Rationalize@{0,11/10,7/5,5/2,1.8};(*{tww[{{0,0},{0,0},{0,0}}],tww[{{0,0},{1,0},{1,0}}],tww[{{0,0},{0,1},{0,1}}],tww[{{0,0},{1,1},{1,1}}],tww[{{0,1},{1,0},{1,1}}]}=weightTr2d*)


{tww[{{0,0},{0,0},{0,0}}],tww[{{0,0},{1,0},{1,0}}],tww[{{0,0},{0,1},{0,1}}],tww[{{0,0},{1,1},{1,1}}],tww[{{0,1},{1,0},{1,1}}]}=weightTr2d;


(* ::Input:: *)
(*matrixform2/@Table[{(showbetter32d/@orderokTrPar2d)[[i]],{N@weightTr2d[[i]]}},{i,5}]*)


(* ::Input:: *)
(*(*tww/@orderokTrPar2d;*)*)


(* ::Input:: *)
(*(*tww/@orderokTrPar*)*)


(* ::Subsubsection::Closed:: *)
(*Squares 2d*)


sqparvec2d=Sort@DeleteDuplicates[Sort/@Abs/@Select[Flatten[Table[{Tuples[Range[-1,1],2][[i]],Tuples[Range[-1,1],2][[j]],Tuples[Range[-1,1],2][[y]],Tuples[Range[-1,1],2][[i]]+Tuples[Range[-1,1],2][[j]]+Tuples[Range[-1,1],2][[y]]},{i,1,9},{j,1,9},{y,1,9}],2],!MemberQ[Abs/@Flatten@#,2|3]&]];


(* ::Input:: *)
(*Length@sqparvec2d*)


sqpar2d=showbetter42d/@Sort/@reverse2/@sqparvec2d;


(* ::Text:: *)
(*New order*)


orderoksqpar2d={{{0,0},{0,0},{0,0},{0,0}},{{0,0},{0,0},{q1,0},{q1,0}},{{0,0},{0,0},{0,q2},{0,q2}},{{0,0},{0,0},{q1,q2},{q1,q2}},
{{0,0},{0,q2},{q1,0},{q1,q2}},
{{q1,0},{q1,0},{q1,0},{q1,0}},{{0,q2},{0,q2},{0,q2},{0,q2}},{{0,q2},{0,q2},{q1,0},{q1,0}},
{{q1,0},{q1,0},{q1,q2},{q1,q2}},{{0,q2},{0,q2},{q1,q2},{q1,q2}},{{q1,q2},{q1,q2},{q1,q2},{q1,q2}}
};


(* ::Input:: *)
(*tableform2/@Table[{orderoksqpar2d[[i]],{}},{i,11}];*)


(* ::Text:: *)
(*Assigning*)


(* ::Input:: *)
(*{{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "0"},*)
(*{"0", "0"},*)
(*{"0", "0"},*)
(*{"0", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "0", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "0"},*)
(*{"0", "0"},*)
(*{"q1", "0"},*)
(*{"q1", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "1.2", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "0"},*)
(*{"0", "0"},*)
(*{"0", "q2"},*)
(*{"0", "q2"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "1.5", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "0"},*)
(*{"0", "0"},*)
(*{"q1", "q2"},*)
(*{"q1", "q2"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.6", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "0"},*)
(*{"0", "q2"},*)
(*{"q1", "0"},*)
(*{"q1", "q2"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "1.9", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"q1", "0"},*)
(*{"q1", "0"},*)
(*{"q1", "0"},*)
(*{"q1", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "1.2", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "q2"},*)
(*{"0", "q2"},*)
(*{"0", "q2"},*)
(*{"0", "q2"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "1.6", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "q2"},*)
(*{"0", "q2"},*)
(*{"q1", "0"},*)
(*{"q1", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"q1", "0"},*)
(*{"q1", "0"},*)
(*{"q1", "q2"},*)
(*{"q1", "q2"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.4", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "q2"},*)
(*{"0", "q2"},*)
(*{"q1", "q2"},*)
(*{"q1", "q2"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.5", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"q1", "q2"},*)
(*{"q1", "q2"},*)
(*{"q1", "q2"},*)
(*{"q1", "q2"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.2", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)}};*)


(* ::Input:: *)
(*Rationalize@Flatten[%[[All,2]]]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Flatten@Table[Position[sqpar2d,orderoksqpar2d[[i]]][[1]],{i,11}]*)
(*sqparvecordOK2d=sqparvec2d[[%]]*)


weightSq2d={0,6/5,3/2,13/5,19/10,6/5,8/5,12/5,5/2,11/5};


(* ::Input:: *)
(*sqww/@sqparvecordOK2d;*)


{sqww[{{0,0},{0,0},{0,0},{0,0}}],sqww[{{0,0},{0,0},{0,1},{0,1}}],sqww[{{0,0},{0,0},{1,0},{1,0}}],sqww[{{0,0},{0,0},{1,1},{1,1}}],sqww[{{0,0},{0,1},{1,0},{1,1}}],sqww[{{0,1},{0,1},{0,1},{0,1}}],sqww[{{1,0},{1,0},{1,0},{1,0}}],sqww[{{0,1},{0,1},{1,0},{1,0}}],sqww[{{0,1},{0,1},{1,1},{1,1}}],sqww[{{1,0},{1,0},{1,1},{1,1}}],sqww[{{1,1},{1,1},{1,1},{1,1}}]}={0,6/5,3/2,13/5,19/10,6/5,8/5,2,12/5,5/2,11/5};


(* ::Input:: *)
(*(*{sqww[{{0,0},{0,0},{0,0},{0,0}}],sqww[{{0,0},{0,0},{0,1},{0,1}}],sqww[{{0,0},{0,0},{1,0},{1,0}}],sqww[{{0,0},{0,0},{1,1},{1,1}}],sqww[{{0,0},{0,1},{1,0},{1,1}}],sqww[{{0,1},{0,1},{0,1},{0,1}}],sqww[{{1,0},{1,0},{1,0},{1,0}}],sqww[{{0,1},{0,1},{1,0},{1,0}}],sqww[{{0,1},{0,1},{1,1},{1,1}}],sqww[{{1,0},{1,0},{1,1},{1,1}}],sqww[{{1,1},{1,1},{1,1},{1,1}}]}={0,6/5,3/2,13/5,19/10,6/5,8/5,2,12/5,5/2,11/5};*)*)


(* ::Subsection::Closed:: *)
(*w 3d*)


(* ::Subsubsection::Closed:: *)
(*Simpler pieces ({gw, bw, sw, v\[Psi]aw, v\[Psi]bw}) 3d*)


basicsequnce={{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1},{1,1,1}};


(* ::Input:: *)
(*Rationalize/@{*)
(*{0,1,1.3,1.6,2.4,2.8,3.3,5},*)
(*{0,1,1.2,1.4,1.6,1.8,2.2,3},*)
(*{0,1,1.2,1.4,1.8,2,2.4,3.2},*)
(*{0,1,1.2,1.4,1.6,1.8,2.2,3},*)
(*{0,1,1.2,1.4,10,10,10,15}}*)


(* ::Input:: *)
(*(*Clear@v\[Psi]bww*)*)


{gww[{0,0,0}],gww[{1,0,0}],gww[{0,1,0}],gww[{0,0,1}],gww[{1,1,0}],gww[{1,0,1}],gww[{0,1,1}],gww[{1,1,1}]}={0,1,1.3,1.6,2.4,2.8,3.3,5};


(* ::Input:: *)
(*{gww/@basicsequnce,bww/@basicsequnce,sww/@basicsequnce,v\[Psi]aww/@basicsequnce,v\[Psi]bww/@basicsequnce};*)


{{bww[{0,0,0}],bww[{1,0,0}],bww[{0,1,0}],bww[{0,0,1}],bww[{1,1,0}],bww[{1,0,1}],bww[{0,1,1}],bww[{1,1,1}]},{sww[{0,0,0}],sww[{1,0,0}],sww[{0,1,0}],sww[{0,0,1}],sww[{1,1,0}],sww[{1,0,1}],sww[{0,1,1}],sww[{1,1,1}]},{v\[Psi]aww[{0,0,0}],v\[Psi]aww[{1,0,0}],v\[Psi]aww[{0,1,0}],v\[Psi]aww[{0,0,1}],v\[Psi]aww[{1,1,0}],v\[Psi]aww[{1,0,1}],v\[Psi]aww[{0,1,1}],v\[Psi]aww[{1,1,1}]},{v\[Psi]bww[{0,0,0}],v\[Psi]bww[{1,0,0}],v\[Psi]bww[{0,1,0}],v\[Psi]bww[{0,0,1}],v\[Psi]bww[{1,1,0}],v\[Psi]bww[{1,0,1}],v\[Psi]bww[{0,1,1}],v\[Psi]bww[{1,1,1}]}}={{0,1,6/5,7/5,8/5,9/5,11/5,3},{0,1,6/5,7/5,9/5,2,12/5,16/5},{0,1,6/5,7/5,8/5,9/5,11/5,3},{0,1,6/5,7/5,10,10,10,15}};


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(*v\[Psi]aww/@{{0,0,0},{1,0,0},{0,0,1},{0,1,0},{1,1,0},{1,0,1},{0,1,1},{1,1,1}}*)*)


(* ::Input:: *)
(*(*{v\[Psi]aww[{0,0,0}],v\[Psi]aww[{1,0,0}],v\[Psi]aww[{0,0,1}],v\[Psi]aww[{0,1,0}],v\[Psi]aww[{1,1,0}],v\[Psi]aww[{1,0,1}],v\[Psi]aww[{0,1,1}],v\[Psi]aww[{1,1,1}]}={0,1,6/5,7/5,8/5,9/5,11/5,3};*)*)


(* ::Subsubsection::Closed:: *)
(*Triangles 3d*)


(* ::Input:: *)
(*Sort@DeleteDuplicates[Sort/@Abs/@Select[Flatten[Table[{Tuples[Range[-1,1],3][[i]],Tuples[Range[-1,1],3][[j]],Tuples[Range[-1,1],3][[i]]+Tuples[Range[-1,1],3][[j]]},{i,1,27},{j,1,27}],1],!MemberQ[Abs/@Flatten@#,2]&]]*)
(*Length@%*)


(* ::Input:: *)
(*(*Sort/@reverse2/@%*)*)


(* ::Text:: *)
(*New order*)


orderokTrPar={{{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{1,0,0},{1,0,0}},{{0,0,0},{0,1,0},{0,1,0}},{{0,0,0},{0,0,1},{0,0,1}},{{0,0,0},{1,1,0},{1,1,0}},{{0,0,0},{1,0,1},{1,0,1}},{{0,0,0},{0,1,1},{0,1,1}},{{0,0,0},{1,1,1},{1,1,1}},{{0,1,0},{1,0,0},{1,1,0}},{{0,0,1},{1,0,0},{1,0,1}},{{0,0,1},{0,1,0},{0,1,1}},{{0,1,1},{1,0,0},{1,1,1}},{{0,1,0},{1,0,1},{1,1,1}},{{0,0,1},{1,1,0},{1,1,1}},{{0,1,1},{1,0,1},{1,1,0}}};


weightTr=Rationalize@{0,11/10,7/5,9/5,5/2,29/10,7/2,51/10,1.8,2,2.4,3.6,3.5,3.4,2.7};


(* ::Input:: *)
(*matrixform2/@Table[{(showbetter3/@orderokTrPar)[[i]],{N@weightTr[[i]]}},{i,15}]*)


(* ::Input:: *)
(*(*tableform2/@showbetter3/@orderokTrPar*)*)


(* ::Input:: *)
(*(*tww/@orderokTrPar*)*)


(* ::Input:: *)
(*(*Rationalize[{0,1,1.3,1.7,2.4,2.8,3.4,5}+0.1]*)*)


(* ::Input:: *)
(**)


{tww[{{0,0,0},{0,0,0},{0,0,0}}],tww[{{0,0,0},{1,0,0},{1,0,0}}],tww[{{0,0,0},{0,1,0},{0,1,0}}],tww[{{0,0,0},{0,0,1},{0,0,1}}],tww[{{0,0,0},{1,1,0},{1,1,0}}],tww[{{0,0,0},{1,0,1},{1,0,1}}],tww[{{0,0,0},{0,1,1},{0,1,1}}],tww[{{0,0,0},{1,1,1},{1,1,1}}],tww[{{0,1,0},{1,0,0},{1,1,0}}],tww[{{0,0,1},{1,0,0},{1,0,1}}],tww[{{0,0,1},{0,1,0},{0,1,1}}],tww[{{0,1,1},{1,0,0},{1,1,1}}],tww[{{0,1,0},{1,0,1},{1,1,1}}],tww[{{0,0,1},{1,1,0},{1,1,1}}],tww[{{0,1,1},{1,0,1},{1,1,0}}]}={0,11/10,7/5,9/5,5/2,29/10,7/2,51/10,9/5,2,12/5,18/5,7/2,17/5,27/10};


(* ::Subsubsection::Closed:: *)
(*Plot ({gw, bw, sw, v\[Psi]aw, v\[Psi]bw, tw}) 3d*)


(* ::Text:: *)
(*Rationalize /@ {*)
(*   {0, 1, 1.3, 1.6, 2.4, 2.8, 3.3, 5},*)
(*   {0, 1, 1.2, 1.4, 1.6, 1.8, 2.2, 3},*)
(*   {0, 1, 1.2, 1.4, 1.8, 2, 2.4, 3.2},*)
(*   weightTr};*)
(*ListPlot[%, Filling -> Axis, PlotLegends -> {"gw", "bw", "sw", "tw"}, ImageSize -> Large]*)


(* ::Subsubsection::Closed:: *)
(*Squares 3d*)


sqparvec=Sort@DeleteDuplicates[Sort/@Abs/@Select[Flatten[Table[{Tuples[Range[-1,1],3][[i]],Tuples[Range[-1,1],3][[j]],Tuples[Range[-1,1],3][[y]],Tuples[Range[-1,1],3][[i]]+Tuples[Range[-1,1],3][[j]]+Tuples[Range[-1,1],3][[y]]},{i,1,27},{j,1,27},{y,1,27}],2],!MemberQ[Abs/@Flatten@#,2|3]&]];


(* ::Input:: *)
(*Length@sqparvec*)


sqpar=showbetter4/@Sort/@reverse2/@sqparvec;


(* ::Text:: *)
(*New order*)


orderoksqpar={{{{0,0,0},{0,0,0},{0,0,0},{0,0,0}},{{0,0,0},{0,0,0},{q1,0,0},{q1,0,0}},{{0,0,0},{0,0,0},{0,q2,0},{0,q2,0}},{{0,0,0},{0,0,0},{0,0,q3},{0,0,q3}},{{0,0,0},{0,0,0},{q1,q2,0},{q1,q2,0}},{{0,0,0},{0,0,0},{q1,0,q3},{q1,0,q3}},{{0,0,0},{0,0,0},{0,q2,q3},{0,q2,q3}},{{0,0,0},{0,0,0},{q1,q2,q3},{q1,q2,q3}},{{0,0,0},{0,q2,0},{q1,0,0},{q1,q2,0}},{{0,0,0},{0,0,q3},{q1,0,0},{q1,0,q3}},{{0,0,0},{0,0,q3},{0,q2,0},{0,q2,q3}},{{0,0,0},{0,q2,q3},{q1,0,0},{q1,q2,q3}},{{0,0,0},{0,q2,0},{q1,0,q3},{q1,q2,q3}},{{0,0,0},{0,0,q3},{q1,q2,0},{q1,q2,q3}},{{0,0,0},{0,q2,q3},{q1,0,q3},{q1,q2,0}}},
{{{q1,0,0},{q1,0,0},{q1,0,0},{q1,0,0}},{{0,q2,0},{0,q2,0},{q1,0,0},{q1,0,0}},{{q1,0,0},{q1,0,0},{q1,q2,0},{q1,q2,0}},{{0,0,q3},{0,0,q3},{q1,0,0},{q1,0,0}},{{q1,0,0},{q1,0,0},{q1,0,q3},{q1,0,q3}},{{0,q2,q3},{0,q2,q3},{q1,0,0},{q1,0,0}},{{q1,0,0},{q1,0,0},{q1,q2,q3},{q1,q2,q3}},
{{0,q2,0},{0,q2,0},{0,q2,0},{0,q2,0}},{{0,q2,0},{0,q2,0},{q1,q2,0},{q1,q2,0}},{{0,0,q3},{0,0,q3},{0,q2,0},{0,q2,0}},{{0,q2,0},{0,q2,0},{q1,0,q3},{q1,0,q3}},{{0,q2,0},{0,q2,0},{0,q2,q3},{0,q2,q3}},{{0,q2,0},{0,q2,0},{q1,q2,q3},{q1,q2,q3}},
{{q1,q2,0},{q1,q2,0},{q1,q2,0},{q1,q2,0}},{{0,0,q3},{0,0,q3},{q1,q2,0},{q1,q2,0}},{{q1,0,q3},{q1,0,q3},{q1,q2,0},{q1,q2,0}},{{0,q2,q3},{0,q2,q3},{q1,q2,0},{q1,q2,0}},{{q1,q2,0},{q1,q2,0},{q1,q2,q3},{q1,q2,q3}},
{{0,0,q3},{0,0,q3},{0,0,q3},{0,0,q3}},{{0,0,q3},{0,0,q3},{q1,0,q3},{q1,0,q3}},{{0,0,q3},{0,0,q3},{0,q2,q3},{0,q2,q3}},{{0,0,q3},{0,0,q3},{q1,q2,q3},{q1,q2,q3}},
{{q1,0,q3},{q1,0,q3},{q1,0,q3},{q1,0,q3}},{{0,q2,q3},{0,q2,q3},{q1,0,q3},{q1,0,q3}},{{q1,0,q3},{q1,0,q3},{q1,q2,q3},{q1,q2,q3}},
{{0,q2,q3},{0,q2,q3},{0,q2,q3},{0,q2,q3}},{{0,q2,q3},{0,q2,q3},{q1,q2,q3},{q1,q2,q3}},{{q1,q2,q3},{q1,q2,q3},{q1,q2,q3},{q1,q2,q3}}},
{{{0,q2,0},{0,q2,q3},{q1,0,0},{q1,0,q3}},{{0,0,q3},{0,q2,q3},{q1,0,0},{q1,q2,0}},{{0,0,q3},{0,q2,0},{q1,0,q3},{q1,q2,0}},
{{q1,0,0},{q1,0,q3},{q1,q2,0},{q1,q2,q3}},{{0,q2,0},{0,q2,q3},{q1,q2,0},{q1,q2,q3}},{{0,0,q3},{0,q2,q3},{q1,0,q3},{q1,q2,q3}},
{{0,0,q3},{0,q2,0},{q1,0,0},{q1,q2,q3}}}
};


(* ::Input:: *)
(*Length/@orderoksqpar*)


(* ::Input:: *)
(*(*Rationalize[{{0,11/10,7/5,9/5,5/2,29/10,7/2,51/10,9/5,2,12/5,18/5,7/2,17/5,27/10}+0.1}]*)*)


tmppp0={0,6/5,3/2,19/10,13/5,3,18/5,26/5,19/10,21/10,5/2,37/10,18/5,7/2,14/5};


tmppp={{{q1,0,0},{q1,0,0},{q1,0,0},{q1,0,0}},{{0,q2,0},{0,q2,0},{q1,0,0},{q1,0,0}},{{q1,0,0},{q1,0,0},{q1,q2,0},{q1,q2,0}},{{0,0,q3},{0,0,q3},{q1,0,0},{q1,0,0}},{{q1,0,0},{q1,0,0},{q1,0,q3},{q1,0,q3}},{{0,q2,q3},{0,q2,q3},{q1,0,0},{q1,0,0}},{{q1,0,0},{q1,0,0},{q1,q2,q3},{q1,q2,q3}},
{{0,q2,0},{0,q2,0},{0,q2,0},{0,q2,0}},{{0,q2,0},{0,q2,0},{q1,q2,0},{q1,q2,0}},{{0,0,q3},{0,0,q3},{0,q2,0},{0,q2,0}},{{0,q2,0},{0,q2,0},{q1,0,q3},{q1,0,q3}},{{0,q2,0},{0,q2,0},{0,q2,q3},{0,q2,q3}},{{0,q2,0},{0,q2,0},{q1,q2,q3},{q1,q2,q3}},
{{q1,q2,0},{q1,q2,0},{q1,q2,0},{q1,q2,0}},{{0,0,q3},{0,0,q3},{q1,q2,0},{q1,q2,0}},{{q1,0,q3},{q1,0,q3},{q1,q2,0},{q1,q2,0}},{{0,q2,q3},{0,q2,q3},{q1,q2,0},{q1,q2,0}},{{q1,q2,0},{q1,q2,0},{q1,q2,q3},{q1,q2,q3}},
{{0,0,q3},{0,0,q3},{0,0,q3},{0,0,q3}},{{0,0,q3},{0,0,q3},{q1,0,q3},{q1,0,q3}},{{0,0,q3},{0,0,q3},{0,q2,q3},{0,q2,q3}},{{0,0,q3},{0,0,q3},{q1,q2,q3},{q1,q2,q3}},
{{q1,0,q3},{q1,0,q3},{q1,0,q3},{q1,0,q3}},{{0,q2,q3},{0,q2,q3},{q1,0,q3},{q1,0,q3}},{{q1,0,q3},{q1,0,q3},{q1,q2,q3},{q1,q2,q3}},
{{0,q2,q3},{0,q2,q3},{0,q2,q3},{0,q2,q3}},{{0,q2,q3},{0,q2,q3},{q1,q2,q3},{q1,q2,q3}},{{q1,q2,q3},{q1,q2,q3},{q1,q2,q3},{q1,q2,q3}}};


tmppp2={{{0,q2,0},{0,q2,q3},{q1,0,0},{q1,0,q3}},{{0,0,q3},{0,q2,q3},{q1,0,0},{q1,q2,0}},{{0,0,q3},{0,q2,0},{q1,0,q3},{q1,q2,0}},
{{q1,0,0},{q1,0,q3},{q1,q2,0},{q1,q2,q3}},{{0,q2,0},{0,q2,q3},{q1,q2,0},{q1,q2,q3}},{{0,0,q3},{0,q2,q3},{q1,0,q3},{q1,q2,q3}},
{{0,0,q3},{0,q2,0},{q1,0,0},{q1,q2,q3}}};


(* ::Input:: *)
(*tableform2/@Table[{tmppp[[i]],{}},{i,28}];*)


(* ::Text:: *)
(*Assigning*)


(* ::Input:: *)
(*{{\!\(\**)
(*TagBox[GridBox[{*)
(*{"q1", "0", "0"},*)
(*{"q1", "0", "0"},*)
(*{"q1", "0", "0"},*)
(*{"q1", "0", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "1.2", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "q2", "0"},*)
(*{"0", "q2", "0"},*)
(*{"q1", "0", "0"},*)
(*{"q1", "0", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"q1", "0", "0"},*)
(*{"q1", "0", "0"},*)
(*{"q1", "q2", "0"},*)
(*{"q1", "q2", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.4", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "0", "q3"},*)
(*{"0", "0", "q3"},*)
(*{"q1", "0", "0"},*)
(*{"q1", "0", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.2", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"q1", "0", "0"},*)
(*{"q1", "0", "0"},*)
(*{"q1", "0", "q3"},*)
(*{"q1", "0", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.5", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "q2", "q3"},*)
(*{"0", "q2", "q3"},*)
(*{"q1", "0", "0"},*)
(*{"q1", "0", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.7", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"q1", "0", "0"},*)
(*{"q1", "0", "0"},*)
(*{"q1", "q2", "q3"},*)
(*{"q1", "q2", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.9", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "q2", "0"},*)
(*{"0", "q2", "0"},*)
(*{"0", "q2", "0"},*)
(*{"0", "q2", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "1.6", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "q2", "0"},*)
(*{"0", "q2", "0"},*)
(*{"q1", "q2", "0"},*)
(*{"q1", "q2", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.5", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "0", "q3"},*)
(*{"0", "0", "q3"},*)
(*{"0", "q2", "0"},*)
(*{"0", "q2", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.4", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "q2", "0"},*)
(*{"0", "q2", "0"},*)
(*{"q1", "0", "q3"},*)
(*{"q1", "0", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.6", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "q2", "0"},*)
(*{"0", "q2", "0"},*)
(*{"0", "q2", "q3"},*)
(*{"0", "q2", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.6", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "q2", "0"},*)
(*{"0", "q2", "0"},*)
(*{"q1", "q2", "q3"},*)
(*{"q1", "q2", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "3", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"q1", "q2", "0"},*)
(*{"q1", "q2", "0"},*)
(*{"q1", "q2", "0"},*)
(*{"q1", "q2", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.2", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "0", "q3"},*)
(*{"0", "0", "q3"},*)
(*{"q1", "q2", "0"},*)
(*{"q1", "q2", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.5", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"q1", "0", "q3"},*)
(*{"q1", "0", "q3"},*)
(*{"q1", "q2", "0"},*)
(*{"q1", "q2", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.7", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "q2", "q3"},*)
(*{"0", "q2", "q3"},*)
(*{"q1", "q2", "0"},*)
(*{"q1", "q2", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.9", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"q1", "q2", "0"},*)
(*{"q1", "q2", "0"},*)
(*{"q1", "q2", "q3"},*)
(*{"q1", "q2", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "3.1", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "0", "q3"},*)
(*{"0", "0", "q3"},*)
(*{"0", "0", "q3"},*)
(*{"0", "0", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "1.8", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "0", "q3"},*)
(*{"0", "0", "q3"},*)
(*{"q1", "0", "q3"},*)
(*{"q1", "0", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.6", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "0", "q3"},*)
(*{"0", "0", "q3"},*)
(*{"0", "q2", "q3"},*)
(*{"0", "q2", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.8", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "0", "q3"},*)
(*{"0", "0", "q3"},*)
(*{"q1", "q2", "q3"},*)
(*{"q1", "q2", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "3.3", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"q1", "0", "q3"},*)
(*{"q1", "0", "q3"},*)
(*{"q1", "0", "q3"},*)
(*{"q1", "0", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.3", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "q2", "q3"},*)
(*{"0", "q2", "q3"},*)
(*{"q1", "0", "q3"},*)
(*{"q1", "0", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "3", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"q1", "0", "q3"},*)
(*{"q1", "0", "q3"},*)
(*{"q1", "q2", "q3"},*)
(*{"q1", "q2", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "3.2", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "q2", "q3"},*)
(*{"0", "q2", "q3"},*)
(*{"0", "q2", "q3"},*)
(*{"0", "q2", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.5", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "q2", "q3"},*)
(*{"0", "q2", "q3"},*)
(*{"q1", "q2", "q3"},*)
(*{"q1", "q2", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "3.3", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"q1", "q2", "q3"},*)
(*{"q1", "q2", "q3"},*)
(*{"q1", "q2", "q3"},*)
(*{"q1", "q2", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "3.4", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)}};*)


(* ::Input:: *)
(*Rationalize@Flatten[%[[All,2]]]*)


(* ::Input:: *)
(*tableform2/@Table[{tmppp2[[i]],{}},{i,Length@tmppp2}];*)


(* ::Input:: *)
(*{{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "q2", "0"},*)
(*{"0", "q2", "q3"},*)
(*{"q1", "0", "0"},*)
(*{"q1", "0", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.7", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "0", "q3"},*)
(*{"0", "q2", "q3"},*)
(*{"q1", "0", "0"},*)
(*{"q1", "q2", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.8", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "0", "q3"},*)
(*{"0", "q2", "0"},*)
(*{"q1", "0", "q3"},*)
(*{"q1", "q2", "0"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "2.9", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"q1", "0", "0"},*)
(*{"q1", "0", "q3"},*)
(*{"q1", "q2", "0"},*)
(*{"q1", "q2", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "3.3", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "q2", "0"},*)
(*{"0", "q2", "q3"},*)
(*{"q1", "q2", "0"},*)
(*{"q1", "q2", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "3.4", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "0", "q3"},*)
(*{"0", "q2", "q3"},*)
(*{"q1", "0", "q3"},*)
(*{"q1", "q2", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "3.5", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)},{\!\(\**)
(*TagBox[GridBox[{*)
(*{"0", "0", "q3"},*)
(*{"0", "q2", "0"},*)
(*{"q1", "0", "0"},*)
(*{"q1", "q2", "q3"}*)
(*},*)
(*GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},*)
(*GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[2.0999999999999996`]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\),\!\(\**)
(*TagBox[*)
(*RowBox[{"{", "3", "}"}],*)
(*Function[BoxForm`e$, TableForm[BoxForm`e$]]]\)}};*)


(* ::Input:: *)
(*Rationalize@Flatten[%[[All,2]]]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Flatten@Table[Position[sqpar,Flatten[orderoksqpar,1][[i]]][[1]],{i,50}]*)
(*sqparvecordOK=sqparvec[[%]]*)


weightSq=Flatten@{{0,6/5,3/2,19/10,13/5,3,18/5,26/5,19/10,21/10,5/2,37/10,18/5,7/2,14/5},{6/5,2,12/5,11/5,5/2,27/10,29/10,8/5,5/2,12/5,13/5,13/5,3,11/5,5/2,27/10,29/10,31/10,9/5,13/5,14/5,33/10,23/10,3,16/5,5/2,33/10,17/5},{27/10,14/5,29/10,33/10,17/5,7/2,3}};


(* ::Input:: *)
(*sqww/@sqparvecordOK;*)


{sqww[{{0,0,0},{0,0,0},{0,0,0},{0,0,0}}],sqww[{{0,0,0},{0,0,0},{0,0,1},{0,0,1}}],sqww[{{0,0,0},{0,0,0},{0,1,0},{0,1,0}}],sqww[{{0,0,0},{0,0,0},{1,0,0},{1,0,0}}],sqww[{{0,0,0},{0,0,0},{0,1,1},{0,1,1}}],sqww[{{0,0,0},{0,0,0},{1,0,1},{1,0,1}}],sqww[{{0,0,0},{0,0,0},{1,1,0},{1,1,0}}],sqww[{{0,0,0},{0,0,0},{1,1,1},{1,1,1}}],sqww[{{0,0,0},{0,0,1},{0,1,0},{0,1,1}}],sqww[{{0,0,0},{0,0,1},{1,0,0},{1,0,1}}],sqww[{{0,0,0},{0,1,0},{1,0,0},{1,1,0}}],sqww[{{0,0,0},{0,0,1},{1,1,0},{1,1,1}}],sqww[{{0,0,0},{0,1,0},{1,0,1},{1,1,1}}],sqww[{{0,0,0},{0,1,1},{1,0,0},{1,1,1}}],sqww[{{0,0,0},{0,1,1},{1,0,1},{1,1,0}}],sqww[{{0,0,1},{0,0,1},{0,0,1},{0,0,1}}],sqww[{{0,0,1},{0,0,1},{0,1,0},{0,1,0}}],sqww[{{0,0,1},{0,0,1},{0,1,1},{0,1,1}}],sqww[{{0,0,1},{0,0,1},{1,0,0},{1,0,0}}],sqww[{{0,0,1},{0,0,1},{1,0,1},{1,0,1}}],sqww[{{0,0,1},{0,0,1},{1,1,0},{1,1,0}}],sqww[{{0,0,1},{0,0,1},{1,1,1},{1,1,1}}],sqww[{{0,1,0},{0,1,0},{0,1,0},{0,1,0}}],sqww[{{0,1,0},{0,1,0},{0,1,1},{0,1,1}}],sqww[{{0,1,0},{0,1,0},{1,0,0},{1,0,0}}],sqww[{{0,1,0},{0,1,0},{1,0,1},{1,0,1}}],sqww[{{0,1,0},{0,1,0},{1,1,0},{1,1,0}}],sqww[{{0,1,0},{0,1,0},{1,1,1},{1,1,1}}],sqww[{{0,1,1},{0,1,1},{0,1,1},{0,1,1}}],sqww[{{0,1,1},{0,1,1},{1,0,0},{1,0,0}}],sqww[{{0,1,1},{0,1,1},{1,0,1},{1,0,1}}],sqww[{{0,1,1},{0,1,1},{1,1,0},{1,1,0}}],sqww[{{0,1,1},{0,1,1},{1,1,1},{1,1,1}}],sqww[{{1,0,0},{1,0,0},{1,0,0},{1,0,0}}],sqww[{{1,0,0},{1,0,0},{1,0,1},{1,0,1}}],sqww[{{1,0,0},{1,0,0},{1,1,0},{1,1,0}}],sqww[{{1,0,0},{1,0,0},{1,1,1},{1,1,1}}],sqww[{{1,0,1},{1,0,1},{1,0,1},{1,0,1}}],sqww[{{1,0,1},{1,0,1},{1,1,0},{1,1,0}}],sqww[{{1,0,1},{1,0,1},{1,1,1},{1,1,1}}],sqww[{{1,1,0},{1,1,0},{1,1,0},{1,1,0}}],sqww[{{1,1,0},{1,1,0},{1,1,1},{1,1,1}}],sqww[{{1,1,1},{1,1,1},{1,1,1},{1,1,1}}],sqww[{{0,0,1},{0,1,0},{1,0,1},{1,1,0}}],sqww[{{0,0,1},{0,1,1},{1,0,0},{1,1,0}}],sqww[{{0,1,0},{0,1,1},{1,0,0},{1,0,1}}],sqww[{{0,0,1},{0,1,1},{1,0,1},{1,1,1}}],sqww[{{0,1,0},{0,1,1},{1,1,0},{1,1,1}}],sqww[{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}],sqww[{{0,0,1},{0,1,0},{1,0,0},{1,1,1}}]}={0,6/5,3/2,19/10,13/5,3,18/5,26/5,19/10,21/10,5/2,37/10,18/5,7/2,14/5,6/5,2,12/5,11/5,5/2,27/10,29/10,8/5,5/2,12/5,13/5,13/5,3,11/5,5/2,27/10,29/10,31/10,9/5,13/5,14/5,33/10,23/10,3,16/5,5/2,33/10,17/5,27/10,14/5,29/10,33/10,17/5,7/2,3};


(* ::Input:: *)
(*(*{sqww[{{0,0,0},{0,0,0},{0,0,0},{0,0,0}}],sqww[{{0,0,0},{0,0,0},{0,0,1},{0,0,1}}],sqww[{{0,0,0},{0,0,0},{0,1,0},{0,1,0}}],sqww[{{0,0,0},{0,0,0},{1,0,0},{1,0,0}}],sqww[{{0,0,0},{0,0,0},{0,1,1},{0,1,1}}],sqww[{{0,0,0},{0,0,0},{1,0,1},{1,0,1}}],sqww[{{0,0,0},{0,0,0},{1,1,0},{1,1,0}}],sqww[{{0,0,0},{0,0,0},{1,1,1},{1,1,1}}],sqww[{{0,0,0},{0,0,1},{0,1,0},{0,1,1}}],sqww[{{0,0,0},{0,0,1},{1,0,0},{1,0,1}}],sqww[{{0,0,0},{0,1,0},{1,0,0},{1,1,0}}],sqww[{{0,0,0},{0,0,1},{1,1,0},{1,1,1}}],sqww[{{0,0,0},{0,1,0},{1,0,1},{1,1,1}}],sqww[{{0,0,0},{0,1,1},{1,0,0},{1,1,1}}],sqww[{{0,0,0},{0,1,1},{1,0,1},{1,1,0}}],sqww[{{0,0,1},{0,0,1},{0,0,1},{0,0,1}}],sqww[{{0,0,1},{0,0,1},{0,1,0},{0,1,0}}],sqww[{{0,0,1},{0,0,1},{0,1,1},{0,1,1}}],sqww[{{0,0,1},{0,0,1},{1,0,0},{1,0,0}}],sqww[{{0,0,1},{0,0,1},{1,0,1},{1,0,1}}],sqww[{{0,0,1},{0,0,1},{1,1,0},{1,1,0}}],sqww[{{0,0,1},{0,0,1},{1,1,1},{1,1,1}}],sqww[{{0,1,0},{0,1,0},{0,1,0},{0,1,0}}],sqww[{{0,1,0},{0,1,0},{0,1,1},{0,1,1}}],sqww[{{0,1,0},{0,1,0},{1,0,0},{1,0,0}}],sqww[{{0,1,0},{0,1,0},{1,0,1},{1,0,1}}],sqww[{{0,1,0},{0,1,0},{1,1,0},{1,1,0}}],sqww[{{0,1,0},{0,1,0},{1,1,1},{1,1,1}}],sqww[{{0,1,1},{0,1,1},{0,1,1},{0,1,1}}],sqww[{{0,1,1},{0,1,1},{1,0,0},{1,0,0}}],sqww[{{0,1,1},{0,1,1},{1,0,1},{1,0,1}}],sqww[{{0,1,1},{0,1,1},{1,1,0},{1,1,0}}],sqww[{{0,1,1},{0,1,1},{1,1,1},{1,1,1}}],sqww[{{1,0,0},{1,0,0},{1,0,0},{1,0,0}}],sqww[{{1,0,0},{1,0,0},{1,0,1},{1,0,1}}],sqww[{{1,0,0},{1,0,0},{1,1,0},{1,1,0}}],sqww[{{1,0,0},{1,0,0},{1,1,1},{1,1,1}}],sqww[{{1,0,1},{1,0,1},{1,0,1},{1,0,1}}],sqww[{{1,0,1},{1,0,1},{1,1,0},{1,1,0}}],sqww[{{1,0,1},{1,0,1},{1,1,1},{1,1,1}}],sqww[{{1,1,0},{1,1,0},{1,1,0},{1,1,0}}],sqww[{{1,1,0},{1,1,0},{1,1,1},{1,1,1}}],sqww[{{1,1,1},{1,1,1},{1,1,1},{1,1,1}}],sqww[{{0,0,1},{0,1,0},{1,0,1},{1,1,0}}],sqww[{{0,0,1},{0,1,1},{1,0,0},{1,1,0}}],sqww[{{0,1,0},{0,1,1},{1,0,0},{1,0,1}}],sqww[{{0,0,1},{0,1,1},{1,0,1},{1,1,1}}],sqww[{{0,1,0},{0,1,1},{1,1,0},{1,1,1}}],sqww[{{1,0,0},{1,0,1},{1,1,0},{1,1,1}}],sqww[{{0,0,1},{0,1,0},{1,0,0},{1,1,1}}]}={0,6/5,3/2,19/10,13/5,3,18/5,26/5,19/10,21/10,5/2,37/10,18/5,7/2,14/5,6/5,2,12/5,11/5,5/2,27/10,29/10,8/5,5/2,12/5,13/5,13/5,3,11/5,5/2,27/10,29/10,31/10,9/5,13/5,14/5,33/10,23/10,3,16/5,5/2,33/10,17/5,27/10,14/5,29/10,33/10,17/5,7/2,3};*)*)


(* ::Input:: *)
(**)


(* ::Text:: *)
(*Plot with Sqw also*)


(* ::Input:: *)
(*Rationalize/@{*)
(*{0,1,1.3,1.7,2.4,2.8,3.4,5},*)
(*{0,1,1.2,1.4,1.6,1.8,2.2,3},*)
(*{0,1,1.2,1.4,1.8,2,2.4,3.2},*)
(*weightTr,weightSq};*)
(*ListPlot[%,Filling->Axis,PlotLegends->{"gw","bw","sw","tw","sqw"},ImageSize->Large]*)


(* ::Subsection::Closed:: *)
(*Substitution*)


subWW={gw-> gww, bw-> bww, sw-> sww,v\[Psi]aw-> v\[Psi]aww,v\[Psi]bw-> v\[Psi]bww, tw-> tww,sqw-> sqww,v\[Kappa]w-> v\[Kappa]ww,v\[Theta]w->v\[Theta]ww ,v\[Theta]2w-> v\[Theta]2ww,v\[Sigma]w->v\[Sigma]ww,v\[Sigma]2w->v\[Sigma]2ww };


(* ::Subsection::Closed:: *)
(*Procedure*)


(* ::Subsubsection::Closed:: *)
(*Procedure 2d*)


bestReparam4ptv02d[o_,grp_]:=Module[{conf12,grpconf,matOktmp,weights,positionbest},
conf12=GtP4ptsdConfig[0,o,grp][[1]];
grpconf=tovector2/@(conf12[[2]]);
matOktmp=matricesTransOKconf2d[DeleteDuplicates@grpconf];
weights=N[Table[Total@Flatten[weightprov[conf12[[1]],Abs/@Table[grpconf[[i]] . matOktmp[[j]],{i,Length@grpconf}]]],{j,Length@matOktmp}]/.subWW];
positionbest=Flatten@Position[weights,Min[weights]];
matOktmp[[positionbest]](*\[LeftDoubleBracket]-1\[RightDoubleBracket]*)
]


reparmom2[m_]:=Module[{nrep},
nrep=m . {q1,q2};
{q1-> nrep[[1]],q2-> nrep[[2]]}
]


(*bestReparam4ptv02dinside[o_,grp_]:=Module[{conf12,grpconf,matOktmp,weights,positionbest},
conf12=GtP4ptsdConfig[0,o,grp][[1]];
grpconf=tovector2/@(conf12[[2]]);
matOktmp=matricesTransOKconf2d[DeleteDuplicates@grpconf];
weights=N[Table[Total@Flatten[weightprov[conf12[[1]],Abs/@Table[grpconf[[i]] . matOktmp[[j]],{i,Length@grpconf}]]],{j,Length@matOktmp}]/.subWW];
positionbest=Flatten@Position[weights,Min[weights]];
Print[weights,positionbest];
matOktmp[[positionbest]](*\[LeftDoubleBracket]-1\[RightDoubleBracket]*)
]*)


bestReparam4ptv0lev42d[o_,grp_]:=Module[{conf12,grpconf,matOktmp,weights,positionbest},
conf12=GtP4ptsdConfiglevel4[0,o,grp][[1]];
grpconf=tovector2/@(conf12[[2]]);
matOktmp=matricesTransOKconf2d[DeleteDuplicates@grpconf];
weights=N[Table[Total@Flatten[weightprovlev4[conf12[[1]],Abs/@Table[grpconf[[i]] . matOktmp[[j]],{i,Length@grpconf}]]],{j,Length@matOktmp}]/.subWW];
positionbest=Flatten@Position[weights,Min[weights]];
matOktmp[[positionbest]](*\[LeftDoubleBracket]-1\[RightDoubleBracket]*)
]


(* ::Subsubsection::Closed:: *)
(*Procedure 3d*)


bestReparam4ptv0[o_,grp_]:=Module[{conf12,grpconf,matOktmp,weights,positionbest},
conf12=GtP4ptsdConfig[0,o,grp][[1]];
grpconf=tovector3/@(conf12[[2]]);
matOktmp=matricesTransOKconf[DeleteDuplicates@grpconf];
weights=N[Table[Total@Flatten[weightprov[conf12[[1]],Abs/@Table[grpconf[[i]] . matOktmp[[j]],{i,Length@grpconf}]]],{j,Length@matOktmp}]/.subWW];
positionbest=Flatten@Position[weights,Min[weights]];
matOktmp[[positionbest]](*\[LeftDoubleBracket]-1\[RightDoubleBracket]*)
]


reparmom3[m_]:=Module[{nrep},
nrep=m . {q1,q2,q3};
{q1-> nrep[[1]],q2-> nrep[[2]],q3-> nrep[[3]]}
]


(* ::Input:: *)
(*(*bestReparam4ptv0[6,{fourpt1PIbstc2W[0,6,standardweights]\[LeftDoubleBracket]order\[CapitalGamma]4o6\[RightDoubleBracket]\[LeftDoubleBracket]95\[RightDoubleBracket]}]*)*)


(* ::Input:: *)
(*(*bestReparam4ptv0lev4[o_,grp_]:=Module[{conf12,grpconf,matOktmp,weights,positionbest},*)
(*conf12=GtP4ptsdConfiglevel4[0,o,grp]\[LeftDoubleBracket]1\[RightDoubleBracket];*)
(*grpconf=tovector3/@(conf12\[LeftDoubleBracket]2\[RightDoubleBracket]);*)
(*matOktmp=matricesTransOKconf[DeleteDuplicates@grpconf];*)
(*weights=N[Table[Total@Flatten[weightprov[conf12\[LeftDoubleBracket]1\[RightDoubleBracket],Abs/@Table[grpconf\[LeftDoubleBracket]i\[RightDoubleBracket].matOktmp\[LeftDoubleBracket]j\[RightDoubleBracket],{i,Length@grpconf}]]],{j,Length@matOktmp}]/.subWW];*)
(*positionbest=Flatten@Position[weights,Min[weights]];*)
(*matOktmp\[LeftDoubleBracket]positionbest\[RightDoubleBracket](*\[LeftDoubleBracket]-1\[RightDoubleBracket]*)*)
(*]*)*)


bestReparam4ptv0lev4[o_,grp_]:=Module[{conf12,grpconf,matOktmp,weights,positionbest},
conf12=GtP4ptsdConfiglevel4[0,o,grp][[1]];
grpconf=tovector3/@(conf12[[2]]);
matOktmp=matricesTransOKconf[DeleteDuplicates@grpconf];
weights=N[Table[Total@Flatten[weightprovlev4[conf12[[1]],Abs/@Table[grpconf[[i]] . matOktmp[[j]],{i,Length@grpconf}]]],{j,Length@matOktmp}]/.subWW];
positionbest=Flatten@Position[weights,Min[weights]];
matOktmp[[positionbest]](*\[LeftDoubleBracket]-1\[RightDoubleBracket]*)
]


(* ::Subsection::Closed:: *)
(*Reassigning*)


(* ::Subsubsection::Closed:: *)
(*l = 2*)


reassignation2llev2o[o_,o6G4l2_]:=Module[{listdiagr2l,ll,tmpgrpform,changes},
listdiagr2l=o6G4l2[[All,1]];
ll=Length@listdiagr2l;
tmpgrpform=fourpt1PIbstc2[0,o][[order\[CapitalGamma]4o[o]]];
changes=Table[bestReparam4ptv02d[o,{tmpgrpform[[listdiagr2l[[i]]]]}],{i,ll}];
Table[{o6G4l2[[j,1]],o6G4l2[[j,2]]/.reparmom2[changes[[j,-1]]]},{j,ll}]
]


reassignation2llev4o[o_,o6G4l2_]:=Module[{listdiagr2l,ll,tmpgrpform,changes},
listdiagr2l=o6G4l2[[All,1]];
ll=Length@listdiagr2l;
tmpgrpform=fourpt1PIbstc4[0,o][[order\[CapitalGamma]4o[o]]];
changes=Table[bestReparam4ptv0lev42d[o,{tmpgrpform[[listdiagr2l[[i]]]]}],{i,ll}];
Table[{o6G4l2[[j,1]],o6G4l2[[j,2]]/.reparmom2[changes[[j,-1]]]},{j,ll}]
]


(* ::Subsubsection::Closed:: *)
(*l = 3*)


reassignation3llev2o[o_,o6G4l3_]:=Module[{listdiagr,ll,tmpgrpform,changes},
listdiagr=o6G4l3[[All,1]];
ll=Length@listdiagr;
tmpgrpform=fourpt1PIbstc2[0,o][[order\[CapitalGamma]4o[o]]];
changes=Table[bestReparam4ptv0[o,{tmpgrpform[[listdiagr[[i]]]]}],{i,ll}];
Table[{o6G4l3[[j,1]],o6G4l3[[j,2]]/.reparmom3[changes[[j,-1]]]},{j,ll}]
]


reassignation3llev2ao[o_,o6G4l3_]:=Module[{listdiagr,ll,tmpgrpform,changes},
listdiagr=o6G4l3[[All,1]];
ll=Length@listdiagr;
tmpgrpform=fourpt1PIbstc2a[0,o][[order\[CapitalGamma]4o[o]]];
changes=Table[bestReparam4ptv0[o,{tmpgrpform[[listdiagr[[i]]]]}],{i,ll}];
Table[{o6G4l3[[j,1]],o6G4l3[[j,2]]/.reparmom3[changes[[j,-1]]]},{j,ll}]
]


reassignation3llev2oW[o_,o6G4l3_,w_]:=Module[{listdiagr,ll,tmpgrpform,changes},
listdiagr=o6G4l3[[All,1]];
ll=Length@listdiagr;
tmpgrpform=fourpt1PIbstc2W[0,o,w][[order\[CapitalGamma]4o[o]]];
changes=Table[bestReparam4ptv0[o,{tmpgrpform[[listdiagr[[i]]]]}],{i,ll}];
Table[{o6G4l3[[j,1]],o6G4l3[[j,2]]/.reparmom3[changes[[j,-1]]]},{j,ll}]
]


(* ::Text:: *)
(*Weights not assigned*)


reassignation3llev4o[o_,o6G4l3_]:=Module[{listdiagr,ll,tmpgrpform,changes},
listdiagr=o6G4l3[[All,1]];
ll=Length@listdiagr;
tmpgrpform=fourpt1PIbstc4[0,o][[order\[CapitalGamma]4o[o]]];
changes=Table[bestReparam4ptv0lev4[o,{tmpgrpform[[listdiagr[[i]]]]}],{i,ll}];
Table[{o6G4l3[[j,1]],o6G4l3[[j,2]]/.reparmom3[changes[[j,-1]]]},{j,ll}]
]


(* ::Subsection::Closed:: *)
(*Generalization to \[CapitalGamma]^(2)*)


(* ::Subsubsection::Closed:: *)
(*Function momenta*)


(* ::Text:: *)
(*Level 3a (or lower)*)


GtP2ptsdConfig[n3_,n4_,twopt1PIod_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,trinvert,\[Psi]avert,vertex\[Psi]apos,\[Psi]bvert,vertex\[Psi]bpos,squarevert,vertexposquad,propagatorspos,deauxOK,vertexTripos,bubvert,vertexposbub,sunvert,vertexpossuns},
diaglist=twopt1PIod;dl=Length[diaglist];Table[deltas={};daux=replaceTad\[Tau]\[Beta][diaglist[[lj]][[2]]];pos=DeleteCases[VertexList[daux],0|o1|o2|o3|o4];If[Length[pos]==1,diaglist[[lj]][[1]] suntad^Count[pos,Subscript[\[Tau], _]] bubbl2tr^Count[pos,Subscript[\[Beta], _]],paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length[pos],idx++,xlist=Select[daux,MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},pos[[idx]]]&];xpos=DeleteDuplicates[Flatten[Table[Position[daux,xlist[[k]]],{k,1,Length[xlist]}]]];xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length[xlist]}];AppendTo[deltas,\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(Length[xpossign]\)]\(xpossign[\([k, 2]\)]\ paux[\([xpossign[\([k, 1]\)]]\)]\)\)==0];];deltasol=Reduce[deltas]/. And->List/. Equal:>Rule;(*Print[deltasol];*)
paux2=DeleteCases[DeleteDuplicates[Flatten[paux/. deltasol/. Plus->List]/. Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];paux3=paux/. deltasol/. Table[If[paux2=={},Null,paux2[[k]]->qauxify[k]],{k,1,Length[paux2]}];(*Print[paux3];*)
deauxOK=daux[[3;;-1]];
propagatorspos=Select[Range[Length@deauxOK],!MemberQ[deauxOK[[#]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4]&];
trinvert=Cases[pos,Subscript[t, _]];
vertexTripos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},trinvert[[i]]]&)]],{i,Length[trinvert]}];
\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Psi]avert[[i]]]&)]],{i,Length[\[Psi]avert]}];
\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Psi]bvert[[i]]]&)]],{i,Length[\[Psi]bvert]}];
squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},squarevert[[i]]]&)]],{i,Length[squarevert]}];
bubvert=Cases[pos,Subscript[b, _]];
vertexposbub=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},bubvert[[i]]]&)]],{i,Length[bubvert]}];
sunvert=Cases[pos,Subscript[s, _]];
vertexpossuns=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},sunvert[[i]]]&)]],{i,Length[sunvert]}];

{{propagatorspos,vertexposbub,vertexpossuns,vertex\[Psi]apos,vertex\[Psi]bpos,vertexTripos,vertexposquad},
paux3[[3;;-1]]/.{pext->0,-pext->0}}],{lj,1,dl}]]


(* ::Text:: *)
(*Level 4*)


GtP2ptsdConfiglevel4[n3_,n4_,twopt1PIod_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,trinvert,\[Psi]avert,vertex\[Psi]apos,\[Psi]bvert,vertex\[Psi]bpos,squarevert,vertexposquad,propagatorspos,deauxOK,vertexTripos,bubvert,vertexposbub,sunvert,vertexpossuns,vertex\[Kappa]pos,vertexpostheta,vertexpostheta2,vertexpossigma,vertexpossigma2,\[Kappa]vert,thetavert,theta2vert,sigmavert,sigma2vert},
diaglist=twopt1PIod;dl=Length[diaglist];Table[deltas={};daux=replaceTad\[Tau]\[Beta][diaglist[[lj]][[2]]];pos=DeleteCases[VertexList[daux],0|o1|o2|o3|o4];If[Length[pos]==1,diaglist[[lj]][[1]] suntad^Count[pos,Subscript[\[Tau], _]] bubbl2tr^Count[pos,Subscript[\[Beta], _]],paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length[pos],idx++,xlist=Select[daux,MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},pos[[idx]]]&];xpos=DeleteDuplicates[Flatten[Table[Position[daux,xlist[[k]]],{k,1,Length[xlist]}]]];xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length[xlist]}];AppendTo[deltas,\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(Length[xpossign]\)]\(xpossign[\([k, 2]\)]\ paux[\([xpossign[\([k, 1]\)]]\)]\)\)==0];];deltasol=Reduce[deltas]/. And->List/. Equal:>Rule;(*Print[deltasol];*)
paux2=DeleteCases[DeleteDuplicates[Flatten[paux/. deltasol/. Plus->List]/. Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];paux3=paux/. deltasol/. Table[If[paux2=={},Null,paux2[[k]]->qauxify[k]],{k,1,Length[paux2]}];(*Print[paux3];*)
deauxOK=daux[[3;;-1]];
propagatorspos=Select[Range[Length@deauxOK],!MemberQ[deauxOK[[#]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4]&];
trinvert=Cases[pos,Subscript[t, _]];
vertexTripos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},trinvert[[i]]]&)]],{i,Length[trinvert]}];
\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Psi]avert[[i]]]&)]],{i,Length[\[Psi]avert]}];
\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Psi]bvert[[i]]]&)]],{i,Length[\[Psi]bvert]}];
\[Kappa]vert=Cases[pos,Subscript[\[Kappa], _]];
vertex\[Kappa]pos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Kappa]vert[[i]]]&)]],{i,Length[\[Kappa]vert]}];
squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},squarevert[[i]]]&)]],{i,Length[squarevert]}];
bubvert=Cases[pos,Subscript[b, _]];
vertexposbub=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},bubvert[[i]]]&)]],{i,Length[bubvert]}];
sunvert=Cases[pos,Subscript[s, _]];
vertexpossuns=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},sunvert[[i]]]&)]],{i,Length[sunvert]}];
thetavert=Cases[pos,Subscript[\[Theta], _]];
vertexpostheta=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},thetavert[[i]]]&)]],{i,Length[thetavert]}];
theta2vert=Cases[pos,Subscript[\[Theta]2, _]];
vertexpostheta2=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},theta2vert[[i]]]&)]],{i,Length[theta2vert]}];
sigmavert=Cases[pos,Subscript[\[Sigma], _]];
vertexpossigma=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},sigmavert[[i]]]&)]],{i,Length[sigmavert]}];
sigma2vert=Cases[pos,Subscript[\[Sigma]2, _]];
vertexpossigma2=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},sigma2vert[[i]]]&)]],{i,Length[sigma2vert]}];

{{propagatorspos,vertexposbub,vertexpossuns,vertex\[Psi]apos,vertex\[Psi]bpos,vertexTripos,vertexposquad,vertex\[Kappa]pos,vertexpostheta,vertexpostheta2,vertexpossigma,vertexpossigma2},
paux3[[3;;-1]]/. {pext->0,-pext->0}}],{lj,1,dl}]]


(* ::Subsubsection::Closed:: *)
(*Procedure 2d*)


bestReparam2ptv02d[o_,grp_]:=Module[{conf12,grpconf,matOktmp,weights,positionbest},
conf12=GtP2ptsdConfig[0,o,grp][[1]];
grpconf=tovector2/@(conf12[[2]]);
matOktmp=matricesTransOKconf2d[DeleteDuplicates@grpconf];
weights=N[Table[Total@Flatten[weightprov[conf12[[1]],Abs/@Table[grpconf[[i]] . matOktmp[[j]],{i,Length@grpconf}]]],{j,Length@matOktmp}]/.subWW];
positionbest=Flatten@Position[weights,Min[weights]];
matOktmp[[positionbest]](*\[LeftDoubleBracket]-1\[RightDoubleBracket]*)
]


bestReparam2ptv0lev42d[o_,grp_]:=Module[{conf12,grpconf,matOktmp,weights,positionbest},
conf12=GtP2ptsdConfiglevel4[0,o,grp][[1]];
grpconf=tovector2/@(conf12[[2]]);
matOktmp=matricesTransOKconf2d[DeleteDuplicates@grpconf];
weights=N[Table[Total@Flatten[weightprovlev4[conf12[[1]],Abs/@Table[grpconf[[i]] . matOktmp[[j]],{i,Length@grpconf}]]],{j,Length@matOktmp}]/.subWW];
positionbest=Flatten@Position[weights,Min[weights]];
matOktmp[[positionbest]](*\[LeftDoubleBracket]-1\[RightDoubleBracket]*)
]


(* ::Subsubsection::Closed:: *)
(*Procedure 3d*)


bestReparam2ptv0[o_,grp_]:=Module[{conf12,grpconf,matOktmp,weights,positionbest},
conf12=GtP2ptsdConfig[0,o,grp][[1]];
grpconf=tovector3/@(conf12[[2]]);
matOktmp=matricesTransOKconf[DeleteDuplicates@grpconf];
weights=N[Table[Total@Flatten[weightprov[conf12[[1]],Abs/@Table[grpconf[[i]] . matOktmp[[j]],{i,Length@grpconf}]]],{j,Length@matOktmp}]/.subWW];
positionbest=Flatten@Position[weights,Min[weights]];
matOktmp[[positionbest]](*\[LeftDoubleBracket]-1\[RightDoubleBracket]*)
]


bestReparam2ptv0lev4[o_,grp_]:=Module[{conf12,grpconf,matOktmp,weights,positionbest},
conf12=GtP2ptsdConfiglevel4[0,o,grp][[1]];
grpconf=tovector3/@(conf12[[2]]);
matOktmp=matricesTransOKconf[DeleteDuplicates@grpconf];
weights=N[Table[Total@Flatten[weightprov[conf12[[1]],Abs/@Table[grpconf[[i]] . matOktmp[[j]],{i,Length@grpconf}]]],{j,Length@matOktmp}]/.subWW];
positionbest=Flatten@Position[weights,Min[weights]];
matOktmp[[positionbest]](*\[LeftDoubleBracket]-1\[RightDoubleBracket]*)
]


(* ::Subsubsection::Closed:: *)
(*Reassigning l = 2*)


reassignationG22llev2o[o_,o6G4l2_]:=Module[{listdiagr2l,ll,tmpgrpform,changes},
listdiagr2l=o6G4l2[[All,1]];
ll=Length@listdiagr2l;
tmpgrpform=twopt1PIbstc2[0,o][[order\[CapitalGamma]2o[o]]];
changes=Table[bestReparam2ptv02d[o,{tmpgrpform[[listdiagr2l[[i]]]]}],{i,ll}];
Table[{o6G4l2[[j,1]],o6G4l2[[j,2]]/.reparmom2[changes[[j,-1]]]},{j,ll}]
]


reassignationG22llev2of[o_,o6G4l2_]:=Module[{listdiagr2l,ll,tmpgrpform,tmpgrpformf,changes},
listdiagr2l=o6G4l2[[All,1]];
ll=Length@listdiagr2l;
tmpgrpform=twopt1PIbstc2[0,o][[order\[CapitalGamma]2o[o]]];
tmpgrpformf=Transpose[{tmpgrpform[[All,1]],Flatten/@(tmpgrpform[[All,2]])}];
changes=Table[bestReparam2ptv02d[o,{tmpgrpformf[[listdiagr2l[[i]]]]}],{i,ll}];
Table[{o6G4l2[[j,1]],o6G4l2[[j,2]]/.reparmom2[changes[[j,-1]]]},{j,ll}]
]


reassignationG22llev4o[o_,o6G4l2_]:=Module[{listdiagr2l,ll,tmpgrpform,changes},
listdiagr2l=o6G4l2[[All,1]];
ll=Length@listdiagr2l;
tmpgrpform=twopt1PIbstc4[0,o][[order\[CapitalGamma]2o[o]]];
changes=Table[bestReparam2ptv0lev42d[o,{tmpgrpform[[listdiagr2l[[i]]]]}],{i,ll}];
Table[{o6G4l2[[j,1]],o6G4l2[[j,2]]/.reparmom2[changes[[j,-1]]]},{j,ll}]
]


(* ::Subsubsection::Closed:: *)
(*Reassigning l = 3*)


reassignationG23llev2o[o_,o6G4l3_]:=Module[{listdiagr,ll,tmpgrpform,changes},
listdiagr=o6G4l3[[All,1]];
ll=Length@listdiagr;
tmpgrpform=twopt1PIbstc2[0,o][[order\[CapitalGamma]2o[o]]];
changes=Table[bestReparam2ptv0[o,{tmpgrpform[[listdiagr[[i]]]]}],{i,ll}];
Table[{o6G4l3[[j,1]],o6G4l3[[j,2]]/.reparmom3[changes[[j,-1]]]},{j,ll}]
]


reassignationG23llev2oW[o_,o6G4l3_,w_]:=Module[{listdiagr,ll,tmpgrpform,changes},
listdiagr=o6G4l3[[All,1]];
ll=Length@listdiagr;
tmpgrpform=twopt1PIbstc2W[0,o,w][[order\[CapitalGamma]2o[o]]];
changes=Table[bestReparam2ptv0[o,{tmpgrpform[[listdiagr[[i]]]]}],{i,ll}];
Table[{o6G4l3[[j,1]],o6G4l3[[j,2]]/.reparmom3[changes[[j,-1]]]},{j,ll}]
]


reassignationG23llev4o[o_,o6G4l3_]:=Module[{listdiagr,ll,tmpgrpform,changes},
listdiagr=o6G4l3[[All,1]];
ll=Length@listdiagr;
tmpgrpform=twopt1PIbstc4[0,o][[order\[CapitalGamma]2o[o]]];
changes=Table[bestReparam2ptv0lev4[o,{tmpgrpform[[listdiagr[[i]]]]}],{i,ll}];
Table[{o6G4l3[[j,1]],o6G4l3[[j,2]]/.reparmom3[changes[[j,-1]]]},{j,ll}]
]


(* ::Subsection::Closed:: *)
(*Generalization to (\[CapitalGamma]^(2))'*)


(* ::Subsubsection::Closed:: *)
(*Function momenta*)


GtP2dptsdConfig[n3_,n4_,twopt1PIod_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,trinvert,\[Psi]avert,vertex\[Psi]apos,\[Psi]bvert,vertex\[Psi]bpos,squarevert,vertexposquad,propagatorspos,deauxOK,vertexTripos,bubvert,vertexposbub,sunvert,vertexpossuns},
diaglist=twopt1PIod;dl=Length[diaglist];Table[deltas={};daux=choicepextnop@replaceTad\[Tau]\[Beta][diaglist[[lj]][[2]]];pos=DeleteCases[VertexList[daux],0|o1|o2|o3|o4];If[Length[pos]==1,diaglist[[lj]][[1]] suntad^Count[pos,Subscript[\[Tau], _]] bubbl2tr^Count[pos,Subscript[\[Beta], _]],paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux-2]];
paux=Insert[paux,pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o1]&][[1]]][[1,1]]];
paux=Insert[paux,-pext,Position[daux,Select[daux,MemberQ[(#/.UndirectedEdge[a_,b_]:>{a,b}),o2]&][[1]]][[1,1]]];
For[idx=1,idx<=Length[pos],idx++,xlist=Select[daux,MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},pos[[idx]]]&];xpos=DeleteDuplicates[Flatten[Table[Position[daux,xlist[[k]]],{k,1,Length[xlist]}]]];xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length[xlist]}];AppendTo[deltas,\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(Length[xpossign]\)]\(xpossign[\([k, 2]\)]\ paux[\([xpossign[\([k, 1]\)]]\)]\)\)==0];];deltasol=Reduce[deltas]/. And->List/. Equal:>Rule;(*Print[deltasol];*)
paux2=DeleteCases[DeleteDuplicates[Flatten[paux/. deltasol/. Plus->List]/. Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];paux3=paux/. deltasol/. Table[If[paux2=={},Null,paux2[[k]]->qauxify[k]],{k,1,Length[paux2]}];(*Print[paux3];*)
deauxOK=daux[[3;;-1]];
propagatorspos=Select[Range[Length@deauxOK],!MemberQ[deauxOK[[#]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4]&];
trinvert=Cases[pos,Subscript[t, _]];
vertexTripos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},trinvert[[i]]]&)]],{i,Length[trinvert]}];
\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Psi]avert[[i]]]&)]],{i,Length[\[Psi]avert]}];
\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Psi]bvert[[i]]]&)]],{i,Length[\[Psi]bvert]}];
squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},squarevert[[i]]]&)]],{i,Length[squarevert]}];
bubvert=Cases[pos,Subscript[b, _]];
vertexposbub=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},bubvert[[i]]]&)]],{i,Length[bubvert]}];
sunvert=Cases[pos,Subscript[s, _]];
vertexpossuns=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},sunvert[[i]]]&)]],{i,Length[sunvert]}];

{{propagatorspos,vertexposbub,vertexpossuns,vertex\[Psi]apos,vertex\[Psi]bpos,vertexTripos,vertexposquad},
paux3[[3;;-1]]/.{pext->0,-pext->0}}],{lj,1,dl}]]


(* ::Subsubsection::Closed:: *)
(*Procedure 2d*)


bestReparam2dptv02d[o_,grp_]:=Module[{conf12,grpconf,matOktmp,weights,positionbest},
conf12=GtP2dptsdConfig[0,o,grp][[1]];
grpconf=tovector2/@(conf12[[2]]);
matOktmp=matricesTransOKconf2d[DeleteDuplicates@grpconf];
weights=N[Table[Total@Flatten[weightprov[conf12[[1]],Abs/@Table[grpconf[[i]] . matOktmp[[j]],{i,Length@grpconf}]]],{j,Length@matOktmp}]/.subWW];
positionbest=Flatten@Position[weights,Min[weights]];
matOktmp[[positionbest]](*\[LeftDoubleBracket]-1\[RightDoubleBracket]*)
]


(* ::Subsubsection::Closed:: *)
(*Procedure 3d*)


bestReparam2dptv0[o_,grp_]:=Module[{conf12,grpconf,matOktmp,weights,positionbest},
conf12=GtP2dptsdConfig[0,o,grp][[1]];
grpconf=tovector3/@(conf12[[2]]);
matOktmp=matricesTransOKconf[DeleteDuplicates@grpconf];
weights=N[Table[Total@Flatten[weightprov[conf12[[1]],Abs/@Table[grpconf[[i]] . matOktmp[[j]],{i,Length@grpconf}]]],{j,Length@matOktmp}]/.subWW];
positionbest=Flatten@Position[weights,Min[weights]];
matOktmp[[positionbest]](*\[LeftDoubleBracket]-1\[RightDoubleBracket]*)
]


(* ::Subsubsection::Closed:: *)
(*Reassigning l = 2*)


reassignationG2d2llev2o[o_,o6G4l2_]:=Module[{listdiagr2l,ll,tmpgrpform,changes},
listdiagr2l=o6G4l2[[All,1]];
ll=Length@listdiagr2l;
tmpgrpform=twopt1PIbstc2[0,o][[order\[CapitalGamma]2o[o]]];
changes=Table[bestReparam2dptv02d[o,{tmpgrpform[[listdiagr2l[[i]]]]}],{i,ll}];
Table[{o6G4l2[[j,1]],o6G4l2[[j,2]]/.reparmom2[changes[[j,-1]]]},{j,ll}]
]


reassignationG2d2llev2oW[o_,o6G4l2_,w_]:=Module[{listdiagr,ll,tmpgrpform,changes},
listdiagr=o6G4l2[[All,1]];
ll=Length@listdiagr;
tmpgrpform=twopt1PIbstc2W[0,o,w][[order\[CapitalGamma]2o[o]]];
changes=Table[bestReparam2dptv02d[o,{tmpgrpform[[listdiagr[[i]]]]}],{i,ll}];
Table[{o6G4l2[[j,1]],o6G4l2[[j,2]]/.reparmom2[changes[[j,-1]]]},{j,ll}]
]


(* ::Subsubsection::Closed:: *)
(*Reassigning l = 3*)


reassignationG2d3llev2o[o_,o6G4l3_]:=Module[{listdiagr,ll,tmpgrpform,changes},
listdiagr=o6G4l3[[All,1]];
ll=Length@listdiagr;
tmpgrpform=twopt1PIbstc2[0,o][[order\[CapitalGamma]2o[o]]];
changes=Table[bestReparam2dptv0[o,{tmpgrpform[[listdiagr[[i]]]]}],{i,ll}];
Table[{o6G4l3[[j,1]],o6G4l3[[j,2]]/.reparmom3[changes[[j,-1]]]},{j,ll}]
]


reassignationG2d3llev2oW[o_,o6G4l3_,w_]:=Module[{listdiagr,ll,tmpgrpform,changes},
listdiagr=o6G4l3[[All,1]];
ll=Length@listdiagr;
tmpgrpform=twopt1PIbstc2W[0,o,w][[order\[CapitalGamma]2o[o]]];
changes=Table[bestReparam2dptv0[o,{tmpgrpform[[listdiagr[[i]]]]}],{i,ll}];
Table[{o6G4l3[[j,1]],o6G4l3[[j,2]]/.reparmom3[changes[[j,-1]]]},{j,ll}]
]


(* ::Subsection::Closed:: *)
(*Generalization to \[CapitalGamma]^(0)*)


(* ::Subsubsection::Closed:: *)
(*Function momenta*)


(* ::Text:: *)
(*Level 3a (or lower)*)


GtP0ptsdConfig[n3_,n4_,zeropt1PIod_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,trinvert,\[Psi]avert,vertex\[Psi]apos,\[Psi]bvert,vertex\[Psi]bpos,squarevert,vertexposquad,propagatorspos,deauxOK,vertexTripos,bubvert,vertexposbub,sunvert,vertexpossuns},
diaglist=zeropt1PIod;dl=Length[diaglist];Table[deltas={};daux=replaceTad\[Tau]\[Beta][diaglist[[lj]][[2]]];pos=DeleteCases[VertexList[daux],0|o1|o2|o3|o4];If[Length[pos]==1,diaglist[[lj]][[1]] suntad^Count[pos,Subscript[\[Tau], _]] bubbl2tr^Count[pos,Subscript[\[Beta], _]],paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux]];
For[idx=1,idx<=Length[pos],idx++,xlist=Select[daux,MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},pos[[idx]]]&];xpos=DeleteDuplicates[Flatten[Table[Position[daux,xlist[[k]]],{k,1,Length[xlist]}]]];xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length[xlist]}];AppendTo[deltas,\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(Length[xpossign]\)]\(xpossign[\([k, 2]\)]\ paux[\([xpossign[\([k, 1]\)]]\)]\)\)==0];];deltasol=Reduce[deltas]/. And->List/. Equal:>Rule;(*Print[deltasol];*)
paux2=DeleteCases[DeleteDuplicates[Flatten[paux/. deltasol/. Plus->List]/. Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];paux3=paux/. deltasol/. Table[If[paux2=={},Null,paux2[[k]]->qauxify[k]],{k,1,Length[paux2]}];(*Print[paux3];*)
deauxOK=daux;
propagatorspos=Select[Range[Length@deauxOK],!MemberQ[deauxOK[[#]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4]&];
trinvert=Cases[pos,Subscript[t, _]];
vertexTripos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},trinvert[[i]]]&)]],{i,Length[trinvert]}];
\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Psi]avert[[i]]]&)]],{i,Length[\[Psi]avert]}];
\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Psi]bvert[[i]]]&)]],{i,Length[\[Psi]bvert]}];
squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},squarevert[[i]]]&)]],{i,Length[squarevert]}];
bubvert=Cases[pos,Subscript[b, _]];
vertexposbub=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},bubvert[[i]]]&)]],{i,Length[bubvert]}];
sunvert=Cases[pos,Subscript[s, _]];
vertexpossuns=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},sunvert[[i]]]&)]],{i,Length[sunvert]}];

{{propagatorspos,vertexposbub,vertexpossuns,vertex\[Psi]apos,vertex\[Psi]bpos,vertexTripos,vertexposquad},
paux3/.{pext->0,-pext->0}}],{lj,1,dl}]]


(* ::Text:: *)
(*Level 4*)


GtP0ptsdConfiglevel4[n3_,n4_,zeropt1PIod_]:=Module[{pos,daux,diaglist,dl,xlist,xpos,xpossign,paux,paux2,paux3,idx,deltas,deltasol,trinvert,\[Psi]avert,vertex\[Psi]apos,\[Psi]bvert,vertex\[Psi]bpos,squarevert,vertexposquad,propagatorspos,deauxOK,vertexTripos,bubvert,vertexposbub,sunvert,vertexpossuns,vertex\[Kappa]pos,vertexpostheta,vertexpostheta2,vertexpossigma,vertexpossigma2,\[Kappa]vert,thetavert,theta2vert,sigmavert,sigma2vert},
diaglist=zeropt1PIod;dl=Length[diaglist];Table[deltas={};daux=replaceTad\[Tau]\[Beta][diaglist[[lj]][[2]]];pos=DeleteCases[VertexList[daux],0|o1|o2|o3|o4];If[Length[pos]==1,diaglist[[lj]][[1]] suntad^Count[pos,Subscript[\[Tau], _]] bubbl2tr^Count[pos,Subscript[\[Beta], _]],paux={p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18}[[1;;Length@daux]];
For[idx=1,idx<=Length[pos],idx++,xlist=Select[daux,MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},pos[[idx]]]&];xpos=DeleteDuplicates[Flatten[Table[Position[daux,xlist[[k]]],{k,1,Length[xlist]}]]];xpossign=Table[{xpos[[k]],(-1)^Which[MemberQ[{daux[[xpos[[k]]]]},pos[[idx]]\[UndirectedEdge]_],1,MemberQ[{daux[[xpos[[k]]]]},_\[UndirectedEdge]pos[[idx]]],0]},{k,1,Length[xlist]}];AppendTo[deltas,\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(k = 1\), \(Length[xpossign]\)]\(xpossign[\([k, 2]\)]\ paux[\([xpossign[\([k, 1]\)]]\)]\)\)==0];];deltasol=Reduce[deltas]/. And->List/. Equal:>Rule;(*Print[deltasol];*)
paux2=DeleteCases[DeleteDuplicates[Flatten[paux/. deltasol/. Plus->List]/. Times[n_,qq_]:>qq/;NumericQ[n]],0|pext];
paux3=paux/. deltasol/. Table[If[paux2=={},Null,paux2[[k]]->qauxify[k]],{k,1,Length[paux2]}];(*Print[paux3];*)
deauxOK=daux;
propagatorspos=Select[Range[Length@deauxOK],!MemberQ[deauxOK[[#]],Subscript[b, _]|Subscript[s, _]|Subscript[t, _]|Subscript[c, _]|Subscript[\[Psi]a, _]|Subscript[\[Psi]b, _]|o1|o2|o3|o4]&];
trinvert=Cases[pos,Subscript[t, _]];
vertexTripos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},trinvert[[i]]]&)]],{i,Length[trinvert]}];
\[Psi]avert=Cases[pos,Subscript[\[Psi]a, _]];
vertex\[Psi]apos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Psi]avert[[i]]]&)]],{i,Length[\[Psi]avert]}];
\[Psi]bvert=Cases[pos,Subscript[\[Psi]b, _]];
vertex\[Psi]bpos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Psi]bvert[[i]]]&)]],{i,Length[\[Psi]bvert]}];
\[Kappa]vert=Cases[pos,Subscript[\[Kappa], _]];
vertex\[Kappa]pos=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},\[Kappa]vert[[i]]]&)]],{i,Length[\[Kappa]vert]}];
squarevert=Cases[pos,Subscript[c, _]];
vertexposquad=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},squarevert[[i]]]&)]],{i,Length[squarevert]}];
bubvert=Cases[pos,Subscript[b, _]];
vertexposbub=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},bubvert[[i]]]&)]],{i,Length[bubvert]}];
sunvert=Cases[pos,Subscript[s, _]];
vertexpossuns=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},sunvert[[i]]]&)]],{i,Length[sunvert]}];
thetavert=Cases[pos,Subscript[\[Theta], _]];
vertexpostheta=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},thetavert[[i]]]&)]],{i,Length[thetavert]}];
theta2vert=Cases[pos,Subscript[\[Theta]2, _]];
vertexpostheta2=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},theta2vert[[i]]]&)]],{i,Length[theta2vert]}];
sigmavert=Cases[pos,Subscript[\[Sigma], _]];
vertexpossigma=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},sigmavert[[i]]]&)]],{i,Length[sigmavert]}];
sigma2vert=Cases[pos,Subscript[\[Sigma]2, _]];
vertexpossigma2=Table[Flatten[Position[deauxOK,_?(MemberQ[#1/. a_\[UndirectedEdge]b_:>{a,b},sigma2vert[[i]]]&)]],{i,Length[sigma2vert]}];

{{propagatorspos,vertexposbub,vertexpossuns,vertex\[Psi]apos,vertex\[Psi]bpos,vertexTripos,vertexposquad,vertex\[Kappa]pos,vertexpostheta,vertexpostheta2,vertexpossigma,vertexpossigma2},
paux3/. {pext->0,-pext->0}}],{lj,1,dl}]]


(* ::Subsubsection::Closed:: *)
(*Procedure 2d*)


bestReparam0ptv02d[o_,grp_]:=Module[{conf12,grpconf,matOktmp,weights,positionbest},
conf12=GtP0ptsdConfig[0,o,grp][[1]];
grpconf=tovector2/@(conf12[[2]]);
matOktmp=matricesTransOKconf2d[DeleteDuplicates@grpconf];
weights=N[Table[Total@Flatten[weightprov[conf12[[1]],Abs/@Table[grpconf[[i]] . matOktmp[[j]],{i,Length@grpconf}]]],{j,Length@matOktmp}]/.subWW];
positionbest=Flatten@Position[weights,Min[weights]];
matOktmp[[positionbest]](*\[LeftDoubleBracket]-1\[RightDoubleBracket]*)
]


bestReparam0ptv0lev42d[o_,grp_]:=Module[{conf12,grpconf,matOktmp,weights,positionbest},
conf12=GtP0ptsdConfiglevel4[0,o,grp][[1]];
grpconf=tovector2/@(conf12[[2]]);
matOktmp=matricesTransOKconf2d[DeleteDuplicates@grpconf];
weights=N[Table[Total@Flatten[weightprovlev4[conf12[[1]],Abs/@Table[grpconf[[i]] . matOktmp[[j]],{i,Length@grpconf}]]],{j,Length@matOktmp}]/.subWW];
positionbest=Flatten@Position[weights,Min[weights]];
matOktmp[[positionbest]](*\[LeftDoubleBracket]-1\[RightDoubleBracket]*)
]


(* ::Subsubsection::Closed:: *)
(*Procedure 3d*)


bestReparam0ptv0[o_,grp_]:=Module[{conf12,grpconf,matOktmp,weights,positionbest},
conf12=GtP0ptsdConfig[0,o,grp][[1]];
grpconf=tovector3/@(conf12[[2]]);
matOktmp=matricesTransOKconf[DeleteDuplicates@grpconf];
weights=N[Table[Total@Flatten[weightprov[conf12[[1]],Abs/@Table[grpconf[[i]] . matOktmp[[j]],{i,Length@grpconf}]]],{j,Length@matOktmp}]/.subWW];
positionbest=Flatten@Position[weights,Min[weights]];
matOktmp[[positionbest]](*\[LeftDoubleBracket]-1\[RightDoubleBracket]*)
]


bestReparam0ptv0lev4[o_,grp_]:=Module[{conf12,grpconf,matOktmp,weights,positionbest},
conf12=GtP0ptsdConfiglevel4[0,o,grp][[1]];
grpconf=tovector3/@(conf12[[2]]);
matOktmp=matricesTransOKconf[DeleteDuplicates@grpconf];
weights=N[Table[Total@Flatten[weightprov[conf12[[1]],Abs/@Table[grpconf[[i]] . matOktmp[[j]],{i,Length@grpconf}]]],{j,Length@matOktmp}]/.subWW];
positionbest=Flatten@Position[weights,Min[weights]];
matOktmp[[positionbest]](*\[LeftDoubleBracket]-1\[RightDoubleBracket]*)
]


(* ::Subsubsection::Closed:: *)
(*Reassigning l = 2*)


reassignationG02llev2o[o_,o6G4l2_]:=Module[{listdiagr2l,ll,tmpgrpform,changes},
listdiagr2l=o6G4l2[[All,1]];
ll=Length@listdiagr2l;
tmpgrpform=zeropt1PIbstc2[0,o][[order\[CapitalGamma]0o[o]]];
changes=Table[bestReparam0ptv02d[o,{tmpgrpform[[listdiagr2l[[i]]]]}],{i,ll}];
Table[{o6G4l2[[j,1]],o6G4l2[[j,2]]/.reparmom2[changes[[j,-1]]]},{j,ll}]
]


reassignationG02llev2of[o_,o6G4l2_]:=Module[{listdiagr2l,ll,tmpgrpform,tmpgrpformf,changes},
listdiagr2l=o6G4l2[[All,1]];
ll=Length@listdiagr2l;
tmpgrpform=zeropt1PIbstc2[0,o][[order\[CapitalGamma]0o[o]]];
tmpgrpformf=Transpose[{tmpgrpform[[All,1]],Flatten/@(tmpgrpform[[All,2]])}];
changes=Table[bestReparam0ptv02d[o,{tmpgrpformf[[listdiagr2l[[i]]]]}],{i,ll}];
Table[{o6G4l2[[j,1]],o6G4l2[[j,2]]/.reparmom2[changes[[j,-1]]]},{j,ll}]
]


reassignationG02llev4o[o_,o6G4l2_]:=Module[{listdiagr2l,ll,tmpgrpform,changes},
listdiagr2l=o6G4l2[[All,1]];
ll=Length@listdiagr2l;
tmpgrpform=zeropt1PIbstc4[0,o][[order\[CapitalGamma]0o[o]]];
changes=Table[bestReparam0ptv0lev42d[o,{tmpgrpform[[listdiagr2l[[i]]]]}],{i,ll}];
Table[{o6G4l2[[j,1]],o6G4l2[[j,2]]/.reparmom2[changes[[j,-1]]]},{j,ll}]
]


(* ::Subsubsection::Closed:: *)
(*Reassigning l = 3*)


reassignationG03llev2o[o_,o6G4l3_]:=Module[{listdiagr,ll,tmpgrpform,changes},
listdiagr=o6G4l3[[All,1]];
ll=Length@listdiagr;
tmpgrpform=zeropt1PIbstc2[0,o][[order\[CapitalGamma]0o[o]]];
changes=Table[bestReparam0ptv0[o,{tmpgrpform[[listdiagr[[i]]]]}],{i,ll}];
Table[{o6G4l3[[j,1]],o6G4l3[[j,2]]/.reparmom3[changes[[j,-1]]]},{j,ll}]
]


reassignationG03llev2oW[o_,o6G4l3_,w_]:=Module[{listdiagr,ll,tmpgrpform,changes},
listdiagr=o6G4l3[[All,1]];
ll=Length@listdiagr;
tmpgrpform=zeropt1PIbstc2W[0,o,w][[order\[CapitalGamma]0o[o]]];
changes=Table[bestReparam0ptv0[o,{tmpgrpform[[listdiagr[[i]]]]}],{i,ll}];
Table[{o6G4l3[[j,1]],o6G4l3[[j,2]]/.reparmom3[changes[[j,-1]]]},{j,ll}]
]


reassignationG03llev4o[o_,o6G4l3_]:=Module[{listdiagr,ll,tmpgrpform,changes},
listdiagr=o6G4l3[[All,1]];
ll=Length@listdiagr;
tmpgrpform=zeropt1PIbstc4[0,o][[order\[CapitalGamma]0o[o]]];
changes=Table[bestReparam0ptv0lev4[o,{tmpgrpform[[listdiagr[[i]]]]}],{i,ll}];
Table[{o6G4l3[[j,1]],o6G4l3[[j,2]]/.reparmom3[changes[[j,-1]]]},{j,ll}]
]


(* ::Subsection::Closed:: *)
(*l = 4 (Here just permutations, no linear combination)*)


(* ::Subsubsection::Closed:: *)
(*General functions*)


tovector4[x_]:={If[MemberQ[x,q1,All],Coefficient[x,q1],0],If[MemberQ[x,q2,All],Coefficient[x,q2],0],If[MemberQ[x,q3,All],Coefficient[x,q3],0],If[MemberQ[x,q4,All],Coefficient[x,q4],0]}


matricesPermOK4d=Permutations@Permutations@{1,0,0,0};


(* ::Text:: *)
(*Length@matricesPermOK4d*)


(* ::Subsubsection::Closed:: *)
(*Simpler pieces ({gw}) 4d*)


(* ::Input:: *)
(*(*Tuples[Range[0,1],4]*)*)


basicsequnce4d={{0,0,0,0},{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},
{1,1,0,0},{1,0,1,0},{1,0,0,1},{0,1,1,0},{0,1,0,1},{0,0,1,1},{1,1,1,0},{1,1,0,1},{1,0,1,1},{0,1,1,1},{1,1,1,1}};


(* ::Input:: *)
(*(*Length@basicsequnce4d*)*)


(* ::Input:: *)
(*(*Table[{basicsequnce4d[[i]].{q1,q2,q3,q4},{}},{i,Length@basicsequnce4d}]*)*)


(* ::Input:: *)
(*{{0,{0}},{q1,{1}},{q2,{1.3}},{q3,{1.6}},{q4,{1.7}},{q1+q2,{2.4}},{q1+q3,{2.8}},{q1+q4,{2.9}},{q2+q3,{3.5}},{q2+q4,{3.6}},{q3+q4,{4}},{q1+q2+q3,{5}},{q1+q2+q4,{5.1}},{q1+q3+q4,{5.5}},{q2+q3+q4,{5.9}},{q1+q2+q3+q4,{7}}};*)


(* ::Text:: *)
(*Rationalize@%[[All, 2, 1]]*)


(* ::Input:: *)
(*(*gww/@basicsequnce4d*)*)


{gww[{0,0,0,0}],gww[{1,0,0,0}],gww[{0,1,0,0}],gww[{0,0,1,0}],gww[{0,0,0,1}],gww[{1,1,0,0}],gww[{1,0,1,0}],gww[{1,0,0,1}],gww[{0,1,1,0}],gww[{0,1,0,1}],gww[{0,0,1,1}],gww[{1,1,1,0}],gww[{1,1,0,1}],gww[{1,0,1,1}],gww[{0,1,1,1}],gww[{1,1,1,1}]}={0,1,13/10,8/5,17/10,12/5,14/5,29/10,7/2,18/5,4,5,51/10,11/2,59/10,7};


(* ::Subsubsection::Closed:: *)
(*Procedure 4d and reassigning permutation l = 4*)


bestPerm4ptv14d[o_,intgr_]:=Module[{sub1all,sub1allbutProp,momProp,momList,momconf,weights,positionbest},
sub1all={Gt[a_]:> 1,bubbl[a_]:> 1,suns[a_]:> 1,trinfnc[a_,b_,c_]:> 1,squarfnc[a_,b_,c_,d_]:> 1,suntad-> 1,bubbl2tr-> 1,analytic\[Psi]a[a_,b_,c_]:> 1};
sub1allbutProp={bubbl[a_]:> 1,suns[a_]:> 1,trinfnc[a_,b_,c_]:> 1,squarfnc[a_,b_,c_,d_]:> 1,suntad-> 1,bubbl2tr-> 1,analytic\[Psi]a[a_,b_,c_]:> 1};
momProp=(intgr/.sub1allbutProp)/(intgr/.sub1all);
momList=(List@@Expand@momProp)/.{Power->Apply[Sequence]@*ConstantArray}/.Gt[a_]:>a;
(*Print[momList];*)
momconf=tovector4/@momList;
weights=N[Table[Total@Flatten[gw/@Abs/@Table[momconf[[i]] . matricesPermOK4d[[j]],{i,Length@momconf}]],{j,Length@matricesPermOK4d}]/.subWW];
(*Print[weights];*)
positionbest=Flatten@Position[weights,Min[weights]];
matricesPermOK4d[[positionbest]][[1]]
]


bestPermv14d[o_,intgr_]:=Module[{sub1all,sub1allbutProp,momProp,momList,momconf,weights,positionbest},
sub1all={Gt[a_]:> 1,bubbl[a_]:> 1,suns[a_]:> 1,trinfnc[a_,b_,c_]:> 1,squarfnc[a_,b_,c_,d_]:> 1,suntad-> 1,bubbl2tr-> 1,analytic\[Psi]a[a_,b_,c_]:> 1};
sub1allbutProp={bubbl[a_]:> 1,suns[a_]:> 1,trinfnc[a_,b_,c_]:> 1,squarfnc[a_,b_,c_,d_]:> 1,suntad-> 1,bubbl2tr-> 1,analytic\[Psi]a[a_,b_,c_]:> 1};
momProp=(intgr/.sub1allbutProp)/(intgr/.sub1all);
momList=(List@@Expand@momProp)/.{Power->Apply[Sequence]@*ConstantArray}/.Gt[a_]:>a;
(*Print[momList];*)
momconf=tovector4/@momList;
weights=N[Table[Total@Flatten[gw/@Abs/@Table[momconf[[i]] . matricesPermOK4d[[j]],{i,Length@momconf}]],{j,Length@matricesPermOK4d}]/.subWW];
(*Print[weights];*)
positionbest=Flatten@Position[weights,Min[weights]];
matricesPermOK4d[[positionbest]][[1]]
]


reparmom4[m_]:=Module[{nrep},
nrep=m . {q1,q2,q3,q4};
{q1-> nrep[[1]],q2-> nrep[[2]],q3-> nrep[[3]],q4-> nrep[[4]]}
]


reassignation4llev2oW[o_,o6G4l3_,w_]:=Module[{listdiagr,ll,tmpgrpform,changes},
listdiagr=o6G4l3[[All,1]];
ll=Length@listdiagr;
tmpgrpform=GtP4ptslevel2Wass[0,o,w][[order\[CapitalGamma]4o[o]]];
changes=Table[bestPerm4ptv14d[o,tmpgrpform[[listdiagr[[i]]]]],{i,ll}];
Table[{o6G4l3[[j,1]],tmpgrpform[[listdiagr[[j]]]]/.reparmom4[changes[[j]]]},{j,ll}]
]


(* ::Text:: *)
(*(*reassignation4llev2oWOK[o_,o6G4l3_]:=Module[{listdiagr,ll,tmpgrpform,changes},*)
(*listdiagr=o6G4l3[[All,1]];*)
(*ll=Length@o6G4l3;*)
(*changes=Table[bestPerm4ptv14d[o,o6G4l3[[i,2]]],{i,ll}];*)
(*Table[{o6G4l3[[j,1]],o6G4l3[[j,2]]/.reparmom4[changes[[j]]]},{j,ll}]*)
(*]*)*)


reassignation4llev2oWOK[o_,o6G4l3_]:=Module[{listdiagr,ll,tmpgrpform,changes},
listdiagr=o6G4l3[[All,1]];
ll=Length@o6G4l3;
changes=Table[bestPermv14d[o,o6G4l3[[i,2]]],{i,ll}];
Table[{o6G4l3[[j,1]],o6G4l3[[j,2]]/.reparmom4[changes[[j]]]},{j,ll}]
]


(* ::Subsection::Closed:: *)
(*l = 5 (Here just permutations, no linear combination)*)


(* ::Subsubsection::Closed:: *)
(*General functions*)


tovector5[x_]:={If[MemberQ[x,q1,All],Coefficient[x,q1],0],If[MemberQ[x,q2,All],Coefficient[x,q2],0],If[MemberQ[x,q3,All],Coefficient[x,q3],0],If[MemberQ[x,q4,All],Coefficient[x,q4],0],If[MemberQ[x,q5,All],Coefficient[x,q5],0]}


matricesPermOK5d=Permutations@Permutations@{1,0,0,0,0};


(* ::Text:: *)
(*Length@matricesPermOK5d*)


(* ::Subsubsection::Closed:: *)
(*Simpler pieces ({gw}) 5d*)


(* ::Input:: *)
(*(*Tuples[Range[0,1],5]*)*)


(* ::Input:: *)
(*reverse2@(Tuples[Range[0,1],5]);*)


basicsequnce5d={{0,0,0,0,0},{1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1},
{1,1,0,0,0},{1,0,1,0,0},{1,0,0,1,0},{1,0,0,0,1},
{0,1,1,0,0},{0,1,0,1,0},{0,1,0,0,1},{0,0,1,1,0},{0,0,1,0,1},{0,0,0,1,1},
{1,1,1,0,0},{1,1,0,1,0},{1,1,0,0,1},{1,0,1,1,0},{1,0,1,0,1},{1,0,0,1,1},{0,1,1,1,0},{0,1,1,0,1},{0,1,0,1,1},{0,0,1,1,1},
{1,1,1,1,0},{1,1,1,0,1},{1,1,0,1,1},{1,0,1,1,1},{0,1,1,1,1},{1,1,1,1,1}};


(* ::Input:: *)
(*Table[{basicsequnce5d[[i]] . {q1,q2,q3,q4,q5},{}},{i,Length@basicsequnce5d}];*)


(* ::Input:: *)
(*{{0,{0}},{q1,{1}},{q2,{1.3}},{q3,{1.6}},{q4,{1.7}},{q5,{1.8}},{q1+q2,{2.4}},{q1+q3,{2.8}},{q1+q4,{2.9}},{q1+q5,{3}},{q2+q3,{3.5}},{q2+q4,{3.6}},{q2+q5,{3.7}},{q3+q4,{4}},{q3+q5,{4.1}},{q4+q5,{4.2}},{q1+q2+q3,{5}},{q1+q2+q4,{5.1}},{q1+q2+q5,{5.2}},{q1+q3+q4,{5.5}},{q1+q3+q5,{5.6}},{q1+q4+q5,{5.7}},{q2+q3+q4,{6}},{q2+q3+q5,{6.1}},{q2+q4+q5,{6.2}},{q3+q4+q5,{6.3}},{q1+q2+q3+q4,{7}},{q1+q2+q3+q5,{7.1}},{q1+q2+q4+q5,{7.2}},{q1+q3+q4+q5,{7.3}},{q2+q3+q4+q5,{7.4}},{q1+q2+q3+q4+q5,{8.5}}};*)


(* ::Text:: *)
(*Rationalize@%[[All, 2, 1]]*)


(* ::Input:: *)
(*(*gww/@basicsequnce5d*)*)


{gww[{0,0,0,0,0}],gww[{1,0,0,0,0}],gww[{0,1,0,0,0}],gww[{0,0,1,0,0}],gww[{0,0,0,1,0}],gww[{0,0,0,0,1}],gww[{1,1,0,0,0}],gww[{1,0,1,0,0}],gww[{1,0,0,1,0}],gww[{1,0,0,0,1}],gww[{0,1,1,0,0}],gww[{0,1,0,1,0}],gww[{0,1,0,0,1}],gww[{0,0,1,1,0}],gww[{0,0,1,0,1}],gww[{0,0,0,1,1}],gww[{1,1,1,0,0}],gww[{1,1,0,1,0}],gww[{1,1,0,0,1}],gww[{1,0,1,1,0}],gww[{1,0,1,0,1}],gww[{1,0,0,1,1}],gww[{0,1,1,1,0}],gww[{0,1,1,0,1}],gww[{0,1,0,1,1}],gww[{0,0,1,1,1}],gww[{1,1,1,1,0}],gww[{1,1,1,0,1}],gww[{1,1,0,1,1}],gww[{1,0,1,1,1}],gww[{0,1,1,1,1}],gww[{1,1,1,1,1}]}={0,1,13/10,8/5,17/10,9/5,12/5,14/5,29/10,3,7/2,18/5,37/10,4,41/10,21/5,5,51/10,26/5,11/2,28/5,57/10,6,61/10,31/5,63/10,7,71/10,36/5,73/10,37/5,17/2};


(* ::Subsubsection::Closed:: *)
(*Procedure 5d and reassigning permutation l = 5*)


bestPermv15d[o_,intgr_]:=Module[{sub1all,sub1allbutProp,momProp,momList,momconf,weights,positionbest},
sub1all={Gt[a_]:> 1,bubbl[a_]:> 1,suns[a_]:> 1,trinfnc[a_,b_,c_]:> 1,squarfnc[a_,b_,c_,d_]:> 1,suntad-> 1,bubbl2tr-> 1,analytic\[Psi]a[a_,b_,c_]:> 1};
sub1allbutProp={bubbl[a_]:> 1,suns[a_]:> 1,trinfnc[a_,b_,c_]:> 1,squarfnc[a_,b_,c_,d_]:> 1,suntad-> 1,bubbl2tr-> 1,analytic\[Psi]a[a_,b_,c_]:> 1};
momProp=(intgr/.sub1allbutProp)/(intgr/.sub1all);
momList=(List@@Expand@momProp)/.{Power->Apply[Sequence]@*ConstantArray}/.Gt[a_]:>a;
(*Print[momList];*)
momconf=tovector5/@momList;
weights=N[Table[Total@Flatten[gw/@Abs/@Table[momconf[[i]] . matricesPermOK5d[[j]],{i,Length@momconf}]],{j,Length@matricesPermOK5d}]/.subWW];
(*Print[weights];*)
positionbest=Flatten@Position[weights,Min[weights]];
matricesPermOK5d[[positionbest]][[1]]
]


reparmom5[m_]:=Module[{nrep},
nrep=m . {q1,q2,q3,q4,q5};
{q1-> nrep[[1]],q2-> nrep[[2]],q3-> nrep[[3]],q4-> nrep[[4]],q5-> nrep[[5]]}
]


(* ::Input:: *)
(*(*reassignation5llev2oW[o_,o6G4l3_,w_]:=Module[{listdiagr,ll,tmpgrpform,changes},*)
(*listdiagr=o6G4l3[[All,1]];*)
(*ll=Length@listdiagr;*)
(*tmpgrpform=GtP4ptslevel2Wass[0,o,w]\[LeftDoubleBracket]order\[CapitalGamma]4o[o]\[RightDoubleBracket];*)
(*changes=Table[bestPermv15d[o,tmpgrpform\[LeftDoubleBracket]listdiagr\[LeftDoubleBracket]i\[RightDoubleBracket]\[RightDoubleBracket]],{i,ll}];*)
(*Table[{o6G4l3[[j,1]],tmpgrpform\[LeftDoubleBracket]listdiagr\[LeftDoubleBracket]j\[RightDoubleBracket]\[RightDoubleBracket]/.reparmom5[changes\[LeftDoubleBracket]j\[RightDoubleBracket]]},{j,ll}]*)
(*]*)*)


(* ::Input:: *)
(*(*reassignation4llev2oWOK[o_,o6G4l3_]:=Module[{listdiagr,ll,tmpgrpform,changes},*)
(*listdiagr=o6G4l3[[All,1]];*)
(*ll=Length@o6G4l3;*)
(*changes=Table[bestPerm4ptv14d[o,o6G4l3\[LeftDoubleBracket]i,2\[RightDoubleBracket]],{i,ll}];*)
(*Table[{o6G4l3[[j,1]],o6G4l3\[LeftDoubleBracket]j,2\[RightDoubleBracket]/.reparmom4[changes\[LeftDoubleBracket]j\[RightDoubleBracket]]},{j,ll}]*)
(*]*)*)


reassignation5llev2oWOK[o_,o6G4l3_]:=Module[{listdiagr,ll,tmpgrpform,changes},
listdiagr=o6G4l3[[All,1]];
ll=Length@o6G4l3;
changes=Table[bestPermv15d[o,o6G4l3[[i,2]]],{i,ll}];
Table[{o6G4l3[[j,1]],o6G4l3[[j,2]]/.reparmom5[changes[[j]]]},{j,ll}]
]


(* ::Subsection::Closed:: *)
(*Pext & Comparison weights (for \[CapitalGamma]^(2)')*)


Table[tuplesPext[i]=Tuples[Range[-1,1],i],{i,5}];


qqupto5={q1,q2,q3,q4,q5};


shifPext[l_,j_]:=Table[qqupto5[[i]]->qqupto5[[i]]+tuplesPext[l][[j,i]]pext,{i,l}]


coef2[x_]:=Coefficient[x,pext]


meqpext[x_]:=MemberQ[x,pext,Infinity]


subwpext={Gt->wGt,bubbl->wbubbl,suns->wsuns,squarfnc->wsquarfnc,trinfnc->wtrinfnc};


{wGt[0],wGt[1],wbubbl[0],wbubbl[1],wsuns[0],wsuns[1],wtrinfnc[0],wtrinfnc[1],wsquarfnc[0],wsquarfnc[1]}={0,1,0,2/5,0,2/5,0,3/5,0,10};


rassPext[l_,int_]:=Module[{listint,typesint,ll,llor,momint,ltuples,momintrep,momintrepflat,repok,countext,typesintw,weights,positionbest,choice},
listint=(List@@(int)/.{Power->Apply[Sequence]@*ConstantArray});
typesint=DeleteCases[Head/@listint,Symbol|Integer|Rational];
{ll,llor}=Length/@{typesint,listint};
momint=listint[[llor-ll+1;;-1]]/.{Gt[a_]:> {a},bubbl[a_]:> {a},suns[a_]:> {a},squarfnc[a_,b_,c_,d_]:> {a,b,c,d},trinfnc[a_,b_,c_]:> {a,b,c}};
ltuples=Length@tuplesPext[l];
momintrep=Table[momint/.shifPext[l,i],{i,ltuples}];
momintrepflat=Abs[coef2/@Flatten/@momintrep];
repok={};
For[i=1,i<=ltuples,i++,
If[MemberQ[momintrepflat[[i]],2,Infinity]||MemberQ[momintrepflat[[i]],3,Infinity]||MemberQ[momintrepflat[[i]],4,Infinity],Null;,AppendTo[repok,i]]
];
countext=Table[meqpext/@momintrep[[repok]][[i]],{i,Length@repok}]/.{True-> 1,False-> 0};
typesintw=typesint/.subwpext;
weights=Total/@Table[Table[typesintw[[i]]@countext[[j,i]],{i,ll}],{j,Length@repok}];
positionbest=Flatten@Position[weights,Min[weights]];
choice=repok[[positionbest]][[-1]];
int/.shifPext[l,choice]
]


rassPextv2[l_,int_]:=Module[{listint,typesint,ll,llor,momint,ltuples,momintrep,momintrepflat,repok,countext,typesintw,typesintwC,weights,weightsC,weightstot,positionbest,choice},
listint=(List@@(int)/.{Power->Apply[Sequence]@*ConstantArray});
typesint=DeleteCases[Head/@listint,Symbol|Integer|Rational];
{ll,llor}=Length/@{typesint,listint};
momint=listint[[llor-ll+1;;-1]]/.{Gt[a_]:> {a},bubbl[a_]:> {a},suns[a_]:> {a},squarfnc[a_,b_,c_,d_]:> {a,b,c,d},trinfnc[a_,b_,c_]:> {a,b,c}};
ltuples=Length@tuplesPext[l];
momintrep=Table[momint/.shifPext[l,i],{i,ltuples}];
momintrepflat=Abs[coef2/@Flatten/@momintrep];
repok={};
For[i=1,i<=ltuples,i++,
If[MemberQ[momintrepflat[[i]],2,Infinity]||MemberQ[momintrepflat[[i]],3,Infinity]||MemberQ[momintrepflat[[i]],4,Infinity],Null;,AppendTo[repok,i]]
];
countext=Table[meqpext/@momintrep[[repok]][[i]],{i,Length@repok}]/.{True-> 1,False-> 0};
typesintw=typesint/.subwpext;
typesintwC=typesint/.subwCpext;
weights=Total/@Table[Table[typesintw[[i]]@countext[[j,i]],{i,ll}],{j,Length@repok}];
weightsC=Total/@Table[Table[typesintwC[[i]]@countext[[j,i]],{i,ll}],{j,Length@repok}];
weightstot=weights+weightsC/10;
positionbest=Flatten@Position[weightstot,Min[weightstot]];
choice=repok[[positionbest]][[-1]];
int/.shifPext[l,choice]
]


(* ::Text:: *)
(*Comparison weights, select the best set*)


subwWpext={Gt->wWGt,bubbl->wWbubbl,suns->wWsuns,squarfnc->wWsquarfnc,trinfnc->wWtrinfnc};


{wWGt[0],wWGt[1],wWbubbl[0],wWbubbl[1],wWsuns[0],wWsuns[1],wWtrinfnc[0],wWtrinfnc[1],wWsquarfnc[0],wWsquarfnc[1]}={0,1,0,2/5,0,2/5,0,3/5,0,10};


subwCpext={Gt->wCGt,bubbl->wCbubbl,suns->wCsuns,squarfnc->wCsquarfnc,trinfnc->wCtrinfnc};


{wCGt[0],wCGt[1],wCbubbl[0],wCbubbl[1],wCsuns[0],wCsuns[1],wCtrinfnc[0],wCtrinfnc[1],wCsquarfnc[0],wCsquarfnc[1]}={0,1,0,1,0,1,0,2,0,10};


weigforWeight[int_]:=Module[{listint,typesint,ll,llor,momint,ltuples,momintrep,momintrepflat,repok,countext,typesintwW,typesintwC,weightsW, weightsC,positionbest,choice},
listint=(List@@(int)/.{Power->Apply[Sequence]@*ConstantArray});
typesint=DeleteCases[Head/@listint,Symbol|Integer|Rational];
{ll,llor}=Length/@{typesint,listint};
momint=listint[[llor-ll+1;;-1]]/.{Gt[a_]:> {a},bubbl[a_]:> {a},suns[a_]:> {a},squarfnc[a_,b_,c_,d_]:> {a,b,c,d},trinfnc[a_,b_,c_]:> {a,b,c}};
(*Print[typesint,momint];*)
countext=(meqpext/@momint)/.{True-> 1,False-> 0};
typesintwW=typesint/.subwWpext;
typesintwC=typesint/.subwCpext;
weightsW=Total@Table[typesintwW[[i]]@countext[[i]],{i,ll}];
weightsC=Total@Table[typesintwC[[i]]@countext[[i]],{i,ll}];
{Count[typesint,Gt],weightsW,weightsC}
]


weigforWeightbest[int_]:=Module[{listint,typesint,ll,llor,momint,ltuples,momintrep,momintrepflat,repok,countext,typesintwW,typesintwC,weightsW, weightsC,positionbest,choice},
listint=(List@@(int)/.{Power->Apply[Sequence]@*ConstantArray});
typesint=DeleteCases[Head/@listint,Symbol|Integer|Rational];
{ll,llor}=Length/@{typesint,listint};
momint=listint[[llor-ll+1;;-1]]/.{Gt[a_]:> {a},bubbl[a_]:> {a},suns[a_]:> {a},squarfnc[a_,b_,c_,d_]:> {a,b,c,d},trinfnc[a_,b_,c_]:> {a,b,c}};
(*Print[typesint,momint];*)
countext=(meqpext/@momint)/.{True-> 1,False-> 0};
typesintwW=typesint/.subwWpext;
typesintwC=typesint/.subwCpext;
weightsW=Total@Table[typesintwW[[i]]@countext[[i]],{i,ll}];
weightsC=Total@Table[typesintwC[[i]]@countext[[i]],{i,ll}];
{Count[typesint,Gt]/20+weightsW+weightsC/100}
]


bestWeightsG2d[listlistint_]:=Module[{ll,lll,numbers,wei,i,twei,posw,j,poswf,intfor,jj},
ll=Length@listlistint;
lll=Length@(listlistint[[1]]);
(*Print[{ll,lll}];*)
numbers=listlistint[[1,All,1]];
wei={};
For[i=1,i<=ll,i++,AppendTo[wei,Flatten[weigforWeightbest/@listlistint[[i,All,2]]]]];
twei=Transpose[wei];
(*Print[wei,twei];*)
posw={};
For[j=1,j<= lll,j++,AppendTo[posw,Position[twei[[j]],Min[twei[[j]]]][[1]]]];
poswf=Flatten@posw;
intfor={};
(*Print[poswf];*)
For[jj=1,jj<= lll,jj++,AppendTo[intfor,listlistint[[poswf[[jj]],jj]]]];
intfor
]


bestWeightsG2dprintchoice[listlistint_]:=Module[{ll,lll,numbers,wei,i,twei,posw,j,poswf,intfor,jj},
ll=Length@listlistint;
lll=Length@(listlistint[[1]]);
(*Print[{ll,lll}];*)
numbers=listlistint[[1,All,1]];
wei={};
For[i=1,i<=ll,i++,AppendTo[wei,Flatten[weigforWeightbest/@listlistint[[i,All,2]]]]];
twei=Transpose[wei];
(*Print[wei,twei];*)
posw={};
For[j=1,j<= lll,j++,AppendTo[posw,Position[twei[[j]],Min[twei[[j]]]][[1]]]];
poswf=Flatten@posw;
intfor={};
(*Print[poswf];*)
For[jj=1,jj<= lll,jj++,AppendTo[intfor,listlistint[[poswf[[jj]],jj]]]];
intfor
]


(* ::Section::Closed:: *)
(*Derivation*)


(* ::Subsection::Closed:: *)
(*Functions*)


(* ::Text:: *)
(*Substitutions from derivatives of the blocks to their values and the scalar products*)


derPropbubbl={
	Derivative[2][Gt][a_]:>1/3 (Gt[a]^2-4Gt[a]^3),
	Derivative[2][bubbl][a_]:>-(1/(6\[Pi]))(1/(4+vec3dsq[a])^2),
	Derivative[1][Gt][a_]*Derivative[1][Gt][b_]:>2vec3dprod[a,b] (Gt[a]^2 Gt[b]^2)/3,
	Derivative[1][Gt][a_]^n_/;(n>= 2&&IntegerQ[n/2]):>2vec3dsq[a] Gt[a]^(2n)/3,
	Derivative[1][Gt][a_]*Derivative[1][bubbl][b_]:>-(vec3dprod[a,b]/Sqrt[vec3dsq[b]]) * (Gt[a]^2 bublder1[b])/3,
	Derivative[1][bubbl][a_]*Derivative[1][bubbl][b_]:>vec3dprod[a,b]/(Sqrt[vec3dsq[a]]*Sqrt[vec3dsq[b]]) * (bublder1[a]bublder1[b])/6,
	Derivative[1][bubbl][a_]^n_/;(n>= 2&&IntegerQ[n/2]):>bublder1[a]^n/6,
	Derivative[1][Gt][a_]*Derivative[1][suns][b_]:>-(vec3dprod[a,b]/Sqrt[vec3dsq[b]]) * Gt[a]^2/3 sunsder1[b],
	Derivative[2][suns][a_]:>-(1/(96 \[Pi]^2))(1/(9+vec3dsq[a])),
	Derivative[1][bubbl][a_]*Derivative[1][suns][b_]:>vec3dprod[a,b]/(Sqrt[vec3dsq[a]]*Sqrt[vec3dsq[b]]) * bublder1[a]/6 * sunsder1[b],
	Derivative[1][suns][a_]*Derivative[1][suns][b_]:>vec3dprod[a,b]/(Sqrt[vec3dsq[a]]*Sqrt[vec3dsq[b]]) * (sunsder1[a]*sunsder1[b])/6};


derTriangl1={
	Derivative[1,0,0][trinfnc][a_,b_,c_]*Derivative[1][Gt][d_]:> -2vec3dprod[a,d] (tder1[a,b,c]Gt[d]^2)/3,
	Derivative[0,1,0][trinfnc][a_,b_,c_]*Derivative[1][Gt][d_]:>-2 vec3dprod[b,d] (tder1[b,c,a]Gt[d]^2)/3,
	Derivative[0,0,1][trinfnc][a_,b_,c_]*Derivative[1][Gt][d_]:> -2vec3dprod[c,d] (tder1[c,a,b]Gt[d]^2)/3,
	Derivative[1,0,0][trinfnc][a_,b_,c_]*Derivative[1][bubbl][d_]:> vec3dprod[a,d]/Sqrt[vec3dsq[d]] (tder1[a,b,c]bublder1[d])/3,
	Derivative[0,1,0][trinfnc][a_,b_,c_]*Derivative[1][bubbl][d_]:> vec3dprod[b,d]/Sqrt[vec3dsq[d]] (tder1[b,c,a]bublder1[d])/3,
	Derivative[0,0,1][trinfnc][a_,b_,c_]*Derivative[1][bubbl][d_]:>vec3dprod[c,d]/Sqrt[vec3dsq[d]] (tder1[c,a,b]bublder1[d])/3,
	Derivative[1,0,0][trinfnc][a_,b_,c_]*Derivative[1][suns][d_]:> vec3dprod[a,d] tder1[a,b,c]/(3Sqrt[vec3dsq[d]]) sunsder1[d],
	Derivative[0,1,0][trinfnc][a_,b_,c_]*Derivative[1][suns][d_]:> vec3dprod[b,d] tder1[b,c,a]/(3Sqrt[vec3dsq[d]]) sunsder1[d],
	Derivative[0,0,1][trinfnc][a_,b_,c_]*Derivative[1][suns][d_]:> vec3dprod[c,d] tder1[c,a,b]/(3Sqrt[vec3dsq[d]]) sunsder1[d]};


derTriangl2={
	Derivative[2,0,0][trinfnc][a_,b_,c_]:> tder1[a,b,c]+(2vec3dsq[a] tder2[a,b,c])/3,
	Derivative[0,2,0][trinfnc][a_,b_,c_]:>tder1[b,c,a]+(2vec3dsq[b] tder2[b,c,a])/3,
	Derivative[0,0,2][trinfnc][a_,b_,c_]:>tder1[c,a,b]+ (2vec3dsq[c] tder2[c,a,b])/3,
	Derivative[1,1,0][trinfnc][a_,b_,c_]:>( 2 vec3dprod[a,b] tder11[a,b,c])/3,
	Derivative[1,0,1][trinfnc][a_,b_,c_]:>( 2 vec3dprod[a,c] tder11[c,a,b])/3,
	Derivative[0,1,1][trinfnc][a_,b_,c_]:>( 2 vec3dprod[b,c] tder11[b,c,a])/3,
	Derivative[1,0,0][trinfnc][a_,b_,c_]*Derivative[1,0,0][trinfnc][d_,e_,f_]:> 2vec3dprod[a,d] (tder1[a,b,c]tder1[d,e,f])/3,
	Derivative[1,0,0][trinfnc][a_,b_,c_]*Derivative[0,1,0][trinfnc][d_,e_,f_]:> 2vec3dprod[a,e] (tder1[a,b,c]tder1[e,f,d])/3,
	Derivative[1,0,0][trinfnc][a_,b_,c_]*Derivative[0,0,1][trinfnc][d_,e_,f_]:> 2vec3dprod[a,f] (tder1[a,b,c]tder1[f,d,e])/3,
	Derivative[0,1,0][trinfnc][a_,b_,c_]*Derivative[0,1,0][trinfnc][d_,e_,f_]:> 2vec3dprod[b,e] (tder1[b,c,a]tder1[e,f,d])/3,
	Derivative[0,1,0][trinfnc][a_,b_,c_]*Derivative[0,0,1][trinfnc][d_,e_,f_]:> 2vec3dprod[b,f] (tder1[b,c,a]tder1[f,d,e])/3,
	Derivative[0,0,1][trinfnc][a_,b_,c_]*Derivative[0,0,1][trinfnc][d_,e_,f_]:> 2vec3dprod[c,f] (tder1[c,a,b]tder1[f,d,e])/3};


(* ::Text:: *)
(*deriveintegrandv13: from function to derive respect to pext to the double derivative*)


deriveintegrandv13[exp_]:=Module[{integrand,symbolic,subsymbolic,subsymbolic2,result},
integrand=exp/.Gt[a_]/;MemberQ[a,-pext]:>Gt[-a]/. bubbl[a_]/;MemberQ[a,-pext]:> bubbl[-a];
symbolic=D[integrand,{pext,2}];
subsymbolic=Expand[symbolic]/.derPropbubbl; subsymbolic2=subsymbolic/.derTriangl1/.derTriangl2;
result=subsymbolic2/.pext->0/. Subscript[0, n_]:>0
]


deriveintegrandv13r[exp_]:=Refine@deriveintegrandv13@exp


deriveintegrandv13s[exp_]:=Simplify@deriveintegrandv13@exp


deriveintegrandv13sc[exp_]:=Simplify[deriveintegrandv13@exp,TimeConstraint->5]


bublder1[a_]:=(1/(\[Pi] Sqrt[vec3dsq[a]] (8+2vec3dsq[a]))-bubbl[a]/Sqrt[vec3dsq[a]])


sunsder1[x_]:=-(1/(32 \[Pi]^2))(2/(Sqrt[vec3dsq[x]] (1+vec3dsq[x]/9))+(2Sqrt[vec3dsq[x]])/(9 (1+vec3dsq[x]/9))-(6 ArcTan[Sqrt[vec3dsq[x]]/3])/vec3dsq[x])


(* ::Section::Closed:: *)
(*Substitutions of the analytical blocks (from Integrands with q-momenta to Integrands ready-to-be-integrated)*)


(* ::Subsection::Closed:: *)
(*Vector and their products*)


(* ::Text:: *)
(*Scalar product generalization*)


vec2d[x_]:=Module[{p=Expand@x,pli},If[Head@p===Symbol,p/.z_:>{Subscript[z, 0],Subscript[z, 1]},If[Head@p===Times&&NumericQ[p/.Times[num_,vec_]:>num],p/.Times[num_,vec_]:>Times[num,vec2d[vec]],If[Head@p===Plus,pli=List@@p;\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(n = 1\), \(Length@pli\)]\(vec2d[pli[\([n]\)]]\)\)]]]]
vec2dsq[p_]:=Module[{arg,tmp,tmpsq},arg=Piecewise[{{{p},Head@p===Symbol||Head@p==Times},{List@@Expand[p],Head@Expand[p]===Plus}}];If[Length@arg>1,tmp=Plus@@(vec2d/@(Expand[p]/.Plus->List));tmpsq=tmp . tmp,If[Length@arg==1,tmpsq=vec2d[p] . vec2d[p]]];Simplify@Expand[tmpsq/.{Subscript[a_, 0]:>Sqrt[Subscript[a, 2]]Cos[Subscript[a, 3]],Subscript[a_, 1]:>Sqrt[Subscript[a, 2]]Sin[Subscript[a, 3]]}]]


vec2dsqcart[p_]:=Module[{arg,tmp,tmpsq},arg=Piecewise[{{{p},Head@p===Symbol||Head@p==Times},{List@@Expand[p],Head@Expand[p]===Plus}}];If[Length@arg>1,tmp=Plus@@(vec2d/@(Expand[p]/.Plus->List));tmpsq=tmp . tmp,If[Length@arg==1,tmpsq=vec2d[p] . vec2d[p]]];tmpsq]


vec3d[x_]:=Module[{p=Expand@x,pli},If[Head@p===Symbol,p/.z_:>{Subscript[z, 0],Subscript[z, 1],Subscript[z, 2]},If[Head@p===Times&&NumericQ[p/.Times[num_,vec_]:>num],p/.Times[num_,vec_]:>Times[num,vec3d[vec]],If[Head@p===Plus,pli=List@@p;\!\(
\*SubsuperscriptBox[\(\[Sum]\), \(n = 1\), \(Length@pli\)]\(vec3d[pli[\([n]\)]]\)\)]]]]


vec3dsq[p_]:=Module[{arg,tmp,tmpsq},arg=Piecewise[{{{p},Head@p===Symbol||Head@p==Times},{List@@Expand[p],Head@Expand[p]===Plus}}];If[Length@arg>1,tmp=Plus@@(vec3d/@(Expand[p]/.Plus->List));tmpsq=tmp . tmp,If[Length@arg==1,tmpsq=vec3d[p] . vec3d[p]]];If[p===0,0,Simplify@Expand[tmpsq/.{Subscript[a_, 0]:>Subscript[a, 3]Sin[Subscript[a, 4]]Cos[Subscript[a, 5]],Subscript[a_, 1]:>Subscript[a, 3]Sin[Subscript[a, 4]]Sin[Subscript[a, 5]],Subscript[a_, 2]:>Subscript[a, 3]Cos[Subscript[a, 4]]}]]]


vec3dnorm[p_]:=Module[{arg,tmp,tmpsq},arg=Piecewise[{{{p},Head@p===Symbol||Head@p==Times},{List@@Expand[p],Head@Expand[p]===Plus}}];If[Length@arg>1,tmp=Plus@@(vec3d/@(Expand[p]/.Plus->List));tmpsq=Sqrt[tmp . tmp],If[Length@arg==1,tmpsq=Sqrt[vec3d[p] . vec3d[p]]]];If[p===0,0,Simplify@Expand[tmpsq/.{Subscript[a_, 0]:>Subscript[a, 3]Sin[Subscript[a, 4]]Cos[Subscript[a, 5]],Subscript[a_, 1]:>Subscript[a, 3]Sin[Subscript[a, 4]]Sin[Subscript[a, 5]],Subscript[a_, 2]:>Subscript[a, 3]Cos[Subscript[a, 4]]}]]]


vec3dsqcart[p_]:=Module[{arg,tmp,tmpsq},arg=Piecewise[{{{p},Head@p===Symbol||Head@p==Times},{List@@Expand[p],Head@Expand[p]===Plus}}];If[Length@arg>1,tmp=Plus@@(vec3d/@(Expand[p]/.Plus->List));tmpsq=tmp . tmp,If[Length@arg==1,tmpsq=vec3d[p] . vec3d[p]]];tmpsq]


vec3dsqcart[p_]:=Module[{arg,tmp,tmpsq},arg=Piecewise[{{{p},Head@p===Symbol||Head@p==Times},{List@@Expand[p],Head@Expand[p]===Plus}}];If[Length@arg>1,tmp=Plus@@(vec3d/@(Expand[p]/.Plus->List));tmpsq=tmp . tmp,If[Length@arg==1,tmpsq=vec3d[p] . vec3d[p]]];tmpsq]


vec3dprod[p_,q_]:=Module[{argp,tmpp,tmpsq,argq,tmpq},argp=Piecewise[{{{p},Head@p===Symbol||Head@p==Times},{List@@Expand[p],Head@Expand[p]===Plus}}];argq=Piecewise[{{{q},Head@q===Symbol||Head@q==Times},{List@@Expand[q],Head@Expand[q]===Plus}}];
If[Length@argp>1,tmpp=Plus@@(vec3d/@(Expand[p]/.Plus->List)),If[Length@argp==1,tmpp=vec3d[p]]];
If[Length@argq>1,tmpq=Plus@@(vec3d/@(Expand[q]/.Plus->List)),If[Length@argq==1,tmpq=vec3d[q]]];
tmpsq=tmpp . tmpq;If[p===0||q===0,0,Simplify@Expand[tmpsq/.{Subscript[a_, 0]:>Subscript[a, 3]Sin[Subscript[a, 4]]Cos[Subscript[a, 5]],Subscript[a_, 1]:>Subscript[a, 3]Sin[Subscript[a, 4]]Sin[Subscript[a, 5]],Subscript[a_, 2]:>Subscript[a, 3]Cos[Subscript[a, 4]]}]]]


(* ::Input:: *)
(**)


modmine[x_]:=Sqrt[\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i\), \(3\)]
\*SuperscriptBox[\(x[\([i]\)]\), \(2\)]\)]


modmineS[x_]:=Simplify@Sqrt[\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i\), \(3\)]
\*SuperscriptBox[\(x[\([i]\)]\), \(2\)]\)]


vec3dcross[p_,q_]:=Module[{argp,tmpp,tmpsq,argq,tmpq},argp=Piecewise[{{{p},Head@p===Symbol||Head@p==Times},{List@@Expand[p],Head@Expand[p]===Plus}}];argq=Piecewise[{{{q},Head@q===Symbol||Head@q==Times},{List@@Expand[q],Head@Expand[q]===Plus}}];
If[Length@argp>1,tmpp=Plus@@(vec3d/@(Expand[p]/.Plus->List)),If[Length@argp==1,tmpp=vec3d[p]]];
If[Length@argq>1,tmpq=Plus@@(vec3d/@(Expand[q]/.Plus->List)),If[Length@argq==1,tmpq=vec3d[q]]];
tmpsq=tmpp\[Cross]tmpq;
If[p===0||q===0,0,Simplify@Expand[tmpsq/.{Subscript[a_, 0]:>Subscript[a, 3]Sin[Subscript[a, 4]]Cos[Subscript[a, 5]],Subscript[a_, 1]:>Subscript[a, 3]Sin[Subscript[a, 4]]Sin[Subscript[a, 5]],Subscript[a_, 2]:>Subscript[a, 3]Cos[Subscript[a, 4]]}]]]


vec3dscalcros[aa_,bb_,cc_]:=Module[{},
Simplify[(vec3d[aa]/.{Subscript[a_, 0]:>Subscript[a, 3]Sin[Subscript[a, 4]]Cos[Subscript[a, 5]],Subscript[a_, 1]:>Subscript[a, 3]Sin[Subscript[a, 4]]Sin[Subscript[a, 5]],Subscript[a_, 2]:>Subscript[a, 3]Cos[Subscript[a, 4]]}) . vec3dcross[bb,cc]]
]


(* ::Text:: *)
(*Example*)


(* ::Input:: *)
(*(*Refine[squarfnc[q3,-q1,q2,q1-q2-q3]/.squarfnc[a_,b_,c_,d_]\[RuleDelayed]  sq@@(Join[vec3dnorm/@{a,b,c,d,a+b,b+c},{vec3dprod[a,a+b]/(vec3dnorm[a]vec3dnorm[a+b]),vec3dprod[c,b+c]/(vec3dnorm[c]vec3dnorm[b+c]),(vec3dcross[a,b]).(vec3dcross[b,c])/((modmine@vec3dcross[a,b])(modmine@vec3dcross[b,c])),vec3dprod[b,a+b],vec3dprod[b,c],vec3dprod[a,b],vec3dprod[b,c+b],(vec3dscalcros[a,b,c])^2}]),assum]/. Subscript[pext1, n_]\[RuleDelayed]0/. Subscript[pext2, n_]\[RuleDelayed]0/. Subscript[pext3, n_]\[RuleDelayed]0/. Subscript[0, n_]\[RuleDelayed]0/.Subscript[q1, 4]\[Rule]0/.Subscript[q1, 5]\[Rule]0/.Subscript[q2, 5]\[Rule]0*)*)


(* ::Subsection::Closed:: *)
(*Substitutions, external momenta identically 0 and analytics (from Integrands with q-momenta to Integrands ready-to-be-integrated)*)


(* ::Subsubsection::Closed:: *)
(*Momenta at zero*)


pextat0={pext->0, pext1->0, pext2->0, pext3->0};


complTad\[Tau]={suntad->-(Log[4/3]/(128 \[Pi]^3)(*+(cren/(8\[Pi]))*))};


complTad\[Beta]={bubbl2tr->1/(12288 \[Pi]^2)};


bub0={bubbl[a_]/;(a===0)->1/(8\[Pi])};


trinfnc0={ trinfnc[a_,b_,c_]/;(a===0&&b=!=0):> 1/(8\[Pi](vec3dsq[b]+4)), trinfnc[a_,b_,c_]/;(b===0&&a=!=0):> 1/(8\[Pi](vec3dsq[a]+4)),trinfnc[a_,b_,c_]/;(c===0&&a=!=0):> 1/(8\[Pi](vec3dsq[a]+4)),
trinfnc[a_,b_,c_]/;(a===b===c===0):> 1/(32\[Pi])};


trinfncD0={ tder1[a_,b_,c_]/;(a===0&&b=!=0):>(-6-vec3dsq[b])/(96 \[Pi] (4+vec3dsq[b])^2) , tder1[a_,b_,c_]/;(b===0&&a=!=0):>-(1/(16 \[Pi] (4+vec3dsq[a])^2)),tder1[a_,b_,c_]/;(c===0&&a=!=0):>-(1/(16 \[Pi] (4+vec3dsq[a])^2)), tder1[a_,b_,c_]/;(a===b===c===0):>-(1/(256 \[Pi])),
 tder2[a_,b_,c_]/;(a===0&&b=!=0):>(80+30 vec3dsq[b]+3(vec3dsq[b])^2)/(960 \[Pi] (4+vec3dsq[b])^3) , tder2[a_,b_,c_]/;(b===0&&a=!=0):>1/(12 \[Pi] (4+vec3dsq[a])^3), tder2[a_,b_,c_]/;(c===0&&a=!=0):>1/(12 \[Pi] (4+vec3dsq[a])^3), tder2[a_,b_,c_]/;(a===b===c===0):>1/(768 \[Pi]),
 tder11[a_,b_,c_]/;(a===0&&b=!=0):>(8+vec3dsq[b])/(192 \[Pi] (4+vec3dsq[b])^3) , tder11[a_,b_,c_]/;(b===0&&a=!=0):>(8+vec3dsq[a])/(192 \[Pi] (4+vec3dsq[a])^3), tder11[a_,b_,c_]/;(c===0&&a=!=0):>1/(24 \[Pi] (4+vec3dsq[b])^3), tder11[a_,b_,c_]/;(a===b===c===0):>1/(1536 \[Pi])};


squarfnc0={squarfnc[a_,b_,c_,d_]/;(c===0&&a=!=0&&b=!=0&&d=!=0):>sq0mom[Sqrt[vec3dsq[d]],Sqrt[vec3dsq[a]],Sqrt[vec3dsq[b]]],squarfnc[a_,b_,c_,d_]/;(b===0&&a=!=0&&c=!=0&&d=!=0):>sq0mom[Sqrt[vec3dsq[c]],Sqrt[vec3dsq[d]],Sqrt[vec3dsq[a]]],squarfnc[a_,b_,c_,d_]/;(a===0&&c=!=0&&b=!=0&&d=!=0):>sq0mom[Sqrt[vec3dsq[b]],Sqrt[vec3dsq[c]],Sqrt[vec3dsq[d]]],squarfnc[a_,b_,c_,d_]/;(d===0&&a=!=0&&b=!=0&&c=!=0):>sq0mom[Sqrt[vec3dsq[a]],Sqrt[vec3dsq[b]],Sqrt[vec3dsq[c]]],

squarfnc[a_,b_,c_,d_]/;(a===-c&&a=!=0&&b=!=0&&c=!=0&&d=!=0):>sqmom1212[Sqrt[vec3dsq[a]],Sqrt[vec3dsq[b]],Sqrt[vec3dsq[a+b]]],squarfnc[a_,b_,c_,d_]/;(a===-b&&a=!=0&&b=!=0&&c=!=0&&d=!=0):>sqmom1122[Sqrt[vec3dsq[a]],Sqrt[vec3dsq[c]],Sqrt[vec3dsq[a+c]]],
squarfnc[a_,b_,c_,d_]/;(c===-b&&a=!=0&&b=!=0&&c=!=0&&d=!=0):>sqmom1122[Sqrt[vec3dsq[b]],Sqrt[vec3dsq[d]],Sqrt[vec3dsq[b+d]]],

squarfnc[a_,b_,c_,d_]/;(a===b===0&&c=!=0):> (vec3dsq[c]+8)/(32\[Pi] (vec3dsq[c]+4)^2),squarfnc[a_,b_,c_,d_]/;(b===c===0&&a=!=0):> (vec3dsq[a]+8)/(32\[Pi] (vec3dsq[a]+4)^2),squarfnc[a_,b_,c_,d_]/;(a===c===0&&b=!=0):> 1/(4\[Pi] (vec3dsq[b]+4)^2),squarfnc[a_,b_,c_,d_]/;(a===d===0&&b=!=0):> (vec3dsq[b]+8)/(32\[Pi] (vec3dsq[b]+4)^2),squarfnc[a_,b_,c_,d_]/;(b===d===0&&c=!=0):>1/(4\[Pi] (vec3dsq[c]+4)^2),squarfnc[a_,b_,c_,d_]/;(d===c===0&&a=!=0):> (vec3dsq[a]+8)/(32\[Pi] (vec3dsq[a]+4)^2),
squarfnc[a_,b_,c_,d_]/;(a===b===c===d===0):> 1/(64\[Pi])};


\[Psi]afuncsub={\[Psi]afnc[a_,b_,c_]/;(a===0&&b=!=0):>ArcTan[Sqrt[vec3dsq[b]]/3]/(32 \[Pi]^2 Sqrt[vec3dsq[b]]) ,\[Psi]afnc[a_,b_,c_]/;(b===0&&a=!=0):>ArcTan[Sqrt[vec3dsq[a]]/3]/(32 \[Pi]^2 Sqrt[vec3dsq[a]]),\[Psi]afnc[a_,b_,c_]/;(c===0&&a=!=0):>ArcTan[Sqrt[vec3dsq[a]]/3]/(32 \[Pi]^2 Sqrt[vec3dsq[a]]),
\[Psi]afnc[a_,b_,c_]/;(a===b===c===0):>1/(96 \[Pi]^2)};


momExt0SimplAnalytic[int_]:=int/.pextat0/.complTad\[Tau]/. complTad\[Beta]/.bub0/.trinfncD0/.trinfnc0/.squarfnc0/.\[Psi]afuncsub


(* ::Subsubsection::Closed:: *)
(*Propagator and bubbles (and sunsets)*)


PropBubbl[int_]:=int/.Gt[p_]:> 1/(vec3dsq[p]+1)/.bubbl[x_]:> ArcTan[ Sqrt[vec3dsq[x]]/2]/(4\[Pi] Sqrt[vec3dsq[x]])


PropBubblSun[int_]:=int/.Gt[p_]:> 1/(vec3dsq[p]+1)/.bubbl[x_]:> ArcTan[ Sqrt[vec3dsq[x]]/2]/(4\[Pi] Sqrt[vec3dsq[x]])/.suns[x_]/;x=!=0:>-(1/(32 \[Pi]^2)) ((6 ArcTan[Sqrt[vec3dsq[x]]/3])/Sqrt[vec3dsq[x]]+Log[1+vec3dsq[x]/9]-2)/.suns[x_]/;x===0:>0


(* ::Subsubsection::Closed:: *)
(*Numerical substitutions*)


\[Psi]bfuncsubNUM={\[Psi]bfnc[a_,b_,c_]/;(a===0&&b=!=0):>\[Psi]bfuncNum[Sqrt[vec3dsq[b]]],\[Psi]bfnc[a_,b_,c_]/;(b===0&&a=!=0):>\[Psi]bfuncNum[Sqrt[vec3dsq[a]]],\[Psi]bfnc[a_,b_,c_]/;(c===0&&a=!=0):>\[Psi]bfuncNum[Sqrt[vec3dsq[a]]],
\[Psi]bfnc[a_,b_,c_]/;(a===b===c===0):>1/(96 \[Pi]^2),\[Kappa]fnc[a_,b_,c_,d_]/;(a===b===c===d===0):>1/(576 \[Pi]^2)};


funcsubNUM={\[Psi]bfnc[a_,b_,c_]/;(a===0&&b=!=0):>\[Psi]bfuncNum[Sqrt[vec3dsq[b]]],\[Psi]bfnc[a_,b_,c_]/;(b===0&&a=!=0):>\[Psi]bfuncNum[Sqrt[vec3dsq[a]]],\[Psi]bfnc[a_,b_,c_]/;(c===0&&a=!=0):>\[Psi]bfuncNum[Sqrt[vec3dsq[a]]],
\[Psi]bfnc[a_,b_,c_]/;(a===b===c===0):>1/(96 \[Pi]^2),\[Kappa]fnc[a_,b_,c_,d_]/;(b===d===0&&a==-c&&a=!=0):>\[Kappa]fncNum[Sqrt[vec3dsq[a]]],\[Kappa]fnc[a_,b_,c_,d_]/;(a===c===0&&b==-d&&b=!=0):>\[Kappa]fncNum[Sqrt[vec3dsq[b]]],\[Kappa]fnc[a_,b_,c_,d_]/;(a===b===c===d===0):>1/(576 \[Pi]^2),t2b[a_]:>t2bNum[Sqrt[vec3dsq[a]]],c1s[a_]:>c1sNum[Sqrt[vec3dsq[a]]],c1s1b[a_]:>c1s1bNum[Sqrt[vec3dsq[a]]],c3b[a_]:>c3bNum[Sqrt[vec3dsq[a]]]};


PropNum0[int_]:=int/.{t2b[x_]/;x===0:>t2b0,c3b[x_]/;x===0:>c3b0,c1s1b[x_]/;x===0:>c1s1b0,c1s[x_]/;x===0:>c1s0}


(* ::Subsubsection::Closed:: *)
(*Numerical values at momentum zero*)


(* ::Input:: *)
(*t2b0=(*N[NIntegrate[Simplify[ArcTan[p/2]^2/(32 (1+p^2) \[Pi]^4)],{p,0,\[Infinity]},WorkingPrecision\[Rule]70,AccuracyGoal\[Rule]30,PrecisionGoal\[Rule]30],30];*)0.00025873486451138085253348322719701228818103419772821970883946957897545422779`50.*)


(* ::Input:: *)
(*c3b0=(*N[NIntegrate[Simplify[ArcTan[x/2]^3/(128 \[Pi]^5 (x+x^3))],{x,0,\[Infinity]},WorkingPrecision\[Rule]70,AccuracyGoal\[Rule]30,PrecisionGoal\[Rule]30],30];*)5.17740191981638141989610945325097609798863848392241187897266200052629111`50.*^-6*)


(* ::Input:: *)
(*c1s1b0=(*N[NIntegrate[Simplify[-((ArcTan[x/2] (6 ArcTan[x/3]+x (-2+Log[1/9 (9+x^2)])))/(256 \[Pi]^5 (1+x^2)^2))],{x,0,\[Infinity]},WorkingPrecision\[Rule]70,AccuracyGoal\[Rule]30,PrecisionGoal\[Rule]30],30];*)-7.4301230310515429595160772331202367411613205717899609039235647713635093`50.*^-7*)


(* ::Input:: *)
(*c1s0=(1-8 Log[2]+Log[81])/(2048 \[Pi]^3);*)


(* ::Input:: *)
(*\[Kappa]0000=1/(576 \[Pi]^2);*)


(* ::Subsubsection::Closed:: *)
(*Triangles for Mathematica*)


dettrian123[k1_,k2_,k3_]:=(k1+k2-k3) (k1-k2+k3) (-k1+k2+k3) (k1+k2+k3)


trian123[k1_,k2_,k3_]:=1/(4\[Pi]) ArcTan[Sqrt[k1^2 k2^2 k3^2+dettrian123[k1,k2,k3]]/(8+ k1^2+k2^2+k3^2)]/Sqrt[k1^2 k2^2 k3^2+dettrian123[k1,k2,k3]]


PropTrian\[Psi]bMath[int_]:=int/.trinfnc[p_,q_,r_]:> t123[Sqrt[vec3dsq[p]],Sqrt[vec3dsq[q]],Sqrt[vec3dsq[r]]]/.\[Psi]bfuncsubNUM


PropTrianNumMath[int_]:=int/.trinfnc[p_,q_,r_]:> t123[Sqrt[vec3dsq[p]],Sqrt[vec3dsq[q]],Sqrt[vec3dsq[r]]]/.funcsubNUM


(* ::Input:: *)
(*(*subTrianDer={tder1[a_,b_,c_]\[RuleDelayed]1/2trider1@@(Sqrt/@vec3dsq/@{a,b,c}),tder2[a_,b_,c_]\[RuleDelayed]1/4trider2@@(Sqrt/@vec3dsq/@{a,b,c}),tder11[a_,b_,c_]\[RuleDelayed]1/4trider11@@(Sqrt/@vec3dsq/@{a,b,c})};*)*)


subTrianDer={tder1[a_,b_,c_]:>1/2 triDer1@@(Sqrt/@vec3dsq/@{a,b,c}),tder2[a_,b_,c_]:>1/4 triDer2@@(Sqrt/@vec3dsq/@{a,b,c}),tder11[a_,b_,c_]:>1/4 triDer11@@(Sqrt/@vec3dsq/@{a,b,c})};


PropTrianDer[int_]:=int/.subTrianDer


(* ::Input:: *)
(*(*functionder2[intg_]:=intg/.tder1[a_,b_,c_]\[RuleDelayed]1/4trider1@@(Sqrt/@vec3dsq/@{a,b,c})/.tder2[a_,b_,c_]\[RuleDelayed]1/4trider2@@(Sqrt/@vec3dsq/@{a,b,c})/.tder11[a_,b_,c_]\[RuleDelayed]1/4trider11@@(Sqrt/@vec3dsq/@{a,b,c})*)*)


explicitfunc={t123[p_,q_,r_]:> trian123[p,q,r],
sq0mom[p_,q_,r_]:> sq0mom123[p,q,r],sqmom1212[p_,q_,r_]:> squarmom1212[p,q,r],sqmom1122[p_,q_,r_]:> squarmom1122[p,q,r],
triDer1[p_,q_,r_]:> trider1[p,q,r],triDer2[p_,q_,r_]:> trider2[p,q,r],triDer11[p_,q_,r_]:> trider11[p,q,r]};


(* ::Input:: *)
(*(*PropTrianSq0Mathsub[int_]:=int/.t123[p_,q_,r_]\[RuleDelayed] trian123[p,q,r]/.sq0mom[p_,q_,r_]\[RuleDelayed] sq0mom123[p,q,r]/.sqmom1212[p_,q_,r_]\[RuleDelayed] squarmom1212[p,q,r]/.sqmom1122[p_,q_,r_]\[RuleDelayed] squarmom1122[p,q,r]*)*)


PropTrianSq0Mathsub[int_]:=int/.explicitfunc


(* ::Subsubsection::Closed:: *)
(*Soft limits of squares for Mathematica*)


sq0mom123[k1_,k2_,k3_]:=(k2^2 (k1^2-k2^2+k3^2)trian123[k1,k2,k3]+(-2 k1^4+2 k2^2 k3^2-2 k3^4+k1^2 (4 k3^2+k2^2 (2+k3^2)))/(4 \[Pi](4+k1^2) (4+k3^2)))/(2(k1^2 k2^2 k3^2+dettrian123[k1,k2,k3]))


squarmom1212[k1s_,k2s_,kxs_]:=(-(ArcTan[Sqrt[k1s^4 (-1+2 k2s^2)-(k2s^2-kxs^2)^2+k1s^2 (2 k2s^4+2 kxs^2-k2s^2 (-2+kxs^2))]/(8+3 k1s^2+3 k2s^2-kxs^2)]/Sqrt[k1s^4 (-1+2 k2s^2)-(k2s^2-kxs^2)^2+k1s^2 (2 k2s^4+2 kxs^2-k2s^2 (-2+kxs^2))])+4\[Pi] trian123[k1s,k2s,kxs])/(2(k1s^2+k2s^2-kxs^2) \[Pi])


squarmom1122[k1s_,k3s_,kys_]:=1/2 ((2 k1s^4+2 k3s^2 (k3s^2-kys^2)-k1s^2 (2 kys^2+k3s^2 (4+kys^2)))/(4\[Pi] (4+k1s^2) (4+k3s^2) )-kys^2 (k1s^2+k3s^2-kys^2) trian123[k1s,k3s,kys])/(k1s^4+(k3s^2-kys^2)^2-k1s^2 (2 kys^2+k3s^2 (2+kys^2)))


(* ::Input:: *)
(*(*FullSimplify[sq0mom123[k1,k2,k3]-squarmom1122[k1,k3,k2]]*)*)


(* ::Subsubsection::Closed:: *)
(*Derivatives of triangles for Mathematica*)


trider1[k1_,k2_,k3_]:=(( k1^2-k2^2-k3^2)/(4\[Pi](4+k1^2))-(2 k1^2-2 k3^2-k2^2 (2+k3^2))trian123[k1,k2,k3])/(k1^4+(k2^2-k3^2)^2-k1^2 (2 k3^2+k2^2 (2+k3^2)))


trider11[k1_,k2_,k3_]:=1/(4\[Pi] ((k1^2 k2^2 k3^2+dettrian123[k1,k2,k3])^2) ) (((-k1^2+k2^2+k3^2)  (2 (k2^2-k3^2)-k1^2 (2+k3^2)))/(4+k1^2)+((k1^2-k2^2+k3^2)(2 k1^2-2 k3^2-k2^2 (2+k3^2)))/(4+k2^2)+1/((4+k1^2) (4+k2^2) (8+k1^2+k2^2+k3^2) ) (2 k1^8-k1^6 k2^2 (8+k3^2)+6 k1^4 (8 k3^2+4 k2^2 k3^2+k2^4 (2+k3^2))+k1^2 (24 k2^4 k3^2-32 k3^2 (-4+k3^2)-k2^6 (8+k3^2)+k2^2 k3^2 (96-8 k3^2+k3^4))+2 (k2^8+24 k2^4 k3^2-16 k2^2 k3^2 (-4+k3^2)-k3^4 (64+8 k3^2+k3^4)))-3 (-2 k2^2+2 k3^2+k1^2 (2+k3^2)) (2 k1^2-2 k3^2-k2^2 (2+k3^2))4\[Pi] trian123[k1,k2,k3])-(2 (2+k3^2)trian123[k1,k2,k3]-(2((k1^2-k2^2)^2-k3^4))/(4\[Pi](4+k1^2) (4+k2^2)(8+k1^2+k2^2+k3^2)))/(k1^2 k2^2 k3^2+dettrian123[k1,k2,k3])


trider2[k1_,k2_,k3_]:=1/(4 \[Pi] ((k1^2 k2^2 k3^2+dettrian123[k1,k2,k3])^2) ) ((4 k1^2-2(2 k3^2+k2^2 (2+k3^2)))((-k1^2+k2^2+k3^2)/(4+k1^2)+(2 k1^2-2 k3^2-k2^2 (2+k3^2))4\[Pi] trian123[k1,k2,k3])-((k1^2-k2^2-k3^2) (2 k1^2-2 k3^2-k2^2 (2+k3^2)))/(4+k1^2)+(-2 k1^2+2 k3^2+k2^2 (2+k3^2))^2 4\[Pi] trian123[k1,k2,k3])-1/(4 \[Pi](k1^2 k2^2 k3^2+dettrian123[k1,k2,k3]) ) (2/(4+k1^2)+(2 (-k1^2+k2^2+k3^2))/(4+k1^2)^2-16\[Pi] trian123[k1,k2,k3])


(* ::Subsection::Closed:: *)
(*Assumptions and simplifications*)


(*qaux={q1,q2,q3,q4,q5,q6,q7,q8,q9};*)


assum={Subscript[q1, 3]>0&&Subscript[q2, 3]>0&&Subscript[q3, 3]>0&&Subscript[q4, 3]>0&&Subscript[q5, 3]>0&&Subscript[q6, 3]>0&&Subscript[q7, 3]>0&&Subscript[q8, 3]>0};


simplroots[intgrand_]:=Refine[intgrand/. Subscript[0, n_]:>0/.Subscript[q1, 4]->0/.Subscript[q1, 5]->0/.Subscript[q2, 5]->0,assum]


simplrootsS[intgrand_]:=Simplify[intgrand/. Subscript[0, n_]:>0/.Subscript[q1, 4]->0/.Subscript[q1, 5]->0/.Subscript[q2, 5]->0,assum]


(* ::Section::Closed:: *)
(*Composition*)


(* ::Subsection::Closed:: *)
(*Composition low l for conversion*)


prepMathl1[xx_]:=(\!\(
\*SubsuperscriptBox[\(q1\), \(3\), \(2\)]/\((2\ \((\[Pi]^2)\)\ )\)\) xx)/.Subscript[q1, 3]->x


prepMathl2[xx_]:=(xx (\!\(
\*SubsuperscriptBox[\(q1\), \(3\), \(2\)] 
\*SubsuperscriptBox[\(q2\), \(3\), \(2\)]\))/(8 \[Pi]^4) Sin[Subscript[q2, 4]])/.{Subscript[q2, 4]-> \[Theta],Subscript[q1, 3]-> x ,Subscript[q2, 3]-> y}


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(*fromqtoreadyStep1[int_]:=simplroots@PropTrianNumMath@PropBubblSun@momExt0SimplAnalytic@int*)*)


(* ::Input:: *)
(*(*fromqtoreadyStep1a[int_]:=simplroots@PropTrian\[Psi]bMath@PropBubblSun@momExt0SimplAnalytic@int*)*)


fromqtoreadyStep1[int_]:=simplroots@PropTrianDer@PropTrianNumMath@PropNum0@PropBubblSun@momExt0SimplAnalytic@int


fromqtoreadyStep2[int_]:=PropTrianSq0Mathsub@int


(* ::Input:: *)
(**)


(* ::Input:: *)
(*(*fromqtoReadyl0[int_]:=PropNum0@fromqtoreadyStep1a@int*)*)


(* ::Input:: *)
(*(*fromqtoReadyl1[int_]:=prepMathl1@PropNum0@fromqtoreadyStep1@int*)*)


fromqtoReadyl0[int_]:=fromqtoreadyStep1@int


fromqtoReadyl1[int_]:=prepMathl1@fromqtoreadyStep1@int


fromqtoMiddlel2[int_]:=prepMathl2@fromqtoreadyStep1@int


fromqtoReadyl2[int_]:=Simplify@PropTrianSq0Mathsub@prepMathl2@fromqtoreadyStep1@int


fromqtoReadyl2ref[int_]:=Refine@PropTrianSq0Mathsub@prepMathl2@fromqtoreadyStep1@int


(* ::Subsection::Closed:: *)
(*Main functions for interface*)


olFunc[o_,l_,distribution_,objectdistibuted_]:=Transpose[{distribution[[1+l]],objectdistibuted[[distribution[[1+l]]]]}]


applylistdouble[list2_,func_]:=Transpose[{list2[[All,1]],func/@(list2[[All,2]])}]


multiplylistdouble[list2_,fact_]:=Transpose[{list2[[All,1]],fact*(list2[[All,2]])}]


olFuncgiven[o_,objectdistibuted_,selected_]:=Transpose[{selected,objectdistibuted[[selected]]}]


deleterepljoin[original_,bettornot_,suppl_]:=Sort@Join[Select[original,And@@Table[#[[1]]=!=Flatten[bettornot][[i]],{i,Length[Flatten[bettornot]]}]&],suppl]


deletereplacecombine[original_,bettornot_,suppl_]:=Sort@Join[Select[original,And@@Table[#[[1]]=!=Flatten[bettornot][[i]],{i,Length[Flatten[bettornot]]}]&],Select[suppl,Or@@Table[#[[1]]==Flatten[bettornot][[i]],{i,Length[Flatten[bettornot]]}]&]]


(* ::Subsection::Closed:: *)
(*Dealing with numerical vertices*)


posNum[x_]:=Module[{thosewithnum},
thosewithnum=Select[x,MemberQ[#,Subscript[\[Psi]b, _]|Subscript[\[Kappa], _]|Subscript[\[Sigma], _]|Subscript[\[Sigma]2, _]|Subscript[\[Theta], _]|Subscript[\[Theta]2, _]|Subscript[\[Gamma], _]|Subscript[\[Gamma]\[Kappa], _]|Subscript[\[Gamma]c, _],Infinity]&];
Flatten[Table[Position[x,thosewithnum[[i]]],{i,Length@thosewithnum}],1]
]


pos\[Psi]a[x_]:=Module[{thosewithnum},
thosewithnum=Select[x,MemberQ[#,Subscript[\[Psi]a, _],Infinity]&];
Flatten[Table[Position[x,thosewithnum[[i]]],{i,Length@thosewithnum}],1]
]


(* ::Input:: *)
(*(*Transpose[{l\[CapitalGamma]4level4poslo[o]\[LeftDoubleBracket]1+l\[RightDoubleBracket],(fourpt1PIbstc4[0,o]\[LeftDoubleBracket]order\[CapitalGamma]4o[o]\[RightDoubleBracket]\[LeftDoubleBracket]l\[CapitalGamma]4level4poslo[o]\[LeftDoubleBracket]1+l\[RightDoubleBracket]\[RightDoubleBracket])}];*)
(*o6G4l2Edges=Table[{%\[LeftDoubleBracket]i,1\[RightDoubleBracket],%\[LeftDoubleBracket]i,2,2\[RightDoubleBracket]},{i,Length@%}]*)*)


(* ::Input:: *)
(*(*posnum@o6G4l2Edges*)*)


countwithVert[edges_,vert_]:=If[Length@vert==0,Count[edges,x_/;MemberQ[x,Subscript[vert, _],Infinity]],Count[edges,x_/;MemberQ[x,Alternatives@@Table[Subscript[vert[[i]], _],{i,Length@vert}],Infinity]]]


distrNumVert[edges_,vertlist_]:=Join[{{Length@edges}},Table[countwithVert[edges,vertlist[[i]]],{i,Length@vertlist}]]


(* ::Input:: *)
(*(*distrNumVert[o6G4l2Edges,{b,c,\[Psi]b,\[Kappa],{\[Sigma],\[Sigma]2,\[Theta], \[Theta]2,\[Gamma],\[Gamma]\[Kappa],\[Gamma]c}}]*)*)


(* ::Input:: *)
(*(*l=1;Transpose[{l\[CapitalGamma]4level4poslo[o]\[LeftDoubleBracket]1+l\[RightDoubleBracket],(fourpt1PIbstc4[0,o]\[LeftDoubleBracket]order\[CapitalGamma]4o[o]\[RightDoubleBracket]\[LeftDoubleBracket]l\[CapitalGamma]4level4poslo[o]\[LeftDoubleBracket]1+l\[RightDoubleBracket]\[RightDoubleBracket])}];*)
(*edgesG4l1o[o]=Table[{%\[LeftDoubleBracket]i,1\[RightDoubleBracket],%\[LeftDoubleBracket]i,2,2\[RightDoubleBracket]},{i,Length@%}]*)*)


(* ::Text:: *)
(*Given the list of diagrams with numerical pieces (provided by posNum) the function comparisonG4Tabl2v3 creates two list of those that without numerical substitutions have the same number of effective loop (so to be replace, with substitutions at level 2) and those with one or more so to bee computed with the tabulations*)


comparisonlevG4ol[o_,l_,edge_]:=Module[{pos2l,posOK,loopsub2,loopsub4,comapar,notbet,bet},
pos2l=(posNum@edge)[[All,1]];
posOK=intG4ol[o,l][[pos2l,1]];
loopsub2=l\[CapitalGamma]4level2o[o];
loopsub4=l\[CapitalGamma]4level4o[o];
comapar=Table[{loopsub2[[posOK[[i]]]],loopsub4[[posOK[[i]]]]},{i,Length@posOK}];
{notbet,bet}={{},{}};

For[j=1,j<=Length@comapar,j++,If[comapar[[j]]== {l,l},AppendTo[notbet,posOK[[j]]],AppendTo[bet,posOK[[j]]]]];
{notbet,bet}
]


comparisonlev\[Psi]aG4ol[o_,l_,edge_]:=Module[{pos2l,posOK,posus,loopsub2,loopsub4,comapar,notbet,bet},
pos2l=(pos\[Psi]a@edge)[[All,1]];
posOK=intG4ol[o,l][[pos2l,1]];
(*posus=order\[CapitalGamma]2o[o]\[LeftDoubleBracket]posOK\[RightDoubleBracket];*)
loopsub2=l\[CapitalGamma]4level2o[o];
loopsub4=l\[CapitalGamma]4level4o[o];
comapar=Table[{loopsub2[[posOK[[i]]]],loopsub4[[posOK[[i]]]]},{i,Length@posOK}];
{notbet,bet}={{},{}};

For[j=1,j<=Length@comapar,j++,If[comapar[[j]]== {l,l},AppendTo[notbet,posOK[[j]]],AppendTo[bet,posOK[[j]]]]];
{notbet,bet}
]


comparisonlevG2ol[o_,l_,edge_]:=Module[{pos2l,posOK,loopsub2,loopsub4,comapar,notbet,bet},
pos2l=(posNum@edge)[[All,1]];
posOK=intG2ol[o,l][[pos2l,1]];
loopsub2=l\[CapitalGamma]2level2o[o];
loopsub4=l\[CapitalGamma]2level4o[o];
comapar=Table[{loopsub2[[posOK[[i]]]],loopsub4[[posOK[[i]]]]},{i,Length@posOK}];
{notbet,bet}={{},{}};

For[j=1,j<=Length@comapar,j++,If[comapar[[j]]== {l,l},AppendTo[notbet,posOK[[j]]],AppendTo[bet,posOK[[j]]]]];
{notbet,bet}
]


comparisonlev\[Psi]aG2ol[o_,l_,edge_]:=Module[{pos2l,posOK,posus,loopsub2,loopsub4,comapar,notbet,bet},
pos2l=(pos\[Psi]a@edge)[[All,1]];
posOK=intG2ol[o,l][[pos2l,1]];
(*posus=order\[CapitalGamma]2o[o]\[LeftDoubleBracket]posOK\[RightDoubleBracket];*)
loopsub2=l\[CapitalGamma]2level2o[o];
loopsub4=l\[CapitalGamma]2level4o[o];
comapar=Table[{loopsub2[[posOK[[i]]]],loopsub4[[posOK[[i]]]]},{i,Length@posOK}];
{notbet,bet}={{},{}};

For[j=1,j<=Length@comapar,j++,If[comapar[[j]]== {l,l},AppendTo[notbet,posOK[[j]]],AppendTo[bet,posOK[[j]]]]];
{notbet,bet}
]


(* ::Section::Closed:: *)
(*Numerical pieces: Tabulations (Just comment)*)


(* ::Text:: *)
(*Import of the tabulations*)


(* ::Input:: *)
(*pathTabulations="SubTAB/subdiagram_tables";*)


(* ::Subsection::Closed:: *)
(*Functions for tabulation*)


(* ::Subsubsection::Closed:: *)
(*Conformal mapping*)


(* ::Input:: *)
(*conformalmap[a_]:=(4 #)/(a (1-#)^2)&;*)


(* ::Input:: *)
(*conformalmapinv[a_]:=(Sqrt[1+a #]-1)/(Sqrt[1+a #]+1)&;*)


(* ::Subsubsection::Closed:: *)
(*Asymptotic behavior*)


(* ::Input:: *)
(*asmodel[a_,b_][q_]:=a/q^2+b/q^3*)


(* ::Input:: *)
(*asympF[f_,OptionsPattern[{prec->100,q0->10^6}]]:=Module[{q0val=OptionValue[q0],aslist,bfit,asmodelb,asfit,coeffasymp,precval=OptionValue[prec],F},*)
(*F[q_]:=NIntegrate[f[q],{p,0,\[Infinity]},WorkingPrecision->10*precval,PrecisionGoal->precval,MaxRecursion->100];*)
(*aslist=ParallelTable[{q0val*10^n,F[q0val*10^n]},{n,0,4,1/4}];*)
(*bfit=Solve[asmodel[a,b][q]==aslist[[1,2]]/. q->aslist[[1,1]],b];*)
(*asmodelb[a_][q2_]:=asmodel[a,b][q2]/.bfit;*)
(**)
(*asfit=NonlinearModelFit[aslist[[2;;-1]],asmodelb[a][q],a,q];*)
(*If[1-asfit["AdjustedRSquared"]>0.01,Print["CHECK Asymptotic Behavior: fitted \!\(\*SuperscriptBox[\(R\), \(2\)]\) = ",asfit["AdjustedRSquared"]];*)
(*Print[aslist]];*)
(*coeffasymp=({a,b}/.bfit/.asfit["BestFitParameters"])[[1]]*)
(*]*)


(* ::Text:: *)
(*Better fit*)


(* ::Input:: *)
(*getexponent[xylist_]:=Mean@Table[-Log[xylist[[k+1,2]]/xylist[[k,2]]]/Log[xylist[[k+1,1]]/xylist[[k,1]]],{k,1,Length@xylist-1}]*)


(* ::Input:: *)
(*asympfit[xylist_]:=Module[{\[Alpha],c},\[Alpha]=getexponent[xylist];c=xylist[[1,2]]*xylist[[1,1]]^\[Alpha];{\[Alpha],c}]*)


(* ::Input:: *)
(*fasymp[parlist_][q_]:=parlist[[2]]/q^parlist[[1]]*)


(* ::Input:: *)
(*fasymp01[parlist_,OptionsPattern[{map->conformalmap[1]}]][x_]:=fasymp[parlist][OptionValue[map][x]]*)


(* ::Input:: *)
(*(*tabulateF01asymp[f_,OptionsPattern[{prec\[Rule]50,npt\[Rule]1000,map\[Rule]conformalmap[1]}]]:=Module[{eps=1/OptionValue[npt],c0p,Fx0,F,xpts,fpts,ymp,precval=OptionValue[prec],rmap=OptionValue[map]},*)
(*c0p=Simplify[SeriesCoefficient[f[q],{q,0,0}],p>0];Fx0=SetPrecision[NIntegrate[c0p,{p,0,\[Infinity]},WorkingPrecision\[Rule]5*precval,PrecisionGoal\[Rule]precval,MaxRecursion\[Rule]100],precval];*)
(*xpts=Join[Table[x,{x,eps,1-eps,eps}],Table[x,{x,1-4eps/5,1-eps/5,eps/5}]];*)
(*F[x_]:=SetPrecision[NIntegrate[f[rmap[x]],{p,0,\[Infinity]},WorkingPrecision\[Rule]10*precval,PrecisionGoal\[Rule]precval,MaxRecursion\[Rule]200],precval];*)
(*fpts=ParallelTable[{xpts\[LeftDoubleBracket]i\[RightDoubleBracket],F[xpts\[LeftDoubleBracket]i\[RightDoubleBracket]]},{i,1,Length@xpts}];*)
(*ymp=asympfit[MapAt[rmap,fpts\[LeftDoubleBracket]-5;;-2\[RightDoubleBracket],{All,1}]];*)
(*{ToString[rmap[x],InputForm],Join[{{0,Fx0}},fpts,{{1,0}}],ymp}*)
(*]*)*)


(* ::Input:: *)
(*(*tabulateF01asymp2[f_,OptionsPattern[{prec\[Rule]50,npt\[Rule]1000,map\[Rule]conformalmap[1]}]]:=Module[{eps=1/OptionValue[npt],c0p,Fx0,F,xpts,fpts,ymp,precval=OptionValue[prec],rmap=OptionValue[map]},*)
(*c0p=Limit[f[p,\[Theta]][q],q->0,Assumptions\[Rule](p>0&&0<\[Theta]<\[Pi])];*)
(*Fx0=Integrate[c0p,{p,0,\[Infinity]},{\[Theta],0,\[Pi]}];*)
(*xpts=Join[Table[x,{x,eps,1-eps,eps}],Table[x,{x,1-4eps/5,1-eps/5,eps/5}]];*)
(*F[x_]:=SetPrecision[NIntegrate[f[p,\[Theta]][rmap[x]],{p,0,\[Infinity]},{\[Theta],0,\[Pi]},WorkingPrecision\[Rule]3*precval,PrecisionGoal\[Rule]precval,MaxRecursion\[Rule]200,Method->{"GlobalAdaptive","SingularityHandler"->"DuffyCoordinates","MaxErrorIncreases"->3000,"SingularityDepth"->1}],precval];*)
(*fpts=ParallelTable[{xpts\[LeftDoubleBracket]i\[RightDoubleBracket],F[xpts\[LeftDoubleBracket]i\[RightDoubleBracket]]},{i,1,Length@xpts}];*)
(*ymp=asympfit[MapAt[rmap,fpts\[LeftDoubleBracket]-5;;-2\[RightDoubleBracket],{All,1}]];*)
(*{ToString[rmap[x],InputForm],Join[{{0,Fx0}},fpts,{{1,0}}],ymp}*)
(*]*)*)


(* ::Subsubsection::Closed:: *)
(*Measure*)


(* ::Input:: *)
(*meas[x_]=TB2map[x]^2 TB2map'[x]//Simplify*)


(* ::Subsection::Closed:: *)
(*TB1)  Triangle bubble*)


(* ::Input:: *)
(*GraphPlot[{Subscript[o, 1]\[UndirectedEdge]x,x\[UndirectedEdge]y,z\[UndirectedEdge]y,y\[UndirectedEdge]Subscript[o, 2],x\[UndirectedEdge]z,z \[UndirectedEdge]y}]*)


(* ::Subsubsection::Closed:: *)
(*Interpolation functions*)


(* ::Input:: *)
(*ClearAll[TB1,TB1map,TB1c]*)
(*TB1=Import[pathTabulations<>"/TB1.dat"];*)
(*TB1map[x_]=ToExpression@TB1[[1,1]];*)
(*TB1c=ToExpression@TB1[[2;;-1]];*)


(* ::Input:: *)
(*ClearAll[TB1interp,TB1approx]*)
(*TB1interp=Interpolation[TB1c,Method->"Hermite",InterpolationOrder->10];*)
(*TB1approx[x_?NumericQ/;x<=TB1c[[-5,1]]]:=TB1interp[x]*)
(*TB1approx[x_?NumericQ/;x>TB1c[[-5,1]]]:=fasymp01[asympfit[MapAt[TB1map,TB1c[[-5;;-2]],{All,1}]],map->TB1map][x];*)


(* ::Subsection::Closed:: *)
(*TB2) Triangle bubble bubble*)


(* ::Input:: *)
(*GraphPlot[{Subscript[o, 1]\[UndirectedEdge]x,x\[UndirectedEdge]y,z\[UndirectedEdge]y,x\[UndirectedEdge]z,x\[UndirectedEdge]z,z \[UndirectedEdge]y,y\[UndirectedEdge]Subscript[o, 2]}]*)


(* ::Subsubsection::Closed:: *)
(*Interpolation functions*)


(* ::Input:: *)
(*ClearAll[TB2,TB2map,TB2c]*)
(*TB2=Import[pathTabulations<>"/TB2.dat"];*)
(*TB2map[x_]=ToExpression@TB2[[1,1]];*)
(*TB2c=ToExpression@TB2[[2;;-1]];*)


(* ::Input:: *)
(*ClearAll[TB2interp,TB2approx]*)
(*TB2interp=Interpolation[TB2c,Method->"Hermite",InterpolationOrder->10];*)
(*TB2approx[x_?NumericQ/;x<=TB2c[[-5,1]]]:=TB2interp[x]*)
(*TB2approx[x_?NumericQ/;x>TB2c[[-5,1]]]:=fasymp01[asympfit[MapAt[TB2map,TB2c[[-5;;-2]],{All,1}]],map->TB2map][x];*)


(* ::Subsection::Closed:: *)
(*QSB0) Boxed sunset*)


(* ::Input:: *)
(*GraphPlot[{Subscript[o, 1]\[UndirectedEdge]x,y\[UndirectedEdge]Subscript[o, 2],x\[UndirectedEdge]y,x\[UndirectedEdge]z,z \[UndirectedEdge]w,z \[UndirectedEdge]w,z \[UndirectedEdge]w,w\[UndirectedEdge]y}]*)


(* ::Subsubsection::Closed:: *)
(*Interpolation functions*)


(* ::Input:: *)
(*ClearAll[QSB0,QSB0map,QSB0c]*)
(*QSB0=Import[pathTabulations<>"/QSB0.dat"];*)
(*QSB0map[x_]=ToExpression@QSB0[[1,1]];*)
(*QSB0c=ToExpression@QSB0[[2;;-1]];*)


(* ::Input:: *)
(*ClearAll[QSB0interp,QSB0approx]*)
(*QSB0interp=Interpolation[QSB0c,Method->"Hermite",InterpolationOrder->10];*)
(*QSB0approx[x_?NumericQ/;x<=QSB0c[[-5,1]]]:=QSB0interp[x]*)
(*QSB0approx[x_?NumericQ/;x>QSB0c[[-5,1]]]:=fasymp01[asympfit[MapAt[QSB0map,QSB0c[[-5;;-2]],{All,1}]],map->QSB0map][x];*)


(* ::Subsection::Closed:: *)
(*QSB1) Boxed sunset bubble*)


(* ::Input:: *)
(*GraphPlot[{Subscript[o, 1]\[UndirectedEdge]x,y\[UndirectedEdge]Subscript[o, 2],x\[UndirectedEdge]y,x\[UndirectedEdge]y,x\[UndirectedEdge]z,z \[UndirectedEdge]w,z \[UndirectedEdge]w,z \[UndirectedEdge]w,w\[UndirectedEdge]y}]*)


(* ::Subsubsection::Closed:: *)
(*Interpolation functions*)


(* ::Input:: *)
(*ClearAll[QSB1,QSB1map,QSB1c]*)
(*QSB1=Import[pathTabulations<>"/QSB1.dat"];*)
(*QSB1map[x_]=ToExpression@QSB1[[1,1]];*)
(*QSB1c=ToExpression@QSB1[[2;;-1]];*)


(* ::Input:: *)
(*ClearAll[QSB1interp,QSB1approx]*)
(*QSB1interp=Interpolation[QSB1c,Method->"Hermite",InterpolationOrder->10];*)
(*QSB1approx[x_?NumericQ/;x<=QSB1c[[-5,1]]]:=QSB1interp[x]*)
(*QSB1approx[x_?NumericQ/;x>QSB1c[[-5,1]]]:=fasymp01[asympfit[MapAt[QSB1map,QSB1c[[-5;;-2]],{All,1}]],map->QSB1map][x];*)


(* ::Subsection::Closed:: *)
(*QB3) Boxed 3 bubbles*)


(* ::Input:: *)
(*GraphPlot[{Subscript[o, 1]\[UndirectedEdge]x,y\[UndirectedEdge]Subscript[o, 2],x\[UndirectedEdge]y,x\[UndirectedEdge]z,x\[UndirectedEdge]z,z \[UndirectedEdge]w,z \[UndirectedEdge]w,w\[UndirectedEdge]y,w\[UndirectedEdge]y}]*)


(* ::Subsubsection::Closed:: *)
(*Interpolation functions*)


(* ::Input:: *)
(*ClearAll[QB3,QB3map,QB3c]*)
(*QB3=Import[pathTabulations<>"/QB3.dat"];*)
(*QB3map[x_]=ToExpression@QB3[[1,1]];*)
(*QB3c=ToExpression@QB3[[2;;-1]];*)


(* ::Input:: *)
(*ClearAll[QB3interp,QB3approx]*)
(*QB3interp=Interpolation[QB3c,Method->"Hermite",InterpolationOrder->10];*)
(*QB3approx[x_?NumericQ/;x<=QB3c[[-5,1]]]:=QB3interp[x]*)
(*QB3approx[x_?NumericQ/;x>QB3c[[-5,1]]]:=fasymp01[asympfit[MapAt[QB3map,QB3c[[-5;;-2]],{All,1}]],map->QB3map][x];*)


(* ::Subsection::Closed:: *)
(*TT) Triangle triangle (Kite)*)


(* ::Input:: *)
(*GraphPlot[{Subscript[o, 1]\[UndirectedEdge]x,x\[UndirectedEdge]y,x\[UndirectedEdge]z,y\[UndirectedEdge]z,y\[UndirectedEdge]w,z\[UndirectedEdge]w,w\[UndirectedEdge]Subscript[o, 2]}]*)


(* ::Subsubsection::Closed:: *)
(*Interpolation functions*)


(* ::Input:: *)
(*ClearAll[TTKF,TTKFmap,TTKFc]*)
(*TTKF=Import[pathTabulations<>"/TTK.dat"];*)
(*TTKFmap[x_]=ToExpression@TTKF[[1,1]];*)
(*TTKFc=ToExpression@TTKF[[2;;-1]];*)


(* ::Input:: *)
(*ClearAll[TTKFinterp,TTKFapprox]*)
(*TTKFinterp=Interpolation[TTKFc,Method->"Hermite",InterpolationOrder->10];*)
(*TTKFapprox[x_?NumericQ/;x<=TTKFc[[-5,1]]]:=TTKFinterp[x]*)
(*TTKFapprox[x_?NumericQ/;x>TTKFc[[-5,1]]]:=fasymp01[asympfit[MapAt[TTKFmap,TTKFc[[-5;;-2]],{All,1}]],map->TTKFmap][x];*)


(* ::Subsection::Closed:: *)
(**)


(* ::Text:: *)
(*Write the new integrands*)


(* ::Subsection::Closed:: *)
(*Tabulation change of variable*)


(* ::Text:: *)
(*Map for the variable of the numeric function from {0, \[Infinity]} to {0, 1}*)


(* ::Input:: *)
(*numMap[x_]:=(16 x)/(1-x)^2*)


(* ::Text:: *)
(*Getting ready the integrand for the integration*)


(* ::Input:: *)
(*numSub={\[Kappa]fncNum[y_]:>TTKFapprox[x],(\[Kappa]fncNum[z_])^y_:>(TTKFapprox[x])^y,\[Psi]bfuncNum[y_]:>TB1approx[x],(\[Psi]bfuncNum[z_])^y_:>(TB1approx[x])^y,c1sNum[y_]:>QSB0approx[x],(c1sNum[z_])^y_:>(QSB0approx[x])^y,t2bNum[y_]:>TB2approx[x],(t2bNum[z_])^y_:>(TB2approx[x])^y,c1s1bNum[y_]:>QSB1approx[x],(c1s1bNum[z_])^y_:>(QSB1approx[x])^y,c3bNum[y_]:>QB3approx[x],(c3bNum[z_])^y_:>(QB3approx[x])^y};*)


(* ::Input:: *)
(*numSubApply[integrand_]:=Simplify@If[MemberQ[integrand,\[Kappa]fncNum[x_]|(\[Kappa]fncNum[x_])^y_|\[Psi]bfuncNum[x_]|(\[Psi]bfuncNum[z_])^y_|c1sNum[y_]|(c1sNum[z_])^y_|c1s1bNum[y_]|(c1s1bNum[z_])^y_|c3bNum[y_]|(c3bNum[z_])^y_|t2bNum[y_]|(t2bNum[z_])^y_],(Simplify[integrand/x^2]/.x->numMap[x]/.numSub) meas[x],*)
(*integrand]*)


(* ::Input:: *)
(*(*numSubApply@((16384 \[Pi] ArcTan[x/2]^2 \[Kappa]fncNum[x])/(1+x^2))*)*)


(* ::Input:: *)
(*(*numSubApply@((65536 \[Pi]^2 x ArcTan[x/2] t2bNum[x])/(1+x^2)^3)*)*)


(* ::Subsection::Closed:: *)
(*Tabulations at l = 2*)


(* ::Input:: *)
(*specialNum=(\[Kappa]fncNum[x_]|(\[Kappa]fncNum[x_])^y_|\[Psi]bfuncNum[x_]|(\[Psi]bfuncNum[x_])^y_|c1sNum[x_]|(c1sNum[x_])^y_|c1s1bNum[x_]|(c1s1bNum[x_])^y_|c3bNum[x_]|(c3bNum[x_])^y_|t2bNum[x_]|(t2bNum[x_])^y_);*)


(* ::Input:: *)
(*varMapfun[int_]:=Cases[int,specialNum]*)


(* ::Input:: *)
(*varMap[int_]:=Cases[int,specialNum]/.specialNum:>x*)


(* ::Input:: *)
(*from2ltomorel[ooG4l_]:=Module[{listvar,posnoko},*)
(*listvar=varMap/@(ooG4l[[All,2]]);*)
(*posnoko=Complement[Range[1,Length@listvar],Flatten@Position[listvar,{}|{x}|{y}|{x,y}|{y,x}]];*)
(*Table[ooG4l[[posnoko[[i]],1]],{i,Length@posnoko}]*)
(*]*)


(* ::Input:: *)
(*from2ltomorelv2[ooG4l_]:=Module[{listvar,listvarv2,posnoko},*)
(*listvar=varMap/@(ooG4l[[All,2]]);*)
(*listvarv2=Sort/@DeleteDuplicates/@listvar;*)
(*posnoko=Complement[Range[1,Length@listvar],Flatten@Position[listvarv2,{}|{x}|{y}|{x,y}]];*)
(*Table[ooG4l[[posnoko[[i]],1]],{i,Length@posnoko}]*)
(*]*)


(* ::Input:: *)
(*numSubx={\[Kappa]fncNum[y_]:>TTKFapprox[x],(\[Kappa]fncNum[z_])^y_:>(TTKFapprox[x])^y,\[Psi]bfuncNum[y_]:>TB1approx[x],(\[Psi]bfuncNum[z_])^y_:>(TB1approx[x])^y,c1sNum[y_]:>QSB0approx[x],(c1sNum[z_])^y_:>(QSB0approx[x])^y,t2bNum[y_]:>TB2approx[x],(t2bNum[z_])^y_:>(TB2approx[x])^y,c1s1bNum[y_]:>QSB1approx[x],(c1s1bNum[z_])^y_:>(QSB1approx[x])^y,c3bNum[y_]:>QB3approx[x],(c3bNum[z_])^y_:>(QB3approx[x])^y};*)


(* ::Input:: *)
(*numSuby={\[Kappa]fncNum[yy_]:>TTKFapprox[y],(\[Kappa]fncNum[zz_])^yy_:>(TTKFapprox[y])^yy,\[Psi]bfuncNum[yy_]:>TB1approx[y],(\[Psi]bfuncNum[zz_])^yy_:>(TB1approx[y])^yy,c1sNum[yy_]:>QSB0approx[y],(c1sNum[zz_])^yy_:>(QSB0approx[y])^yy,t2bNum[yy_]:>TB2approx[y],(t2bNum[zz_])^yy_:>(TB2approx[y])^yy,c1s1bNum[yy_]:>QSB1approx[y],(c1s1bNum[zz_])^yy_:>(QSB1approx[y])^yy,c3bNum[yy_]:>QB3approx[y],(c3bNum[zz_])^yy_:>(QB3approx[y])^yy};*)


(* ::Input:: *)
(*numSubvar[xxx_]:={\[Kappa]fncNum[yy_]:>TTKFapprox[xxx],(\[Kappa]fncNum[zz_])^yy_:>(TTKFapprox[xxx])^yy,\[Psi]bfuncNum[yy_]:>TB1approx[xxx],(\[Psi]bfuncNum[zz_])^yy_:>(TB1approx[xxx])^yy,c1sNum[yy_]:>QSB0approx[xxx],(c1sNum[zz_])^yy_:>(QSB0approx[xxx])^yy,t2bNum[yy_]:>TB2approx[xxx],(t2bNum[zz_])^yy_:>(TB2approx[xxx])^yy,c1s1bNum[yy_]:>QSB1approx[xxx],(c1s1bNum[zz_])^yy_:>(QSB1approx[xxx])^yy,c3bNum[yy_]:>QB3approx[xxx],(c3bNum[zz_])^yy_:>(QB3approx[xxx])^yy};*)


(* ::Input:: *)
(*numSubvarfromx[xxx_]:={\[Kappa]fncNum[x]:>  TTKFapprox[xxx],(\[Kappa]fncNum[x])^yy_:>(TTKFapprox[xxx])^yy,\[Psi]bfuncNum[x]:>TB1approx[xxx],(\[Psi]bfuncNum[x])^yy_:>(TB1approx[xxx])^yy,c1sNum[x]:>QSB0approx[xxx],(c1sNum[x])^yy_:>(QSB0approx[xxx])^yy,t2bNum[x]:>TB2approx[xxx],(t2bNum[x])^yy_:>(TB2approx[xxx])^yy,c1s1bNum[x]:>QSB1approx[xxx],(c1s1bNum[x])^yy_:>(QSB1approx[xxx])^yy,c3bNum[x]:>QB3approx[xxx],(c3bNum[x])^yy_:>(QB3approx[xxx])^yy};*)


(* ::Input:: *)
(*numSubvarfromy[xxx_]:={\[Kappa]fncNum[y]:>  TTKFapprox[xxx],(\[Kappa]fncNum[y])^yy_:>(TTKFapprox[xxx])^yy,\[Psi]bfuncNum[y]:>TB1approx[xxx],(\[Psi]bfuncNum[y])^yy_:>(TB1approx[xxx])^yy,c1sNum[y]:>QSB0approx[xxx],(c1sNum[y])^yy_:>(QSB0approx[xxx])^yy,t2bNum[y]:>TB2approx[xxx],(t2bNum[y])^yy_:>(TB2approx[xxx])^yy,c1s1bNum[y]:>QSB1approx[xxx],(c1s1bNum[y])^yy_:>(QSB1approx[xxx])^yy,c3bNum[y]:>QB3approx[xxx],(c3bNum[y])^yy_:>(QB3approx[xxx])^yy};*)


(* ::Input:: *)
(*numSubl2v1[ooG4l_]:=Module[{ll,listvar,posnoko},*)
(*ll=Length@ooG4l;*)
(*listvar=varMap/@(ooG4l[[All,2]]);*)
(*Table[Simplify@Which[listvar[[i]]=={},{ooG4l[[i,1]],ooG4l[[i,2]],0},listvar[[i]]=== {x},{ooG4l[[i,1]],(Simplify[ooG4l[[i,2]]/x^2]/.x->numMap[p]/.numSubvar[p]) meas[p],1},listvar[[i]]=== {y},{ooG4l[[i,1]],(Simplify[ooG4l[[i,2]]/y^2]/.y->numMap[r]/.numSubvar[r]) meas[r],2},listvar[[i]]=== {x,y}||listvar[[i]]=== {y,x},{ooG4l[[i,1]],(Simplify[ooG4l[[i,2]]/(x^2 y^2)]/. numSubvarfromx[p]/. numSubvarfromy[r]/.x->numMap[p]/.y->numMap[r]) meas[p]meas[r],3}],{i,ll}]*)
(*]*)


(* ::Input:: *)
(*numSubl2v1nos[ooG4l_]:=Module[{ll,listvar,posnoko},*)
(*ll=Length@ooG4l;*)
(*listvar=varMap/@(ooG4l[[All,2]]);*)
(*Table[Which[listvar[[i]]=={},{ooG4l[[i,1]],ooG4l[[i,2]],0},listvar[[i]]=== {x},{ooG4l[[i,1]],(Simplify[ooG4l[[i,2]]/x^2]/.x->numMap[p]/.numSubvar[p]) meas[p],1},listvar[[i]]=== {y},{ooG4l[[i,1]],(Simplify[ooG4l[[i,2]]/y^2]/.y->numMap[r]/.numSubvar[r]) meas[r],2},listvar[[i]]=== {x,y}||listvar[[i]]=== {y,x},{ooG4l[[i,1]],(Simplify[ooG4l[[i,2]]/(x^2 y^2)]/. numSubvarfromx[p]/. numSubvarfromy[r]/.x->numMap[p]/.y->numMap[r]) meas[p]meas[r],3}],{i,ll}]*)
(*]*)


(* ::Input:: *)
(*nintduffymaxerrP[integrand_,maxerr_]:=AbsoluteTiming@NIntegrate[integrand,{p,0,1},{y,0,\[Infinity]},{\[Theta],0,\[Pi]},WorkingPrecision->20,Method->{"GlobalAdaptive","SingularityHandler"->"DuffyCoordinates","MaxErrorIncreases"->maxerr,"SingularityDepth"->1,"SymbolicProcessing"->0}]*)


(* ::Input:: *)
(*nintduffymaxerrR[integrand_,maxerr_]:=AbsoluteTiming@NIntegrate[integrand,{x,0,\[Infinity]},{r,0,1},{\[Theta],0,\[Pi]},WorkingPrecision->20,Method->{"GlobalAdaptive","SingularityHandler"->"DuffyCoordinates","MaxErrorIncreases"->maxerr,"SingularityDepth"->1,"SymbolicProcessing"->0}]*)


(* ::Input:: *)
(*nintduffymaxerrPR[integrand_,maxerr_]:=AbsoluteTiming@NIntegrate[integrand,{p,0,1},{r,0,1},{\[Theta],0,\[Pi]},WorkingPrecision->20,Method->{"GlobalAdaptive","SingularityHandler"->"DuffyCoordinates","MaxErrorIncreases"->maxerr,"SingularityDepth"->1,"SymbolicProcessing"->0}]*)


(* ::Chapter:: *)
(*Integration and Export for integration (from Integrands to Numbers or Files)*)


(* ::Section::Closed:: *)
(*Integration*)


(* ::Subsection::Closed:: *)
(*l = 1 (both with just analytic pieces and with numerical ones)*)


nint1l[integrand_]:=NIntegrate[Simplify[integrand],{x,0,\[Infinity]},WorkingPrecision->80,AccuracyGoal->30,PrecisionGoal->30]


nint1lNum[integrand_]:=NIntegrate[Simplify[integrand],{x,0,1},WorkingPrecision->80,AccuracyGoal->30,PrecisionGoal->30]


nint1lnumSubApply[integrand_]:=Simplify@If[MemberQ[integrand,\[Kappa]fncNum[x_]|(\[Kappa]fncNum[x_])^y_|\[Psi]bfuncNum[x_]|(\[Psi]bfuncNum[z_])^y_|c1sNum[y_]|(c1sNum[z_])^y_|c1s1bNum[y_]|(c1s1bNum[z_])^y_|c3bNum[y_]|(c3bNum[z_])^y_|t2bNum[y_]|(t2bNum[z_])^y_],nint1lNum[(Simplify[integrand/x^2]/.x->numMap[x]/.numSub) meas[x]],
nint1l[integrand]]


nint1lnumSubApplyAround[integrand_]:=Simplify@If[MemberQ[integrand,\[Kappa]fncNum[x_]|(\[Kappa]fncNum[x_])^y_|\[Psi]bfuncNum[x_]|(\[Psi]bfuncNum[z_])^y_|c1sNum[y_]|(c1sNum[z_])^y_|c1s1bNum[y_]|(c1s1bNum[z_])^y_|c3bNum[y_]|(c3bNum[z_])^y_|t2bNum[y_]|(t2bNum[z_])^y_],Around[nint1lNum[(Simplify[integrand/x^2]/.x->numMap[x]/.numSub) meas[x]],10^-18],
Around[nint1l[integrand],10^-25]]


(* ::Input:: *)
(*(*l=1;o6G4l1=Transpose[{l\[CapitalGamma]4level4poslo[o]\[LeftDoubleBracket]1+l\[RightDoubleBracket],(integrglevel4o[o]\[LeftDoubleBracket]l\[CapitalGamma]4level4poslo[o]\[LeftDoubleBracket]1+l\[RightDoubleBracket]\[RightDoubleBracket]nickelNorm[o]/o!)}]*)*)


(* ::Input:: *)
(*(*o6G4l1OKMathst1=Table[{o6G4l1\[LeftDoubleBracket]i,1\[RightDoubleBracket],(prepMathl1@PropNum0@simplroots@PropTrianNumMath@PropBubblSun@momExt0SimplAnalytic@o6G4l1\[LeftDoubleBracket]i,2\[RightDoubleBracket])},{i,Length@o6G4l1}]*)*)


(* ::Input:: *)
(*(*o6G4l1Values=Transpose[{o6G4l1OKMathst1\[LeftDoubleBracket]All,1\[RightDoubleBracket],nint1lnumSubApplyAround/@(o6G4l1OKMathst1\[LeftDoubleBracket]All,2\[RightDoubleBracket])}]*)*)


(* ::Subsection::Closed:: *)
(*l = 2*)


nintduffymaxerr[integrand_,maxerr_]:=AbsoluteTiming@NIntegrate[integrand,{x,0,\[Infinity]},{y,0,\[Infinity]},{\[Theta],0,\[Pi]},WorkingPrecision->20,Method->{"GlobalAdaptive","SingularityHandler"->"DuffyCoordinates","MaxErrorIncreases"->maxerr,"SingularityDepth"->1,"SymbolicProcessing"->0}]


(* ::Section::Closed:: *)
(*l = 2 related content*)


(* ::Subsection:: *)
(*Subs*)


(* ::Subsection::Closed:: *)
(*Assign error to collected one*)


(* ::Input:: *)
(*err2l={x_/;x<10^-13-> 10^-12,y_/; 10^-12>y>10^-13-> 10^-11,y_/; 10^-11>y>10^-12-> 10^-10,y_/;5 10^-10>y>10^-11-> 10^-9,y_/;5 10^-9>y>5 10^-10-> 10^-8,y_/;5 10^-8>y>5 10^-9-> 10^-7};*)


err2lUPDATED={x_/;x<10^-17-> 10^-15,y_/; 10^-14>y>10^-17-> 10^-13,y_/; 10^-12>y>10^-14-> 10^-12,y_/; 10^-12>y>10^-13-> 10^-11,y_/; 10^-11>y>10^-12-> 10^-10,y_/;5 10^-10>y>10^-11-> 10^-9,y_/;5 10^-9>y>5 10^-10-> 10^-8,y_/;5 10^-8>y>5 10^-9-> 10^-7,y_/;5 10^-7>y>5 10^-8-> 10^-6};


(* ::Input:: *)
(*(*Abs[resultso6G4l2me20000-resultso6G4l2me15000]/.err2l;*)
(*resultso6G4l2me20000Ar=Table[Around[resultso6G4l2me20000\[LeftDoubleBracket]i\[RightDoubleBracket],%\[LeftDoubleBracket]i\[RightDoubleBracket]],{i,52}]*)*)


(* ::Input:: *)
(*(*o6G4l2values2e5=(resultso6G4l2me20000Ar/.x_/;x==resultso6G4l2me20000Ar[[13]]\[Rule] o6G4d47replace)*)*)


(* ::Subsection::Closed:: *)
(*Import and export for l = 2 (Just comment)*)


(* ::Subsubsection::Closed:: *)
(*Export*)


(* ::Input:: *)
(*(*textwlv8[intmp_,i_, maxerr_,version_]:=Module[{},*)
(*"\"This is the integration of diagram" <>ToString[i]<>"\"\n\n*)
(*path = \"../SubTAB/subdiagram_tables\";\n\n*)
(*i = " <>ToString[i]<>";\n*)
(*maxerr = " <>ToString[maxerr]<>";\n*)
(*numberintegrand = " <>ToString[intmp\[LeftDoubleBracket]i,1\[RightDoubleBracket]]<>";\n\n*)
(*integrationspace = " <>ToString[intmp\[LeftDoubleBracket]i,3\[RightDoubleBracket]]<>";\n\n*)
(*If[integrationspace=!=0, *)
(*getexponent[xylist_]:= Mean @Table[-Log[xylist[[k+1,2]]/xylist[[k,2]]]/Log[xylist[[k+1,1]]/(xylist[[k,1]])],{k,1,Length@xylist-1}]; \n*)
(*asympfit[xylist_]:=Module[{alpha,c},alpha = getexponent[xylist]; c = xylist[[1,2]]*xylist[[1,1]]^alpha;{alpha,c}];\n*)
(*fasymp01mine[{a_,b_},x_]:= 16^(-a) b (x/(1-x)^2)^(-a);\n\n*)
(**)
(*TB1=Import[path<>\"/TB1.dat\"];*)
(*TB1map[x_]=ToExpression@TB1[[1,1]];*)
(*TB1c=ToExpression@TB1\[LeftDoubleBracket]2;;-1\[RightDoubleBracket];*)
(*TB1interp=Interpolation[TB1c,Method\[Rule]\"Hermite\",InterpolationOrder\[Rule]10];*)
(*TB1approx[x_?NumericQ/;x\[LessEqual]TB1c[[-5,1]]]:=TB1interp[x];*)
(*TB1approx[x_?NumericQ/;x>TB1c[[-5,1]]]:=fasymp01mine[asympfit[MapAt[TB1map,TB1c[[-5;;-2]],{All,1}]],x];\n*)
(*TB2=Import[path<>\"/TB2.dat\"];*)
(*TB2map[x_]=ToExpression@TB2[[1,1]];*)
(*TB2c=ToExpression@TB2[[2;;-1]];*)
(*TB2interp=Interpolation[TB2c,Method\[Rule]\"Hermite\",InterpolationOrder\[Rule]10];*)
(*TB2approx[x_?NumericQ/;x\[LessEqual]TB2c[[-5,1]]]:=TB2interp[x];*)
(*TB2approx[x_?NumericQ/;x>TB2c[[-5,1]]]:=fasymp01mine[asympfit[MapAt[TB2map,TB2c[[-5;;-2]],{All,1}]],x];\n*)
(*QSB0=Import[path<>\"/QSB0.dat\"];*)
(*QSB0map[x_]=ToExpression@QSB0[[1,1]];*)
(*QSB0c=ToExpression@QSB0\[LeftDoubleBracket]2;;-1\[RightDoubleBracket];*)
(*QSB0interp=Interpolation[QSB0c,Method\[Rule]\"Hermite\",InterpolationOrder\[Rule]10];*)
(*QSB0approx[x_?NumericQ/;x\[LessEqual]QSB0c\[LeftDoubleBracket]-5,1\[RightDoubleBracket]]:=QSB0interp[x];*)
(*QSB0approx[x_?NumericQ/;x>QSB0c\[LeftDoubleBracket]-5,1\[RightDoubleBracket]]:=fasymp01mine[asympfit[MapAt[QSB0map,QSB0c[[-5;;-2]],{All,1}]],x];\n*)
(*QSB1=Import[path<>\"/QSB1.dat\"];*)
(*QSB1map[x_]=ToExpression@QSB1[[1,1]];*)
(*QSB1c=ToExpression@QSB1[[2;;-1]];*)
(*QSB1interp=Interpolation[QSB1c,Method\[Rule]\"Hermite\",InterpolationOrder\[Rule]10];*)
(*QSB1approx[x_?NumericQ/;x\[LessEqual]QSB1c[[-5,1]]]:=QSB1interp[x];*)
(*QSB1approx[x_?NumericQ/;x>QSB1c[[-5,1]]]:=fasymp01mine[asympfit[MapAt[QSB1map,QSB1c[[-5;;-2]],{All,1}]],x];\n*)
(*QB3=Import[path<>\"/QB3.dat\"];*)
(*QB3map[x_]=ToExpression@QB3[[1,1]];*)
(*QB3c=ToExpression@QB3[[2;;-1]];*)
(*QB3interp=Interpolation[QB3c,Method\[Rule]\"Hermite\",InterpolationOrder\[Rule]10];*)
(*QB3approx[x_?NumericQ/;x\[LessEqual]QB3c[[-5,1]]]:=QB3interp[x];*)
(*QB3approx[x_?NumericQ/;x>QB3c[[-5,1]]]:=fasymp01mine[asympfit[MapAt[QB3map,QB3c[[-5;;-2]],{All,1}]],x];\n*)
(*TTKF=Import[path<>\"/TTK.dat\"];*)
(*TTKFmap[x_]=ToExpression@TTKF[[1,1]];*)
(*TTKFc=ToExpression@TTKF[[2;;-1]];*)
(*TTKFinterp=Interpolation[TTKFc,Method\[Rule]\"Hermite\",InterpolationOrder\[Rule]10];*)
(*TTKFapprox[x_?NumericQ/;x\[LessEqual]TTKFc[[-5,1]]]:=TTKFinterp[x];*)
(*TTKFapprox[x_?NumericQ/;x>TTKFc[[-5,1]]]:=fasymp01mine[asympfit[MapAt[TTKFmap,TTKFc[[-5;;-2]],{All,1}]],x];\n*)
(*];\n\n*)
(*integrand = " <>ToString[intmp\[LeftDoubleBracket]i,2\[RightDoubleBracket],InputForm]<>";\n\n*)
(*nintduffymaxerr[integrand_,maxerr_]:=AbsoluteTiming@NIntegrate[integrand,{x,0,\[Infinity]},{y,0,\[Infinity]},{\[Theta],0,\[Pi]},WorkingPrecision\[Rule]20,Method\[Rule]{\"GlobalAdaptive\",\"SingularityHandler\"\[Rule]\"DuffyCoordinates\",\"MaxErrorIncreases\"\[Rule]maxerr,\"SingularityDepth\"\[Rule]1,\"SymbolicProcessing\"\[Rule]0}]\n*)
(*nintduffymaxerrP[integrand_,maxerr_]:=AbsoluteTiming@NIntegrate[integrand,{p,0,1},{y,0,\[Infinity]},{\[Theta],0,\[Pi]},WorkingPrecision\[Rule]20,Method\[Rule]{\"GlobalAdaptive\",\"SingularityHandler\"\[Rule]\"DuffyCoordinates\",\"MaxErrorIncreases\"\[Rule]maxerr,\"SingularityDepth\"\[Rule]1,\"SymbolicProcessing\"\[Rule]0}]\n*)
(*nintduffymaxerrR[integrand_,maxerr_]:=AbsoluteTiming@NIntegrate[integrand,{x,0,\[Infinity]},{r,0,1},{\[Theta],0,\[Pi]},WorkingPrecision\[Rule]20,Method\[Rule]{\"GlobalAdaptive\",\"SingularityHandler\"\[Rule]\"DuffyCoordinates\",\"MaxErrorIncreases\"\[Rule]maxerr,\"SingularityDepth\"\[Rule]1,\"SymbolicProcessing\"\[Rule]0}]\n*)
(*nintduffymaxerrPR[integrand_,maxerr_]:=AbsoluteTiming@NIntegrate[integrand,{p,0,1},{r,0,1},{\[Theta],0,\[Pi]},WorkingPrecision\[Rule]20,Method\[Rule]{\"GlobalAdaptive\",\"SingularityHandler\"\[Rule]\"DuffyCoordinates\",\"MaxErrorIncreases\"\[Rule]maxerr,\"SingularityDepth\"\[Rule]1,\"SymbolicProcessing\"\[Rule]0}]\n\n*)
(*nintduffychosen[integrationspace_,integrand_,maxerr_]:=Which[integrationspace==0, nintduffymaxerr[integrand,maxerr], integrationspace==1, nintduffymaxerrP[integrand,maxerr], integrationspace==2, nintduffymaxerrR[integrand,maxerr], integrationspace==3, nintduffymaxerrPR[integrand,maxerr]]\n\n*)
(*result = nintduffychosen[integrationspace,integrand,maxerr];\n*)
(*out = {i, numberintegrand, maxerr, result[[1]], result[[2]]};\n\n*)
(*Export[\"./result_v"<>ToString[version]<>"_diag"<>ToString[i]<>"_me"<>ToString[maxerr]<>".txt\", out,\"Text\"]*)
(*"]*)*)


(* ::Input:: *)
(*textwlv8Ul[intmp_,i_, maxerr_,version_]:=Module[{},*)
(*"\"This is the integration of diagram" <>ToString[i]<>"\"\n\n*)
(*path = \"../subdiagram_tables\";\n\n*)
(*i = " <>ToString[i]<>";\n*)
(*maxerr = " <>ToString[maxerr]<>";\n*)
(*numberintegrand = " <>ToString[intmp[[i,1]]]<>";\n\n*)
(*integrationspace = " <>ToString[intmp[[i,3]]]<>";\n\n*)
(*If[integrationspace=!=0, *)
(*getexponent[xylist_]:= Mean @Table[-Log[xylist[[k+1,2]]/xylist[[k,2]]]/Log[xylist[[k+1,1]]/(xylist[[k,1]])],{k,1,Length@xylist-1}]; \n*)
(*asympfit[xylist_]:=Module[{alpha,c},alpha = getexponent[xylist]; c = xylist[[1,2]]*xylist[[1,1]]^alpha;{alpha,c}];\n*)
(*fasymp01mine[{a_,b_},x_]:= 16^(-a) b (x/(1-x)^2)^(-a);\n\n*)
(**)
(*TB1=Import[path<>\"/TB1.dat\"];*)
(*TB1map[x_]=ToExpression@TB1[[1,1]];*)
(*TB1c=ToExpression@TB1\[LeftDoubleBracket]2;;-1\[RightDoubleBracket];*)
(*TB1interp=Interpolation[TB1c,Method\[Rule]\"Hermite\",InterpolationOrder\[Rule]10];*)
(*TB1approx[x_?NumericQ/;x\[LessEqual]TB1c[[-5,1]]]:=TB1interp[x];*)
(*TB1approx[x_?NumericQ/;x>TB1c[[-5,1]]]:=fasymp01mine[asympfit[MapAt[TB1map,TB1c[[-5;;-2]],{All,1}]],x];\n*)
(*TB2=Import[path<>\"/TB2.dat\"];*)
(*TB2map[x_]=ToExpression@TB2[[1,1]];*)
(*TB2c=ToExpression@TB2[[2;;-1]];*)
(*TB2interp=Interpolation[TB2c,Method\[Rule]\"Hermite\",InterpolationOrder\[Rule]10];*)
(*TB2approx[x_?NumericQ/;x\[LessEqual]TB2c[[-5,1]]]:=TB2interp[x];*)
(*TB2approx[x_?NumericQ/;x>TB2c[[-5,1]]]:=fasymp01mine[asympfit[MapAt[TB2map,TB2c[[-5;;-2]],{All,1}]],x];\n*)
(*QSB0=Import[path<>\"/QSB0.dat\"];*)
(*QSB0map[x_]=ToExpression@QSB0[[1,1]];*)
(*QSB0c=ToExpression@QSB0\[LeftDoubleBracket]2;;-1\[RightDoubleBracket];*)
(*QSB0interp=Interpolation[QSB0c,Method\[Rule]\"Hermite\",InterpolationOrder\[Rule]10];*)
(*QSB0approx[x_?NumericQ/;x\[LessEqual]QSB0c\[LeftDoubleBracket]-5,1\[RightDoubleBracket]]:=QSB0interp[x];*)
(*QSB0approx[x_?NumericQ/;x>QSB0c\[LeftDoubleBracket]-5,1\[RightDoubleBracket]]:=fasymp01mine[asympfit[MapAt[QSB0map,QSB0c[[-5;;-2]],{All,1}]],x];\n*)
(*QSB1=Import[path<>\"/QSB1.dat\"];*)
(*QSB1map[x_]=ToExpression@QSB1[[1,1]];*)
(*QSB1c=ToExpression@QSB1[[2;;-1]];*)
(*QSB1interp=Interpolation[QSB1c,Method\[Rule]\"Hermite\",InterpolationOrder\[Rule]10];*)
(*QSB1approx[x_?NumericQ/;x\[LessEqual]QSB1c[[-5,1]]]:=QSB1interp[x];*)
(*QSB1approx[x_?NumericQ/;x>QSB1c[[-5,1]]]:=fasymp01mine[asympfit[MapAt[QSB1map,QSB1c[[-5;;-2]],{All,1}]],x];\n*)
(*QB3=Import[path<>\"/QB3.dat\"];*)
(*QB3map[x_]=ToExpression@QB3[[1,1]];*)
(*QB3c=ToExpression@QB3[[2;;-1]];*)
(*QB3interp=Interpolation[QB3c,Method\[Rule]\"Hermite\",InterpolationOrder\[Rule]10];*)
(*QB3approx[x_?NumericQ/;x\[LessEqual]QB3c[[-5,1]]]:=QB3interp[x];*)
(*QB3approx[x_?NumericQ/;x>QB3c[[-5,1]]]:=fasymp01mine[asympfit[MapAt[QB3map,QB3c[[-5;;-2]],{All,1}]],x];\n*)
(*TTKF=Import[path<>\"/TTK.dat\"];*)
(*TTKFmap[x_]=ToExpression@TTKF[[1,1]];*)
(*TTKFc=ToExpression@TTKF[[2;;-1]];*)
(*TTKFinterp=Interpolation[TTKFc,Method\[Rule]\"Hermite\",InterpolationOrder\[Rule]10];*)
(*TTKFapprox[x_?NumericQ/;x\[LessEqual]TTKFc[[-5,1]]]:=TTKFinterp[x];*)
(*TTKFapprox[x_?NumericQ/;x>TTKFc[[-5,1]]]:=fasymp01mine[asympfit[MapAt[TTKFmap,TTKFc[[-5;;-2]],{All,1}]],x];\n*)
(*];\n\n*)
(*integrand = " <>ToString[intmp[[i,2]],InputForm]<>";\n\n*)
(*nintduffymaxerr[integrand_,maxerr_]:=AbsoluteTiming@NIntegrate[integrand,{x,0,\[Infinity]},{y,0,\[Infinity]},{\[Theta],0,\[Pi]},WorkingPrecision\[Rule]20,Method\[Rule]{\"GlobalAdaptive\",\"SingularityHandler\"\[Rule]\"DuffyCoordinates\",\"MaxErrorIncreases\"\[Rule]maxerr,\"SingularityDepth\"\[Rule]1,\"SymbolicProcessing\"\[Rule]0}]\n*)
(*nintduffymaxerrP[integrand_,maxerr_]:=AbsoluteTiming@NIntegrate[integrand,{p,0,1},{y,0,\[Infinity]},{\[Theta],0,\[Pi]},WorkingPrecision\[Rule]20,Method\[Rule]{\"GlobalAdaptive\",\"SingularityHandler\"\[Rule]\"DuffyCoordinates\",\"MaxErrorIncreases\"\[Rule]maxerr,\"SingularityDepth\"\[Rule]1,\"SymbolicProcessing\"\[Rule]0}]\n*)
(*nintduffymaxerrR[integrand_,maxerr_]:=AbsoluteTiming@NIntegrate[integrand,{x,0,\[Infinity]},{r,0,1},{\[Theta],0,\[Pi]},WorkingPrecision\[Rule]20,Method\[Rule]{\"GlobalAdaptive\",\"SingularityHandler\"\[Rule]\"DuffyCoordinates\",\"MaxErrorIncreases\"\[Rule]maxerr,\"SingularityDepth\"\[Rule]1,\"SymbolicProcessing\"\[Rule]0}]\n*)
(*nintduffymaxerrPR[integrand_,maxerr_]:=AbsoluteTiming@NIntegrate[integrand,{p,0,1},{r,0,1},{\[Theta],0,\[Pi]},WorkingPrecision\[Rule]20,Method\[Rule]{\"GlobalAdaptive\",\"SingularityHandler\"\[Rule]\"DuffyCoordinates\",\"MaxErrorIncreases\"\[Rule]maxerr,\"SingularityDepth\"\[Rule]1,\"SymbolicProcessing\"\[Rule]0}]\n\n*)
(*nintduffychosen[integrationspace_,integrand_,maxerr_]:=Which[integrationspace==0, nintduffymaxerr[integrand,maxerr], integrationspace==1, nintduffymaxerrP[integrand,maxerr], integrationspace==2, nintduffymaxerrR[integrand,maxerr], integrationspace==3, nintduffymaxerrPR[integrand,maxerr]]\n\n*)
(*result = nintduffychosen[integrationspace,integrand,maxerr];\n*)
(*out = {i, numberintegrand, maxerr, result[[1]], result[[2]]};\n\n*)
(*Export[\"./result_v"<>ToString[version]<>"_diag"<>ToString[i]<>"_me"<>ToString[maxerr]<>".txt\", out,\"Text\"]*)
(*"]*)


(* ::Input:: *)
(*(*exptr[integr_,version_,maxerr_]:=Module[{dir1},*)
(*dir1="o7G4l2v"<>ToString[version];*)
(*CreateDirectory[dir1];*)
(*Table[Export[dir1<>"/intgr_v"<>ToString[version]<>"_d"<>ToString[i]<>"_me"<>ToString[maxerr]<>".wl",textwlv8[integr,i,maxerr,version],"Text"],{i,Length@integr}]*)
(*]*)*)


(* ::Input:: *)
(*exptr[integr_,version_,maxerr_]:=Module[{dir1},*)
(*dir1="o7G4l2v"<>ToString[version];*)
(*CreateDirectory[dir1];*)
(*Table[Export[dir1<>"/intgr_v"<>ToString[version]<>"_d"<>ToString[i]<>"_me"<>ToString[maxerr]<>".wl",textwlv8Ul[integr,i,maxerr,version],"Text"],{i,Length@integr}]*)
(*]*)


(* ::Input:: *)
(*(*exptr[o7G4l2OKMathst2num,3,15000]*)*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*mineslurmv2[version_,start_,end_,maxerr_]:=Table["echo \"RUNNING intgr_v"<>ToString[version]<>"_d"<>ToString[i]<>"_me"<>ToString[maxerr]<>".wl\"\nwolfram -script intgr_v"<>ToString[version]<>"_d"<>ToString[i]<>"_me"<>ToString[maxerr]<>".wl &\n",{i,start,end}]*)


(* ::Input:: *)
(*slurmscriptMathUlyssesOK[version_,start_,end_,maxerr_,mem_,time_,part_,name_]:=Module[{strjob},*)
(*strjob=mineslurmv2[version,start,end,maxerr];*)
(*"#!/bin/bash*)
(**)
(**)
(*#SBATCH--job-name="<>name <>"*)
(*#*)
(*#SBATCH--ntasks=1*)
(*#SBATCH--ntasks-per-node="<>ToString[end-start+1]<>"*)
(*#SBATCH--cpus-per-task=1*)
(*#*)
(*#SBATCH--mem="<>mem <>"*)
(*#*)
(*#SBATCH--partition="<>part <>"*)
(*#SBATCH--time="<>time <>":00:00*)
(*#SBATCH--output=%x.o%j*)
(*#SBATCH--error=%x.e%j*)
(**)
(*module load mathematica/12.1.1*)
(**)
(**)
(*"<>strjob <>"*)
(**)
(*wait"] *)


(* ::Input:: *)
(*slurmMath[position_,version_,totend_,maxerr_,mem_,time_,part_,name_]:=Module[{qr,nscripts,jj,st,en},*)
(*qr=QuotientRemainder[totend,40];*)
(*nscripts=qr[[1]]+If[qr[[2]]>0,1,0];*)
(*For[jj=0, jj<nscripts,jj++,*)
(*{st,en}={40 jj+1,Min[40 (jj+1),totend]};*)
(*Export["./"<>position <>"/job"<>ToString[jj]<>".slurm",slurmscriptMathUlyssesOK[version,st,en,maxerr,mem,time,part,name],"Text"]]*)
(*]*)


(* ::Input:: *)
(*slurmMathv2[position_,version_,totend_,maxerr_,mem_,time_,part_,name_]:=Module[{qr,nscripts,jj,st,en},*)
(*qr=QuotientRemainder[totend,40];*)
(*nscripts=qr[[1]]+If[qr[[2]]>0,1,0];*)
(*For[jj=0, jj<nscripts,jj++,*)
(*{st,en}={40 jj+1,Min[40 (jj+1),totend]};*)
(*Export["./"<>position <>"/job"<>ToString[jj]<>".slurm",slurmscriptMathUlyssesOK[version,st,en,maxerr,mem,time,part,StringJoin[name,ToString@jj]],"Text"]]*)
(*]*)


(* ::Input:: *)
(*(*slurmMath["o7G4l2v1",1,298,500,"20000mb","01","regular1","Same"]*)*)


(* ::Input:: *)
(*(*slurmMath["o7G4l2v2",2,298,10000,"20000mb","03","regular1","Same"]*)*)


(* ::Input:: *)
(*(*slurmMathv2["o7G4l2v3",3,298,15000,"20000mb","04","regular1","Same"]*)*)


(* ::Input:: *)
(*(*slurmMath["o7G4l2v4",4,298,20000,"20000mb","05","regular1","Same"]*)*)


(* ::Input:: *)
(*(*slurmMathv2["o7G4l2v5",5,298,30000,"20000mb","08","regular1","Same"]*)*)


(* ::Input:: *)
(**)


(* ::Subsubsection::Closed:: *)
(*Import*)


(* ::Input:: *)
(*(*import[version_,i_,maxerr_]:=Import["Ulysses/ResultsMathematica/o7G4v"<>ToString[version]<>"/result_v"<>ToString[version]<>"_diag"<>ToString[i]<>"_me"<>ToString[maxerr]<>".txt","Table"]*)*)


(* ::Input:: *)
(*(*Table[import[1,i,500],{i,298}];*)*)


(* ::Input:: *)
(*(*resultso7G4l2me500=Transpose@ToExpression[{%\[LeftDoubleBracket]All,2,1\[RightDoubleBracket],%\[LeftDoubleBracket]All,5,1\[RightDoubleBracket]}]*)*)


(* ::Title:: *)
(*Visible Functions*)


(* ::Section::Closed:: *)
(*Aliases*)


(* ::Subsubsection::Closed:: *)
(*O(N) and cubic factors*)


\[CapitalNu]:=NComponents[]


\[CapitalChi]:=XCubicRatio[]


(* ::Subsubsection::Closed:: *)
(*Subdiagrams*)


\[ScriptQ][i_]:=Momentum[i]
\[ScriptQ][i_,sub_]:=Momentum[i,sub]
Momentum/:Subscript[Momentum[i_],sub_]:=Momentum[i,sub]


\[ScriptP][i_]:=ExternalMomentum[i]
\[ScriptP][i_,sub_]:=ExternalMomentum[i,sub]
ExternalMomentum/:Subscript[ExternalMomentum[i_],sub_]:=ExternalMomentum[i,sub]


\[ScriptCapitalG][q_]:=Propagator[q]


\[ScriptCapitalB][q_]:=BubbleSubdiagram[q]


\[ScriptCapitalS][q_]:=SunsetSubdiagram[q]


\[ScriptCapitalT][q1_,q2_,q3_]:=TriangleSubdiagram[q1,q2,q3]


\[ScriptCapitalQ][q1_,q2_,q3_,q4_]:=SquareSubdiagram[q1,q2,q3,q4]


\[ScriptT]\[ScriptCapitalS]:=TadSunsetSubdiagram[]


\[ScriptT]\[ScriptCapitalT]\[ScriptCapitalB]:=TadTriangleBubblesSubdiagram[]


(* ::Section::Closed:: *)
(*StandardForm for symbols*)


(* ::Subsubsection::Closed:: *)
(*O(N) and cubic factors*)


NComponents/:MakeBoxes[NComponents[],StandardForm]:=RowBox[{"\[CapitalNu]"}]


XCubicRatio/:MakeBoxes[XCubicRatio[],StandardForm]:=RowBox[{"\[CapitalChi]"}]


(* ::Subsubsection::Closed:: *)
(*Subdiagrams*)


Propagator/:MakeBoxes[Propagator[qs_],StandardForm]:=RowBox[{"\[ScriptCapitalG]","[",ToBoxes[qs,StandardForm],"]"}]


BubbleSubdiagram/:MakeBoxes[BubbleSubdiagram[qs_],StandardForm]:=RowBox[{"\[ScriptCapitalB]","[",ToBoxes[qs,StandardForm],"]"}]


SunsetSubdiagram/:MakeBoxes[SunsetSubdiagram[qs_],StandardForm]:=RowBox[{"\[ScriptCapitalS]","[",ToBoxes[qs,StandardForm],"]"}]


TriangleSubdiagram/:MakeBoxes[TriangleSubdiagram[q1_,q2_,q3_],StandardForm]:=RowBox[{"\[ScriptCapitalT]","[",ToBoxes[q1,StandardForm],",",ToBoxes[q2,StandardForm],",",ToBoxes[q3,StandardForm],"]"}]


SquareSubdiagram/:MakeBoxes[SquareSubdiagram[q1_,q2_,q3_,q4_],StandardForm]:=RowBox[{"\[ScriptCapitalQ]","[",ToBoxes[q1,StandardForm],",",ToBoxes[q2,StandardForm],",",ToBoxes[q3,StandardForm],",",ToBoxes[q4,StandardForm],"]"}]


TadSunsetSubdiagram/:MakeBoxes[TadSunsetSubdiagram[],StandardForm]:=RowBox[{"\[ScriptT]\[ScriptCapitalS]"}]


TadTriangleBubblesSubdiagram/:MakeBoxes[TadTriangleBubblesSubdiagram[],StandardForm]:=RowBox[{"\[ScriptT]\[ScriptCapitalT]\[ScriptCapitalB]"}]


(* ::Subsubsection:: *)
(*Momenta representation*)


(*internalMomentumList = Map[ToString,qaux]*)
Momentum/:MakeBoxes[Momentum[n_Integer/;n>0],StandardForm]:=RowBox[{StringJoin["\[ScriptQ][",ToString[n],"]"]}]
Momentum/:MakeBoxes[Momentum[n_Integer/;n>0,sub_],StandardForm]:=RowBox[{SubscriptBox[StringJoin["\[ScriptQ][",ToString[n],"]"],ToBoxes@sub]}]


(*externalMomentumList = Map[ToString,pextlist]*)
ExternalMomentum/:MakeBoxes[ExternalMomentum[n_Integer/;n>0],StandardForm]:=RowBox[{StringJoin["\[ScriptP][",ToString[n],"]"]}]
ExternalMomentum/:MakeBoxes[ExternalMomentum[n_Integer/;n>0,sub_],StandardForm]:=RowBox[{SubscriptBox[StringJoin["\[ScriptP][",ToString[n],"]"],ToBoxes@sub]}]


(* ::Section:: *)
(*Functions*)


(* ::Subsection:: *)
(*Labels*)


ImportedIndicesQ[n_,o3_,o4_]:=Module[{indicesHead},
indicesHead=Head@nickelIndices[n,o3,o4];
If[indicesHead===nickelIndices,False,True]
]


ImportIndices[n_,o3_,o4_]:=(nickelIndices[n,o3,o4]=importLables1PI[
		FileNameJoin[{PacletManager`PacletResource["GSberveglieri/Phi4tools","EDGES"],ToString[n]<>"pt"<>ToString[o3]<>ToString[o4]<>".txt"}]
	]);


NickelIndex[
	n_/;IntegerQ[n] && NonNegative[n] && EvenQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&EvenQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	d_/;IntegerQ[d] && Positive[d]]:=Module[{},
If[n=!=0&&n=!=2&&n=!=4,Message[NickelIndex::nout],
If[n=!=4 && o3==0 && o4==1 && d==1, Message[NickelIndex::tadpole],
If[n==4 && o3==2 && o4==0, Message[NickelIndex::dtoobig,d,0],
If[o4>9||(o4==9&&o3>0),Message[NickelIndex::nvert],
If[ImportedIndicesQ[n,o3,o4],Null;,ImportIndices[n,o3,o4];];
If[d>Length[nickelIndices[n,o3,o4]],Message[NickelIndex::dtoobig,d,Length[nickelIndices[n,o3,o4]]],
nickelIndices[n,o3,o4][[d]]]]]]]
]


NickelIndex[
	n_/;IntegerQ[n] && NonNegative[n] && EvenQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&EvenQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4]]:=Module[{},
If[n=!=0&&n=!=2&&n=!=4,Message[NickelIndex::nout],
If[n=!=4 && o3==0 && o4==1, Message[NickelIndex::tadpole],
If[n==4 && o3==2 && o4==0, {},
If[o4>9||(o4==9&&o3>0),Message[NickelIndex::nvert],
If[ImportedIndicesQ[n,o3,o4],Null;,ImportIndices[n,o3,o4];];
nickelIndices[n,o3,o4]]]]]
]


NickelIndex::nout = "The number of external legs is not 0, 2 or 4. This(These) diagram(s) is(are) not in this package."


NickelIndex::tadpole = "Within this renomalization scheme the diagrams with tadpoles have been removed because they are identically zero."


(*NickelIndex::nvert = "The number of quartic vertices is greater than 9 or equal 9 but there are cubic vertices as well. This(These) diagram(s) is(are) not in this package."*)


NickelIndex::nvert = "Diagram(s) not in this package."


NickelIndex::dtoobig = "Argument `1` is too big, take d\[LessEqual]`2`."


NickelIndex[
	n_/;IntegerQ[n] && NonNegative[n] && OddQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&& OddQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	d_/;IntegerQ[d] && Positive[d]]:=Module[{},
Message[NickelIndex::nout]
]


NickelIndex[
	n_/;IntegerQ[n] && NonNegative[n] && OddQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&& OddQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4]]:=Module[{},
Message[NickelIndex::nout]
]


(* ::Subsection:: *)
(*Symmetry factors*)


renameCubicFactors = {n->NComponents[], x->XCubicRatio[]};


SymmetryFactorDiagram[
	n_/;IntegerQ[n] && NonNegative[n] && EvenQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&EvenQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	d_/;IntegerQ[d] && Positive[d],
	OptionsPattern[]]:=
Module[{subs=OptionValue["Tensor"]},
	If[o3==0 && o4==0, Message[SymmetryFactorDiagram::invalid],
	If[n=!=0&&n=!=2&&n=!=4,
		Message[SymmetryFactorDiagram::nout],
		If[o3=!=0,
			Return[Missing["Not computed"]],
			If[n=!=4 && o3==0 && o4==1,
				Message[SymmetryFactorDiagram::tadpole],
				(*If[n==0 && o3==0 && o4>1 && o4<4 && d==1,
					Return[Missing["Not provided by the current version of the package"]],*)
					If[n==0 && o3==0 && o4>1 && o4<4 && d>1,
						Message[SymmetryFactorDiagram::dtoobig,d,1],
						If[o4==9,
						Return[Missing["Not computed"]],
						If[o4>=10,
							Message[SymmetryFactorDiagram::nvert],
							If[d>Length@symfactorGno[n,o4],
								Message[SymmetryFactorDiagram::dtoobig,d,Length[symfactorGno[n,o4]]],
								Which[
									subs=="O(N)", If[n==0||n==2, Return[symfactorGno[n,o4][[d]]/.x-> 0], Return[symfactorGno[n,o4][[d,1;;-2]]/.x-> 0]],
									subs== "Cubic", Return[symfactorGno[n,o4][[d]]],
									And@@{subs=!="O(N)",subs=!="Cubic"},
									Message[SymmetryFactorDiagram::tensor];
									If[n==0||n==2, Return[symfactorGno[n,o4][[d]]/.x-> 0], Return[symfactorGno[n,o4][[d,1;;-2]]/.x-> 0]]
								]
							]
						]
						]
					]
				(*]*)
			]
		]
	]
	]
]/.renameCubicFactors


SymmetryFactorDiagram::invalid = "Invalid arguments."


SymmetryFactorDiagram::nout = "The number of external legs is not 0, 2 or 4. This(These) diagram(s) is(are) not in this package."


SymmetryFactorDiagram::tensor = "The options for \"Tensor\" are: \"O(N)\" and \"Cubic\". Defaulting to \"Tensor\" -> \"O(N)\"."


SymmetryFactorDiagram::tadpole = "Within this renomalization scheme the diagrams with tadpoles have been removed because they are identically zero."


SymmetryFactorDiagram::nvert = "Diagram(s) not in this package."


SymmetryFactorDiagram::dtoobig = "Argument `1` is too big, take d\[LessEqual]`2`."


SymmetryFactorDiagram[
	n_/;IntegerQ[n] && NonNegative[n] && EvenQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&EvenQ[o3], 
	o4_/;IntegerQ[o4] && NonNegative[o4],
	OptionsPattern[]]:=
Module[{subs=OptionValue["Tensor"]},
	If[o3==0 && o4==0, Message[SymmetryFactorDiagram::invalid],
	If[n=!=0&&n=!=2&&n=!=4,
		Message[SymmetryFactorDiagram::nout],
		If[o3=!=0,
			Return[Missing["Not computed"]],
			If[n=!=4 && o3==0 && o4==1,
				Message[SymmetryFactorDiagram::tadpole],
				(*If[n==0 && o3==0 && o4>1 && o4<4,
					Return[Missing["Not provided by the current version of the package"]],*)
					If[o4==9,
						Return[Missing["Not computed"]],
						If[o4>=10,
							Message[SymmetryFactorDiagram::nvert],
						Which[
							subs=="O(N)", If[n==0||n==2, Return[symfactorGno[n,o4]/.x-> 0], Return[symfactorGno[n,o4][[All,1;;-2]]/.x-> 0]],
							subs== "Cubic", Return[symfactorGno[n,o4]],
							And@@{subs=!="O(N)",subs=!="Cubic"},
							Message[SymmetryFactorDiagram::tensor];
							If[n==0||n==2, Return[symfactorGno[n,o4]/.x-> 0], Return[symfactorGno[n,o4][[All,1;;-2]]/.x-> 0]]
						]
					]
					]
				(*]*)
			]
		]
	]
	]
]/.renameCubicFactors


SymmetryFactorDiagram[
	n_/;IntegerQ[n] && NonNegative[n] && OddQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&& OddQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	d_/;IntegerQ[d] && Positive[d],
	OptionsPattern[]]:=
Module[{subs=OptionValue["Tensor"]},
		Message[SymmetryFactorDiagram::nout]]


SymmetryFactorDiagram[
	n_/;IntegerQ[n] && NonNegative[n] && OddQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&& OddQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	OptionsPattern[]]:=
Module[{subs=OptionValue["Tensor"]},
		Message[SymmetryFactorDiagram::nout]]


(* ::Subsection:: *)
(*Graphs*)


ImportedEgdesQ[n_, o3_, o4_] :=
    Module[{edgesHead},
        edgesHead =
            Head @
                Which[
                    n == 0, zeropt1PI[o3, o4],
                    n == 2, twopt1PI[o3, o4],
                    n == 4, fourpt1PI[o3, o4]
                ];
        If[edgesHead === zeropt1PI || edgesHead === twopt1PI || edgesHead === fourpt1PI,
            False,
            True
        ]
    ]


ImportEdges[n_,o3_,o4_]:=Which[
	n==0,zeropt1PI[o3,o4]=importEdges1PI[
		FileNameJoin[{PacletManager`PacletResource["GSberveglieri/Phi4tools","EDGES"],ToString[n]<>"pt"<>ToString[o3]<>ToString[o4]<>".txt"}]
	];,
	n==2,twopt1PI[o3,o4]=importEdges1PI[
		FileNameJoin[{PacletManager`PacletResource["GSberveglieri/Phi4tools","EDGES"],ToString[n]<>"pt"<>ToString[o3]<>ToString[o4]<>".txt"}]
	];,
	n==4,fourpt1PI[o3,o4]=importEdges1PI[
		FileNameJoin[{PacletManager`PacletResource["GSberveglieri/Phi4tools","EDGES"],ToString[n]<>"pt"<>ToString[o3]<>ToString[o4]<>".txt"}]
	];]


VisualizeDiagram[
	n_/;IntegerQ[n] && NonNegative[n]&&EvenQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&EvenQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	d_/;IntegerQ[d] && Positive[d],
	OptionsPattern["Substitutions"->"Nothing"]]:=
Module[
	{subs=OptionValue["Substitutions"],tmpEdgesall,tmpEdges,tmpsub2},
	If[o3==0 && o4==0, Message[VisualizeDiagram::invalid],
	If[n=!=0&&n=!=2&&n=!=4,Message[VisualizeDiagram::nout],
	If[(o3>0 && o4+o3>8)||(o3==0 && o4>9),Message[VisualizeDiagram::nvert],
	If[n=!=4 && o3==0 && o4==1 && d==1, Message[VisualizeDiagram::tadpole],
	If[n==4 && o3==2 && o4==0 && d==1, Null,
	If[Not@ImportedEgdesQ[n,o3,o4],ImportEdges[n,o3,o4]];
	tmpEdgesall=Which[n==0,zeropt1PI[o3,o4],n==2,twopt1PI[o3,o4],n==4,fourpt1PI[o3,o4]];
	If[
		d>Length@tmpEdgesall,
		Message[VisualizeDiagram::dtoobig,d,Length@tmpEdgesall],
		tmpEdges=renameExternalPoints[tmpEdgesall[[d,2]]];
		Which[
			subs== "Nothing", visualizeGraphSpring@Flatten@tmpEdges,
			subs== "Sunsets", If[n==0 && o3==0 && o4==2 && d==1,visualizeGraph@{0\[UndirectedEdge]1,0\[UndirectedEdge]\!\(\*SubscriptBox[\("\<s\>"\), \({0, 1}\)]\),\!\(\*SubscriptBox[\("\<s\>"\), \({0, 1}\)]\)\[UndirectedEdge]1},visualizeGraphSpring@Flatten@subsjusts@tmpEdges],
			subs== "Analytics",If[n==0 && o3==0 && o4==2 && d==1,visualizeGraph@{0\[UndirectedEdge]1,0\[UndirectedEdge]\!\(\*SubscriptBox[\("\<s\>"\), \({0, 1}\)]\),\!\(\*SubscriptBox[\("\<s\>"\), \({0, 1}\)]\)\[UndirectedEdge]1},
				tmpsub2=visualizeGraph@replaceTad\[Tau]\[Beta]@subsBubTrianSquares@tmpEdges;
				If[
					loopnumber@tmpsub2<=2,
					tmpsub2,
					visualizeGraphSpring@replaceTad\[Tau]\[Beta]@subsBubTrianSquaresW[tmpEdges,invertedweights]
				]],
			And@@{subs=!="Nothing",subs=!="Sunsets",subs=!="Analytics"},Message[VisualizeDiagram::subs];visualizeGraphSpring@Flatten@tmpEdges]
	]]]]]]
]


VisualizeDiagram[
	n_/;IntegerQ[n] && NonNegative[n]&&EvenQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&EvenQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	OptionsPattern["Substitutions"->"Nothing"]]:=
Module[
	{subs=OptionValue["Substitutions"],sub2,sub2invW,effloops},
	If[o3==0 && o4==0, Message[VisualizeDiagram::invalid],
	If[n=!=0&&n=!=2&&n=!=4,Message[VisualizeDiagram::nout],
	If[(o3>0 && o4+o3>8)||(o3==0 && o4>9), Message[VisualizeDiagram::nvert],
	If[n=!=4 && o3==0 && o4==1, Message[VisualizeDiagram::tadpole],
	If[n==4 && o3==2 && o4==0, {},
	If[Not[ImportedEgdesQ[n,o3,o4]],ImportEdges[n,o3,o4]];
	Which[
		subs== "Nothing",visualizersp[n,o3,o4],
		subs== "Sunsets",If[n==0 && o3==0 && o4==2,{visualizeGraph@{0\[UndirectedEdge]1,0\[UndirectedEdge]\!\(\*SubscriptBox[\("\<s\>"\), \({0, 1}\)]\),\!\(\*SubscriptBox[\("\<s\>"\), \({0, 1}\)]\)\[UndirectedEdge]1}},visualizerspS[n,o3,o4]],
		subs== "Analytics",If[n==0 && o3==0 && o4==2,{visualizeGraph@{0\[UndirectedEdge]1,0\[UndirectedEdge]\!\(\*SubscriptBox[\("\<s\>"\), \({0, 1}\)]\),\!\(\*SubscriptBox[\("\<s\>"\), \({0, 1}\)]\)\[UndirectedEdge]1}},
			sub2=visualizerspSub2[n,o3,o4];
			sub2invW=visualizerspSub2W[n,o3,o4,invertedweights];
			effloops=loopnumber/@sub2;
			Table[If[effloops[[i]]>2,sub2invW[[i]],sub2[[i]]],{i,Length@sub2}]],
	And@@{subs=!="Nothing",subs=!="Sunsets",subs=!="Analytics"},Message[VisualizeDiagram::subs];visualizersp[n,o3,o4]]
	]]]]]
]


VisualizeDiagram::invalid = "Invalid arguments."


VisualizeDiagram::nout = "The number of external legs is not 0, 2 or 4. This(These) diagram(s) is(are) not in this package."


VisualizeDiagram::tadpole = "Within this renomalization scheme the diagrams with tadpoles have been removed because they are identically zero."


VisualizeDiagram::subs = "The options for \"Substitutions\" are: \"Nothing\", \"Sunsets\", and \"Analytics\". Defaulting to \"Tensor\" -> \"Nothing\"."


VisualizeDiagram::nvert = "Diagram(s) not in this package."


VisualizeDiagram::dtoobig = "Argument `1` is too big, take d\[LessEqual]`2`."


VisualizeDiagram[
	n_/;IntegerQ[n] && NonNegative[n]&&OddQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&OddQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	d_/;IntegerQ[d] && Positive[d],
	OptionsPattern["Substitutions"->"Nothing"]]:=
Module[
	{subs=OptionValue["Substitutions"],tmpEdgesall,tmpEdges,tmpsub2},
	Message[VisualizeDiagram::nout]
]


VisualizeDiagram[
	n_/;IntegerQ[n] && NonNegative[n]&&OddQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&OddQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	OptionsPattern["Substitutions"->"Nothing"]]:=
Module[
	{subs=OptionValue["Substitutions"],tmpEdgesall,tmpEdges,tmpsub2},
	Message[VisualizeDiagram::nout]
]


(* ::Subsection:: *)
(*Finding*)


(* ::Text:: *)
(*From label or graph find {n, {o3, o4}}*)


fromlabelCountnptorderB[string_]:=Module[{nn,vv,ll,o3},
nn=StringCount[string,"e"];
vv=StringCount[string,"|"];
ll=StringLength[string];
o3=nn+6vv-2ll;
(*Print[{nn,vv,ll}];*)
{nn,{o3,vv-o3}}
]


fromGraphCountnptorderB[graph_]:=Module[{vvnn,nn,vv,ll,o3},
vvnn=VertexList@graph;
ll=EdgeCount@graph;
nn=Count[VertexDegree@graph,1];

vv=Length@vvnn-nn;

o3=nn+4vv-2ll;
(*Print[{nn,vv,ll}];*)
{nn,{o3,vv-o3}}
]


(* ::Text:: *)
(*From graph find the isomorphic one in our lists (given {n, {o3, o4}}, in case we have not saved this list it is imported)*)


FinddGraphedges[grpp_,listt_]:=Module[{grp,list,nodplgrp,nodpllist,ll,matchsimpleList,trues,j,posTrue,lltrue,rulesIso,matchList,posMultTrue},
grp=EdgeList@grpp;
list=EdgeList/@listt;
nodplgrp=DeleteDuplicates@grp;
nodpllist=DeleteDuplicates/@list;
ll=Length@list;
matchsimpleList=Table[IsomorphicGraphQ[Graph[nodpllist[[i]]],Graph@nodplgrp],{i,ll}];
trues=Count[matchsimpleList,True];
If [trues==0,"Diagram not present",
posTrue=Flatten@Position[matchsimpleList,True];
(*Print[posTrue];*)
lltrue=Length@posTrue;
rulesIso=Table[Normal[FindGraphIsomorphism[Graph@nodplgrp,Graph[nodpllist[[posTrue[[j]]]]]]][[1]],{j,lltrue}];
(*Print[rulesIso];*)
(*Print[{Sort[Sort/@(grp/.rulesIso\[LeftDoubleBracket]1\[RightDoubleBracket])],Sort[Sort/@list\[LeftDoubleBracket]posTrue\[LeftDoubleBracket]1\[RightDoubleBracket]\[RightDoubleBracket]]}];*)
matchList=Table[Sort[Sort/@(grp/.rulesIso[[j]])]==Sort[Sort/@list[[posTrue[[j]]]]],{j,lltrue}];
posMultTrue=Flatten@Position[matchList,True];
posTrue[[posMultTrue]]
]
]


(* ::Text:: *)
(*From graph find {n, {o3, o4}, {d}}*)


(*FindGraph[graph_?GraphQ]:=Module[{no3o4,n,o3,o4,tmpEdgesall,d},
no3o4=fromGraphCountnptorderB[graph];
n=no3o4[[1]];{o3,o4}=no3o4[[2]];
If[n=!=0&&n=!=2&&n=!=4,"This diagram is not in this package",
If[OddQ[o3],"These diagrams do not exist",
If[o4>9||(o4==9&&o3>0),"This diagram is not in this package",
If[ImportedEgdesQ[n,o3,o4],Null;,ImportEdges[n,o3,o4];];
tmpEdgesall=Which[n==0,zeropt1PI[o3,o4],n==2,twopt1PI[o3,o4],n==4,fourpt1PI[o3,o4]][[All,2]];
d=FinddGraphedges[graph,tmpEdgesall];
{n,{o3,o4},d}
]]]
]*)


fromGraphCountv2[graph_]:=Module[{vvnn,nn,vv,ll,vdegrees,test,o3,o4},
vvnn=VertexList@graph;
ll=EdgeCount@graph;
vdegrees=SortBy[Tally@VertexDegree@graph,First];
nn=Count[VertexDegree@graph,1];
If[DeleteCases[vdegrees,{1,_}|{3,_}|{4,_}]=!={},test=0,test=1];
o3=Count[VertexDegree@graph,3];
o4=Count[VertexDegree@graph,4];
{test,{nn,{o3,o4}}}
]


FindGraph[graph_?GraphQ]:=Module[{assess,no3o4,n,o3,o4,tmpEdgesall,d},
assess=fromGraphCountv2[graph];
If[assess[[1]]==0,Message[InformationDiagram::typevert];,
no3o4=assess[[2]];
n=no3o4[[1]];{o3,o4}=no3o4[[2]];
If[n=!=0&&n=!=2&&n=!=4,Message[InformationDiagram::nout];,
If[(o3>0 && o4+o3>8)||(o3==0 && o4>9),Message[InformationDiagram::nvert];,
If[ImportedEgdesQ[n,o3,o4],Null;,ImportEdges[n,o3,o4];];
tmpEdgesall=Which[n==0,zeropt1PI[o3,o4],n==2,twopt1PI[o3,o4],n==4,fourpt1PI[o3,o4]][[All,2]];
d=FinddGraphedges[graph,tmpEdgesall];
{n,{o3,o4},d}
]]]
]


FindGraph::typevert = "There are vertices that are neither cubic nor quartic. This diagram is not in this package."


FindGraph::nout = "The number of external legs is not 0, 2 or 4. This diagram is not in this package."


FindGraph::nvert = "Diagram(s) not in this package."


(* ::Text:: *)
(*Draw diagram from Nickel index (also not minimal)*)


attachto[v1_,listforv1_]:=Map[UndirectedEdge[v1,#]&,listforv1]


DrawGraph[string_]:=Module[{str=string,connections,lcon,connectionsExpr,posext,lposext,test1,test2,ii=1,edgelist},
connections=StringSplit[str,"|",All][[1;;-2]];
lcon=Length@connections;
connectionsExpr=Table[Table[StringPart[connections[[j]],i],{i,StringLength[connections[[j]]]}],{j,lcon}];
posext=Position[connectionsExpr,"e",All];
lposext=Length@posext;
test1=Max[ToExpression/@DeleteCases[Flatten[connectionsExpr],"e"]];
test2=Min/@ToExpression/@DeleteCases[DeleteCases[connectionsExpr,"e",All],{}];
If[test1>Length@connectionsExpr-1,Message[DrawGraph::nstndrdstring];,
If[Or@@Table[test2[[jj]]<(jj-1),{jj,Length@test2}],Message[DrawGraph::nstndrdstring];,
While[ii<=lposext&&ii<10,
connectionsExpr=ReplacePart[connectionsExpr,posext[[ii]]->"e["<>ToString[ii]<>"]"];
ii++];
edgelist=Flatten@Table[attachto[ToString[i-1],connectionsExpr[[i]]],{i,lcon}];
visualizeGraphsimpleSpring[edgelist]]]
]


DrawGraph::nstndrdstring = "The input is not a standard string."


(* ::Text:: *)
(*From Label find {n, {o3, o4}, {d}}*)


(*FindIndex[string_?StringQ]:=Module[{no3o4,n,o3,o4,tmpEdgesall,d},
no3o4=fromlabelCountnptorderB[string];
n=no3o4[[1]];{o3,o4}=no3o4[[2]];
If[n=!=0&&n=!=2&&n=!=4,"This diagram is not in this package",
If[OddQ[o3],"These diagrams do not exist",
If[o4>9||(o4==9&&o3>0),"This diagram is not in this package",
If[ImportedIndicesQ[n,o3,o4],Null;,ImportIndices[n,o3,o4];];
(*tmpEdgesall=Which[n==0,zeropt1PI[o3,o4],n==2,twopt1PI[o3,o4],n==4,fourpt1PI[o3,o4]][[All,2]];*)
d=Position[nickelIndices[n,o3,o4],string][[1]];
{n,{o3,o4},d}
]]]
]*)


FindIndex[string_?StringQ]:=Module[{no3o4,n,o3,o4,tmpEdgesall,dall,d,out},
no3o4=fromlabelCountnptorderB[string];
n=no3o4[[1]];{o3,o4}=no3o4[[2]];
If[n=!=0&&n=!=2&&n=!=4,out=0;{out},
If[!(NonNegative[o3] && EvenQ[o3] && NonNegative[o4] && (o4<=8 || o4==9&&o3==0)),out=0;{out},
If[ImportedIndicesQ[n,o3,o4],Null;,ImportIndices[n,o3,o4];];
dall=Position[nickelIndices[n,o3,o4],string];
If[dall=={},out=0;{out},out=1;d=dall[[1]];
{n,{o3,o4},d}
]]]
]


(* ::Subsection:: *)
(*Integrands*)


changeNameMinus={bubbl[-x_]:> bubble[x],bubbl[-x_-y_]:> bubble[x+y],suns[-x_]:> sunset[x],Gt[-x_]:> prop[x],Gt[-x_-y_]:> prop[x+y]};


changeName={
	bubbl[x_] :> bubble[x],
	suns[x_] :> sunset[x],
	Gt[x_] :> prop[x],
	trinfnc[x_,y_,z_] :> triangle[x,y,z],
	squarfnc[x_,y_,z_,w_] :> square[x,y,z,w],
	suntad -> tadSunset,
	bubbl2tr -> tadTrianBub,
	pext -> p1,
	-pext -> - p1};


renameSymbols = {
	prop[x_] :> Propagator[x],
	bubble[x_] :> BubbleSubdiagram[x],
	sunset[x_] :> SunsetSubdiagram[x],
	triangle[x_,y_,z_] :> TriangleSubdiagram[x,y,z],
	square[x_,y_,z_,w_] :> SquareSubdiagram[x,y,z,w],
	tadSunset -> TadSunsetSubdiagram[],
	tadTrianBub -> TadTriangleBubblesSubdiagram[]
};


(*renameMomenta = Table[ToExpression[internalMomentumList[[k]]]->Momentum[k],{k,1,Length@internalMomentumList}]*)


(* ::Text:: *)
(*Extract the internal momentum variables "qi" from an expression*)


(*qvariables[expr_]:=Select[Variables@Level[expr,{-1}],StringStartsQ[ToString[#],"q"] || StringStartsQ[ToString[#],StringJoin[ToString@Context[qaux],"q"]] &]*)


(* ::Text:: *)
(*Substitute the internal momentum variables "qi" in an expression with Momentum[j] with minimal j *)


(*substituteQVariables[expr_] :=
    Module[{momvars = qvariables[expr], n, replaceList},
        n = Length[momvars];
        replaceList = Table[momvars[[i]] -> Momentum[i], {i, 1, n}];
        expr /. replaceList
    ]*)


substituteQVariables[expr_] :=
    Module[{replaceList},
        replaceList = Table[qaux[[i]] -> Momentum[i], {i, 1, Length@qaux}];
        expr /. replaceList
    ]


(* ::Text:: *)
(*Extract the external momentum variables  from an expression*)


(*pvariables[expr_]:=Select[Variables@Level[expr,{-1}], StringStartsQ[ToString[#],StringJoin[ToString@Context[pextlist],"p"]] &]*)


(* ::Text:: *)
(*Substitute the internal momentum variables "p+" in an expression with ExternalMomentum[j] with minimal j *)


(*substitutePVariables[expr_] :=
    Module[{momvars = pvariables[expr], n, replaceList},
        n = Length[momvars];
        replaceList = Table[momvars[[i]] -> ExternalMomentum[i], {i, 1, n}];
        expr /. replaceList
    ]*)


substitutePVariables[expr_] :=
    Module[{replaceList},
        replaceList = Table[pextlist[[i]] -> ExternalMomentum[i], {i, 1, Length@pextlist}];
        expr /. replaceList
    ]


(*renameExternalMomenta = Table[ToExpression[externalMomentumList[[k]]]->ExternalMomentum[k],{k,1,Length@externalMomentumList}]*)


ImportedIntegrandsQ[n_, o3_, o4_] := Module[{integrandsHead},
	If[o3 > 0,
		Return[False],
		integrandsHead = Head[integrand1PIno[n, o4]];
		If[ SameQ[integrandsHead, integrand1PIno],
			If[
				IntegrandFileExistsQ[n,o4],
				integrand1PIno[n, o4]=importtxtIntegrands[IntegrandFileName[n,o4]]; Return[True],
				Return[False]
			],
			Return[True]
		]
	]
]


IntegrandDiagram[
	n_/;IntegerQ[n] && NonNegative[n] && EvenQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&EvenQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	d_/;IntegerQ[d] && Positive[d],
	OptionsPattern[]]:=
Module[{subs=OptionValue["Substitutions"],tmpEdgesall,tmpFactor,tmpEdges,tmpfunc,tmpfuncnocomptad,tmpsub2},
	If[o3==0 && o4==0, Message[IntegrandDiagram::invalid],
	If[n=!=0&&n=!=2&&n=!=4,
		Message[IntegrandDiagram::nout],
		If[OptionValue["ExternalMomentum"]===True,
			If[n==2,
				IntegrandDiagramD[n,o3,o4,d,"Substitutions"->OptionValue["Substitutions"]],
				Message[IntegrandDiagram::extmom];
			],
			If[OptionValue["ExternalMomentum"]===False,Null;,Message[IntegrandDiagram::optextmom]];
			If[(o3>0 && o4+o3>8)||(o3==0 && o4>9),
				Message[IntegrandDiagram::nvert],
				If[n=!=4 && o3==0 && o4==1 && d==1,
					Message[IntegrandDiagram::tadpole],
					If[ImportedEgdesQ[n,o3,o4],
						Null;,
						ImportEdges[n,o3,o4];
					];
					tmpEdgesall=Which[
						n==0,zeropt1PI[o3,o4],
						n==2,twopt1PI[o3,o4],
						n==4,fourpt1PI[o3,o4]
					];
					If[d>Length@tmpEdgesall,
						Message[IntegrandDiagram::dtoobig,d,Length@tmpEdgesall],
						tmpFactor=tmpEdgesall[[d,1]];
						tmpEdges=tmpEdgesall[[d,2]];
						tmpfunc=Which[
							n==0, GtP0ptslevelgiven[o3,o4,{#}]&,
							n==2, GtP2ptslevelgiven[o3,o4,{#}]&,
							n==4, GtP4ptslevelgiven[o3,o4,{#}]&
						];
						tmpfuncnocomptad=Which[
							n==0, GtP0ptslevelgivennocompltad[o3,o4,{#}]&,
							n==2, GtP2ptslevelgivennocompltad[o3,o4,{#}]&,
							n==4, GtP4ptslevelgivennocompltad[o3,o4,{#}]&
						];
						Which[
							subs== "Nothing", tmpfunc[{tmpFactor,tmpEdges}][[1]]//substituteQVariables,
							subs== "Sunsets", If[n==0 && o3==0 && o4==2 && d==1,substituteQVariables[(1/48)prop[q1]sunset[q1]],
													tmpfuncnocomptad[{tmpFactor,subsjusts@tmpEdges}][[1]]//substituteQVariables],
							subs== "Analytics",If[n==0 && o3==0 && o4==2 && d==1,substituteQVariables[(1/48)prop[q1]sunset[q1]],
							substituteQVariables@If[ImportedIntegrandsQ[n,o3,o4],
								integrand1PIno[n,o4][[d,-1]],
								tmpsub2=visualizeGraph@replaceTad\[Tau]\[Beta]@subsBubTrianSquares@tmpEdges;
								If[loopnumber@tmpsub2<=2,
									(tmpfunc[{tmpFactor,subsBubTrianSquares@tmpEdges}])[[1]],
									tmpfunc[{tmpFactor,subsBubTrianSquaresW[tmpEdges,invertedweights]}][[1]]
								]
							]],
							And@@{subs=!="Nothing",subs=!="Sunsets",subs=!="Analytics"},
							Message[IntegrandDiagram::subs];
							tmpfunc[{tmpFactor,tmpEdges}][[1]]//substituteQVariables
						]
					]
				]
			]
		]/.changeNameMinus/.changeName // substitutePVariables
	]]/.renameSymbols
]


IntegrandDiagram[
	n_/;IntegerQ[n] && NonNegative[n]&&EvenQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&EvenQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],OptionsPattern[]]:=
Module[{subs=OptionValue["Substitutions"],tmpEdgesfull,tmpEdgesfullS,tmpgrph,effloops,tmpfunc,tmpfuncnocomptad,tmpsub2,tmpfuncSub2,tmpfuncSub2W},
	Map[substitutePVariables,
	If[o3==0 && o4==0, Message[IntegrandDiagram::invalid],
	If[n=!=0&&n=!=2&&n=!=4,
		Message[IntegrandDiagram::nout],
		If[OptionValue["ExternalMomentum"]===True,
			If[n==2,
				IntegrandDiagramsD[n,o3,o4,"Substitutions"->OptionValue["Substitutions"]],
				Message[IntegrandDiagram::extmom];
			],
			If[OptionValue["ExternalMomentum"]===False,Null;,Message[IntegrandDiagram::optextmom]];
			If[(o3>0 && o4+o3>8)||(o3==0 && o4>9),
				Message[IntegrandDiagram::nvert],
				If[n=!=4 && o3==0 && o4==1,
					Message[IntegrandDiagram::tadpole],
					If[ImportedEgdesQ[n,o3,o4],
						Null;,
						ImportEdges[n,o3,o4];
					];
					tmpEdgesfull=Which[
						n==0, zeropt1PI[o3,o4],
						n==2, twopt1PI[o3,o4],
						n==4, fourpt1PI[o3,o4]
					];
					tmpfunc=Which[
						n==0, GtP0ptslevelgiven[o3,o4,#]&,
						n==2, GtP2ptslevelgiven[o3,o4,#]&,
						n==4, GtP4ptslevelgiven[o3,o4,#]&
					];
					tmpfuncnocomptad=Which[
						n==0, GtP0ptslevelgivennocompltad[o3,o4,#]&,
						n==2, GtP2ptslevelgivennocompltad[o3,o4,#]&,
						n==4, GtP4ptslevelgivennocompltad[o3,o4,#]&
					];
					Which[
						subs=="Nothing", Map[substituteQVariables,tmpfunc[tmpEdgesfull]],
						subs== "Sunsets",If[n==0 && o3==0 && o4==2,{substituteQVariables[(1/48)prop[q1]sunset[q1]]},
										tmpEdgesfullS=Which[n==0,zeropt1PIs[o3,o4],n==2,twopt1PIs[o3,o4],n==4,fourpt1PIs[o3,o4]];
										Map[substituteQVariables,tmpfuncnocomptad[tmpEdgesfullS]]],
						subs== "Analytics",If[n==0 && o3==0 && o4==2,{substituteQVariables[(1/48)prop[q1]sunset[q1]]},
						Map[substituteQVariables,
						If[ImportedIntegrandsQ[n,o3,o4],
							integrand1PIno[n,o4][[All,-1]],
							tmpgrph=visualizerSub2[n,o3,o4];
							effloops=loopnumber/@visualizerSub2[n,o3,o4];
							tmpfuncSub2=Which[n==0,GtP0ptslevel2[o3,o4],n==2,GtP2ptslevel2[o3,o4],n==4,GtP4ptslevel2[o3,o4]];
							tmpfuncSub2W=Which[n==0,GtP0ptslevel2W[o3,o4,invertedweights],n==2,GtP2ptslevel2W[o3,o4,invertedweights],n==4,GtP4ptslevel2W[o3,o4,invertedweights]];
							Table[If[effloops[[i]]<=2,tmpfuncSub2[[i]],tmpfuncSub2W[[i]]],{i,Length@tmpgrph}]]	
						]],
						And@@{subs=!="Nothing",subs=!="Sunsets",subs=!="Analytics"},
							Message[IntegrandDiagram::subs];
							Map[substituteQVariables,tmpfunc[tmpEdgesfull]]
					]
				]
			]
		]/.changeNameMinus/.changeName
	]]/.renameSymbols
	]
]


IntegrandDiagram::invalid = "Invalid arguments."


IntegrandDiagram::nout = "The number of external legs is not 0, 2 or 4. This(These) diagram(s) is(are) not in this package."


(*IntegrandDiagram::odd = "This(These) diagram(s) does(do) not exist."*)


IntegrandDiagram::extmom = "ExternalMomentum option is only available for n=2."


IntegrandDiagram::tadpole = "Within this renomalization scheme the diagrams with tadpoles have been removed because they are identically zero."


IntegrandDiagram::subs = "The options for \"Substitutions\" are: \"Nothing\", \"Sunsets\", and \"Analytics\". Defaulting to \"Substitutions\" -> \"Nothing\"."


IntegrandDiagram::optextmom = "The options for \"ExternalMomentum\" are: True or False. Defaulting to \"ExternalMomentum\" -> False."


IntegrandDiagram::nvert = "Diagram(s) not in this package."


IntegrandDiagram::dtoobig = "Argument `1` is too big, take d\[LessEqual]`2`."


IntegrandDiagram[
	n_/;IntegerQ[n] && NonNegative[n] && OddQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&OddQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	d_/;IntegerQ[d] && Positive[d],
	OptionsPattern[]]:=
Module[{subs=OptionValue["Substitutions"],tmpEdgesall,tmpFactor,tmpEdges,tmpfunc,tmpfuncnocomptad,tmpsub2},
		Message[IntegrandDiagram::nout]
]


IntegrandDiagram[
	n_/;IntegerQ[n] && NonNegative[n] && OddQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&OddQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	OptionsPattern[]]:=
Module[{subs=OptionValue["Substitutions"],tmpEdgesall,tmpFactor,tmpEdges,tmpfunc,tmpfuncnocomptad,tmpsub2},
		Message[IntegrandDiagram::nout]
]


(* ::Subsection:: *)
(*Integrands \[CapitalGamma]^(2)'*)


IntegrandDiagramsD[
	n_/;IntegerQ[n] && NonNegative[n]&&EvenQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	OptionsPattern["Substitutions"->"Nothing"]]:=
Module[
	{subs=OptionValue["Substitutions"],tmpEdgesfull,tmpEdgesfullS,tmpgrph,effloops,tmpfunc,tmpsub2,tmpfuncSub2,tmpfuncSub2W},
	Map[substitutePVariables, 
	If[n=!= 2,"These integrands are just for \!\(\*SuperscriptBox[\(\[CapitalGamma]\), \((2)\)]\)'",
	If[OddQ[o3]||o3<0||o4<0,"This diagram does not exist",
	If[(o3>0 && o4+o3>8)||(o3==0 && o4>9), Message[IntegrandDiagram::nvert],
	If[n=!=4 && o3==0 && o4==1, Message[IntegrandDiagram::tadpole],
	If[ImportedEgdesQ[n,o3,o4],Null;,ImportEdges[n,o3,o4];];
	tmpEdgesfull=twopt1PI[o3,o4];
	(*tmpfunc=GtP2ptslevelgivender[o3,o4,#]&;*)
	
	Which[
		subs== "Nothing", Map[substituteQVariables, GtP2ptslevelgivender[o3,o4,tmpEdgesfull]],
		subs== "Sunsets", tmpEdgesfullS=twopt1PIs[o3,o4]; Map[substituteQVariables, GtP2ptslevelgivendernocompltad[o3,o4,tmpEdgesfullS]],
		subs== "Analytics",
			tmpgrph=visualizerSub2[n,o3,o4];
			effloops=loopnumber/@visualizerSub2[n,o3,o4];
			tmpfuncSub2=GtP2ptslevel2v1der[o3,o4];
			Map[substituteQVariables,Table[rassPext[effloops[[i]],tmpfuncSub2[[i]]],{i,Length@tmpEdgesfull}]],
		And@@{subs=!="Nothing",subs=!="Sunsets",subs=!="Analytics"},
			Message[IntegrandDiagram::subs];
			Map[substituteQVariables, GtP2ptslevelgivender[o3,o4,tmpEdgesfull]]
	]
]]]]/.changeNameMinus/.changeName/.changeName
]
]


IntegrandDiagramD[
	n_/;IntegerQ[n] && NonNegative[n]&&EvenQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	d_/;IntegerQ[d] && Positive[d],
	OptionsPattern["Substitutions"->"Nothing"]]:=
Module[
	{subs=OptionValue["Substitutions"],tmpFactor,tmpEdges,tmpgrph,effloops,tmpfunc,tmpsub2,tmpfuncSub2,tmpfuncSub2W},
	If[n=!= 2,"These integrand are just for \!\(\*SuperscriptBox[\(\[CapitalGamma]\), \((2)\)]\)'",
	If[OddQ[o3]||o3<0||o4<0,"These diagrams do not exist",
	If[(o3>0 && o4+o3>8)||(o3==0 && o4>9),Message[IntegrandDiagram::nvert],
	If[n=!=4 && o3==0 && o4==1 && d==1,Message[IntegrandDiagram::tadpole],
	If[ImportedEgdesQ[n,o3,o4],Null;,ImportEdges[n,o3,o4];];
	If[d>Length@twopt1PI[o3,o4],Message[IntegrandDiagram::dtoobig,d,Length@twopt1PI[o3,o4]],
	tmpFactor=twopt1PI[o3,o4][[d,1]];
	tmpEdges=twopt1PI[o3,o4][[d,2]];

	Which[
		subs== "Nothing",
			GtP2ptslevelgivender[o3,o4,{{tmpFactor,tmpEdges}}][[1]]//substituteQVariables,
		subs== "Sunsets",
			GtP2ptslevelgivendernocompltad[o3,o4,{{tmpFactor,subsjusts@tmpEdges}}][[1]]//substituteQVariables,
		subs== "Analytics",
			tmpsub2=visualizeGraph@replaceTad\[Tau]\[Beta]@subsBubTrianSquares@tmpEdges;
			tmpfuncSub2=GtP2ptslevel2v1der[o3,o4];
			substituteQVariables[ 
				rassPext[loopnumber@tmpsub2,GtP2ptslevelgivender[o3,o4,{{tmpFactor,subsBubTrianSquares@tmpEdges}}][[1]]]],
		And@@{subs=!="Nothing",subs=!="Sunsets",subs=!="Analytics"},
			Message[IntegrandDiagram::subs];
			GtP2ptslevelgivender[o3,o4,{{tmpFactor,tmpEdges}}][[1]]//substituteQVariables
	]
]]]]]/.changeNameMinus/.changeName/.changeName // substitutePVariables
]


(* ::Subsection:: *)
(*Values*)


ValueDiagram[
	n_/;IntegerQ[n] && NonNegative[n] && EvenQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&EvenQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	d_/;IntegerQ[d] && Positive[d],
	OptionsPattern[]]:=
Module[{allvalues},
	If[o3==0 && o4==0, Message[ValueDiagram::invalid],
	If[n=!=0&&n=!=2&&n=!=4,
		Message[ValueDiagram::nout],
		If[OptionValue["Derivative"],
			If[n==2,
				ValueDiagramD[n,o3,o4,d],
				Message[ValueDiagram::extmom]
			],
			If[o4>8||o3=!=0,
				Return[Missing["Not computed"]],
				If[n=!=4 && o3==0 && o4==1 && d==1,
					Message[ValueDiagram::tadpole],
					If[n==0 && o3==0 && o4>1 && o4<4 && d==1, 0,
						If[n==0 && o3==0 && o4>1 && o4<4 && d>1,
							Message[ValueDiagram::dtoobig,d,1],
							allvalues=resultGIno[ToString@n,o4];
							If[d>Length@allvalues,
								Message[ValueDiagram::dtoobig,d,Length@allvalues],
								Return[allvalues[[d,2]]]
							]
						]
					]
				]
			]
		]
	]
	]
]


ValueDiagramD[
	n_/;IntegerQ[n] && NonNegative[n] && EvenQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	d_/;IntegerQ[d] && Positive[d]]:=
Module[{allvalues},
	If[OddQ[o3]||o3<0||o4<0,
		Return[Missing["This diagram does not exist"]],
		If[n=!= 2,
			Return[Missing["These values are just for \!\(\*SuperscriptBox[\(\[CapitalGamma]\), \((2)\)]\)'"]],
			If[o4>8||o3=!=0,
				Return[Missing["Not computed"]],
				allvalues=resultGIno["2D",o4];
				If[d>Length@allvalues,
					Message[ValueDiagram::dtoobig,d,Length@allvalues],
					Return[allvalues[[d,2]]]
				]
			]
		]
	]
]


ValueDiagram[
	n_/;IntegerQ[n] && NonNegative[n] && EvenQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&EvenQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	OptionsPattern[]]:=
Module[{allvalues},
	If[o3==0 && o4==0, Message[ValueDiagram::invalid],
		If[n=!=0&&n=!=2&&n=!=4,
			Message[ValueDiagram::nout],
			If[Not[MemberQ[{True,False},OptionValue["Derivative"]]],
				Message[ValueDiagram::deriv],
				If[OptionValue["Derivative"],
					If[n==2,
						ValueDiagramsD[n,o3,o4],
						Message[ValueDiagram::extmom]
					],
					If[o4>8||o3=!=0,
						Return[Missing["Not computed"]],
						If[n=!=4 && o3==0 && o4==1,
							Message[ValueDiagram::tadpole],
							If[n==0 && o3==0 && o4>1 && o4<4, 
								{0},
								allvalues=resultGIno[ToString@n,o4][[All,2]];
								Return[allvalues]
							]
						]
					]
				]
			]
		]
	]
]


ValueDiagramsD[
	n_/;IntegerQ[n] && NonNegative[n] && EvenQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4]]:=
Module[{allvalues},
	If[OddQ[o3]||o3<0||o4<0,
		Return[Missing["These diagrams do not exist"]],
		If[n=!=2,
			Return[Missing["These values are just for \!\(\*SuperscriptBox[\(\[CapitalGamma]\), \((2)\)]\)'"]],
			If[o4>8||o3=!=0,
				Return[Missing["Not computed"]],
				allvalues=resultGIno["2D",o4][[All,2]];
				Return[allvalues]
			]
		]
	]
]


ValueDiagram::invalid = "Invalid arguments."


ValueDiagram::nout = "The number of external legs is not 0, 2 or 4. This(These) diagram(s) is(are) not in this package."


ValueDiagram::extmom = "Derivative option only available for n=2."


ValueDiagram::dtoobig = "Argument `1` is too big, take d\[LessEqual]`2`."


ValueDiagram::tadpole = "Within this renomalization scheme the diagrams with tadpoles have been removed because they are identically zero."


ValueDiagram::nvert = "Diagram(s) not in this package."


ValueDiagram::deriv = "Available option values for \"Derivative\" are: True, False."


ValueDiagram[
	n_/;IntegerQ[n] && NonNegative[n] && OddQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&OddQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	d_/;IntegerQ[d] && Positive[d],
	OptionsPattern[]]:=
Module[{allvalues},
		Message[ValueDiagram::nout]
]


ValueDiagram[
	n_/;IntegerQ[n] && NonNegative[n] && OddQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&OddQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	OptionsPattern[]]:=
Module[{allvalues},
		Message[ValueDiagram::nout]
]


(* ::Subsection:: *)
(*Composition*)


InformationDiagramV0[n_/;IntegerQ[n] && NonNegative[n]&&EvenQ[n],o3_/;IntegerQ[o3] && NonNegative[o3],o4_/;IntegerQ[o4] && NonNegative[o4],d_/;IntegerQ[d] && Positive[d],OptionsPattern["Substitutions"->"Nothing"]]:=Module[{subs=OptionValue["Substitutions"],nicklab,visdiag,intdiag},
nicklab=NickelIndex[n,o3,o4,d];
visdiag=VisualizeDiagram[n,o3,o4,d,"Substitutions"->subs];
intdiag=IntegrandDiagram[n,o3,o4,d,"Substitutions"->subs];
{{n,{o3,o4},{d}},nicklab,visdiag,intdiag}
]


(* ::Input:: *)
(*(*InformationDiagramIntegrand[n_,o3_,o4_,d_,OptionsPattern["Substitutions"->"Nothing"]]:=Module[{subs=OptionValue["Substitutions"],nicklab,visdiag,intdiag},*)
(*nicklab=NickelIndex[n,o3,o4,d];*)
(*visdiag=VisualizeDiagram[n,o3,o4,d,"Substitutions"->subs];*)
(*intdiag=IntegrandDiagram[n,o3,o4,d,"Substitutions"->subs];*)
(*If[subs==="Nothing",{{n,{o3,o4},{d}},nicklab,visdiag,intdiag}, {{n,{o3,o4},{d}},nicklab,VisualizeDiagram[n,o3,o4,d],visdiag,intdiag}]*)
(*]*)*)


InformationDiagramIntegrand[n_/;IntegerQ[n] && NonNegative[n]&&EvenQ[n],o3_/;IntegerQ[o3] && NonNegative[o3],o4_/;IntegerQ[o4] && NonNegative[o4],d_/;IntegerQ[d] && Positive[d],OptionsPattern["Substitutions"->"Nothing"]]:=Module[{subs=OptionValue["Substitutions"],nicklab,visdiag,intdiag},
If[OddQ[o3],"This diagram does not exist",
If[o4>9||(o4==9&&o3>0),"This diagram is not in this package",
If[d>Length@NickelIndex[n,o3,o4],"This diagram does not exist, take d\[LessEqual]"<>ToString[Length@NickelIndex[n,o3,o4]],
nicklab=NickelIndex[n,o3,o4,d];
visdiag=VisualizeDiagram[n,o3,o4,d,"Substitutions"->subs];
intdiag=IntegrandDiagram[n,o3,o4,d,"Substitutions"->subs];
If[subs==="Nothing",{{n,{o3,o4},{d}},nicklab,visdiag,intdiag}, {{n,{o3,o4},{d}},nicklab,VisualizeDiagram[n,o3,o4,d],visdiag,intdiag}]
]]]
]


(* ::Input:: *)
(*(*InformationDiagram[n_,o3_,o4_,d_]:=Module[{nicklab,visdiag,valdiag},*)
(*nicklab=NickelIndex[n,o3,o4,d];*)
(*visdiag=VisualizeDiagram[n,o3,o4,d];*)
(*valdiag=ValueDiagram[n,o3,o4,d];*)
(*If[o4>8||o3=!=0,{{n,{o3,o4},{d}},nicklab,visdiag},{{n,{o3,o4},{d}},nicklab,visdiag,valdiag}]*)
(*]*)*)


(* ::Input:: *)
(*(*InformationDiagramAll[n_,o3_,o4_,d_,OptionsPattern["Substitutions"->"Nothing"]]:=Module[{subs=OptionValue["Substitutions"],nicklab,visdiag,intdiag,valdiag},*)
(*nicklab=NickelIndex[n,o3,o4,d];*)
(*visdiag=VisualizeDiagram[n,o3,o4,d,"Substitutions"->subs];*)
(*intdiag=IntegrandDiagram[n,o3,o4,d,"Substitutions"->subs];*)
(*valdiag=ValueDiagram[n,o3,o4,d];*)
(*If[subs==="Nothing",If[o4>8||o3=!=0,{{n,{o3,o4},{d}},nicklab,visdiag,intdiag},{{n,{o3,o4},{d}},nicklab,visdiag,intdiag,valdiag}], If[o4>8||o3=!=0,{{n,{o3,o4},{d}},nicklab,VisualizeDiagram[n,o3,o4,d],visdiag,intdiag},{{n,{o3,o4},{d}},nicklab,VisualizeDiagram[n,o3,o4,d],visdiag,intdiag,valdiag}]]*)
(*]*)*)
(**)
(*(*InformationDiagram[*)
(*	n_/;IntegerQ[n] && NonNegative[n] && EvenQ[n],*)
(*	o3_/;IntegerQ[o3] && NonNegative[o3],*)
(*	o4_/;IntegerQ[o4] && NonNegative[o4],*)
(*	d_/;IntegerQ[d] && Positive[d],*)
(*	OptionsPattern[]]:=*)
(*Module[{subs=OptionValue["Substitutions"],substens=OptionValue["Tensor"],nicklab,factdiag,visdiag,intdiag,valdiag,ans},*)
(*	If[OddQ[o3], "This diagram does not exist",*)
(*		If[o4>9||(o4>=9 && o3>0),*)
(*			"This diagram is not in this package",*)
(*			If[d>Length@NickelIndex[n,o3,o4],*)
(*				"This diagram does not exist, take d\[LessEqual]"<>ToString[Length@NickelIndex[n,o3,o4]],*)
(*				nicklab=NickelIndex[n,o3,o4,d];*)
(*				visdiag=VisualizeDiagram[n,o3,o4,d,"Substitutions"->subs];*)
(*				intdiag=IntegrandDiagram[n,o3,o4,d,"Substitutions"->subs];*)
(*				valdiag=ValueDiagram[n,o3,o4,d];*)
(*				ans={{n,{o3,o4},{d}},nicklab};*)
(*				If[o4<9 && o3==0 && substens=!="Identity",*)
(*					factdiag=SymmetryFactorDiagram[n,o3,o4,d,"Tensor"->substens];*)
(*					AppendTo[ans,factdiag];*)
(*				];*)
(*				AppendTo[ans,VisualizeDiagram[n,o3,o4,d]];*)
(*				If[subs=!="Nothing",*)
(*					AppendTo[ans,visdiag];*)
(*				];*)
(*				If[OptionValue["ShowIntegrand"],*)
(*					AppendTo[ans,intdiag];*)
(*				];*)
(*				If[o4<9 && o3==0,*)
(*					AppendTo[ans,valdiag];*)
(*				];*)
(*				Return[ans]*)
(*			]*)
(*		]*)
(*	]*)
(*]*)*)


InformationDiagram[
	n_/;IntegerQ[n] && NonNegative[n] && EvenQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&EvenQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	d_/;IntegerQ[d] && Positive[d],
	OptionsPattern[]]:=
Module[{subs=OptionValue["Substitutions"],substens=OptionValue["Tensor"],nicklab,factdiag,intdiag,valdiag,ans},
	If[n=!=0&&n=!=2&&n=!=4,
		Message[InformationDiagram::nout],
		If[n=!=4 && o3==0 && o4==1,
			Message[InformationDiagram::tadpole],
			If[(o3>0 && o4+o3>8)||(o3==0 && o4>9),
				Message[InformationDiagram::nvert],
				If[d>Length@NickelIndex[n,o3,o4],
					Message[InformationDiagram::dtoobig,d,Length@NickelIndex[n,o3,o4]],
					nicklab=NickelIndex[n,o3,o4,d];
					(*ans=<|"Element"->ToString[{n,{o3,o4},{d}}],"Index"->nicklab|>;*)
					ans=<|"External Legs"-> n, "Cubic Vertices"->o3, "Quartic Vertices"->o4, "List Number"->d,"Nickel Index"->nicklab|>;
					(*If[o4<9 && o3==0,*)
					Which[
						substens==="Identity", AppendTo[ans,"Symmetry Factor"->1],
						substens==="O(N)"||substens==="Cubic",
						factdiag=SymmetryFactorDiagram[n,o3,o4,d,"Tensor"->substens];
						AppendTo[ans,"Symmetry Factor"->factdiag],
						True, None
					];
					AppendTo[ans,"Diagram"->VisualizeDiagram[n,o3,o4,d]];
					If[Not[MemberQ[{"Nothing","Sunsets","Analytics"},subs]],
						subs="Nothing"; Message[InformationDiagram::subs]];
					If[subs==="Sunsets"||subs==="Analytics",
						AppendTo[ans,"Simplified Diagram"->VisualizeDiagram[n,o3,o4,d,"Substitutions"->subs]];
					];
					If[OptionValue["ShowIntegrand"],
						intdiag=IntegrandDiagram[n,o3,o4,d,"Substitutions"->subs];
						AppendTo[ans,"Integrand"->intdiag];
					];
					(*If[o4<9 && o3==0,*)
					valdiag=ValueDiagram[n,o3,o4,d];
					AppendTo[ans,"Value"->valdiag](*,AppendTo[ans,"Value"->Missing[]]]*);
					Return[Dataset[ans]]
				]
			]
		]
	]
]


(*InformationDiagram[x_Graph]:=InformationDiagram@@(Quiet[Flatten[FindGraph[x]],{Flatten::normal}])*)


InformationDiagram[x_Graph,OptionsPattern[]]:=Module[{subs=OptionValue["Substitutions"],substens=OptionValue["Tensor"]},
InformationDiagram@@(Quiet[Join[(Quiet[Flatten[FindGraph[x]],{Flatten::normal}]),{"Substitutions"->subs,"Tensor"->substens,"ShowIntegrand"->OptionValue["ShowIntegrand"]}],{Join::heads}])
]


(*InformationDiagram[x_String]:=Module[{look,tmp},
look=FindIndex[x];
tmp=If[look=={0},FindGraph[DrawGraph[x]],look];
InformationDiagram@@(Quiet[Flatten[tmp],{Flatten::normal}])
]*)


InformationDiagram[x_String,OptionsPattern[]]:=Module[{look,tmp,subs=OptionValue["Substitutions"],substens=OptionValue["Tensor"]},
look=FindIndex[x];
tmp=If[look=={0},FindGraph[DrawGraph[x]],look];
InformationDiagram@@(Quiet[Join[(Quiet[Flatten[tmp],{Flatten::normal}]),{"Substitutions"->subs,"Tensor"->substens,"ShowIntegrand"->OptionValue["ShowIntegrand"]}],{Join::heads}])
]


(* ::Input:: *)
(*(*InformationDiagramAll[n_,o3_,o4_,d_,OptionsPattern[{"Substitutions"->"Nothing","Tensor"->"O(N)"}]]:=Module[{subs=OptionValue["Substitutions"],substens=OptionValue["Tensor"],nicklab,factdiag,visdiag,intdiag,valdiag},*)
(*nicklab=NickelIndex[n,o3,o4,d];*)
(*factdiag=SymmetryFactorDiagram[n,o3,o4,d,"Tensor"->substens];*)
(*visdiag=VisualizeDiagram[n,o3,o4,d,"Substitutions"->subs];*)
(*intdiag=IntegrandDiagram[n,o3,o4,d,"Substitutions"->subs];*)
(*valdiag=ValueDiagram[n,o3,o4,d];*)
(*If[subs==="Nothing",If[o4>8||o3=!=0,{{n,{o3,o4},{d}},nicklab,visdiag,intdiag},{{n,{o3,o4},{d}},nicklab,visdiag,factdiag,intdiag,valdiag}],If[o4>8||o3=!=0,{{n,{o3,o4},{d}},nicklab,VisualizeDiagram[n,o3,o4,d],visdiag,intdiag},{{n,{o3,o4},{d}},nicklab,VisualizeDiagram[n,o3,o4,d],visdiag,factdiag,intdiag,valdiag}]]*)
(*]*)*)


(*InformationDiagramAll[
	n_/;IntegerQ[n] && NonNegative[n] && EvenQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	d_/;IntegerQ[d] && Positive[d],
	OptionsPattern[]]:=
Module[{subs=OptionValue["Substitutions"],substens=OptionValue["Tensor"],nicklab,factdiag,visdiag,intdiag,valdiag},
	If[OddQ[o3],
		"This diagram does not exist",
		If[o4>9||(o4==9&&o3>0),
			"This diagram is not in this package",
			If[d>Length@NickelIndex[n,o3,o4],
				"This diagram does not exist, take d\[LessEqual]"<>ToString[Length@NickelIndex[n,o3,o4]],
				nicklab=NickelIndex[n,o3,o4,d];
				factdiag=SymmetryFactorDiagram[n,o3,o4,d,"Tensor"->substens];
				visdiag=VisualizeDiagram[n,o3,o4,d,"Substitutions"->subs];
				intdiag=IntegrandDiagram[n,o3,o4,d,"Substitutions"->subs];
				valdiag=ValueDiagram[n,o3,o4,d];
				If[subs==="Nothing",
					If[o4>8||o3=!=0,
						{{n,{o3,o4},{d}},nicklab,visdiag,intdiag},
						{{n,{o3,o4},{d}},nicklab,visdiag,factdiag,intdiag,valdiag}
					],
					If[o4>8||o3=!=0,
						{{n,{o3,o4},{d}},nicklab,VisualizeDiagram[n,o3,o4,d],visdiag,intdiag},
						{{n,{o3,o4},{d}},nicklab,VisualizeDiagram[n,o3,o4,d],visdiag,factdiag,intdiag,valdiag}
					]
				]
			]
		]
	]
]*)


InformationDiagram::nout = "The number of external legs is not 0, 2 or 4. This(These) diagram(s) is(are) not in this package."


InformationDiagram::nvert = "Diagram(s) not in this package."


InformationDiagram::dtoobig = "Argument `1` is too big, take d\[LessEqual]`2`."


InformationDiagram::typevert = "There are vertices that are neither cubic nor quartic. This diagram is not in this package."


InformationDiagram::subs = "The options for \"Substitutions\" are: \"Nothing\", \"Sunsets\", and \"Analytics\". Defaulting to \"Substitutions\" -> \"Nothing\"."


InformationDiagram::tadpole = "Within this renomalization scheme the diagrams with tadpoles have been removed because they are identically zero."


InformationDiagram[
	n_/;IntegerQ[n] && NonNegative[n] && OddQ[n],
	o3_/;IntegerQ[o3] && NonNegative[o3]&&OddQ[o3],
	o4_/;IntegerQ[o4] && NonNegative[o4],
	d_/;IntegerQ[d] && Positive[d],
	OptionsPattern[]]:=
Module[{subs=OptionValue["Substitutions"],substens=OptionValue["Tensor"],nicklab,factdiag,intdiag,valdiag,ans},
		Message[InformationDiagram::nout]
]


(* ::Section:: *)
(*Last substitution*)


(* ::Subsection:: *)
(*Values of analytic functions and substitutions*)


(* ::Subsubsection::Closed:: *)
(*Analytic Simple pieces*)


Prop[x_]:=1/(x^2+1)


BubbleV[x_]:=If[x=!=0 ,ArcTan[x/2]/(4\[Pi] x),1/(8 \[Pi])]


SunsetV[x_]:=If[x=!=0 ,-(1/(32 \[Pi]^2)) ((6 ArcTan[x/3])/x+Log[1+x^2/9]-2),0]


tadSunsetValue=-(Log[4/3]/(128 \[Pi]^3));tadTrianBubValue=1/(12288 \[Pi]^2);


AreaTriangleSquaredtimes4[k1_,k2_,k3_]:=(k1+k2-k3) (k1-k2+k3) (-k1+k2+k3) (k1+k2+k3)


TriangleV[k1_,k2_,k3_]:=If[k1=!=0&&k2=!=0&&k3=!=0,1/(4\[Pi]) ArcTan[Sqrt[k1^2 k2^2 k3^2+AreaTriangleSquaredtimes4[k1,k2,k3]]/(8+ k1^2+k2^2+k3^2)]/Sqrt[k1^2 k2^2 k3^2+AreaTriangleSquaredtimes4[k1,k2,k3]],If[k1===0,1/(8\[Pi] (k2^2+4)),1/(8\[Pi] (k1^2+4))]]


(* ::Subsubsection::Closed:: *)
(*Vectors*)


vectorNorm3d[p_]:=Module[{arg,tmp,tmpsq},arg=Piecewise[{{{p},Head@p===Symbol||Head@p==Times},{List@@Expand[p],Head@Expand[p]===Plus}}];If[Length@arg>1,tmp=Plus@@(vec3d/@(Expand[p]/.Plus->List));tmpsq=Sqrt[tmp . tmp],If[Length@arg==1,tmpsq=Sqrt[vec3d[p] . vec3d[p]]]];If[p===0,0,Simplify@Expand[tmpsq/.{Subscript[a_, 0]:>\!\(\*SubscriptBox[\(a\), \("\<\[Rho]\>"\)]\)Sin[\!\(\*SubscriptBox[\(a\), \("\<\[Theta]\>"\)]\)]Cos[\!\(\*SubscriptBox[\(a\), \("\<\[Phi]\>"\)]\)],Subscript[a_, 1]:>\!\(\*SubscriptBox[\(a\), \("\<\[Rho]\>"\)]\)Sin[\!\(\*SubscriptBox[\(a\), \("\<\[Phi]\>"\)]\)]Sin[\!\(\*SubscriptBox[\(a\), \("\<\[Theta]\>"\)]\)],Subscript[a_, 2]:>\!\(\*SubscriptBox[\(a\), \("\<\[Rho]\>"\)]\)Cos[\!\(\*SubscriptBox[\(a\), \("\<\[Theta]\>"\)]\)]}]]]


vectornameComponents={Subscript[x_, 3]:>\!\(\*SubscriptBox[\(x\), \("\<\[Rho]\>"\)]\),Subscript[x_, 4]:>\!\(\*SubscriptBox[\(x\), \("\<\[Theta]\>"\)]\),Subscript[x_, 5]:>\!\(\*SubscriptBox[\(x\), \("\<\[Phi]\>"\)]\)};


(* ::Subsubsection:: *)
(*Applications*)


substitutionAnalytical={prop[x_]:> Prop@vectorNorm3d@x,bubble[x_]:> BubbleV@vectorNorm3d@x,sunset[x_]:> SunsetV@vectorNorm3d@x,triangle[x_,y_,z_]:> TriangleV@@vectorNorm3d/@{x,y,z},tadSunset-> tadSunsetValue,tadTrianBub-> tadTrianBubValue};


(* ::Text:: *)
(**)


substitutionSquare={square[x_,y_,z_,w_]:> squareV@@{x,y,z,w}};


substitutionVectors={prop[x_]:> propag@vectorNorm3d@x,bubble[x_]:> bubbleV@vectorNorm3d@x,sunset[x_]:> sunsetV@vectorNorm3d@x,triangle[x_,y_,z_]:> triangleV@@vectorNorm3d/@{x,y,z},
squareNo0[x_,y_,z_,w_,q_,r_]:>squareNo0V@@vectorNorm3d/@{x,y,z,w,q,r},square0mom123[x_,y_,z_]:> square0mom123V@@vectorNorm3d/@{x,y,z},squaremom1212[x_,y_,z_]:> squaremom1212V@@vectorNorm3d/@{x,y,z},squaremom1122[x_,y_,z_]:> squaremom1122V@@vectorNorm3d/@{x,y,z},squaremom0011[x_]:> squaremom0011V[vectorNorm3d@x],squaremom0101[x_]:> squaremom0101V[vectorNorm3d@x],triangleD1[x_,y_,z_]:> triangleD1V@@vectorNorm3d/@{x,y,z},triangleD11[x_,y_,z_]:> triangleD11V@@vectorNorm3d/@{x,y,z},triangleD2[x_,y_,z_]:> triangleD2V@@vectorNorm3d/@{x,y,z}};


substitutionFunctions={propag[x_]:> Prop@x,bubbleV[x_]:> BubbleV@x,sunsetV[x_]:> SunsetV@x,triangleV[x_,y_,z_]:> TriangleV@@{x,y,z},tadSunset-> tadSunsetValue,tadTrianBub-> tadTrianBubValue,
squareNo0V[x_,y_,z_,w_,q_,r_]:> SquareNo0V[x,y,z,w,q,r],square0mom123V[x_,y_,z_]:> Square0mom123V[x,y,z],squaremom1212V[x_,y_,z_]:> Squaremom1212V[x,y,z],squaremom1122V[x_,y_,z_]:> Squaremom1122V[x,y,z],squaremom0011V[x_]:>Squaremom0011V[x],squaremom0101V[x_]:>Squaremom0101V[x],triangleD1V[x_,y_,z_]:> TriangleD1V[x,y,z],triangleD11V[x_,y_,z_]:> TriangleD11V[x,y,z],triangleD2V[x_,y_,z_]:> TriangleD2V[x,y,z]};


(* ::Text:: *)
(**)


positivityRadialComponents={\!\(\*SubscriptBox[\(q1\), \("\<\[Rho]\>"\)]\)>0&&\!\(\*SubscriptBox[\(q2\), \("\<\[Rho]\>"\)]\)>0&&\!\(\*SubscriptBox[\(q3\), \("\<\[Rho]\>"\)]\)>0&&\!\(\*SubscriptBox[\(q4\), \("\<\[Rho]\>"\)]\)>0&&\!\(\*SubscriptBox[\(q5\), \("\<\[Rho]\>"\)]\)>0&&\!\(\*SubscriptBox[\(q6\), \("\<\[Rho]\>"\)]\)>0&&\!\(\*SubscriptBox[\(q7\), \("\<\[Rho]\>"\)]\)>0&&\!\(\*SubscriptBox[\(q8\), \("\<\[Rho]\>"\)]\)>0&&\!\(\*SubscriptBox[\(q9\), \("\<\[Rho]\>"\)]\)>0&&\!\(\*SubscriptBox[\(q10\), \("\<\[Rho]\>"\)]\)>0};


refinewithSymmetries[intgrand_]:=Refine[intgrand/. Subscript[0, n_]:>0/.\!\(\*SubscriptBox[\(q1\), \("\<\[Theta]\>"\)]\)->0/.\!\(\*SubscriptBox[\(q1\), \("\<\[Phi]\>"\)]\)->0/.\!\(\*SubscriptBox[\(q2\), \("\<\[Phi]\>"\)]\)->0,positivityRadialComponents]


simplifywithSymmetries[intgrand_]:=Simplify[intgrand/. Subscript[0, n_]:>0/.\!\(\*SubscriptBox[\(q1\), \("\<\[Theta]\>"\)]\)->0/.\!\(\*SubscriptBox[\(q1\), \("\<\[Phi]\>"\)]\)->0/.\!\(\*SubscriptBox[\(q2\), \("\<\[Phi]\>"\)]\)->0,positivityRadialComponents]


(*zeroExtMomenta[integrand_]:=Module[
	{pvars=pvariables[integrand]},
	ReplaceAll[integrand, Flatten@Table[{px->0, Subscript[px, a_]:>0},{px,pvars}]]]*)


(* ::Subsubsection::Closed:: *)
(*Square*)


squareV[k1_,k2_,k3_,k4_]:=If[k1=!=0&&k2=!=0&&k3=!=0&&k4=!=0&&k1=!=-k2&&k2=!=-k3&&k1=!=-k3,squareNo0[k1,k2,k3,k4,k1+k2,k2+k3],
Which[k1===0&&k2=!=0&&k3=!=0&&k4=!=0,square0mom123[k2,k3,k4],
k2===0&&k1=!=0&&k3=!=0&&k4=!=0,square0mom123[k3,k4,k1],
k3===0&&k1=!=0&&k2=!=0&&k4=!=0,square0mom123[k4,k1,k2],
k4===0&&k1=!=0&&k2=!=0&&k3=!=0,square0mom123[k1,k2,k3],
k1===-k3&&k1=!=0&&k2=!=0&&k3=!=0&&k4=!=0,squaremom1212[k1,k2,k1+k2],
k1===-k2&&k1=!=0&&k2=!=0&&k3=!=0&&k4=!=0,squaremom1122[k1,k3,k1+k3],
k3===-k2&&k1=!=0&&k2=!=0&&k3=!=0&&k4=!=0,squaremom1122[k2,k4,k2+k4],
(k1===k2===0&&k3=!=0)||(k1===k4===0&&k3=!=0),squaremom0011 [k3],
(k2===k3===0&&k1=!=0)||(k3===k4===0&&k1=!=0),squaremom0011 [k1],
(k2===k4===0&&k3=!=0),squaremom0101 [k3],
(k1===k3===0&&k2=!=0),squaremom0101 [k2],
k1===0&&k2===0&&k3===0&&k4===0,1/(64\[Pi])]
]


Square0mom123V[k1_,k2_,k3_]:=(k2^2 (k1^2-k2^2+k3^2)TriangleV[k1,k2,k3]+(-2 k1^4+2 k2^2 k3^2-2 k3^4+k1^2 (4 k3^2+k2^2 (2+k3^2)))/(4 \[Pi](4+k1^2) (4+k3^2)))/(2(k1^2 k2^2 k3^2+AreaTriangleSquaredtimes4[k1,k2,k3]))


Squaremom1212V[k1s_,k2s_,kxs_]:=(-(ArcTan[Sqrt[k1s^4 (-1+2 k2s^2)-(k2s^2-kxs^2)^2+k1s^2 (2 k2s^4+2 kxs^2-k2s^2 (-2+kxs^2))]/(8+3 k1s^2+3 k2s^2-kxs^2)]/Sqrt[k1s^4 (-1+2 k2s^2)-(k2s^2-kxs^2)^2+k1s^2 (2 k2s^4+2 kxs^2-k2s^2 (-2+kxs^2))])+4\[Pi] TriangleV[k1s,k2s,kxs])/(2(k1s^2+k2s^2-kxs^2) \[Pi])


Squaremom1122V[k1s_,k3s_,kys_]:=1/2 ((2 k1s^4+2 k3s^2 (k3s^2-kys^2)-k1s^2 (2 kys^2+k3s^2 (4+kys^2)))/(4\[Pi] (4+k1s^2) (4+k3s^2) )-kys^2 (k1s^2+k3s^2-kys^2) TriangleV[k1s,k3s,kys])/(k1s^4+(k3s^2-kys^2)^2-k1s^2 (2 kys^2+k3s^2 (2+kys^2)))


Squaremom0011V[k1_]:= (k1^2+8)/(32\[Pi] (k1^2+4)^2)


Squaremom0101V[k1_]:=1/(4\[Pi] (k1^2+4)^2)


(* ::Text:: *)
(*Square in case we are in any of the soft limits*)


D4Matrix[k1_,k2_,k3_,k4_,kx_,ky_]:={{1,1+k1/2,1+kx/2,1+k4/2},{1+k1/2,1,1+k2/2,1+ky/2},{1+kx/2,1+k2/2,1,1+k3/2},{1+k4/2,1+ky/2,1+k3/2,1}}


D4[k1_,k2_,k3_,k4_,kx_,ky_]:=Det@({{1,1+k1/2,1+kx/2,1+k4/2},{1+k1/2,1,1+k2/2,1+ky/2},{1+kx/2,1+k2/2,1,1+k3/2},{1+k4/2,1+ky/2,1+k3/2,1}})


F4Matrix[l_,k1_,k2_,k3_,k4_,kx_,ky_]:=Table[Table[D4Matrix[k1,k2,k3,k4,kx,ky][[i,j]](1-KroneckerDelta[j,l])+KroneckerDelta[j,l],{j,4}],{i,4}]


F4[l_,k1_,k2_,k3_,k4_,kx_,ky_]:=Det@Table[Table[D4Matrix[k1,k2,k3,k4,kx,ky][[i,j]](1-KroneckerDelta[j,l])+KroneckerDelta[j,l],{j,4}],{i,4}]


square1234xy=(fmomenta4[1,k1,k2,k3,k4,kx,ky]Simplify@trian123[Sqrt[k2],Sqrt[k3],Sqrt[ky]]+fmomenta4[2,k1,k2,k3,k4,kx,ky]Simplify@trian123[Sqrt[kx],Sqrt[k3],Sqrt[k4]]+fmomenta4[3,k1,k2,k3,k4,kx,ky]Simplify@trian123[Sqrt[k1],Sqrt[ky],Sqrt[k4]]+fmomenta4[4,k1,k2,k3,k4,kx,ky]Simplify@trian123[Sqrt[k1],Sqrt[k2],Sqrt[kx]])/delta4[k1,k2,k3,k4,kx,ky];


SquareNo0V[k1_,k2_,k3_,k4_,kx_,ky_]:=(F4[1,k1^2,k2^2,k3^2,k4^2,kx^2,ky^2]TriangleV[k2,k3,ky]+F4[2,k1^2,k2^2,k3^2,k4^2,kx^2,ky^2]TriangleV[kx,k3,k4]+F4[3,k1^2,k2^2,k3^2,k4^2,kx^2,ky^2]TriangleV[k1,ky,k4]+F4[4,k1^2,k2^2,k3^2,k4^2,kx^2,ky^2]TriangleV[k1,k2,kx])/D4[k1^2,k2^2,k3^2,k4^2,kx^2,ky^2]


(* ::Subsubsection::Closed:: *)
(*TrianglesD*)


TriangleD1V[k1_,k2_,k3_]:=1/2 If[k1=!=0&&k2=!=0&&k3=!=0,(( k1^2-k2^2-k3^2)/(4\[Pi](4+k1^2))-(2 k1^2-2 k3^2-k2^2 (2+k3^2))TriangleV[k1,k2,k3])/(k1^4+(k2^2-k3^2)^2-k1^2 (2 k3^2+k2^2 (2+k3^2))),If[k1===0,(-6-k2^2)/(96 \[Pi] (4+k2^2)^2),-(1/(16 \[Pi] (4+k1^2)^2))]]


TriangleD2V[k1_,k2_,k3_]:=1/4 If[k1=!=0&&k2=!=0&&k3=!=0,1/(4 \[Pi] ((k1^2 k2^2 k3^2+AreaTriangleSquaredtimes4[k1,k2,k3])^2) ) ((4 k1^2-2(2 k3^2+k2^2 (2+k3^2)))((-k1^2+k2^2+k3^2)/(4+k1^2)+(2 k1^2-2 k3^2-k2^2 (2+k3^2))4\[Pi] TriangleV[k1,k2,k3])-((k1^2-k2^2-k3^2) (2 k1^2-2 k3^2-k2^2 (2+k3^2)))/(4+k1^2)+(-2 k1^2+2 k3^2+k2^2 (2+k3^2))^2 4\[Pi] TriangleV[k1,k2,k3])-1/(4 \[Pi](k1^2 k2^2 k3^2+AreaTriangleSquaredtimes4[k1,k2,k3]) ) (2/(4+k1^2)+(2 (-k1^2+k2^2+k3^2))/(4+k1^2)^2-16\[Pi] TriangleV[k1,k2,k3]),If[k1===0,(80+30 k2^2+3(k2^2)^2)/(960 \[Pi] (4+k2^2)^3),-(1/(12 \[Pi] (4+k1^2)^3))]]


TriangleD11V[k1_,k2_,k3_]:=1/4 If[k1=!=0&&k2=!=0&&k3=!=0,1/(4\[Pi] ((k1^2 k2^2 k3^2+AreaTriangleSquaredtimes4[k1,k2,k3])^2) ) (((-k1^2+k2^2+k3^2)  (2 (k2^2-k3^2)-k1^2 (2+k3^2)))/(4+k1^2)+((k1^2-k2^2+k3^2)(2 k1^2-2 k3^2-k2^2 (2+k3^2)))/(4+k2^2)+1/((4+k1^2) (4+k2^2) (8+k1^2+k2^2+k3^2) ) (2 k1^8-k1^6 k2^2 (8+k3^2)+6 k1^4 (8 k3^2+4 k2^2 k3^2+k2^4 (2+k3^2))+k1^2 (24 k2^4 k3^2-32 k3^2 (-4+k3^2)-k2^6 (8+k3^2)+k2^2 k3^2 (96-8 k3^2+k3^4))+2 (k2^8+24 k2^4 k3^2-16 k2^2 k3^2 (-4+k3^2)-k3^4 (64+8 k3^2+k3^4)))-3 (-2 k2^2+2 k3^2+k1^2 (2+k3^2)) (2 k1^2-2 k3^2-k2^2 (2+k3^2))4\[Pi] TriangleV[k1,k2,k3])-(2 (2+k3^2)TriangleV[k1,k2,k3]-(2((k1^2-k2^2)^2-k3^4))/(4\[Pi](4+k1^2) (4+k2^2)(8+k1^2+k2^2+k3^2)))/(k1^2 k2^2 k3^2+AreaTriangleSquaredtimes4[k1,k2,k3]),Which[k1===0&&k2=!=0,(8+k2^2)/(192 \[Pi] (4+k2^2)^3),k2===0&&k1=!=0,(8+k1^2)/(192 \[Pi] (4+k1^2)^3),k3===0&&k1=!=0,1/(24 \[Pi] (4+k2^2)^3),k1===0&&k2===0&&k3===0,1/(1536 \[Pi])]]


(* ::Subsection:: *)
(*Derivation for presentation*)


(* ::Subsubsection::Closed:: *)
(*Substitutions*)


(* ::Text:: *)
(*Substitutions from derivatives of the blocks to their values and the scalar products*)


derPropagatorsBubbles={
	Derivative[2][prop][a_]:>1/3 (prop[a]^2-4prop[a]^3),
	Derivative[2][bubble][a_]:>-(1/(6\[Pi]))(1/(4+vec3dsq[a])^2),
	Derivative[1][prop][a_]*Derivative[1][prop][b_]:>2vec3dprod[a,b] (prop[a]^2 prop[b]^2)/3,
	Derivative[1][prop][a_]^n_/;(n>= 2&&IntegerQ[n/2]):>2vec3dsq[a]*prop[a]^(2n)/3,
	Derivative[1][prop][a_]*Derivative[1][bubble][b_]:>-(vec3dprod[a,b]/Sqrt[vec3dsq[b]]) (prop[a]^2 bubbleDer1[b])/3,
	Derivative[1][bubble][a_]*Derivative[1][bubble][b_]:>vec3dprod[a,b]/(Sqrt[vec3dsq[a]]*Sqrt[vec3dsq[b]]) (bubbleDer1[a]bubbleDer1[b])/6,
	Derivative[1][bubble][a_]^n_/;(n>= 2&&IntegerQ[n/2]):>bubbleDer1[a]^n/6,
	Derivative[1][prop][a_]*Derivative[1][sunset][b_]:>-(vec3dprod[a,b]/Sqrt[vec3dsq[b]]) prop[a]^2/3 * sunsetDer1[b],
	Derivative[2][sunset][a_]:>-(1/(96 \[Pi]^2))(1/(9+vec3dsq[a])),
	Derivative[1][bubble][a_]*Derivative[1][sunset][b_]:>vec3dprod[a,b]/(Sqrt[vec3dsq[a]]*Sqrt[vec3dsq[b]])*bubbleDer1[a]/6 sunsetDer1[b]};


derTriangles1={
	Derivative[1,0,0][triangle][a_,b_,c_]*Derivative[1][prop][d_]:> -2vec3dprod[a,d] (triangleD1[a,b,c]prop[d]^2)/3,
	Derivative[0,1,0][triangle][a_,b_,c_]*Derivative[1][prop][d_]:>-2 vec3dprod[b,d] (triangleD1[b,c,a]prop[d]^2)/3,
	Derivative[0,0,1][triangle][a_,b_,c_]*Derivative[1][prop][d_]:> -2vec3dprod[c,d] (triangleD1[c,a,b]prop[d]^2)/3,
	Derivative[1,0,0][triangle][a_,b_,c_]*Derivative[1][bubble][d_]:> vec3dprod[a,d]/Sqrt[vec3dsq[d]] (triangleD1[a,b,c]*bubbleDer1[d])/3,
	Derivative[0,1,0][triangle][a_,b_,c_]*Derivative[1][bubble][d_]:> vec3dprod[b,d]/Sqrt[vec3dsq[d]] (triangleD1[b,c,a]bubbleDer1[d])/3,
	Derivative[0,0,1][triangle][a_,b_,c_]*Derivative[1][bubble][d_]:>vec3dprod[c,d]/Sqrt[vec3dsq[d]] (triangleD1[c,a,b]bubbleDer1[d])/3,
	Derivative[1,0,0][triangle][a_,b_,c_]*Derivative[1][sunset][d_]:> vec3dprod[a,d]*triangleD1[a,b,c]/(3Sqrt[vec3dsq[d]]) * sunsetDer1[d],
	Derivative[0,1,0][triangle][a_,b_,c_]*Derivative[1][sunset][d_]:> vec3dprod[b,d]*triangleD1[b,c,a]/(3Sqrt[vec3dsq[d]]) * sunsetDer1[d],
	Derivative[0,0,1][triangle][a_,b_,c_]*Derivative[1][sunset][d_]:> vec3dprod[c,d]*triangleD1[c,a,b]/(3Sqrt[vec3dsq[d]]) * sunsetDer1[d]};


derTriangles2={
	Derivative[2,0,0][triangle][a_,b_,c_]:> triangleD1[a,b,c]+(2vec3dsq[a]*triangleD2[a,b,c])/3,
	Derivative[0,2,0][triangle][a_,b_,c_]:>triangleD1[b,c,a]+(2vec3dsq[b]*triangleD2[b,c,a])/3,
	Derivative[0,0,2][triangle][a_,b_,c_]:>triangleD1[c,a,b]+ (2vec3dsq[c]*triangleD2[c,a,b])/3,
	Derivative[1,1,0][triangle][a_,b_,c_]:>( 2 vec3dprod[a,b]*triangleD11[a,b,c])/3,
	Derivative[1,0,1][triangle][a_,b_,c_]:>( 2 vec3dprod[a,c]*triangleD11[c,a,b])/3,
	Derivative[0,1,1][triangle][a_,b_,c_]:>( 2 vec3dprod[b,c]*triangleD11[b,c,a])/3,
	Derivative[1,0,0][triangle][a_,b_,c_]*Derivative[1,0,0][triangle][d_,e_,f_]:> 2vec3dprod[a,d] (triangleD1[a,b,c]*triangleD1[d,e,f])/3,
	Derivative[1,0,0][triangle][a_,b_,c_]*Derivative[0,1,0][triangle][d_,e_,f_]:> 2vec3dprod[a,e] (triangleD1[a,b,c]*triangleD1[e,f,d])/3,
	Derivative[1,0,0][triangle][a_,b_,c_]*Derivative[0,0,1][triangle][d_,e_,f_]:> 2vec3dprod[a,f] (triangleD1[a,b,c]*triangleD1[f,d,e])/3,
	Derivative[0,1,0][triangle][a_,b_,c_]*Derivative[0,1,0][triangle][d_,e_,f_]:> 2vec3dprod[b,e] (triangleD1[b,c,a]*triangleD1[e,f,d])/3,
	Derivative[0,1,0][triangle][a_,b_,c_]*Derivative[0,0,1][triangle][d_,e_,f_]:> 2vec3dprod[b,f] (triangleD1[b,c,a]*triangleD1[f,d,e])/3,
	Derivative[0,0,1][triangle][a_,b_,c_]*Derivative[0,0,1][triangle][d_,e_,f_]:> 2vec3dprod[c,f] (triangleD1[c,a,b]*triangleD1[f,d,e])/3};


(* ::Subsubsection:: *)
(*Functions*)


(* ::Text:: *)
(*Deriveintegrand: from function to derive respect to p to the double derivative*)


Deriveintegrand[exp_]:=Module[{integrand,symbolic,subsymbolic,subsymbolic2,result},
integrand=exp/.prop[a_]/;MemberQ[a,-p]:>prop[-a]/. bubble[a_]/;MemberQ[a,-p]:> bubble[-a];
symbolic=D[integrand,{p,2}];(*Print[symbolic];*)
subsymbolic=Expand[symbolic]/.derPropagatorsBubbles; subsymbolic2=subsymbolic/.derPropagatorsBubbles/.derTriangles1/.derTriangles2;
result=subsymbolic2/.p->0/. Subscript[0, n_]:>0
]


DeriveintegrandR[exp_]:=Refine@Deriveintegrand@exp


DeriveintegrandS[exp_]:=Simplify@Deriveintegrand@exp


DeriveintegrandSC[exp_]:=Simplify[Deriveintegrand@exp,TimeConstraint->5]


bubbleDer1[a_]:=(1/(\[Pi] Sqrt[vec3dsq[a]] (8+2vec3dsq[a]))-bubble[a]/Sqrt[vec3dsq[a]])


sunsetDer1[x_]:=-(1/(32 \[Pi]^2))(2/(Sqrt[vec3dsq[x]] (1+vec3dsq[x]/9))+(2Sqrt[vec3dsq[x]])/(9 (1+vec3dsq[x]/9))-(6 ArcTan[Sqrt[vec3dsq[x]]/3])/vec3dsq[x])


(* ::Subsection:: *)
(*Functions*)


(*inverseRenameMomenta = Table[Momentum[k]->ToExpression[internalMomentumList[[k]]],{k,1,Length@internalMomentumList}]*)


inverseRenameMomenta = { Momentum[i_] :> qauxify[i] }


(*inverseRenameExternalMomenta = Table[ExternalMomentum[k]->ToExpression[externalMomentumList[[k]]],{k,1,Length@externalMomentumList}]*)


inverseRenameExternalMomenta = { ExternalMomentum[i_] :> pextify[i] }


inverseRenameSubdiagrams = {
	Propagator -> prop,
	BubbleSubdiagram -> bubble,
	SunsetSubdiagram -> sunset,
	TriangleSubdiagram -> triangle,
	SquareSubdiagram -> square,
	TadSunsetSubdiagram[] :> tadSunset,
	TadTriangleBubblesSubdiagram[] :> tadTrianBub
}


renameExtMomentum = { p1 -> p }


(*renameMomentaExplicit = Join[
  Table[Subscript[ToExpression[internalMomentumList[[k]]], "\[Rho]"] -> Momentum[k, "\[Rho]"], {k, 1, Length @ internalMomentumList}],
  Table[Subscript[ToExpression[internalMomentumList[[k]]], "\[Theta]"] -> Momentum[k, "\[Theta]"], {k, 1, Length @ internalMomentumList}],
  Table[Subscript[ToExpression[internalMomentumList[[k]]], "\[Phi]"] -> Momentum[k, "\[Phi]"], {k, 1, Length @ internalMomentumList}]]*)


(*subscriptedMomentaR = { Subscript[Momentum[i_],sub_] :> Momentum[i,sub] }*)


(*renameExternalMomentaExplicit = Join[
	Table[Subscript[ToExpression[externalMomentumList[[k]]],"\[Rho]"] -> ExternalMomentum[k,"\[Rho]"],{k,1,Length@externalMomentumList}],
	Table[Subscript[ToExpression[externalMomentumList[[k]]],"\[Theta]"] -> ExternalMomentum[k,"\[Theta]"],{k,1,Length@externalMomentumList}],
	Table[Subscript[ToExpression[externalMomentumList[[k]]],"\[Phi]"] -> ExternalMomentum[k,"\[Phi]"],{k,1,Length@externalMomentumList}]
]*)


(*subscriptedExtMomentaR = { Subscript[ExternalMomentum[i_],sub_] :> ExternalMomentum[i,sub] }*)


WriteExplicit[integrandfunc_, OptionsPattern[]] :=
  Module[{simp = OptionValue["Simplification"],integrand = ReplaceAll[integrandfunc, Join[inverseRenameMomenta, inverseRenameExternalMomenta,inverseRenameSubdiagrams]]},
    substitutePVariables[substituteQVariables[
    Which[
        simp === "Refine",
          refinewithSymmetries[refinewithSymmetries[integrand /. substitutionSquare
             /. substitutionVectors] /. substitutionFunctions]
        ,
        simp === "Simplify",
          simplifywithSymmetries[simplifywithSymmetries[integrand /. 
            substitutionSquare /. substitutionVectors] /. substitutionFunctions]
        ,
        And@@{simp=!="Refine",simp=!="Simplify"},
			Message[WriteExplicit::simp];
          simplifywithSymmetries[simplifywithSymmetries[integrand /. 
            substitutionSquare /. substitutionVectors] /. substitutionFunctions]
      ]
      ]]
  ]


WriteExplicit::simp = "The options for \"Simplification\" are: \"Refine\" and \"Simplify\". Defaulting to \"Simplification\" -> \"Refine\"."


(* ::Text:: *)
(**)


DeriveAndWriteExplicit::dsquare = "Derivative of the square subdiagram has not been implemented"


DeriveAndWriteExplicit[integrandfunc_] :=
    Module[{der, integrand = ReplaceAll[integrandfunc, Join[inverseRenameMomenta,inverseRenameExternalMomenta, inverseRenameSubdiagrams]]/.renameExtMomentum},
        der = DeriveintegrandSC[integrand];
        substitutePVariables[substituteQVariables[
        If[MemberQ[der, Derivative[2, 0, 0, 0][square][_, _, _, _] |
          Derivative[0, 2, 0, 0][square][_, _, _, _] | Derivative[0, 0, 2, 0][square][_, _, _, _] | Derivative[0, 0, 0, 2][square][_, _, _, _], All],
                    Message[DeriveAndWriteExplicit::dsquare],
                    refinewithSymmetries[refinewithSymmetries[der /. vectornameComponents /. substitutionSquare /. substitutionVectors] /. substitutionFunctions]
        ]
        ]]
    ]


(*CountLoops[intgrand_]:=Module[{int=intgrand,i=1,l},
While[MemberQ[int,Momentum[i],All]&&i<10,i++];
l=i-1
]*)


MomVars[integrand_]:=
Select[
	Variables[Level[integrand,{-2}]],
	Head[#] == Momentum &
]


CountLoops[integrand_]:=
Length@DeleteDuplicates[
	ReplaceAll[
		MomVars[integrand],
		Momentum[i_,___]:>i
	]
]


(* ::Subtitle:: *)
(*End package*)


(* ::Subsection:: *)
(*Set Attributes to functions*)


ToExpression[
   Names["GSberveglieri`Phi4tools`*"], 
   InputForm, 
   Function[sym,SetAttributes[sym, {ReadProtected(*,Locked*)}],HoldFirst] 
]


End[]; (* End `Private` *)


EndPackage[];
