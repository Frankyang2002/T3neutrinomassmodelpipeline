args=Rest[$ScriptCommandLine];
If[Length[args]!=2,Exit[2]];
payload=Import[args[[1]],"RawJSON"];
If[!AssociationQ[payload] || !KeyExistsQ[payload,"report_betas"],Exit[3]];
raw=payload["report_betas"];
If[!AssociationQ[raw],Exit[4]];
make[expression_]:=Module[{expr,terms,rows,fac,c,structure},
  expr=ToExpression[expression,InputForm];
  terms=If[Head[Expand[expr]]===Plus,List@@Expand[expr],{Expand[expr]}];
  rows=Map[Function[term,
     fac=FactorTermsList[term];
     c=First[fac];structure=Times@@Rest[fac];
     If[!NumericQ[c] || !FreeQ[c,_Symbol],c=1;structure=term];
     <|"structure"->ToString[InputForm[structure]],
       "structure_latex"->ToString[TeXForm[structure]],
       "coefficient"->ToString[InputForm[c]],
       "coefficient_latex"->ToString[TeXForm[c]]|>
  ],terms];
  rows
];
Export[args[[2]],<|"terms"->Association@KeyValueMap[
    Function[{name,expression},name->make[expression]],raw]|>,"RawJSON"];
