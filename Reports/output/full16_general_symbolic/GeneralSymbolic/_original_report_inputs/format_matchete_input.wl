args=Rest[$ScriptCommandLine];
If[Length[args]!=2,Exit[2]];
input=args[[1]]; output=args[[2]];
text=Import[input,"Text"];
expr=ToExpression[text,InputForm];
parts=If[Head[expr]===Plus,List@@expr,{expr}];
tex=Map[ToString[TeXForm[#]]&,parts];
Export[output,<|"terms"->tex|>,"RawJSON"];
