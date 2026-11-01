args=Rest[$ScriptCommandLine];
If[Length[args]!=2,Exit[2]];
input=args[[1]]; output=args[[2]];
payload=Import[input,"RawJSON"];
If[!AssociationQ[payload] || !KeyExistsQ[payload,"report_betas"],Exit[3]];
raw=payload["report_betas"];
If[!AssociationQ[raw],Exit[4]];
tex=Association@KeyValueMap[Function[{name,expression},
    name->ToString[TeXForm[ToExpression[expression,InputForm]]]],raw];
Export[output,<|"report_beta_latex"->tex|>,"RawJSON"];
