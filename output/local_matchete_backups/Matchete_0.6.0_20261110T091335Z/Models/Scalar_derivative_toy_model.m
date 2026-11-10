(* ::Package:: *)

(* ::Title:: *)
(*Scalar toy model with derivative interactions*)


(* ::Subtitle:: *)
(*Toy model with EFT derivative interaction*)


(* ::Section:: *)
(*Lagrangian*)


DefineField[\[Phi], Scalar, SelfConjugate-> True, Mass-> {Heavy, M}]
DefineField[\[Psi], Scalar, SelfConjugate-> True, Mass-> {Light, m}]


DefineCoupling[g, SelfConjugate->True]


FreeLag[\[Phi], \[Psi]]+ g[]* CD[a, \[Phi][]]^2 *\[Psi][]
