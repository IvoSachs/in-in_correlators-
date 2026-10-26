(* ::Package:: *)

If[!MemberQ[$Packages, "CosmologicalBaikov`"],
  Get[FileNameJoin[{DirectoryName[$InputFileName], "CosmologicalBaikov.wl"}]]
];

BeginPackage["CosmologicalBaikov011`", {"CosmologicalBaikov`"}];

J011::usage = "J011[s,w,v] denotes Ibar^{0 s w, v 0 0; 0} in the D2-D3 (I011) sector.";
MasterBasis011::usage = "MasterBasis011[] returns the five-master basis for the I011 sector.";
ScalingIdentity011::usage = "ScalingIdentity011[s,w,v] returns the Lambda-power lowering identity.";
ParameterDerivativeA011::usage = "ParameterDerivativeA011[s,w,v] returns the raw a derivative before IBP reduction.";
IBPIdentity011::usage = "IBPIdentity011[type,s,w,v] returns a D2-D3-D4 vector IBP; type is pp, pq, qp, or qq.";
TangentialIBPIdentity011::usage = "TangentialIBPIdentity011[s,w,v] returns the Lambda-preserving tangential IBP.";
TangentialCheck011::usage = "TangentialCheck011[] verifies that the tangential vector field annihilates Lambda and has zero divergence.";
ReductionData011::usage = "ReductionData011[] constructs the staged Hermite/IBP reduction data for I011.";
ReductionRules011::usage = "ReductionRules011[] returns the I011 reduction rules.";
ReduceIntegral011::usage = "ReduceIntegral011[expr] applies the staged I011 reduction rules.";
DifferentialMatrixA011::usage = "DifferentialMatrixA011[] returns A011 in partial_a m011=A011.m011.";
GenericityConditions011::usage = "GenericityConditions011[] gives the elementary nonvanishing conditions assumed by this reduction.";
NumericalCrossCheck011::usage = "NumericalCrossCheck011[] evaluates A011 at the exact regression point (a,b,c,X2,X3,eps)=(5,-1,1,1,2,1/7).";
CutGramPolynomial011::usage = "CutGramPolynomial011[z] returns the D2=D3=0 Gram/Baikov polynomial for the I011 sector.";
LeadingSingularityLocations011::usage = "LeadingSingularityLocations011[] returns the two D4 values on the I011 maximal cut.";
LeadingSingularityBasisMatrix011::usage = "LeadingSingularityBasisMatrix011[] maps {I011,000;0,I011,-100;0} to the two branch-adapted combinations.";
CutLambda011::usage = "CutLambda011[z,aVariable] returns Lambda on D2=D3=0 at D4=z.";

Begin["`Private`"];

ClearAll[J011, z2, z3, z4];

xSP = z2 - x2sq;
ySP = z3 - x3sq;
zSP = z4 - cK/2;

ClearAll[polyToIntegrals011, lambdaRaise011, polyTimesLambdaRaise011];

polyToIntegrals011[poly_, {s_, w_, v_}] := Module[{terms},
  terms = CoefficientRules[Expand[poly], {z2, z3, z4}];
  Total[(Last[#] J011[
      s - First[#][[1]], w - First[#][[2]], v - First[#][[3]]
    ]) & /@ terms]
];

lambdaRaise011[s_, w_, v_] := (
  -2 (1/2 - 2 eps - s - w - v) J011[s, w, v]
  - 2 s x2sq J011[s + 1, w, v]
  - 2 w x3sq J011[s, w + 1, v]
  - v cK J011[s, w, v + 1]
)/Lambda0[];

polyTimesLambdaRaise011[poly_, {s_, w_, v_}] := Module[{terms},
  terms = CoefficientRules[Expand[poly], {z2, z3, z4}];
  Total[(Last[#] lambdaRaise011[
      s - First[#][[1]], w - First[#][[2]], v - First[#][[3]]
    ]) & /@ terms]
];

$contractions011 = <|
  "pp" -> <|
    "div" -> 1 - 2 eps,
    "D2" -> 2 xSP, "D3" -> 0, "D4" -> zSP,
    "L" -> cK zSP - 2 bK xSP
  |>,
  "pq" -> <|
    "div" -> 0,
    "D2" -> 2 zSP, "D3" -> 0, "D4" -> ySP,
    "L" -> cK ySP - 2 bK zSP
  |>,
  "qp" -> <|
    "div" -> 0,
    "D2" -> 0, "D3" -> 2 zSP, "D4" -> xSP,
    "L" -> cK xSP - 2 aK zSP
  |>,
  "qq" -> <|
    "div" -> 1 - 2 eps,
    "D2" -> 0, "D3" -> 2 ySP, "D4" -> zSP,
    "L" -> cK zSP - 2 aK ySP
  |>
|>;

IBPIdentity011::type = "Unknown IBP type `1`. Use pp, pq, qp, or qq.";
IBPIdentity011[type_String, s_Integer, w_Integer, v_Integer] := Module[
  {q, lhs},
  If[!KeyExistsQ[$contractions011, type],
    Message[IBPIdentity011::type, type]; Return[$Failed]
  ];
  q = $contractions011[type];
  lhs = Together[
    q["div"] J011[s, w, v]
    - s polyToIntegrals011[q["D2"], {s + 1, w, v}]
    - w polyToIntegrals011[q["D3"], {s, w + 1, v}]
    - v polyToIntegrals011[q["D4"], {s, w, v + 1}]
    - 1/2 polyTimesLambdaRaise011[q["L"], {s, w, v}]
  ];
  lhs == 0
];

ScalingIdentity011[s_Integer, w_Integer, v_Integer] := (
  Lambda0[] J011[s, w, v, 1] ==
    -2 (1/2 - 2 eps - s - w - v) J011[s, w, v]
    - 2 s x2sq J011[s + 1, w, v]
    - 2 w x3sq J011[s, w + 1, v]
    - v cK J011[s, w, v + 1]
);

(* Since D6=p3^2+b=D3+b-x3sq, differentiating in a gives
   1/2 (D3+b-x3sq) Lambda^(-3/2). *)
ParameterDerivativeA011[s_Integer, w_Integer, v_Integer] := Together[
  1/2 (
    lambdaRaise011[s, w - 1, v]
    + (bK - x3sq) lambdaRaise011[s, w, v]
  )
];

TangentialIBPIdentity011[s_Integer, w_Integer, v_Integer] := Module[{lhs},
  lhs = (
    cK First@IBPIdentity011["pp", s, w, v]
    - 2 aK First@IBPIdentity011["pq", s, w, v]
    + 2 bK First@IBPIdentity011["qp", s, w, v]
    - cK First@IBPIdentity011["qq", s, w, v]
  );
  Together[lhs] == 0
];

TangentialCheck011[] := Module[{xp, yp, zp, lambda, tp, tq},
  lambda = cK zp - bK xp - aK yp + Lambda0[];
  tp = cK (2 xp D[lambda, xp] + zp D[lambda, zp])
    - 2 aK (2 zp D[lambda, xp] + yp D[lambda, zp]);
  tq = 2 bK (xp D[lambda, zp] + 2 zp D[lambda, yp])
    - cK (zp D[lambda, zp] + 2 yp D[lambda, yp]);
  {Expand[tp + tq], Expand[(1 - 2 eps) cK - (1 - 2 eps) cK]}
];

MasterBasis011[] := {
  J011[0, 0, 0],
  J011[1, 0, 0],
  J011[0, 1, 0],
  J011[1, 1, 0],
  J011[1, 1, -1]
};

ClearAll[buildReductionData011];
buildReductionData011[] := Module[
  {
    m, s0, s1, t0, t1, subRules,
    mixed0, eq0, ca0, sol0,
    mixed1, eq1, ca1, sol1, residual1,
    d4, d5, solMat0, solMat1,
    row1, row2, row3, row4, row5,
    rawMatrix, allRules
  },

  m = MasterBasis011[];

  s0 = (4 eps - 1)/(2 x2sq (aK - x2sq));
  s1 = (aK (1 + 2 eps) - x2sq (1 + 4 eps))/
    (2 x2sq (aK - x2sq));
  t0 = (4 eps - 1)/(2 x3sq (bK - x3sq));
  t1 = (bK (1 + 2 eps) - x3sq (1 + 4 eps))/
    (2 x3sq (bK - x3sq));

  subRules = {
    J011[2, 0, 0] -> s0 m[[1]] + s1 m[[2]],
    J011[0, 2, 0] -> t0 m[[1]] + t1 m[[3]],
    J011[1, 0, -1] ->
      cK/(2 aK) (m[[1]] + (aK - x2sq) m[[2]]),
    J011[2, 0, -1] ->
      cK/(2 aK) (m[[2]] + (aK - x2sq) (s0 m[[1]] + s1 m[[2]])),
    J011[0, 1, -1] ->
      cK/(2 bK) (m[[1]] + (bK - x3sq) m[[3]]),
    J011[0, 2, -1] ->
      cK/(2 bK) (m[[3]] + (bK - x3sq) (t0 m[[1]] + t1 m[[3]]))
  };

  mixed0 = {
    J011[2, 1, -1], J011[2, 1, 0],
    J011[1, 2, 0], J011[1, 2, -1]
  };
  eq0 = Numerator[Together[
      First[IBPIdentity011[#, 1, 1, 0]] /. subRules
    ]] & /@ Keys[$contractions011];
  ca0 = CoefficientArrays[eq0, mixed0];
  sol0 = Thread[mixed0 -> LinearSolve[ca0[[2]], -ca0[[1]]]];

  mixed1 = {J011[1, 1, -2], J011[2, 1, -2], J011[1, 2, -2]};
  eq1 = (
    First[IBPIdentity011[#, 1, 1, -1]] /. subRules /. sol0
  ) & /@ Keys[$contractions011];
  ca1 = CoefficientArrays[eq1, mixed1];
  sol1 = Thread[
    mixed1 -> LinearSolve[ca1[[2, 1 ;; 3]], -ca1[[1, 1 ;; 3]]]
  ];
  residual1 = eq1[[4]] /. sol1;

  (* The three lower rows can be reduced directly and are compact. *)
  row1 = bK eps/Lambda0[] m[[1]];
  row2 = (
    (4 eps - 1)/(2 aK (aK - x2sq)) m[[1]]
    + (-1/(2 aK) + eps (
        bK/Lambda0[]
        + (aK - 2 x2sq)/(aK (aK - x2sq))
      )) m[[2]]
  );
  row3 = bK eps/Lambda0[] m[[3]];

  d4 = ParameterDerivativeA011[1, 1, 0] /. subRules;
  d5 = ParameterDerivativeA011[1, 1, -1] /. subRules /. sol0;
  solMat0 = Table[
    Coefficient[mixed0[[i]] /. sol0, m[[j]]],
    {i, 4}, {j, 5}
  ];
  solMat1 = Table[
    Coefficient[mixed1[[i]] /. sol1, m[[j]]],
    {i, 3}, {j, 5}
  ];
  row4 = (
    Table[Coefficient[d4, m[[j]]], {j, 5}]
    + Table[Coefficient[d4, mixed0[[i]]], {i, 4}] . solMat0
  );
  row5 = (
    Table[Coefficient[d5, m[[j]]], {j, 5}]
    + Table[Coefficient[d5, mixed1[[i]]], {i, 3}] . solMat1
  );

  rawMatrix = Join[
    {Table[Coefficient[row1, m[[j]]], {j, 5}]},
    {Table[Coefficient[row2, m[[j]]], {j, 5}]},
    {Table[Coefficient[row3, m[[j]]], {j, 5}]},
    {row4},
    {row5}
  ];
  allRules = Join[subRules, sol0, sol1];

  <|
    "Masters" -> m,
    "SubtopologyRules" -> subRules,
    "MixedV0Rules" -> sol0,
    "MixedVMinus1Rules" -> sol1,
    "MixedMatrices" -> {ca0[[2]], ca1[[2, 1 ;; 3]]},
    "FourthIBPResidual" -> residual1,
    "Rules" -> allRules,
    "RawMatrix" -> rawMatrix,
    "HermiteCoefficients" -> <|
      "s0" -> s0, "s1" -> s1,
      "t0" -> t0, "t1" -> t1
    |>
  |>
];

ReductionData011[] := ReductionData011[] = buildReductionData011[];
ReductionRules011[] := ReductionData011[]["Rules"];

ReduceIntegral011[expr_] := FixedPoint[
  Expand[# /. ReductionRules011[]] &,
  expr
];

Options[DifferentialMatrixA011] = {"SimplifyFunction" -> Identity};
DifferentialMatrixA011[OptionsPattern[]] := Module[{f, matrix},
  f = OptionValue["SimplifyFunction"];
  matrix = ReductionData011[]["RawMatrix"];
  Map[f, matrix, {2}]
];

GenericityConditions011[] := And[
  Lambda0[] != 0,
  aK != 0, bK != 0, cK != 0,
  x2sq != 0, x3sq != 0,
  aK - x2sq != 0,
  bK - x3sq != 0,
  4 eps - 1 != 0
];

NumericalCrossCheck011[] := Together[
  DifferentialMatrixA011[] /. {
    aK -> 5, bK -> -1, cK -> 1,
    x2sq -> 1, x3sq -> 2, eps -> 1/7
  }
];

CutGramPolynomial011[z_] :=
  x2sq x3sq - (z - cK/2)^2;

LeadingSingularityLocations011[] := {
  cK/2 + Sqrt[x2sq x3sq],
  cK/2 - Sqrt[x2sq x3sq]
};

LeadingSingularityBasisMatrix011[] := Module[
  {z = LeadingSingularityLocations011[]},
  {{-z[[1]], 1}, {-z[[2]], 1}}
];

CutLambda011[z_, aVariable_] := (
  cK (z - cK/2)
  + bK x2sq
  + aVariable x3sq
  + cK^2/4
  - aVariable bK
);

End[];
EndPackage[];
