(* ::Package:: *)

If[!MemberQ[$Packages, "CosmologicalBaikov`"],
  Get[FileNameJoin[{DirectoryName[$InputFileName],
    "CosmologicalBaikov.wl"}]]
];
If[!MemberQ[$Packages, "CosmologicalBaikovFlatTriangle`"],
  Get[FileNameJoin[{DirectoryName[$InputFileName],
    "CosmologicalBaikovFlatTriangle.wl"}]]
];
If[!MemberQ[$Packages, "CosmologicalBaikov111`"],
  Get[FileNameJoin[{DirectoryName[$InputFileName],
    "CosmologicalBaikov111.wl"}]]
];

BeginPackage[
  "CosmologicalBaikovFiniteNested`",
  {"CosmologicalBaikovFlatTriangle`"}
];

NMI::usage = "NMI[weight,r,p,q,v] is an inert nested outer master label. r,p,q are line choices (0,1,None) and v=0 or -1 labels the D4 numerator in a two-line sector.";
OneW::usage = "OneW labels the rational inner weight.";
LogYW::usage = "LogYW labels Ly=log(D10/D6).";
LogXW::usage = "LogXW labels Lx=log(D10/D5).";
TriangleW::usage = "TriangleW labels the finite four-dimensional massless triangle T111.";

D1o::usage = "D1o=(p3-p2)^2+x0sq is the original R-line denominator.";
D2o::usage = "D2o=p2^2+x2sq is the original P-line denominator.";
D3o::usage = "D3o=p3^2+x3sq is the Q-line denominator.";
D4o::usage = "D4o=p2 p3+c/2=Q2.Q3.";
D5o::usage = "D5o=p2^2+a=Q2^2.";
D6o::usage = "D6o=p3^2+b=Q3^2.";
D10o::usage = "D10o=(p2-p3)^2+a+b-c=(Q2-Q3)^2.";
LambdaOuter::usage = "LambdaOuter=D4o^2-D5o D6o is the flat Gram determinant.";

RationalOuterBasisA::usage = "RationalOuterBasisA[] returns the twelve ordinary-vacuum masters needed by the a connection after resolving the parallel pairs {D1,D10} and {D2,D5}.";
DirectImageBasis::usage = "DirectImageBasis[weight] returns four ordinary-vacuum masters for LogYW or LogXW, and the eleven square-root-twisted masters for TriangleW.";
FiniteNestedBasisA::usage = "FiniteNestedBasisA[] returns the 31-component finite nested basis ordered as {rational,LogY,LogX,triangle}.";
FiniteNestedBlockDimensionsA::usage = "FiniteNestedBlockDimensionsA[] returns {12,4,4,11}.";
FiniteNestedBlockRangesA::usage = "FiniteNestedBlockRangesA[] gives the four index ranges in the 31-component basis.";

DirectFlatWeightConnectionA::usage = "DirectFlatWeightConnectionA[] gives the exact integrand-level partial_a connection of {1,Ly,Lx,T111}.";
FiniteEquation31SourceA::usage = "FiniteEquation31SourceA[] returns the bubble part of equation (31) in the finite difference basis, including the overall factor 1/2.";
BubblePoleCancellationA::usage = "BubblePoleCancellationA[] verifies that the three unrenormalized bubble coefficients sum to zero.";
ParallelDenominatorIdentitiesA::usage = "ParallelDenominatorIdentitiesA[] returns the identities relating D10 to D1 and D5 to D2.";

OrdinaryVacuumMassConnectionsD1::usage = "OrdinaryVacuumMassConnectionsD1[mR,mP,mQ] returns the three d=1 mass-derivative connections in the basis {RP,RQ,PQ,RPQ}.";
RationalConnectionA::usage = "RationalConnectionA[] returns the explicit 12 by 12 epsilon=0 rational connection.";
LogYToRationalConnectionA::usage = "LogYToRationalConnectionA[] returns the explicit 4 by 12 partial_a map from the Ly masters to rational masters.";
LogXToRationalConnectionA::usage = "LogXToRationalConnectionA[] returns the explicit 4 by 12 partial_a map from the Lx masters to rational masters.";
TriangleHomogeneousConnectionA::usage = "TriangleHomogeneousConnectionA[] returns the eleven-dimensional eps=0 homogeneous triangle block, equal to the earlier square-root-twisted connection.";
FiniteConnectionTemplateA::usage = "FiniteConnectionTemplateA[] returns the 31 by 31 lower-block-triangular connection with all blocks explicit except the triangle-to-lower-weight source blocks.";
FiniteConnectionStrictUpperBlockCheckA::usage = "FiniteConnectionStrictUpperBlockCheckA[] verifies the forced zero blocks above the finite filtration.";
FiniteTriangularConnectionTemplateA::usage = "FiniteTriangularConnectionTemplateA[] returns the epsilon=0 template after the known eleven-dimensional triangle leading-singularity transformation; the remaining transformed source blocks have inert named entries.";
FiniteTriangularStrictUpperCheckA::usage = "FiniteTriangularStrictUpperCheckA[] verifies that every entry strictly above the diagonal of the transformed 31 by 31 template vanishes.";
TransformFiniteConnectionA::usage = "TransformFiniteConnectionA[connection,transform,var] applies a basis change including the derivative gauge term.";
OuterExponentAudit::usage = "OuterExponentAudit[] records why Lambda^j in equation (24) is correct: T111 itself supplies the effective Lambda^(-1/2) factor only in the triangle sector.";
FiniteReductionAssumptionsA::usage = "FiniteReductionAssumptionsA[] lists the nonvanishing quantities assumed by the finite basis and parallel-denominator resolution.";

Begin["`Private`"];

ClearAll[NMI, OneW, LogYW, LogXW, TriangleW];
ClearAll[D1o, D2o, D3o, D4o, D5o, D6o, D10o, LambdaOuter];

(* With the printed Lambda^j measure, the rational sector is an ordinary
   two-loop vacuum family.  Empty and one-line sectors are scaleless in
   dimensional regularization.  Each two-line sector factorizes and each
   three-line family has one sunrise master. *)
RationalOuterBasisA[] := Join[
  Flatten[Table[NMI[OneW, r, p, None, 0],
    {r, 0, 1}, {p, 0, 1}]],
  Table[NMI[OneW, r, None, 0, 0], {r, 0, 1}],
  Table[NMI[OneW, None, p, 0, 0], {p, 0, 1}],
  Flatten[Table[NMI[OneW, r, p, 0, 0],
    {r, 0, 1}, {p, 0, 1}]]
];

ordinaryWeightedBasis[weight_] := {
  NMI[weight, 0, 0, None, 0],
  NMI[weight, 0, None, 0, 0],
  NMI[weight, None, 0, 0, 0],
  NMI[weight, 0, 0, 0, 0]
};

triangleWeightedBasis[] := {
  NMI[weight, None, None, None, 0],
  NMI[weight, None, 0, None, 0],
  NMI[weight, 0, None, None, 0],
  NMI[weight, None, None, 0, 0],
  NMI[weight, 0, 0, None, 0],
  NMI[weight, 0, 0, None, -1],
  NMI[weight, 0, None, 0, 0],
  NMI[weight, 0, None, 0, -1],
  NMI[weight, None, 0, 0, 0],
  NMI[weight, None, 0, 0, -1],
  NMI[weight, 0, 0, 0, 0]
} /. weight -> TriangleW;

DirectImageBasis[LogYW] := ordinaryWeightedBasis[LogYW];
DirectImageBasis[LogXW] := ordinaryWeightedBasis[LogXW];
DirectImageBasis[TriangleW] := triangleWeightedBasis[];

FiniteNestedBasisA[] := Join[
  RationalOuterBasisA[],
  DirectImageBasis[LogYW],
  DirectImageBasis[LogXW],
  DirectImageBasis[TriangleW]
];

FiniteNestedBlockDimensionsA[] := {12, 4, 4, 11};

FiniteNestedBlockRangesA[] := AssociationThread[
  {"Rational", "LogY", "LogX", "Triangle"},
  {Range[1, 12], Range[13, 16], Range[17, 20], Range[21, 31]}
];

DirectFlatWeightConnectionA[] := {
  {0, 0, 0, 0},
  {1/D10o, 0, 0, 0},
  {1/D10o - 1/D5o, 0, 0, 0},
  {
    0,
    (D4o/D5o + (D6o - D4o)/D10o)/(2 LambdaOuter),
    (D4o - D6o)/(2 D10o LambdaOuter),
    D6o/(2 LambdaOuter)
  }
};

FiniteEquation31SourceA[] := 1/2 (
  (D4o/D5o + (D6o - D4o)/D10o) LogYW
  + (D4o - D6o)/D10o LogXW
);

BubblePoleCancellationA[] := Together[Total[{
  -D4o/D5o,
  D4o/D5o + (D6o - D4o)/D10o,
  (D4o - D6o)/D10o
}]];

ParallelDenominatorIdentitiesA[] := {
  D10o == D1o + (CosmologicalBaikov`aK + CosmologicalBaikov`bK
    - CosmologicalBaikov`cK - CosmologicalBaikov`x0sq),
  D5o == D2o + (CosmologicalBaikov`aK
    - CosmologicalBaikov`x2sq)
};

ClearAll[ordinaryDelta, ordinaryVacuumMassConnectionsD1];

ordinaryDelta[mR_, mP_, mQ_] :=
  mR^2 + mP^2 + mQ^2 - 2 mR mP - 2 mR mQ - 2 mP mQ;

OrdinaryVacuumMassConnectionsD1[mR_, mP_, mQ_] := Module[
  {delta, dR, dP, dQ},
  delta = ordinaryDelta[mR, mP, mQ];

  dR = ConstantArray[0, {4, 4}];
  dR[[1, 1]] = -1/(2 mR);
  dR[[2, 2]] = -1/(2 mR);
  dR[[4]] = {
    -(mR - mP + mQ)/(2 mR delta),
    -(mR + mP - mQ)/(2 mR delta),
    1/delta,
    -2 (mR - mP - mQ)/delta
  };

  dP = ConstantArray[0, {4, 4}];
  dP[[1, 1]] = -1/(2 mP);
  dP[[3, 3]] = -1/(2 mP);
  dP[[4]] = {
    (mR - mP - mQ)/(2 mP delta),
    1/delta,
    -(mR + mP - mQ)/(2 mP delta),
    2 (mR - mP + mQ)/delta
  };

  dQ = ConstantArray[0, {4, 4}];
  dQ[[2, 2]] = -1/(2 mQ);
  dQ[[3, 3]] = -1/(2 mQ);
  dQ[[4]] = {
    1/delta,
    (mR - mP - mQ)/(2 mQ delta),
    (-mR + mP - mQ)/(2 mQ delta),
    2 (mR + mP - mQ)/delta
  };

  <|"R" -> dR, "P" -> dP, "Q" -> dQ|>
];

ClearAll[rpIndex, rqIndex, pqIndex, sunriseIndex];
rpIndex[r_Integer, p_Integer] := 1 + 2 r + p;
rqIndex[r_Integer] := 5 + r;
pqIndex[p_Integer] := 7 + p;
sunriseIndex[r_Integer, p_Integer] := 9 + 2 r + p;

RationalConnectionA[] := Module[
  {matrix, rMasses, pMasses, qMass, local, localRow, columns},
  matrix = ConstantArray[0, {12, 12}];
  rMasses = {
    CosmologicalBaikov`x0sq,
    CosmologicalBaikov`aK + CosmologicalBaikov`bK
      - CosmologicalBaikov`cK
  };
  pMasses = {CosmologicalBaikov`x2sq, CosmologicalBaikov`aK};
  qMass = CosmologicalBaikov`x3sq;

  Do[
    matrix[[rpIndex[r, p], rpIndex[r, p]]] =
      -r/(2 rMasses[[r + 1]]) - p/(2 pMasses[[p + 1]]),
    {r, 0, 1}, {p, 0, 1}
  ];
  Do[
    matrix[[rqIndex[r], rqIndex[r]]] = -r/(2 rMasses[[r + 1]]),
    {r, 0, 1}
  ];
  Do[
    matrix[[pqIndex[p], pqIndex[p]]] = -p/(2 pMasses[[p + 1]]),
    {p, 0, 1}
  ];

  Do[
    local = OrdinaryVacuumMassConnectionsD1[
      rMasses[[r + 1]], pMasses[[p + 1]], qMass
    ];
    localRow = r local["R"][[4]] + p local["P"][[4]];
    columns = {
      rpIndex[r, p], rqIndex[r], pqIndex[p], sunriseIndex[r, p]
    };
    matrix[[sunriseIndex[r, p], columns]] = localRow,
    {r, 0, 1}, {p, 0, 1}
  ];
  Map[Together, matrix, {2}]
];

LogYToRationalConnectionA[] := Module[{matrix, deltaR},
  matrix = ConstantArray[0, {4, 12}];
  deltaR = CosmologicalBaikov`aK + CosmologicalBaikov`bK
    - CosmologicalBaikov`cK - CosmologicalBaikov`x0sq;
  matrix[[1, {rpIndex[0, 0], rpIndex[1, 0]}]] = {1, -1}/deltaR;
  matrix[[2, {rqIndex[0], rqIndex[1]}]] = {1, -1}/deltaR;
  matrix[[3, sunriseIndex[1, 0]]] = 1;
  matrix[[4, {sunriseIndex[0, 0], sunriseIndex[1, 0]}]] =
    {1, -1}/deltaR;
  matrix
];

LogXToRationalConnectionA[] := Module[{matrix, deltaR, deltaP},
  matrix = ConstantArray[0, {4, 12}];
  deltaR = CosmologicalBaikov`aK + CosmologicalBaikov`bK
    - CosmologicalBaikov`cK - CosmologicalBaikov`x0sq;
  deltaP = CosmologicalBaikov`aK - CosmologicalBaikov`x2sq;

  matrix[[1, rpIndex[0, 0]]] = 1/deltaR - 1/deltaP;
  matrix[[1, rpIndex[1, 0]]] = -1/deltaR;
  matrix[[1, rpIndex[0, 1]]] = 1/deltaP;

  matrix[[2, {rqIndex[0], rqIndex[1]}]] = {1, -1}/deltaR;
  matrix[[2, sunriseIndex[0, 1]]] = -1;

  matrix[[3, sunriseIndex[1, 0]]] = 1;
  matrix[[3, {pqIndex[0], pqIndex[1]}]] = {-1, 1}/deltaP;

  matrix[[4, sunriseIndex[0, 0]]] = 1/deltaR - 1/deltaP;
  matrix[[4, sunriseIndex[1, 0]]] = -1/deltaR;
  matrix[[4, sunriseIndex[0, 1]]] = 1/deltaP;
  matrix
];

TriangleHomogeneousConnectionA[] :=
  CosmologicalBaikov111`CompactConnectionEpsilonZero[];

ClearAll[ctr, cty, ctx];

FiniteConnectionTemplateA[] := Module[
  {arBlock, yToR, xToR, tToR, tToY, tToX, atBlock, zero},
  arBlock = RationalConnectionA[];
  yToR = LogYToRationalConnectionA[];
  xToR = LogXToRationalConnectionA[];
  tToR = Array[ctr, {11, 12}];
  tToY = Array[cty, {11, 4}];
  tToX = Array[ctx, {11, 4}];
  atBlock = TriangleHomogeneousConnectionA[];
  zero[m_, n_] := ConstantArray[0, {m, n}];
  ArrayFlatten[{
    {arBlock, zero[12, 4], zero[12, 4], zero[12, 11]},
    {yToR, zero[4, 4], zero[4, 4], zero[4, 11]},
    {xToR, zero[4, 4], zero[4, 4], zero[4, 11]},
    {tToR, tToY, tToX, atBlock}
  }]
];

FiniteConnectionStrictUpperBlockCheckA[] := Module[{m, r},
  m = FiniteConnectionTemplateA[];
  r = Values[FiniteNestedBlockRangesA[]];
  And[
    m[[r[[1]], Join[r[[2]], r[[3]], r[[4]]]]] ===
      ConstantArray[0, {12, 19}],
    m[[r[[2]], Join[r[[3]], r[[4]]]]] ===
      ConstantArray[0, {4, 15}],
    m[[r[[3]], r[[4]]]] === ConstantArray[0, {4, 11}],
    m[[r[[2]], r[[3]]]] === ConstantArray[0, {4, 4}],
    m[[r[[3]], r[[2]]]] === ConstantArray[0, {4, 4}]
  ]
];

ClearAll[ctrTri, ctyTri, ctxTri];

FiniteTriangularConnectionTemplateA[] := Module[
  {arBlock, yToR, xToR, tToR, tToY, tToX, atBlock, zero},
  arBlock = RationalConnectionA[];
  yToR = LogYToRationalConnectionA[];
  xToR = LogXToRationalConnectionA[];
  tToR = Array[ctrTri, {11, 12}];
  tToY = Array[ctyTri, {11, 4}];
  tToX = Array[ctxTri, {11, 4}];
  atBlock = CosmologicalBaikov111`CompactTriangularConnectionEpsilonZero[];
  zero[m_, n_] := ConstantArray[0, {m, n}];
  ArrayFlatten[{
    {arBlock, zero[12, 4], zero[12, 4], zero[12, 11]},
    {yToR, zero[4, 4], zero[4, 4], zero[4, 11]},
    {xToR, zero[4, 4], zero[4, 4], zero[4, 11]},
    {tToR, tToY, tToX, atBlock}
  }]
];

FiniteTriangularStrictUpperCheckA[] := Module[{matrix},
  matrix = FiniteTriangularConnectionTemplateA[];
  And @@ Flatten@Table[
    If[j > i, Together[matrix[[i, j]]] === 0, True],
    {i, 31}, {j, 31}
  ]
];

TransformFiniteConnectionA[connection_?MatrixQ, transform_?MatrixQ, var_] :=
  Map[Together, D[transform, var].Inverse[transform]
    + transform.connection.Inverse[transform], {2}];

OuterExponentAudit[] := <|
  "Equation24PrintedExponent" -> 0,
  "Equation24Status" -> "Correct",
  "EffectiveSectorExponents" -> <|
    "Rational" -> 0, "LogY" -> 0, "LogX" -> 0,
    "Triangle" -> 1/2
  |>,
  "Reason" -> "The evaluated four-dimensional triangle T111 already contains its leading-singularity factor Lambda^(-1/2). Therefore no additional half-power belongs in the outer measure of equation (24)."
|>;

FiniteReductionAssumptionsA[] := {
  D5o != 0, D10o != 0, LambdaOuter != 0,
  CosmologicalBaikov`aK - CosmologicalBaikov`x2sq != 0,
  CosmologicalBaikov`aK + CosmologicalBaikov`bK
    - CosmologicalBaikov`cK - CosmologicalBaikov`x0sq != 0
};

End[];
EndPackage[];
