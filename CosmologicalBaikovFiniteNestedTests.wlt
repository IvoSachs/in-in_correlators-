Get[FileNameJoin[{DirectoryName[$TestFileName],
  "CosmologicalBaikovFiniteNested.wl"}]];

VerificationTest[
  Length[CosmologicalBaikovFiniteNested`RationalOuterBasisA[]],
  12,
  TestID -> "finite-rational-basis-count"
]

VerificationTest[
  Length[CosmologicalBaikovFiniteNested`FiniteNestedBasisA[]],
  31,
  TestID -> "finite-nested-basis-count"
]

VerificationTest[
  DuplicateFreeQ[CosmologicalBaikovFiniteNested`FiniteNestedBasisA[]],
  True,
  TestID -> "finite-nested-basis-has-no-duplicates"
]

VerificationTest[
  CosmologicalBaikovFiniteNested`BubblePoleCancellationA[],
  0,
  TestID -> "finite-bubble-common-pole-cancels"
]

VerificationTest[
  Dimensions[CosmologicalBaikovFiniteNested`DirectFlatWeightConnectionA[]],
  {4, 4},
  TestID -> "finite-flat-direct-connection-dimensions"
]

VerificationTest[
  Dimensions[CosmologicalBaikovFiniteNested`FiniteConnectionTemplateA[]],
  {31, 31},
  TestID -> "finite-connection-template-dimensions"
]

VerificationTest[
  CosmologicalBaikovFiniteNested`FiniteConnectionStrictUpperBlockCheckA[],
  True,
  TestID -> "finite-connection-template-filtration"
]

VerificationTest[
  CosmologicalBaikovFiniteNested`OuterExponentAudit[][
    "Equation24Status"],
  "Correct",
  TestID -> "outer-exponent-audit"
]

VerificationTest[
  Dimensions[CosmologicalBaikovFiniteNested`RationalConnectionA[]],
  {12, 12},
  TestID -> "ordinary-rational-connection-dimensions"
]

VerificationTest[
  And @@ Flatten@Table[
    If[j > i,
      CosmologicalBaikovFiniteNested`RationalConnectionA[][[i, j]] === 0,
      True
    ],
    {i, 12}, {j, 12}
  ],
  True,
  TestID -> "ordinary-rational-connection-is-lower-triangular"
]

VerificationTest[
  Module[{mR, mP, mQ, matrices, ar, ap, aq, curvatures},
    matrices = CosmologicalBaikovFiniteNested`OrdinaryVacuumMassConnectionsD1[
      mR, mP, mQ];
    ar = matrices["R"];
    ap = matrices["P"];
    aq = matrices["Q"];
    curvatures = {
      D[ap, mR] - D[ar, mP] + ap.ar - ar.ap,
      D[aq, mR] - D[ar, mQ] + aq.ar - ar.aq,
      D[aq, mP] - D[ap, mQ] + aq.ap - ap.aq
    };
    Map[Together, curvatures, {3}]
  ],
  ConstantArray[0, {3, 4, 4}],
  TestID -> "ordinary-vacuum-mass-connection-zero-curvature"
]

VerificationTest[
  Dimensions[CosmologicalBaikovFiniteNested`LogYToRationalConnectionA[]],
  {4, 12},
  TestID -> "log-y-source-block-dimensions"
]

VerificationTest[
  Dimensions[CosmologicalBaikovFiniteNested`LogXToRationalConnectionA[]],
  {4, 12},
  TestID -> "log-x-source-block-dimensions"
]

VerificationTest[
  Dimensions[CosmologicalBaikovFiniteNested`TriangleHomogeneousConnectionA[]],
  {11, 11},
  TestID -> "triangle-homogeneous-block-dimensions"
]

VerificationTest[
  Dimensions[CosmologicalBaikovFiniteNested`FiniteTriangularConnectionTemplateA[]],
  {31, 31},
  TestID -> "finite-triangular-template-dimensions"
]

VerificationTest[
  CosmologicalBaikovFiniteNested`FiniteTriangularStrictUpperCheckA[],
  True,
  TestID -> "finite-triangular-template-strict-upper-zero"
]
