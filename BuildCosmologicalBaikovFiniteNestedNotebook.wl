Get[FileNameJoin[{
  DirectoryName[$InputFileName],
  "CosmologicalBaikovFiniteNested.wl"
}]];

ClearAll[inputCell];
SetAttributes[inputCell, HoldAllComplete];
inputCell[expr_] := Cell[
  BoxData[ToBoxes[Defer[expr], StandardForm]],
  "Input"
];

notebook = Notebook[
  {
    Cell["Finite bubble-difference basis for the nested triangle", "Title"],
    Cell[
      "This notebook implements the finite convention selected for Section 2. It proves the cancellation of the common bubble divergence, uses Ly=log(D10/D6) and Lx=log(D10/D5), constructs the finite nested outer basis required by the partial_a orbit, and exposes the forced lower-block-triangular connection. Exact derived entries and still-unreduced outer-IBP slots are deliberately kept distinct.",
      "Text"
    ],

    Cell["1. Load the implementation", "Section"],
    inputCell[
      ClearAll["Global`*"];
      notebookDirectory = Quiet@Check[
        NotebookDirectory[], DirectoryName[$InputFileName]
      ];
      Get[FileNameJoin[{
        notebookDirectory,
        "CosmologicalBaikovFiniteNested.wl"
      }]];
    ],

    Cell["2. Why equation (24) has no extra half-power", "Section"],
    Cell[
      "Equation (24) is correct as printed. The evaluated four-dimensional triangle T111 already contains its leading-singularity factor Lambda^(-1/2). Indeed its homogeneous derivative is partial_a T111=(Q3^2/(2 Lambda)) T111+..., exactly the derivative of Lambda^(-1/2). The rational and logarithmic sectors therefore have outer exponent zero, while only the triangle sector has effective exponent 1/2.",
      "Text"
    ],
    inputCell[
      exponentAudit =
        CosmologicalBaikovFiniteNested`OuterExponentAudit[];
      exponentAudit
    ],

    Cell["3. Finite bubble differences", "Section"],
    Cell[
      "Let x=D5=Q2^2, y=D6=Q3^2, z=D4=Q2.Q3 and s=D10=(Q2-Q3)^2. The coefficients of {T011,T101,T110} in the last line of equation (31) sum to zero. Hence the common UV pole and the additive bubble-scheme constant cancel before the p2,p3 integration.",
      "Text"
    ],
    inputCell[
      poleCancellation =
        CosmologicalBaikovFiniteNested`BubblePoleCancellationA[];
      VerificationTest[
        poleCancellation,
        0,
        TestID -> "common-bubble-pole-and-constant-cancel"
      ]
    ],
    Cell[
      "The two finite weights are Ly=T101-T011=log(D10/D6) and Lx=T110-T011=log(D10/D5). The exact finite form of the new term in equation (31) is shown below.",
      "Text"
    ],
    inputCell[
      finiteEquation31Source =
        CosmologicalBaikovFiniteNested`FiniteEquation31SourceA[];
      finiteEquation31Source
    ],

    Cell["4. Exact flat-weight connection", "Section"],
    Cell[
      "In the order {1,Ly,Lx,T111}, partial_a acts as partial_x at fixed p2,p3,b,c. This matrix is exact and lower triangular. It is the local input for the lifted outer IBP system.",
      "Text"
    ],
    inputCell[
      directFlatConnectionA =
        CosmologicalBaikovFiniteNested`DirectFlatWeightConnectionA[];
      MatrixForm[directFlatConnectionA]
    ],

    Cell["5. Why shifted outer denominators are required", "Section"],
    Cell[
      "Differentiating Ly and Lx generates 1/D10 and 1/D5. These denominators are parallel to D1 and D2, respectively. Products are resolved by partial fractions, but both line choices must be retained as distinct mass families unless their mass differences vanish.",
      "Text"
    ],
    inputCell[
      parallelDenominatorIdentities =
        CosmologicalBaikovFiniteNested`ParallelDenominatorIdentitiesA[];
      parallelDenominatorIdentities
    ],

    Cell["6. Finite extended basis for the partial_a orbit", "Section"],
    Cell[
      "The rational part is the union of the four ordinary-vacuum families R_i P_j Q_0, with R_0=D1, R_1=D10, P_0=D2, P_1=D5 and Q_0=D3. Empty and one-line sectors are scaleless in dimensional regularization. There are eight factorized two-line masters and four sunrise masters, hence twelve rational masters. The diagonal logarithmic problem is also an ordinary vacuum family, so each finite logarithm needs only {RP,RQ,PQ,RPQ}. The triangle retains the eleven square-root-twisted masters derived previously. The finite basis therefore has 12+4+4+11=31 components.",
      "Text"
    ],
    inputCell[
      rationalBasis =
        CosmologicalBaikovFiniteNested`RationalOuterBasisA[];
      finiteNestedBasis =
        CosmologicalBaikovFiniteNested`FiniteNestedBasisA[];
      blockDimensions =
        CosmologicalBaikovFiniteNested`FiniteNestedBlockDimensionsA[];
      {
        "rational count" -> Length[rationalBasis],
        "total count" -> Length[finiteNestedBasis],
        "block dimensions" -> blockDimensions,
        "duplicate free" -> DuplicateFreeQ[finiteNestedBasis]
      }
    ],
    inputCell[
      Dataset[AssociationThread[
        Range[Length[finiteNestedBasis]],
        finiteNestedBasis
      ]]
    ],

    Cell["7. Connection and triangular filtration", "Section"],
    Cell[
      "In the order {rational,Ly,Lx,triangle}, the connection is lower block triangular. The 12 by 12 rational block is the d=1 ordinary-vacuum mass connection. Since the original propagator masses x0^2,x2^2,x3^2 are independent of a, the logarithmic diagonal blocks vanish; their derivatives are the explicit rational source maps below. The homogeneous triangle block is exactly the earlier eleven-master epsilon=0 connection. Only the three triangle-to-lower-weight blocks remain as named lifted-IBP slots.",
      "Text"
    ],
    inputCell[
      rationalConnectionA =
        CosmologicalBaikovFiniteNested`RationalConnectionA[];
      logYToRationalA =
        CosmologicalBaikovFiniteNested`LogYToRationalConnectionA[];
      logXToRationalA =
        CosmologicalBaikovFiniteNested`LogXToRationalConnectionA[];
      triangleHomogeneousA =
        CosmologicalBaikovFiniteNested`TriangleHomogeneousConnectionA[];
      {
        MatrixForm[rationalConnectionA],
        MatrixForm[logYToRationalA],
        MatrixForm[logXToRationalA],
        MatrixForm[triangleHomogeneousA]
      }
    ],
    inputCell[
      connectionTemplateA =
        CosmologicalBaikovFiniteNested`FiniteConnectionTemplateA[];
      {
        Dimensions[connectionTemplateA],
        CosmologicalBaikovFiniteNested`FiniteConnectionStrictUpperBlockCheckA[]
      }
    ],
    inputCell[
      blockRanges = Values[
        CosmologicalBaikovFiniteNested`FiniteNestedBlockRangesA[]
      ];
      MatrixForm@Table[
        If[Total[Abs[connectionTemplateA[[blockRanges[[i]],
              blockRanges[[j]]]]]] === 0, 0,
          Symbol["A" <> ToString[i] <> ToString[j]]],
        {i, 4}, {j, 4}
      ]
    ],

    Cell["8. Basis transformations", "Section"],
    Cell[
      "The rational block is already lower triangular: factorized bubbles precede the four sunrise masters. The logarithmic blocks lie below it and have zero diagonal. The only nontrivial leading-singularity transformation is therefore the eleven-dimensional triangle transformation derived earlier. For any transformation T(a), the helper evaluates (partial_a T) T^(-1)+T A T^(-1), including the necessary gauge term.",
      "Text"
    ],
    inputCell[
      triangularTemplateA =
        CosmologicalBaikovFiniteNested`FiniteTriangularConnectionTemplateA[];
      {
        Dimensions[triangularTemplateA],
        CosmologicalBaikovFiniteNested`FiniteTriangularStrictUpperCheckA[]
      }
    ],
    inputCell[
      transformRule = HoldForm[
        triangularA =
          CosmologicalBaikovFiniteNested`TransformFiniteConnectionA[
            fullA, leadingSingularityTransform, aK
          ]
      ];
      transformRule
    ],

    Cell["9. Verification summary", "Section"],
    inputCell[
      TestReport[{
        VerificationTest[Length[rationalBasis], 12,
          TestID -> "rational-sector-count"],
        VerificationTest[Length[finiteNestedBasis], 31,
          TestID -> "finite-total-count"],
        VerificationTest[DuplicateFreeQ[finiteNestedBasis], True,
          TestID -> "basis-is-duplicate-free"],
        VerificationTest[poleCancellation, 0,
          TestID -> "bubble-pole-cancellation"],
        VerificationTest[Dimensions[rationalConnectionA], {12, 12},
          TestID -> "rational-connection-dimensions"],
        VerificationTest[Dimensions[logYToRationalA], {4, 12},
          TestID -> "log-y-source-dimensions"],
        VerificationTest[Dimensions[logXToRationalA], {4, 12},
          TestID -> "log-x-source-dimensions"],
        VerificationTest[Dimensions[triangleHomogeneousA], {11, 11},
          TestID -> "triangle-block-dimensions"],
        VerificationTest[
          CosmologicalBaikovFiniteNested`FiniteTriangularStrictUpperCheckA[],
          True,
          TestID -> "strict-upper-triangular-check"],
        VerificationTest[
          CosmologicalBaikovFiniteNested`FiniteConnectionStrictUpperBlockCheckA[],
          True,
          TestID -> "lower-block-filtration"
        ]
      }]
    ],
    Cell[
      "Scope note. The finite flat connection, the twelve-dimensional rational block, both logarithm-to-rational blocks, and the eleven-dimensional homogeneous triangle block are explicit. The remaining symbols ctr, cty and ctx are the triangle-to-rational, triangle-to-Ly and triangle-to-Lx source blocks. They require a lifted outer IBP reduction of the inhomogeneous terms in the triangle differential equation and are intentionally not presented as already-reduced coefficients.",
      "Text"
    ]
  },
  WindowTitle -> "Cosmological Baikov: finite nested basis",
  StyleDefinitions -> "Default.nb"
];

target = FileNameJoin[{
  DirectoryName[$InputFileName],
  "CosmologicalBaikovFiniteNestedBasis.nb"
}];
Put[notebook, target];
Print[target];
