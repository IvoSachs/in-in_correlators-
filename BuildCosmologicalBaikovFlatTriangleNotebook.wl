(* Build the regulator-independent flat-triangle extension notebook. *)

Get[FileNameJoin[{
  DirectoryName[$InputFileName],
  "CosmologicalBaikovFlatTriangle.wl"
}]];

ClearAll[inputCell];
SetAttributes[inputCell, HoldAllComplete];
inputCell[expr_] := Cell[
  BoxData[ToBoxes[Defer[expr], StandardForm]],
  "Input"
];

notebook = Notebook[
  {
    Cell["Finite flat-triangle extension of the cosmological Baikov system", "Title"],
    Cell[
      "This notebook extracts the part of Section 2 that is independent of the treatment of the UV-divergent flat-space bubbles. It replaces the three bubble subtopologies by two finite bubble differences, constructs the resulting four-component flat connection, proves its integrability, and identifies the lower-block-triangular structure that the complete nested p2,p3 reduction must have.",
      "Text"
    ],

    Cell["1. Load definitions", "Section"],
    inputCell[
      ClearAll["Global`*"];
      notebookDirectory = Quiet@Check[
        NotebookDirectory[],
        DirectoryName[$InputFileName]
      ];
      Get[FileNameJoin[{
        notebookDirectory,
        "CosmologicalBaikovFlatTriangle.wl"
      }]];
    ],
    inputCell[
      invariants = <|
        "x" -> CosmologicalBaikovFlatTriangle`xF,
        "y" -> CosmologicalBaikovFlatTriangle`yF,
        "z" -> CosmologicalBaikovFlatTriangle`zF,
        "s" -> CosmologicalBaikovFlatTriangle`sF,
        "Gram" -> CosmologicalBaikovFlatTriangle`lambdaF
      |>;
      invariants
    ],

    Cell["2. Cancellation of the common bubble divergence", "Section"],
    Cell[
      "Let x=Q2^2, y=Q3^2, z=Q2.Q3 and s=(Q2-Q3)^2. The coefficients of {T011,T101,T110} in the partial_a row of equation (30) add to zero. Since every massless off-shell bubble has the same UV pole, the pole and scheme-dependent additive constant cancel in equation (31).",
      "Text"
    ],
    inputCell[
      bubbleCoefficientsA =
        CosmologicalBaikovFlatTriangle`FlatBubbleCoefficientsA[];
      {
        bubbleCoefficientsA,
        VerificationTest[
          Together[Total[bubbleCoefficientsA]],
          0,
          TestID -> "common-bubble-piece-cancels"
        ]
      }
    ],

    Cell["3. Finite flat-space basis", "Section"],
    Cell[
      "Write B(u)=common-log(u) for a dimensionally regulated bubble. The scheme-independent differences are Ly=B(y)-B(s)=log(s/y) and Lx=B(x)-B(s)=log(s/x). Thus the finite flat basis is {1,Ly,Lx,T111}; three separate bubble masters are unnecessary.",
      "Text"
    ],
    inputCell[
      flatBasis =
        CosmologicalBaikovFlatTriangle`FlatFiniteBasis[];
      flatBasis
    ],

    Cell["4. Flat connection", "Section"],
    Cell[
      "The order {1,Ly,Lx,T111} makes every flat-space connection lower triangular. Here partial_c means varying c=2 y2.y3, hence partial_c=(1/2) partial_z.",
      "Text"
    ],
    inputCell[
      flatConnectionA =
        CosmologicalBaikovFlatTriangle`FlatConnectionA[];
      flatConnectionB =
        CosmologicalBaikovFlatTriangle`FlatConnectionY[];
      flatConnectionC =
        CosmologicalBaikovFlatTriangle`FlatConnectionC[];
      {
        MatrixForm[flatConnectionA],
        MatrixForm[flatConnectionB],
        MatrixForm[flatConnectionC]
      }
    ],
    inputCell[
      curvatureChecks =
        CosmologicalBaikovFlatTriangle`FlatCurvatureChecks[];
      VerificationTest[
        curvatureChecks,
        ConstantArray[0, {3, 4, 4}],
        TestID -> "flat-connection-zero-curvature"
      ]
    ],

    Cell["5. Finite form of the new source in equation (31)", "Section"],
    Cell[
      "After cancellation of the common bubble, the last line of equation (31) contains only Cy Ly+Cx Lx. The expression below includes the overall factor 1/2 appearing in that equation.",
      "Text"
    ],
    inputCell[
      sourceCoefficientsA =
        CosmologicalBaikovFlatTriangle`FlatSourceCoefficientsA[];
      finiteSourceIntegrandA =
        CosmologicalBaikovFlatTriangle`FlatSourceIntegrandA[];
      {sourceCoefficientsA, finiteSourceIntegrandA}
    ],

    Cell["6. Structure of the extended nested connection", "Section"],
    Cell[
      "After the remaining p2,p3 IBP reduction, order the direct-image masters as {rational sector, Ly sector, Lx sector, triangle sector}. Differentiation lowers this filtration, so the complete matrix must have the block pattern shown below. The previous eleven-master matrix is the homogeneous triangle block; the new source maps fill only the blocks below it.",
      "Text"
    ],
    inputCell[
      MatrixForm[
        CosmologicalBaikovFlatTriangle`ExtendedBlockPattern[]
      ]
    ],

    Cell["7. Required convention for the remaining outer IBP reduction", "Section"],
    Cell[
      "Equation (28) writes d^4 L, but T011, T101 and T110 are UV divergent in exactly four dimensions. A unique explicit outer master list and numerical connection therefore require one of two conventions: (i) use d_L=4-2 delta and retain delta until the three bubble sectors have been combined, or (ii) use the finite difference basis {1,Ly,Lx,T111} above from the outset. These lead to different intermediate master lists, although their finite triangle equation agrees. The second convention is smaller and is recommended for an epsilon->0 triangular system.",
      "Text"
    ]
  },
  WindowTitle -> "Cosmological Baikov: flat-triangle extension",
  StyleDefinitions -> "Default.nb"
];

target = FileNameJoin[{
  DirectoryName[$InputFileName],
  "CosmologicalBaikovFlatTriangleExtension.nb"
}];
Put[notebook, target];
Print[target];
