(* Build the executable I011 reduction notebook. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "CosmologicalBaikov.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "CosmologicalBaikov011.wl"}]];

ClearAll[inputCell];
SetAttributes[inputCell, HoldAllComplete];
inputCell[expr_] := Cell[
  BoxData[ToBoxes[Defer[expr], StandardForm]],
  "Input"
];

notebook = Notebook[
  {
    Cell["IBP reduction and differential equations for the I011 pinch", "Title"],
    Cell[
      "This notebook treats the genuinely distinct D2-D3 two-denominator sector. It derives a five-master basis, constructs four vector IBPs in the variables {D2,D3,D4}, reduces the raised powers and D4 numerators, and assembles the exact matrix in partial_a m011=A011 m011.",
      "Text"
    ],

    Cell["1. Load the implementation", "Section"],
    inputCell[
      ClearAll["Global`*"];
      notebookDirectory = Quiet@Check[
        NotebookDirectory[],
        DirectoryName[$InputFileName]
      ];
      Get[FileNameJoin[{notebookDirectory, "CosmologicalBaikov.wl"}]];
      Get[FileNameJoin[{notebookDirectory, "CosmologicalBaikov011.wl"}]];
    ],
    inputCell[
      basis011 = CosmologicalBaikov011`MasterBasis011[];
      basis011
    ],
    Cell[
      "The notation is J011[s,w,v]=Ibar^(0 s w, v 0 0;0). Thus the basis is {I000,000;0, I010,000;0, I001,000;0, I011,000;0, I011,-100;0}. The last integral has one D4 numerator.",
      "Text"
    ],

    Cell["2. Baikov variables and the a derivative", "Section"],
    Cell[
      "For this pinch, p2^2=D2-x2^2, p3^2=D3-x3^2 and p2.p3=D4-c/2. In particular D6=p3^2+b=D3+b-x3^2. Differentiating Lambda^(-1/2) with respect to a inserts D6 Lambda^(-3/2)/2, after which the scaling identity removes the extra inverse Lambda.",
      "Text"
    ],
    inputCell[
      scalingExample =
        CosmologicalBaikov011`ScalingIdentity011[1, 1, 0];
      derivativeExample =
        CosmologicalBaikov011`ParameterDerivativeA011[1, 1, 0];
      {scalingExample, derivativeExample}
    ],

    Cell["3. Vector and tangential IBPs", "Section"],
    Cell[
      "The four generators are p2.partial_p2, p3.partial_p2, p2.partial_p3 and p3.partial_p3. Their Lambda contractions are respectively c z-2 b x, c y-2 b z, c x-2 a z and c z-2 a y. The same tangential combination as in the I110 sector annihilates Lambda.",
      "Text"
    ],
    inputCell[
      ibpsAtTop =
        CosmologicalBaikov011`IBPIdentity011[#, 1, 1, 0] & /@
          {"pp", "pq", "qp", "qq"};
      tangentialCheck = CosmologicalBaikov011`TangentialCheck011[];
      {Length[ibpsAtTop], tangentialCheck}
    ],

    Cell["4. Staged reduction", "Section"],
    Cell[
      "The one-denominator Hermite reductions remove J011[2,0,0] and J011[0,2,0] and their D4-numerator analogues. Four top-sector IBPs then solve for {J211,-1,J211,0,J120,0,J120,-1}. Three independent IBPs at v=-1 remove the remaining D4^2 numerator and raised indices; the fourth is a consistency check.",
      "Text"
    ],
    inputCell[
      reductionData011 = CosmologicalBaikov011`ReductionData011[];
      regressionPoint = {
        CosmologicalBaikov`aK -> 5,
        CosmologicalBaikov`bK -> -1,
        CosmologicalBaikov`cK -> 1,
        CosmologicalBaikov`x2sq -> 1,
        CosmologicalBaikov`x3sq -> 2,
        CosmologicalBaikov`eps -> 1/7
      };
      mixedRanks = MatrixRank[# /. regressionPoint] & /@
        reductionData011["MixedMatrices"];
      fourthResidual = Together[
        reductionData011["FourthIBPResidual"] /. regressionPoint
      ];
      {
        reductionData011["SubtopologyRules"],
        mixedRanks,
        fourthResidual
      }
    ],

    Cell["5. Differential matrix", "Section"],
    inputCell[
      matrix011 = CosmologicalBaikov011`DifferentialMatrixA011[];
      lowerRows011 = Map[Together, matrix011[[1 ;; 3]], {2}];
      {
        Dimensions[matrix011],
        MatrixForm[lowerRows011]
      }
    ],
    Cell[
      "The first master and the D3-only master both obey partial_a I=(b epsilon/lambda0) I. The D2-only row is the same one-denominator equation appearing in the I110 reduction. The exact fourth and fifth rows are retained as rational expressions and can be simplified entrywise with DifferentialMatrixA011[\"SimplifyFunction\"->Together].",
      "Text"
    ],
    inputCell[
      matrixAtRegressionPoint =
        CosmologicalBaikov011`NumericalCrossCheck011[];
      MatrixForm[matrixAtRegressionPoint]
    ],

    Cell["6. Automated checks", "Section"],
    inputCell[
      expectedMatrix = {
        {-4/147, 0, 0, 0, 0},
        {-3/280, -311/2940, 0, 0, 0},
        {0, 0, -4/147, 0, 0},
        {-165/45038, -895/472899, 976/472899, -2526/22519, 264/22519},
        {-72/22519, -1366/472899, -184/157633, 954/22519, -3630/22519}
      };
      {
        VerificationTest[
          matrixAtRegressionPoint,
          expectedMatrix,
          TestID -> "I011-exact-matrix-regression"
        ],
        VerificationTest[
          FreeQ[matrix011, _CosmologicalBaikov011`J011],
          True,
          TestID -> "I011-no-unreduced-integrals"
        ]
      }
    ],
    Cell[
      "An independent global reduction was also performed at the regression point using 400 IBP equations over 223 integrals. Its rational row reduction leaves the same five masters (plus four seed-boundary artifacts) and reproduces all 25 matrix entries shown above.",
      "Text"
    ],

    Cell["7. Usage and genericity", "Section"],
    inputCell[
      {
        CosmologicalBaikov011`GenericityConditions011[],
        Keys[reductionData011]
      }
    ],
    Cell[
      "The displayed rules are generic. In addition to the elementary nonzero factors returned above, the two mixed-sector IBP determinants must not vanish. On a singular locus the master basis should be recomputed rather than obtained by substituting into a formula with a zero denominator.",
      "Text"
    ]
  },
  WindowTitle -> "Cosmological Baikov: I011 reduction",
  StyleDefinitions -> "Default.nb"
];

target = FileNameJoin[{
  DirectoryName[$InputFileName],
  "CosmologicalBaikov011Reduction.nb"
}];
Put[notebook, target];
Print[target];
