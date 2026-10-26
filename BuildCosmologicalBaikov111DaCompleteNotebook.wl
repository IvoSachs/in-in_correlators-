(* Build the complete partial_a I111 reduction and triangularization notebook. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "CosmologicalBaikov111.wl"}]];

ClearAll[inputCell];
SetAttributes[inputCell, HoldAllComplete];
inputCell[expr_] := Cell[
  BoxData[ToBoxes[Defer[expr], StandardForm]],
  "Input"
];

notebook = Notebook[
  {
    Cell["Complete partial_a connection for the I111 sector", "Title"],
    Cell[
      "This notebook combines the I110, I101, and I011 pinch reductions with the algebraic top-sector identity to obtain the complete eleven-master differential equation partial_a m=A_a m. It then constructs an algebraic leading-singularity basis which is lower triangular at epsilon=0 and a three-Riccati basis which is exactly lower triangular at finite epsilon.",
      "Text"
    ],

    Cell["1. Load the reduction packages", "Section"],
    inputCell[
      ClearAll["Global`*"];
      notebookDirectory = Quiet@Check[
        NotebookDirectory[],
        DirectoryName[$InputFileName]
      ];
      Get[FileNameJoin[{notebookDirectory, "CosmologicalBaikov111.wl"}]];
    ],

    Cell["2. Closed eleven-master basis", "Section"],
    Cell[
      "The basis is ordered by topology: the vacuum, the three one-denominator sectors, two masters for each two-denominator pinch, and finally I111. The fourth index of FBI denotes the power of D4; -1 therefore means one D4 numerator.",
      "Text"
    ],
    inputCell[
      fullMasters = CosmologicalBaikov111`FullMasterBasis[];
      masterLabels = {
        "I000", "I010", "I100", "I001",
        "I110", "I110[-D4]",
        "I101", "I101[-D4]",
        "I011", "I011[-D4]",
        "I111"
      };
      Thread[masterLabels -> fullMasters]
    ],

    Cell["3. Obtain the I101 pinch without a new IBP reduction", "Section"],
    Cell[
      "Set q2=p3-p2 and q3=p3. The I101 sector becomes an I011 sector with (a,b,c,x2^2,x3^2) mapped to (a+b-c,b,2b-c,x0^2,x3^2). The auxiliary numerator is D4'=D3-D4+b-x3^2, so the transformed numerator master equals I100+(b-x3^2) I101-I101[-D4].",
      "Text"
    ],
    inputCell[
      Clear[q2, q3, p2, p3];
      originalLambdaInQ =
        CosmologicalBaikov`cK (q3 - q2) q3
        - CosmologicalBaikov`bK (q3 - q2)^2
        - CosmologicalBaikov`aK q3^2
        + CosmologicalBaikov`cK^2/4
        - CosmologicalBaikov`aK CosmologicalBaikov`bK;
      mappedLambda =
        (2 CosmologicalBaikov`bK - CosmologicalBaikov`cK) q2 q3
        - CosmologicalBaikov`bK q2^2
        - (CosmologicalBaikov`aK + CosmologicalBaikov`bK -
            CosmologicalBaikov`cK) q3^2
        + (2 CosmologicalBaikov`bK - CosmologicalBaikov`cK)^2/4
        - (CosmologicalBaikov`aK + CosmologicalBaikov`bK -
            CosmologicalBaikov`cK) CosmologicalBaikov`bK;
      numeratorMapCheck = Expand[
        (p3 - p2) p3 +
          (2 CosmologicalBaikov`bK - CosmologicalBaikov`cK)/2
        - ((p3^2 + CosmologicalBaikov`x3sq)
          - (p2 p3 + CosmologicalBaikov`cK/2)
          + CosmologicalBaikov`bK - CosmologicalBaikov`x3sq)
      ];
      {
        CosmologicalBaikov111`I101ParameterRules[],
        MatrixForm[CosmologicalBaikov111`I101BasisTransform[]],
        Expand[originalLambdaInQ - mappedLambda],
        numeratorMapCheck
      }
    ],
    inputCell[
      connection101 =
        CosmologicalBaikov111`DifferentialMatrixA101[];
      overlapChecks = Map[
        Together,
        CosmologicalBaikov111`OverlapConsistencyChecks[],
        {2}
      ];
      VerificationTest[
        overlapChecks,
        ConstantArray[0, {5, 11}],
        TestID -> "all-pinch-overlaps-agree"
      ]
    ],

    Cell["4. Assemble the complete connection", "Section"],
    Cell[
      "The top-sector identity is the equation immediately following lambda0 (numbered (6) in the supplied PDF). The first displayed form of tilde Lambda0 in the PDF has a c^2/2 typo; the second line, tilde Lambda0=lambda0+b x2^2+a x3^2+c(x0^2-x2^2-x3^2)/2, is algebraically consistent and is used here.",
      "Text"
    ],
    inputCell[
      fullConnection =
        CosmologicalBaikov111`DifferentialMatrixAFull[];
      structuralPattern = Map[
        If[TrueQ[# === 0], 0, 1] &,
        fullConnection,
        {2}
      ];
      {
        Dimensions[fullConnection],
        MatrixForm[structuralPattern]
      }
    ],
    Cell[
      "With e110 selecting the scalar I110 master, the last row is [e110/2+(a-c/2) A110+(b-c/2) A101+(c/2) A011]/tildeLambda0 plus (b-x3^2)/(2 tildeLambda0) on the I111 diagonal. Thus no new IBP master is introduced by partial_a I111.",
      "Text"
    ],
    inputCell[
      reconstructedTopRow = (
        1/2 UnitVector[11, 5]
        + (CosmologicalBaikov`aK - CosmologicalBaikov`cK/2)
            fullConnection[[5]]
        + (CosmologicalBaikov`bK - CosmologicalBaikov`cK/2)
            fullConnection[[7]]
        + CosmologicalBaikov`cK/2 fullConnection[[9]]
      )/CosmologicalBaikov`TildeLambda0[];
      reconstructedTopRow[[11]] =
        (CosmologicalBaikov`bK - CosmologicalBaikov`x3sq)/
          (2 CosmologicalBaikov`TildeLambda0[]);
      regressionRules = {
        CosmologicalBaikov`aK -> 5,
        CosmologicalBaikov`bK -> -1,
        CosmologicalBaikov`cK -> 1,
        CosmologicalBaikov`x0sq -> 4,
        CosmologicalBaikov`x2sq -> 1,
        CosmologicalBaikov`x3sq -> 2,
        CosmologicalBaikov`eps -> 1/7
      };
      VerificationTest[
        Together[
          (fullConnection[[11]] - reconstructedTopRow) /.
            regressionRules
        ],
        ConstantArray[0, 11],
        TestID -> "top-row-reproduces-the-algebraic-I111-identity"
      ]
    ],
    inputCell[
      fullConnectionAtRegressionPoint = Map[
        Together,
        fullConnection /. regressionRules,
        {2}
      ];
      MatrixForm[fullConnectionAtRegressionPoint]
    ],

    Cell["5. Compact epsilon-zero connection", "Section"],
    Cell[
      "The Simplify/Factor observation can be implemented reproducibly by first setting epsilon=0, replacing the exact tildeLambda0 polynomial by an inert atomic abbreviation, and then factoring numerator and denominator. CompactifyConnectionEpsilonZero performs the full construction; it can take several minutes on the inherited expanded I110 expressions. CompactConnectionEpsilonZero loads the precomputed exact result immediately. The abbreviation can always be expanded back before differentiation or numerical evolution.",
      "Text"
    ],
    inputCell[
      fullConnectionS =
        CosmologicalBaikov111`CompactConnectionEpsilonZero[];
      abbreviationRules = {
        CosmologicalBaikov111`TildeLambda0AbbreviationRule[],
        CosmologicalBaikov111`TildeLambda0ExpansionRule[]
      };
      recomputationCommand = HoldForm[
        CosmologicalBaikov111`CompactifyConnectionEpsilonZero[
          fullConnection
        ]
      ];
      compactConnectionStatistics = <|
        "OriginalLeafCount" ->
          LeafCount[fullConnection /. CosmologicalBaikov`eps -> 0],
        "CompactLeafCount" -> LeafCount[fullConnectionS]
      |>;
      {
        abbreviationRules,
        recomputationCommand,
        compactConnectionStatistics,
        MatrixForm[fullConnectionS]
      }
    ],
    inputCell[
      expandedFullConnectionS =
        CosmologicalBaikov111`ExpandTildeLambda0Abbreviation[
          fullConnectionS
        ];
      compactEqualityCheck = DeleteDuplicates[
        Together /@ Flatten[
          (expandedFullConnectionS -
            (fullConnection /. CosmologicalBaikov`eps -> 0)) /.
              DeleteCases[
                regressionRules,
                HoldPattern[CosmologicalBaikov`eps -> _]
              ]
        ]
      ];
      VerificationTest[
        compactEqualityCheck,
        {0},
        TestID -> "compact-full-connection-equals-original"
      ]
    ],

    Cell["6. Leading singularities of all three pinch sectors", "Section"],
    Cell[
      "Each two-denominator sector has two maximal-cut points in the original D4 variable. Since its second master has one D4 numerator, the branch combinations are I[-D4]-zPlus I and I[-D4]-zMinus I.",
      "Text"
    ],
    inputCell[
      cutLocations = <|
        "110" -> CosmologicalBaikov`LeadingSingularityLocations[],
        "101" -> CosmologicalBaikov111`LeadingSingularityLocations101[],
        "011" -> CosmologicalBaikov011`LeadingSingularityLocations011[]
      |>;
      cutBases = <|
        "110" -> CosmologicalBaikov`LeadingSingularityBasisMatrix[],
        "101" -> CosmologicalBaikov111`LeadingSingularityBasisMatrix101[],
        "011" -> CosmologicalBaikov011`LeadingSingularityBasisMatrix011[]
      |>;
      branchMasters = <|
        "110" -> cutBases["110"].fullMasters[[5 ;; 6]],
        "101" -> cutBases["101"].fullMasters[[7 ;; 8]],
        "011" -> cutBases["011"].fullMasters[[9 ;; 10]]
      |>;
      {
        cutLocations,
        branchMasters,
        Determinant /@ cutBases
      }
    ],

    Cell["7. Algebraic lower-triangular form at epsilon=0", "Section"],
    inputCell[
      epsilonZeroRules = {
        CosmologicalBaikov`aK -> 5,
        CosmologicalBaikov`bK -> -1,
        CosmologicalBaikov`cK -> 1,
        CosmologicalBaikov`x0sq -> 4,
        CosmologicalBaikov`x2sq -> 1,
        CosmologicalBaikov`x3sq -> 2
      };
      epsilonZeroTriangularConnection =
        CosmologicalBaikov111`TriangularConnectionEpsilonZero[
          epsilonZeroRules
        ];
      epsilonZeroUpperEntries = DeleteDuplicates[
        Together /@
          CosmologicalBaikov111`StrictUpperEntries[
            epsilonZeroTriangularConnection
          ]
      ];
      {
        MatrixForm[epsilonZeroTriangularConnection],
        VerificationTest[
          epsilonZeroUpperEntries,
          {0},
          TestID -> "complete-epsilon-zero-connection-is-lower-triangular"
        ]
      }
    ],
    Cell[
      "The same rational compaction is effective after the leading-singularity change of basis. The cached symbolic result below is exactly the similarity transform of fullConnectionS and is manifestly lower triangular. Its leaf count is roughly half that of the unsimplified transform of the already compact matrix.",
      "Text"
    ],
    inputCell[
      leadingSingularityTransform =
        CosmologicalBaikov111`LeadingSingularityTransformFull[];
      rawCompactTriangularConnection =
        leadingSingularityTransform . fullConnectionS .
          Inverse[leadingSingularityTransform];
      triangularConnectionS =
        CosmologicalBaikov111`CompactTriangularConnectionEpsilonZero[];
      compactTriangularStatistics = <|
        "BeforeFinalSimplification" ->
          LeafCount[rawCompactTriangularConnection],
        "CompactTriangularLeafCount" -> LeafCount[triangularConnectionS]
      |>;
      compactTriangularUpperEntries = DeleteDuplicates[
        Together /@
          CosmologicalBaikov111`StrictUpperEntries[
            triangularConnectionS
          ]
      ];
      {
        compactTriangularStatistics,
        MatrixForm[triangularConnectionS],
        VerificationTest[
          compactTriangularUpperEntries,
          {0},
          TestID -> "compact-symbolic-connection-is-lower-triangular"
        ]
      }
    ],
    Cell[
      "The diagonal two-master entries have a compact cut interpretation. The first row of each cut basis is supported on zMinus and the second on zPlus, hence the reversal below. Each eigenvalue is -1/2 partial_a log LambdaCut.",
      "Text"
    ],
    inputCell[
      Clear[aa, z];
      cutLambda110[z_, t_] :=
        CosmologicalBaikov`cK (z - CosmologicalBaikov`cK/2)
        + CosmologicalBaikov`bK CosmologicalBaikov`x2sq
        - t (2 z - CosmologicalBaikov`x0sq +
            CosmologicalBaikov`x2sq - CosmologicalBaikov`cK)
        + CosmologicalBaikov`cK^2/4
        - t CosmologicalBaikov`bK;
      cutLambda101[z_, t_] :=
        CosmologicalBaikov`cK (z - CosmologicalBaikov`cK/2)
        - CosmologicalBaikov`bK (2 z - CosmologicalBaikov`x0sq +
            CosmologicalBaikov`x3sq - CosmologicalBaikov`cK)
        + t CosmologicalBaikov`x3sq
        + CosmologicalBaikov`cK^2/4
        - t CosmologicalBaikov`bK;
      cutLambda011[z_, t_] :=
        CosmologicalBaikov011`CutLambda011[z, t];
      cutEigenvalues = <|
        "110" -> (-1/2 D[Log[cutLambda110[#, aa]], aa] & /@
          Reverse[cutLocations["110"]]),
        "101" -> (-1/2 D[Log[cutLambda101[#, aa]], aa] & /@
          Reverse[cutLocations["101"]]),
        "011" -> (-1/2 D[Log[cutLambda011[#, aa]], aa] & /@
          Reverse[cutLocations["011"]])
      |>;
      Map[Simplify, cutEigenvalues, {2}]
    ],

    Cell["8. Exact finite-epsilon triangularization", "Section"],
    Cell[
      "At nonzero epsilon the algebraic cut bases leave O(epsilon) mixing. For each pinch pair {M,K}, define N=K+rho(a) M and order the pair as {N,M}. If its homogeneous block is B, exact triangularity requires rho'=B12 rho^2+(B22-B11) rho-B21. The three sectors are disjoint, so three independent Riccati transformations make the entire eleven-master connection lower triangular, including the final I111 row.",
      "Text"
    ],
    inputCell[
      finiteEpsilonRules = {
        CosmologicalBaikov`bK -> -1,
        CosmologicalBaikov`cK -> 1,
        CosmologicalBaikov`x0sq -> 4,
        CosmologicalBaikov`x2sq -> 1,
        CosmologicalBaikov`x3sq -> 2,
        CosmologicalBaikov`eps -> 1/7
      };
      Clear[rho110, rho101, rho011];
      riccatiEquations =
        CosmologicalBaikov111`ExactTriangularRiccatiEquationsFull[
          {rho110, rho101, rho011},
          CosmologicalBaikov`aK,
          finiteEpsilonRules
        ];
      riccatiInitialBranches = {
        rho110[5] -> -cutLocations["110"],
        rho101[5] -> -cutLocations["101"],
        rho011[5] -> -cutLocations["011"]
      } /. finiteEpsilonRules;
      {riccatiEquations, riccatiInitialBranches}
    ],
    inputCell[
      exactFiniteEpsilonTriangularConnection =
        CosmologicalBaikov111`ExactTriangularConnectionFull[
          {rho110, rho101, rho011},
          CosmologicalBaikov`aK,
          "KinematicRules" -> finiteEpsilonRules
        ];
      finiteEpsilonUpperEntries = DeleteDuplicates[
        Together /@
          CosmologicalBaikov111`StrictUpperEntries[
            exactFiniteEpsilonTriangularConnection
          ]
      ];
      {
        MatrixForm@Map[
          If[TrueQ[Together[#] === 0], 0, 1] &,
          exactFiniteEpsilonTriangularConnection,
          {2}
        ],
        VerificationTest[
          finiteEpsilonUpperEntries,
          {0},
          TestID -> "three-Riccati-bases-give-exact-finite-epsilon-triangularity"
        ]
      }
    ],

    Cell["9. Genericity and use", "Section"],
    inputCell[
      CosmologicalBaikov111`GenericityConditionsFull[]
    ],
    Cell[
      "DifferentialMatrixAFull[] is the exact symbolic connection in the original basis. At a chosen exact kinematic point, TriangularConnectionEpsilonZero[rules] gives the algebraic epsilon=0 triangular matrix. At finite epsilon, keep a free and pass the remaining fixed parameters through the KinematicRules option of ExactTriangularConnectionFull. The leading-singularity branches supply rho(a,0)=-zPlus or -zMinus. Degenerate cuts, vanishing tildeLambda0, or zeros of the denominators listed above require a separate local basis.",
      "Text"
    ]
  },
  WindowTitle -> "Cosmological Baikov: complete partial_a I111 system",
  StyleDefinitions -> "Default.nb"
];

target = FileNameJoin[{
  DirectoryName[$InputFileName],
  "CosmologicalBaikov111DaComplete.nb"
}];
Put[notebook, target];
Print[target];
