(* Build the analytic lower-sector integration notebook without a front end. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "CosmologicalBaikov.wl"}]];

ClearAll[inputCell];
SetAttributes[inputCell, HoldAllComplete];
inputCell[expr_] := Cell[
  BoxData[ToBoxes[Defer[expr], StandardForm]],
  "Input"
];

notebook = Notebook[
  {
    Cell["Integration of the lower 3 x 3 Cosmological-Baikov system", "Title"],
    Cell[
      "Here \"lower block\" means the lower-sector subsystem on the first three masters {I^(000,000;0), I^(010,000;0), I^(100,000;0)}, namely A[[1;;3,1;;3]]. The notebook constructs an analytic fundamental matrix and the boundary-value propagator U(a,a0), verifies them symbolically, and compares with NDSolve.",
      "Text"
    ],

    Cell["1. Load the reduction", "Section"],
    inputCell[
      ClearAll["Global`*"];
      notebookDirectory = Quiet@Check[
        NotebookDirectory[],
        DirectoryName[$InputFileName]
      ];
      Get[FileNameJoin[{notebookDirectory, "CosmologicalBaikov.wl"}]];
    ],
    inputCell[
      masters = CosmologicalBaikov`MasterBasis[];
      connection = CosmologicalBaikov`DifferentialMatrixA[];
      lowerMasters = masters[[1 ;; 3]];
      lowerBlock = Map[Together, connection[[1 ;; 3, 1 ;; 3]], {2}];
      {lowerMasters, MatrixForm[lowerBlock]}
    ],
    Cell[
      "The subsystem is triangular: M1 evolves independently, while M2 and M3 are driven only by M1. M2 and M3 do not mix with one another.",
      "Text"
    ],

    Cell["2. Common quadrature", "Section"],
    Cell[
      "Define lambda(a)=c^2/4-a b and h(a)=a+b-c. After dividing by the homogeneous solutions, both driven equations reduce to the integral of (4 epsilon-1)/2 u^(-1/2-2 epsilon) (u-X)^(epsilon-1). On the real branch u>X>0 it is an incomplete beta function.",
      "Text"
    ],
    inputCell[
      Clear[aa, a0, u, x];
      lam[t_] := CosmologicalBaikov`cK^2/4 -
        t CosmologicalBaikov`bK;
      h[t_] := t + CosmologicalBaikov`bK - CosmologicalBaikov`cK;
      betaPrimitive[u_, x_] :=
        (4 CosmologicalBaikov`eps - 1)/2 *
        x^(-1/2 - CosmologicalBaikov`eps) *
        Beta[
          1 - x/u,
          CosmologicalBaikov`eps,
          CosmologicalBaikov`eps + 1/2
        ];
      primitiveIntegrand[u_, x_] :=
        (4 CosmologicalBaikov`eps - 1)/2 *
        u^(-1/2 - 2 CosmologicalBaikov`eps) *
        (u - x)^(CosmologicalBaikov`eps - 1);
      primitiveCheck = FullSimplify[
        D[betaPrimitive[u, x], u] - primitiveIntegrand[u, x],
        Assumptions -> {
          u > x > 0,
          CosmologicalBaikov`eps > 0
        }
      ];
      primitiveCheck
    ],

    Cell["3. Fundamental matrix", "Section"],
    Cell[
      "The three diagonal homogeneous solutions are g1, g2 and g3 below. A convenient fundamental matrix F has its first column equal to the solution sourced by M1 and its remaining columns equal to the two homogeneous solutions.",
      "Text"
    ],
    inputCell[
      g1[t_] := lam[t]^(-CosmologicalBaikov`eps);
      g2[t_] :=
        lam[t]^(-CosmologicalBaikov`eps) *
        t^(2 CosmologicalBaikov`eps - 1/2) *
        (t - CosmologicalBaikov`x2sq)^(-CosmologicalBaikov`eps);
      g3[t_] :=
        lam[t]^(-CosmologicalBaikov`eps) *
        h[t]^(2 CosmologicalBaikov`eps - 1/2) *
        (h[t] - CosmologicalBaikov`x0sq)^(-CosmologicalBaikov`eps);
      fundamentalMatrix[t_] := {
        {g1[t], 0, 0},
        {
          g2[t] betaPrimitive[t, CosmologicalBaikov`x2sq],
          g2[t],
          0
        },
        {
          g3[t] betaPrimitive[h[t], CosmologicalBaikov`x0sq],
          0,
          g3[t]
        }
      };
      MatrixForm[fundamentalMatrix[aa]]
    ],
    inputCell[
      lowerBlockAA = lowerBlock /. CosmologicalBaikov`aK -> aa;
      branchAssumptions = {
        aa > CosmologicalBaikov`x2sq > 0,
        h[aa] > CosmologicalBaikov`x0sq > 0,
        lam[aa] > 0,
        CosmologicalBaikov`eps > 0
      };
      fundamentalResidual = FullSimplify[
        D[fundamentalMatrix[aa], aa] -
          lowerBlockAA . fundamentalMatrix[aa],
        Assumptions -> branchAssumptions,
        TransformationFunctions -> {Automatic, PowerExpand}
      ];
      VerificationTest[
        fundamentalResidual,
        ConstantArray[0, {3, 3}],
        TestID -> "analytic-fundamental-matrix"
      ]
    ],

    Cell["4. Boundary-value evolution matrix", "Section"],
    Cell[
      "For initial data m(a0)=m0, the solution is m(a)=U(a,a0)m0. Writing U directly in terms of primitive differences avoids irrelevant integration constants and is better behaved when epsilon approaches zero.",
      "Text"
    ],
    inputCell[
      lowerEvolution[t_, t0_] := {
        {g1[t]/g1[t0], 0, 0},
        {
          g2[t]/g1[t0] *
            (betaPrimitive[t, CosmologicalBaikov`x2sq] -
              betaPrimitive[t0, CosmologicalBaikov`x2sq]),
          g2[t]/g2[t0],
          0
        },
        {
          g3[t]/g1[t0] *
            (betaPrimitive[h[t], CosmologicalBaikov`x0sq] -
              betaPrimitive[h[t0], CosmologicalBaikov`x0sq]),
          0,
          g3[t]/g3[t0]
        }
      };
      initialValueCheck = Simplify[lowerEvolution[a0, a0]];
      {MatrixForm[lowerEvolution[aa, a0]], initialValueCheck}
    ],
    inputCell[
      evolutionAssumptions = Join[
        branchAssumptions,
        {
          a0 > CosmologicalBaikov`x2sq,
          h[a0] > CosmologicalBaikov`x0sq,
          lam[a0] > 0
        }
      ];
      evolutionResidual = FullSimplify[
        D[lowerEvolution[aa, a0], aa] -
          lowerBlockAA . lowerEvolution[aa, a0],
        Assumptions -> evolutionAssumptions,
        TransformationFunctions -> {Automatic, PowerExpand}
      ];
      VerificationTest[
        {
          initialValueCheck,
          evolutionResidual
        },
        {
          IdentityMatrix[3],
          ConstantArray[0, {3, 3}]
        },
        TestID -> "boundary-value-evolution-operator"
      ]
    ],
    inputCell[
      Clear[m10, m20, m30];
      initialMasters = {m10, m20, m30};
      lowerSolution[t_, t0_, initial_List] :=
        lowerEvolution[t, t0] . initial;
      lowerSolution[aa, a0, initialMasters]
    ],

    Cell["5. The epsilon=0 limit", "Section"],
    Cell[
      "At epsilon=0 the individual incomplete beta primitive contains an a-independent divergent constant, but it cancels from U(a,a0). A regular primitive for direct epsilon=0 work is ArcTanh[Sqrt[X/u]]/Sqrt[X]. Use its difference at the two endpoints in lowerEvolution.",
      "Text"
    ],
    inputCell[
      epsilonZeroPrimitive[u_, x_] :=
        ArcTanh[Sqrt[x/u]]/Sqrt[x];
      epsilonZeroCheck = FullSimplify[
        D[epsilonZeroPrimitive[u, x], u] -
          (primitiveIntegrand[u, x] /. CosmologicalBaikov`eps -> 0),
        Assumptions -> u > x > 0
      ];
      epsilonZeroCheck
    ],

    Cell["6. Numerical cross-check", "Section"],
    Cell[
      "The following Euclidean example compares the analytic propagator with a direct high-precision NDSolve integration. The final number is the norm of their difference.",
      "Text"
    ],
    inputCell[
      numericalParameters = {
        CosmologicalBaikov`bK -> -1,
        CosmologicalBaikov`cK -> 0,
        CosmologicalBaikov`x0sq -> 1,
        CosmologicalBaikov`x2sq -> 1,
        CosmologicalBaikov`eps -> 1/7
      };
      numericalInitial = {2, -1, 3};
      analyticEndpoint =
        N[lowerEvolution[4, 3] /. numericalParameters, 40] .
          numericalInitial;
      numericalBlock[t_?NumericQ] := N[
        lowerBlock /. numericalParameters /.
          CosmologicalBaikov`aK -> t,
        40
      ];
      numericalSolution = NDSolveValue[
        {
          vector'[t] == numericalBlock[t] . vector[t],
          vector[3] == numericalInitial
        },
        vector,
        {t, 3, 4},
        WorkingPrecision -> 40,
        AccuracyGoal -> 25,
        PrecisionGoal -> 25,
        MaxSteps -> Infinity
      ];
      numericalEndpoint = numericalSolution[4];
      {
        analyticEndpoint,
        numericalEndpoint,
        N[Norm[analyticEndpoint - numericalEndpoint], 10]
      }
    ],

    Cell["Branch and singular-locus remarks", "Section"],
    Cell[
      "The symbolic verification above chooses the real branch a>x2^2>0, h(a)>x0^2>0 and lambda(a)>0. Elsewhere the same formulas define the analytic continuation when all powers and incomplete beta functions are continued consistently along the integration path. Do not evolve through lambda=0, a=0, a=x2^2, h=0 or h=x0^2 without changing local basis or specifying a contour prescription.",
      "Text"
    ]
  },
  WindowTitle -> "Cosmological Baikov: lower 3 x 3 integration",
  StyleDefinitions -> "Default.nb"
];

target = FileNameJoin[{
  DirectoryName[$InputFileName],
  "CosmologicalBaikovLowerBlockIntegration.nb"
}];
Put[notebook, target];
Print[target];
