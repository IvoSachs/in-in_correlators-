(* Build the epsilon=0 integration notebook for the remaining two masters. *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "CosmologicalBaikov.wl"}]];

ClearAll[inputCell];
SetAttributes[inputCell, HoldAllComplete];
inputCell[expr_] := Cell[
  BoxData[ToBoxes[Defer[expr], StandardForm]],
  "Input"
];

notebook = Notebook[
  {
    Cell["Integration of the remaining two Cosmological-Baikov masters", "Title"],
    Cell[
      "This notebook integrates I^(110,000;0) and I^(110,-100;0) at epsilon=0. It first changes to the two leading-singularity combinations, where the homogeneous 2 x 2 block is diagonal, and then applies variation of constants using the already integrated lower 3 x 3 subsystem.",
      "Text"
    ],

    Cell["1. Load and transform the differential system", "Section"],
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
      cutLocations = CosmologicalBaikov`LeadingSingularityLocations[];
      cutBasis = CosmologicalBaikov`LeadingSingularityBasisMatrix[];
      fullCutBasis = ArrayFlatten[{
        {
          IdentityMatrix[3],
          ConstantArray[0, {3, 2}]
        },
        {
          ConstantArray[0, {2, 3}],
          cutBasis
        }
      }];
      epsilonZeroConnection = connection /. CosmologicalBaikov`eps -> 0;
      transformedEpsilonZeroConnection =
        CosmologicalBaikov`TransformConnection[
          epsilonZeroConnection,
          fullCutBasis,
          CosmologicalBaikov`aK
        ];
      lowerBlock0 = epsilonZeroConnection[[1 ;; 3, 1 ;; 3]];
      sourceBlock0 =
        transformedEpsilonZeroConnection[[4 ;; 5, 1 ;; 3]];
      transformedTopBlock0 =
        transformedEpsilonZeroConnection[[4 ;; 5, 4 ;; 5]];
      {
        cutBasis . masters[[4 ;; 5]],
        Dimensions[sourceBlock0]
      }
    ],
    Cell[
      "The transformed masters are NPlus=M5-zPlus M4 and NMinus=M5-zMinus M4. NPlus is supported on the zMinus cut and NMinus on the zPlus cut; this reversal is important in the homogeneous factors below.",
      "Text"
    ],

    Cell["2. Homogeneous top-sector solutions", "Section"],
    Cell[
      "On D1=D2=0, with z=D4, the Baikov polynomial appearing in the leading singularity is LambdaCut(z,a). At epsilon=0 each branch evolves as LambdaCut^(-1/2).",
      "Text"
    ],
    inputCell[
      Clear[aa, a0, ss];
      cutLambda[z_, t_] :=
        CosmologicalBaikov`cK (z - CosmologicalBaikov`cK/2) +
        CosmologicalBaikov`bK CosmologicalBaikov`x2sq -
        t (2 z - CosmologicalBaikov`x0sq +
          CosmologicalBaikov`x2sq - CosmologicalBaikov`cK) +
        CosmologicalBaikov`cK^2/4 -
        t CosmologicalBaikov`bK;
      supportedCutLocations = Reverse[cutLocations];
      topHomogeneous[t_] := DiagonalMatrix[
        cutLambda[#, t]^(-1/2) & /@ supportedCutLocations
      ];
      topHomogeneousRatio[t_, t0_] :=
        topHomogeneous[t] . Inverse[topHomogeneous[t0]];
      {
        supportedCutLocations,
        MatrixForm[topHomogeneous[aa]]
      }
    ],
    inputCell[
      homogeneousCheckPoint = {
        CosmologicalBaikov`aK -> 5,
        CosmologicalBaikov`bK -> -1,
        CosmologicalBaikov`cK -> 1,
        CosmologicalBaikov`x0sq -> 2,
        CosmologicalBaikov`x2sq -> 1
      };
      topBlockAtCheckPoint = Map[
        Together,
        transformedTopBlock0 /. homogeneousCheckPoint,
        {2}
      ];
      predictedTopBlockAtCheckPoint = Map[
        Together,
        D[topHomogeneous[aa], aa] .
          Inverse[topHomogeneous[aa]] /.
          homogeneousCheckPoint /.
          aa -> 5,
        {2}
      ];
      VerificationTest[
        FullSimplify[
          topBlockAtCheckPoint - predictedTopBlockAtCheckPoint
        ],
        ConstantArray[0, {2, 2}],
        TestID -> "cut-leading-singularities-give-top-homogeneous-factors"
      ]
    ],

    Cell["3. Epsilon-zero lower-sector propagator", "Section"],
    Cell[
      "At epsilon=0, M1 is constant. The other two lower masters use the regular primitive phi0(u,X)=ArcTanh[Sqrt[X/u]]/Sqrt[X].",
      "Text"
    ],
    inputCell[
      h[t_] := t + CosmologicalBaikov`bK - CosmologicalBaikov`cK;
      phi0[u_, x_] := ArcTanh[Sqrt[x/u]]/Sqrt[x];
      lowerEvolution0[t_, t0_] := {
        {1, 0, 0},
        {
          t^(-1/2) *
            (phi0[t, CosmologicalBaikov`x2sq] -
              phi0[t0, CosmologicalBaikov`x2sq]),
          Sqrt[t0/t],
          0
        },
        {
          h[t]^(-1/2) *
            (phi0[h[t], CosmologicalBaikov`x0sq] -
              phi0[h[t0], CosmologicalBaikov`x0sq]),
          0,
          Sqrt[h[t0]/h[t]]
        }
      };
      lowerEvolutionCheck = FullSimplify[
        D[lowerEvolution0[aa, a0], aa] -
          (lowerBlock0 /. CosmologicalBaikov`aK -> aa) .
            lowerEvolution0[aa, a0],
        Assumptions -> {
          aa > CosmologicalBaikov`x2sq > 0,
          a0 > CosmologicalBaikov`x2sq,
          h[aa] > CosmologicalBaikov`x0sq > 0,
          h[a0] > CosmologicalBaikov`x0sq
        },
        TransformationFunctions -> {Automatic, PowerExpand}
      ];
      VerificationTest[
        {
          Simplify[lowerEvolution0[a0, a0]],
          lowerEvolutionCheck
        },
        {
          IdentityMatrix[3],
          ConstantArray[0, {3, 3}]
        },
        TestID -> "epsilon-zero-lower-propagator"
      ]
    ],

    Cell["4. Exact variation-of-constants quadrature", "Section"],
    Cell[
      "Let ell0 be the three lower-master values and mTop0={M4(a0),M5(a0)}. In the leading-singularity basis n=P mTop, the exact epsilon=0 solution is n(a)=H(a)[H(a0)^(-1)P mTop0 + Integral(H(s)^(-1) S(s) ULower(s,a0) ell0,ds)]. The original two masters are recovered with P^(-1).",
      "Text"
    ],
    inputCell[
      topSourceKernel[s_, t0_] :=
        Inverse[topHomogeneous[s]] .
        (sourceBlock0 /. CosmologicalBaikov`aK -> s) .
        lowerEvolution0[s, t0];
      remainingTopFormal[t_, t0_, lowerInitial_List, topInitial_List] :=
        Inverse[cutBasis] . topHomogeneous[t] . (
          Inverse[topHomogeneous[t0]] . cutBasis . topInitial +
          Inactive[Integrate][
            topSourceKernel[ss, t0] . lowerInitial,
            {ss, t0, t}
          ]
        );
      HoldForm[
        {M4[aa], M5[aa]} ==
          remainingTopFormal[
            aa,
            a0,
            {M10, M20, M30},
            {M40, M50}
          ]
      ]
    ],

    Cell["5. Numerical evaluator for arbitrary generic kinematics", "Section"],
    Cell[
      "The function below evaluates only the explicit 2 x 3 quadrature. Before NIntegrate it simplifies the source matrix after substituting the requested kinematics, which is much faster than repeatedly evaluating the unsimplified IBP expressions.",
      "Text"
    ],
    inputCell[
      Clear[remainingMastersEpsilonZero];
      remainingMastersEpsilonZero[
        t1_?NumericQ,
        t0_?NumericQ,
        lowerInitial_List,
        topInitial_List,
        kinematicRules_List,
        workingPrecision_Integer : 40
      ] := Module[
        {
          integrationVariable,
          pMatrix,
          sourceExpression,
          kernelExpression,
          integralMatrix,
          nEndpoint,
          lowerEndpoint,
          topEndpoint
        },
        pMatrix = cutBasis /. kinematicRules;
        sourceExpression = Map[
          Together,
          sourceBlock0 /. kinematicRules /.
            CosmologicalBaikov`aK -> integrationVariable,
          {2}
        ];
        kernelExpression = Map[
          Together,
          (
            Inverse[topHomogeneous[integrationVariable]] .
            sourceExpression .
            lowerEvolution0[integrationVariable, t0]
          ) /. kinematicRules,
          {2}
        ];
        integralMatrix = Table[
          NIntegrate[
            kernelExpression[[i, j]],
            {integrationVariable, t0, t1},
            WorkingPrecision -> workingPrecision,
            AccuracyGoal -> Floor[workingPrecision/2] + 5,
            PrecisionGoal -> Floor[workingPrecision/2] + 5,
            MaxRecursion -> 20
          ],
          {i, 2},
          {j, 3}
        ];
        nEndpoint =
          (topHomogeneous[t1] /. kinematicRules) . (
            Inverse[topHomogeneous[t0] /. kinematicRules] .
              pMatrix . topInitial +
            integralMatrix . lowerInitial
          );
        topEndpoint = N[Inverse[pMatrix] . nEndpoint, workingPrecision];
        lowerEndpoint = N[
          (lowerEvolution0[t1, t0] /. kinematicRules) . lowerInitial,
          workingPrecision
        ];
        <|
          "LowerMasters" -> lowerEndpoint,
          "TopMasters" -> topEndpoint,
          "AllMasters" -> Join[lowerEndpoint, topEndpoint]
        |>
      ];
      "remainingMastersEpsilonZero is defined"
    ],

    Cell["6. Independent numerical verification", "Section"],
    inputCell[
      numericalKinematics = {
        CosmologicalBaikov`bK -> -1,
        CosmologicalBaikov`cK -> 1,
        CosmologicalBaikov`x0sq -> 2,
        CosmologicalBaikov`x2sq -> 1
      };
      lowerInitialValues = {2, -1, 3};
      topInitialValues = {1, -2};
      quadratureSolution = remainingMastersEpsilonZero[
        6,
        5,
        lowerInitialValues,
        topInitialValues,
        numericalKinematics,
        40
      ];
      quadratureSolution
    ],
    inputCell[
      fullInitialValues = Join[lowerInitialValues, topInitialValues];
      simplifiedNumericalConnection = Map[
        Together,
        epsilonZeroConnection /. numericalKinematics /.
          CosmologicalBaikov`aK -> tt,
        {2}
      ];
      directSolution = NDSolveValue[
        {
          masterVector'[tt] ==
            simplifiedNumericalConnection . masterVector[tt],
          masterVector[5] == fullInitialValues
        },
        masterVector,
        {tt, 5, 6},
        WorkingPrecision -> 40,
        AccuracyGoal -> 25,
        PrecisionGoal -> 25,
        MaxSteps -> Infinity
      ];
      directEndpoint = directSolution[6];
      comparisonError = N[
        Norm[quadratureSolution["AllMasters"] - directEndpoint],
        10
      ];
      {
        quadratureSolution["AllMasters"],
        directEndpoint,
        comparisonError
      }
    ],

    Cell["Existence and branch qualifications", "Section"],
    Cell[
      "The coefficient matrix and its evolution operator have a finite epsilon->0 limit at generic kinematics. This does not by itself prove that a particular physical integral has a finite limit: its boundary constants may still contain 1/epsilon poles. The formulas require nonzero cutLambda values and a nondegenerate cut basis x0^2 x2^2 != 0. Powers, square roots and ArcTanh must be continued consistently if the integration contour leaves the displayed Euclidean branch or crosses a singular locus.",
      "Text"
    ]
  },
  WindowTitle -> "Cosmological Baikov: remaining masters at epsilon=0",
  StyleDefinitions -> "Default.nb"
];

target = FileNameJoin[{
  DirectoryName[$InputFileName],
  "CosmologicalBaikovRemainingMastersEpsilonZero.nb"
}];
Put[notebook, target];
Print[target];
