(* Build the executable notebook without requiring a running front end. *)

(* Load the package while constructing boxes so exported symbols are written
   with the correct context into the notebook. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "CosmologicalBaikov.wl"}]];

ClearAll[inputCell];
SetAttributes[inputCell, HoldAllComplete];
inputCell[expr_] := Cell[
  BoxData[ToBoxes[Defer[expr], StandardForm]],
  "Input"
];

notebook = Notebook[
  {
    Cell["Leading-singularity basis for the I110 differential system", "Title"],
    Cell[
      "This notebook constructs the two maximal-cut leading-singularity branches of the top sector {I^(110,000;0), I^(110,-100;0)}, changes to branch-adapted linear combinations, and checks the differential equation. It also explains and implements the extra step needed for exact triangularity at finite epsilon.",
      "Text"
    ],

    Cell["1. Load the reduction and extract the top-sector block", "Section"],
    Cell[
      "The package orders the masters as {I^(000,000;0), I^(010,000;0), I^(100,000;0), I^(110,000;0), I^(110,-100;0)}. The last two rows and columns form the homogeneous top-sector block B.",
      "Text"
    ],
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
      topBlock = connection[[4 ;; 5, 4 ;; 5]];
      {masters, Dimensions[connection]}
    ],

    Cell["2. Maximal cut and leading singularities", "Section"],
    Cell[
      "Let z=D4. In terms of Baikov variables, p2^2=D2-x2^2, p2.p3=D4-c/2, and p3^2=D1-D2+2 D4-x0^2+x2^2-c. On D1=D2=0 their Gram determinant has two zeros zPlus and zMinus. These are the two one-dimensional leading-singularity branches.",
      "Text"
    ],
    inputCell[
      Clear[d1, d2, d4];
      p2Squared = d2 - CosmologicalBaikov`x2sq;
      p2Dotp3 = d4 - CosmologicalBaikov`cK/2;
      p3Squared = d1 - d2 + 2 d4 - CosmologicalBaikov`x0sq +
        CosmologicalBaikov`x2sq - CosmologicalBaikov`cK;
      gramPolynomial = Expand[p2Squared p3Squared - p2Dotp3^2];
      cutGram = Factor[gramPolynomial /. {d1 -> 0, d2 -> 0}];
      cutLocations = CosmologicalBaikov`LeadingSingularityLocations[];
      {
        cutGram,
        cutLocations,
        Simplify[CosmologicalBaikov`CutGramPolynomial[#] & /@ cutLocations]
      }
    ],
    Cell[
      "M5=I^(110,-100;0) has one factor of D4 in its numerator relative to M4=I^(110,000;0). Therefore LS_i(M5)/LS_i(M4)=z_i. The combinations M5-zPlus M4 and M5-zMinus M4 each kill one cut branch.",
      "Text"
    ],
    inputCell[
      cutBasis = CosmologicalBaikov`LeadingSingularityBasisMatrix[];
      newTopMasters = cutBasis . masters[[4 ;; 5]];
      leadingSingularityRows = {
        {1, cutLocations[[1]]},
        {1, cutLocations[[2]]}
      };
      branchSupport = Simplify[cutBasis . Transpose[leadingSingularityRows]];
      {
        newTopMasters,
        branchSupport,
        Factor[Det[cutBasis]]
      }
    ],

    Cell["3. Transform the differential equation", "Section"],
    Cell[
      "For n=T m and partial_a m=A m, the transformed connection is Anew=(partial_a T) T^(-1)+T A T^(-1). The derivative term must be retained for any a-dependent basis. Here the cut basis is independent of a, so this term vanishes.",
      "Text"
    ],
    inputCell[
      zero32 = ConstantArray[0, {3, 2}];
      zero23 = ConstantArray[0, {2, 3}];
      fullCutBasis = ArrayFlatten[
        {{IdentityMatrix[3], zero32}, {zero23, cutBasis}}
      ];
      transformedConnection = CosmologicalBaikov`TransformConnection[
        connection, fullCutBasis, CosmologicalBaikov`aK
      ];
      transformedTopBlock = transformedConnection[[4 ;; 5, 4 ;; 5]];
    ],
    Cell[
      "At epsilon=0 the two branch-adapted masters decouple, so the 2 by 2 block is diagonal (and the full five-master matrix is block-lower-triangular). At finite epsilon, the remaining off-diagonal entries are proportional to epsilon; leading singularities alone do not make the exact finite-epsilon connection triangular.",
      "Text"
    ],
    inputCell[
      sampleKinematics = {
        CosmologicalBaikov`aK -> 2, CosmologicalBaikov`bK -> 3,
        CosmologicalBaikov`cK -> 1, CosmologicalBaikov`x0sq -> 5,
        CosmologicalBaikov`x2sq -> 7
      };
      sampleCutBlock = Map[
        Factor,
        Map[Together, transformedTopBlock /. sampleKinematics, {2}],
        {2}
      ];
      sampleOffDiagonal = {
        sampleCutBlock[[1, 2]], sampleCutBlock[[2, 1]]
      };
      {
        MatrixForm[sampleCutBlock /. CosmologicalBaikov`eps -> 0],
        Together[sampleOffDiagonal /. CosmologicalBaikov`eps -> 0],
        Factor /@ Together[sampleOffDiagonal/CosmologicalBaikov`eps]
      }
    ],
    inputCell[
      VerificationTest[
        Together[sampleOffDiagonal /. CosmologicalBaikov`eps -> 0],
        {0, 0},
        TestID -> "cut-basis-diagonal-at-epsilon-zero"
      ]
    ],

    Cell["4. Exact finite-epsilon triangularization", "Section"],
    Cell[
      "To triangularize exactly, define N=M5+rho(a) M4 and order the pair as {N,M4}. Vanishing of the M4 coefficient in partial_a N gives a Riccati equation. Its two epsilon=0 branches are selected by rho(a,0)=-zPlus or -zMinus. This is the missing gauge-transformation step if equation (22) is meant at finite epsilon.",
      "Text"
    ],
    inputCell[
      Clear[rho, aa, b11, b12, b21, b22];
      symbolicBlock = {{b11, b12}, {b21, b22}};
      rhoBasis = CosmologicalBaikov`TriangularizingBasisMatrix[rho, aa];
      riccatiEquation = CosmologicalBaikov`TriangularizingRiccatiEquation[
        symbolicBlock, rho, aa
      ];
      rhoConnection = CosmologicalBaikov`TransformConnection[
        symbolicBlock, rhoBasis, aa
      ];
      triangularityCheck = Simplify[
        rhoConnection[[1, 2]] /.
          First@Solve[riccatiEquation, rho'[aa]]
      ];
      {riccatiEquation, MatrixForm[rhoConnection], triangularityCheck}
    ],
    inputCell[
      actualBlock = topBlock /. CosmologicalBaikov`aK -> aa;
      actualRiccatiEquation =
        CosmologicalBaikov`TriangularizingRiccatiEquation[
        actualBlock, rho, aa
      ];
      initialDataPlus = rho[2] ==
        (-cutLocations[[1]] /. CosmologicalBaikov`aK -> 2);
      {actualRiccatiEquation, initialDataPlus}
    ],

    Cell["5. Numerical finite-epsilon example", "Section"],
    Cell[
      "This optional cell solves the exact Riccati equation at a generic point with epsilon=1/7. The residual shown at the end is the upper-right entry of the transformed top block in the ordered basis {N,M4}; it should be zero to the requested working precision.",
      "Text"
    ],
    inputCell[
      numericalParameters = {
        CosmologicalBaikov`bK -> 3, CosmologicalBaikov`cK -> 1,
        CosmologicalBaikov`x0sq -> 5, CosmologicalBaikov`x2sq -> 7,
        CosmologicalBaikov`eps -> 1/7
      };
      numericalBlock = Map[
        Together,
        topBlock /. numericalParameters /. CosmologicalBaikov`aK -> aa,
        {2}
      ];
      rhoInitial = N[
        -cutLocations[[1]] /. numericalParameters,
        30
      ];
      rhoSolution = NDSolveValue[
        {
          rhoN'[aa] ==
            rhoN[aa]^2 numericalBlock[[1, 2]] +
            rhoN[aa] (numericalBlock[[2, 2]] - numericalBlock[[1, 1]]) -
            numericalBlock[[2, 1]],
          rhoN[2] == rhoInitial
        },
        rhoN,
        {aa, 3/2, 5/2},
        WorkingPrecision -> 30
      ];
      blockAtTwo = numericalBlock /. aa -> 2;
      triangularResidualAtTwo = N[
        rhoSolution'[2] + blockAtTwo[[2, 1]] +
          rhoSolution[2] (blockAtTwo[[1, 1]] - blockAtTwo[[2, 2]]) -
          rhoSolution[2]^2 blockAtTwo[[1, 2]],
        12
      ];
      {rhoInitial, rhoSolution[2], triangularResidualAtTwo}
    ],

    Cell["Interpretation", "Section"],
    Cell[
      "The leading-singularity basis is algebraic and is the natural epsilon=0 choice. It makes the top-sector block diagonal at epsilon=0 and leaves only O(epsilon) mixing. If exact triangularity for nonzero epsilon is required, rho(a) must satisfy the Riccati equation; simply diagonalizing B(a) pointwise is insufficient because it omits (partial_a T) T^(-1). Singular kinematic loci x0^2 x2^2=0 merge the two cut branches and require a different basis.",
      "Text"
    ]
  },
  WindowTitle -> "Cosmological Baikov: leading-singularity basis",
  StyleDefinitions -> "Default.nb"
];

target = FileNameJoin[{
  DirectoryName[$InputFileName],
  "CosmologicalBaikovLeadingSingularities.nb"
}];
Put[notebook, target];
Print[target];
