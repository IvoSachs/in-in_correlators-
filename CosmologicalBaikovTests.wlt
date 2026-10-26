If[!MemberQ[$Packages, "CosmologicalBaikov`"],
  candidateDirectories = DeleteDuplicates@Select[
    {
      Quiet@Check[DirectoryName[$InputFileName], Nothing],
      Directory[],
      FileNameJoin[{Directory[], "outputs"}]
    },
    StringQ
  ];
  packageFile = SelectFirst[
    FileNameJoin[{#, "CosmologicalBaikov.wl"}] & /@
      candidateDirectories,
    FileExistsQ,
    $Failed
  ];
  Get[packageFile];
];

VerificationTest[
  CosmologicalBaikov`TangentialCheck[],
  {0, 0},
  TestID -> "equation-12-is-tangent-and-divergence-free"
]

VerificationTest[
  Module[{p, q, lhs, rhs},
    lhs = -(p - q)^2 CosmologicalBaikov`cK/2
      - p^2 (CosmologicalBaikov`bK - CosmologicalBaikov`cK/2)
      - q^2 (CosmologicalBaikov`aK - CosmologicalBaikov`cK/2)
      + CosmologicalBaikov`cK^2/4
      - CosmologicalBaikov`aK CosmologicalBaikov`bK;
    rhs = CosmologicalBaikov`cK p q
      - CosmologicalBaikov`bK p^2
      - CosmologicalBaikov`aK q^2
      + CosmologicalBaikov`Lambda0[];
    Expand[lhs - rhs]
  ],
  0,
  TestID -> "equations-3-and-4-compact-form"
]

VerificationTest[
  Module[{p, q, d1, d2, d3, lambda},
    d1 = (p - q)^2 + CosmologicalBaikov`x0sq;
    d2 = p^2 + CosmologicalBaikov`x2sq;
    d3 = q^2 + CosmologicalBaikov`x3sq;
    lambda = CosmologicalBaikov`cK p q
      - CosmologicalBaikov`bK p^2
      - CosmologicalBaikov`aK q^2
      + CosmologicalBaikov`Lambda0[];
    Expand[
      lambda - (-CosmologicalBaikov`cK d1/2
        + (CosmologicalBaikov`cK/2 - CosmologicalBaikov`bK) d2
        + (CosmologicalBaikov`cK/2 - CosmologicalBaikov`aK) d3
        + CosmologicalBaikov`TildeLambda0[])
    ]
  ],
  0,
  TestID -> "corrected-equation-7-decomposition"
]

VerificationTest[
  CosmologicalBaikov`MasterBasis[],
  {
    CosmologicalBaikov`CBI[0, 0, 0],
    CosmologicalBaikov`CBI[0, 1, 0],
    CosmologicalBaikov`CBI[1, 0, 0],
    CosmologicalBaikov`CBI[1, 1, 0],
    CosmologicalBaikov`CBI[1, 1, -1]
  },
  TestID -> "equation-11-basis-order"
]

VerificationTest[
  CosmologicalBaikov`NumericalCrossCheck[],
  {
    {-12/161, 0, 0, 0, 0},
    {3/140, -493/3220, 0, 0, 0},
    {3/56, 0, 19/1288, 0, 0},
    {52/6445, 7288/3112935, -4940/207529, 1884/9023, 1100/27069},
    {155/72184, -34022/622587, 301745/4980696, 5975/27069, -4348/27069}
  },
  TestID -> "five-by-five-matrix-generic-rational-point"
]

VerificationTest[
  Together[
    CosmologicalBaikov`ReductionData[]["FourthIBPResidual"] /. {
      CosmologicalBaikov`aK -> 2,
      CosmologicalBaikov`bK -> 3,
      CosmologicalBaikov`cK -> 1,
      CosmologicalBaikov`x0sq -> 5,
      CosmologicalBaikov`x2sq -> 7,
      CosmologicalBaikov`eps -> 1/7
    }
  ],
  0,
  TestID -> "unused-fourth-v-minus-one-ibp-is-consistent"
]

VerificationTest[
  FreeQ[CosmologicalBaikov`DifferentialMatrixA[], _CosmologicalBaikov`CBI],
  True,
  TestID -> "matrix-contains-no-unreduced-integrals"
]

VerificationTest[
  Simplify[
    CosmologicalBaikov`CutGramPolynomial[#] & /@
      CosmologicalBaikov`LeadingSingularityLocations[]
  ],
  {0, 0},
  TestID -> "leading-singularity-points-are-cut-gram-roots"
]

VerificationTest[
  Module[{block, p, rules},
    block = CosmologicalBaikov`DifferentialMatrixA[][[4 ;; 5, 4 ;; 5]];
    p = CosmologicalBaikov`LeadingSingularityBasisMatrix[];
    rules = {
      CosmologicalBaikov`aK -> 2, CosmologicalBaikov`bK -> 3,
      CosmologicalBaikov`cK -> 1, CosmologicalBaikov`x0sq -> 5,
      CosmologicalBaikov`x2sq -> 7, CosmologicalBaikov`eps -> 0
    };
    With[{transformed = Map[Together, p.block.Inverse[p] /. rules, {2}]},
      {transformed[[1, 2]], transformed[[2, 1]]}
    ]
  ],
  {0, 0},
  TestID -> "leading-singularity-basis-diagonalizes-epsilon-zero-top-block"
]

VerificationTest[
  Module[{b11, b12, b21, b22, block, transform, equation, transformed},
    block = {{b11, b12}, {b21, b22}};
    transform = CosmologicalBaikov`TriangularizingBasisMatrix[rho, aa];
    equation = CosmologicalBaikov`TriangularizingRiccatiEquation[block, rho, aa];
    transformed = CosmologicalBaikov`TransformConnection[block, transform, aa];
    Simplify[transformed[[1, 2]] /. First@Solve[equation, rho'[aa]]]
  ],
  0,
  TestID -> "riccati-equation-makes-top-block-lower-triangular"
]
