packageDirectory = DirectoryName[$InputFileName];
Get[FileNameJoin[{packageDirectory, "CosmologicalBaikov111.wl"}]];

VerificationTest[
  CosmologicalBaikov111`FullMasterBasis[],
  {
    CosmologicalBaikov111`FBI[0, 0, 0, 0],
    CosmologicalBaikov111`FBI[0, 1, 0, 0],
    CosmologicalBaikov111`FBI[1, 0, 0, 0],
    CosmologicalBaikov111`FBI[0, 0, 1, 0],
    CosmologicalBaikov111`FBI[1, 1, 0, 0],
    CosmologicalBaikov111`FBI[1, 1, 0, -1],
    CosmologicalBaikov111`FBI[1, 0, 1, 0],
    CosmologicalBaikov111`FBI[1, 0, 1, -1],
    CosmologicalBaikov111`FBI[0, 1, 1, 0],
    CosmologicalBaikov111`FBI[0, 1, 1, -1],
    CosmologicalBaikov111`FBI[1, 1, 1, 0]
  },
  TestID -> "full-eleven-master-basis"
]

VerificationTest[
  Module[{q2, q3, original, transformed},
    original = CosmologicalBaikov`cK (q3 - q2) q3
      - CosmologicalBaikov`bK (q3 - q2)^2
      - CosmologicalBaikov`aK q3^2
      + CosmologicalBaikov`cK^2/4
      - CosmologicalBaikov`aK CosmologicalBaikov`bK;
    transformed = (2 CosmologicalBaikov`bK - CosmologicalBaikov`cK) q2 q3
      - CosmologicalBaikov`bK q2^2
      - (CosmologicalBaikov`aK + CosmologicalBaikov`bK -
          CosmologicalBaikov`cK) q3^2
      + (2 CosmologicalBaikov`bK - CosmologicalBaikov`cK)^2/4
      - (CosmologicalBaikov`aK + CosmologicalBaikov`bK -
          CosmologicalBaikov`cK) CosmologicalBaikov`bK;
    Expand[original - transformed]
  ],
  0,
  TestID -> "I101-momentum-map-preserves-Lambda"
]

VerificationTest[
  Module[{p2, p3, d4Prime, rhs},
    d4Prime = (p3 - p2) p3 +
      (2 CosmologicalBaikov`bK - CosmologicalBaikov`cK)/2;
    rhs = (p3^2 + CosmologicalBaikov`x3sq)
      - (p2 p3 + CosmologicalBaikov`cK/2)
      + CosmologicalBaikov`bK - CosmologicalBaikov`x3sq;
    Expand[d4Prime - rhs]
  ],
  0,
  TestID -> "I101-auxiliary-numerator-map"
]

VerificationTest[
  Map[Together, CosmologicalBaikov111`OverlapConsistencyChecks[], {2}],
  ConstantArray[0, {5, 11}],
  TestID -> "pinch-connections-agree-on-overlaps"
]

VerificationTest[
  Dimensions[CosmologicalBaikov111`DifferentialMatrixAFull[]],
  {11, 11},
  TestID -> "full-connection-dimensions"
]

VerificationTest[
  FreeQ[
    CosmologicalBaikov111`DifferentialMatrixAFull[],
    _CosmologicalBaikov111`FBI
  ],
  True,
  TestID -> "full-connection-has-no-unreduced-integrals"
]

VerificationTest[
  Simplify[
    CosmologicalBaikov111`CutGramPolynomial101[#] & /@
      CosmologicalBaikov111`LeadingSingularityLocations101[]
  ],
  {0, 0},
  TestID -> "I101-leading-singularity-locations"
]

VerificationTest[
  Module[{connection, expected, rules},
    connection = CosmologicalBaikov111`DifferentialMatrixAFull[];
    expected = (
      1/2 UnitVector[11, 5]
      + (CosmologicalBaikov`aK - CosmologicalBaikov`cK/2)
          connection[[5]]
      + (CosmologicalBaikov`bK - CosmologicalBaikov`cK/2)
          connection[[7]]
      + CosmologicalBaikov`cK/2 connection[[9]]
    )/CosmologicalBaikov`TildeLambda0[];
    expected[[11]] =
      (CosmologicalBaikov`bK - CosmologicalBaikov`x3sq)/
        (2 CosmologicalBaikov`TildeLambda0[]);
    rules = {
      CosmologicalBaikov`aK -> 5,
      CosmologicalBaikov`bK -> -1,
      CosmologicalBaikov`cK -> 1,
      CosmologicalBaikov`x0sq -> 4,
      CosmologicalBaikov`x2sq -> 1,
      CosmologicalBaikov`x3sq -> 2,
      CosmologicalBaikov`eps -> 1/7
    };
    Together[(connection[[11]] - expected) /. rules]
  ],
  ConstantArray[0, 11],
  TestID -> "top-sector-row-is-equation-6-reduced"
]

VerificationTest[
  Module[{rules, triangular},
    rules = {
      CosmologicalBaikov`aK -> 5,
      CosmologicalBaikov`bK -> -1,
      CosmologicalBaikov`cK -> 1,
      CosmologicalBaikov`x0sq -> 4,
      CosmologicalBaikov`x2sq -> 1,
      CosmologicalBaikov`x3sq -> 2
    };
    triangular =
      CosmologicalBaikov111`TriangularConnectionEpsilonZero[rules];
    DeleteDuplicates[
      Together /@ CosmologicalBaikov111`StrictUpperEntries[triangular]
    ]
  ],
  {0},
  TestID -> "leading-singularity-basis-triangularizes-full-epsilon-zero-system"
]

VerificationTest[
  Module[{rules, triangular},
    rules = {
      CosmologicalBaikov`bK -> -1,
      CosmologicalBaikov`cK -> 1,
      CosmologicalBaikov`x0sq -> 4,
      CosmologicalBaikov`x2sq -> 1,
      CosmologicalBaikov`x3sq -> 2,
      CosmologicalBaikov`eps -> 1/7
    };
    triangular = CosmologicalBaikov111`ExactTriangularConnectionFull[
      {rho110, rho101, rho011},
      CosmologicalBaikov`aK,
      "KinematicRules" -> rules
    ];
    DeleteDuplicates[
      Together /@ CosmologicalBaikov111`StrictUpperEntries[triangular]
    ]
  ],
  {0},
  TestID -> "three-Riccati-transformations-triangularize-finite-epsilon-system"
]

VerificationTest[
  Module[{rules, compact, original},
    rules = {
      CosmologicalBaikov`aK -> 5,
      CosmologicalBaikov`bK -> -1,
      CosmologicalBaikov`cK -> 1,
      CosmologicalBaikov`x0sq -> 4,
      CosmologicalBaikov`x2sq -> 1,
      CosmologicalBaikov`x3sq -> 2
    };
    compact = CosmologicalBaikov111`ExpandTildeLambda0Abbreviation[
      CosmologicalBaikov111`CompactConnectionEpsilonZero[]
    ];
    original = CosmologicalBaikov111`DifferentialMatrixAFull[] /.
      CosmologicalBaikov`eps -> 0;
    DeleteDuplicates[Together /@ Flatten[(compact - original) /. rules]]
  ],
  {0},
  TestID -> "compact-epsilon-zero-connection-equals-original"
]

VerificationTest[
  Module[{rules, compact, triangular, transform},
    rules = {
      CosmologicalBaikov`aK -> 5,
      CosmologicalBaikov`bK -> -1,
      CosmologicalBaikov`cK -> 1,
      CosmologicalBaikov`x0sq -> 4,
      CosmologicalBaikov`x2sq -> 1,
      CosmologicalBaikov`x3sq -> 2
    };
    compact = CosmologicalBaikov111`CompactConnectionEpsilonZero[];
    triangular =
      CosmologicalBaikov111`CompactTriangularConnectionEpsilonZero[];
    transform =
      CosmologicalBaikov111`LeadingSingularityTransformFull[];
    DeleteDuplicates[
      Together /@ Flatten[
        (triangular - transform.compact.Inverse[transform]) /. rules
      ]
    ]
  ],
  {0},
  TestID -> "compact-triangular-connection-is-the-similarity-transform"
]

VerificationTest[
  DeleteDuplicates[
    Together /@
      CosmologicalBaikov111`StrictUpperEntries[
        CosmologicalBaikov111`CompactTriangularConnectionEpsilonZero[]
      ]
  ],
  {0},
  TestID -> "compact-symbolic-epsilon-zero-connection-is-lower-triangular"
]
