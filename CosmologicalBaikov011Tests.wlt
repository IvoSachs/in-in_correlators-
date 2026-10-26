packageDirectory = DirectoryName[$InputFileName];
Get[FileNameJoin[{packageDirectory, "CosmologicalBaikov.wl"}]];
Get[FileNameJoin[{packageDirectory, "CosmologicalBaikov011.wl"}]];

VerificationTest[
  CosmologicalBaikov011`MasterBasis011[],
  {
    CosmologicalBaikov011`J011[0, 0, 0],
    CosmologicalBaikov011`J011[1, 0, 0],
    CosmologicalBaikov011`J011[0, 1, 0],
    CosmologicalBaikov011`J011[1, 1, 0],
    CosmologicalBaikov011`J011[1, 1, -1]
  },
  TestID -> "I011-five-master-basis"
]

VerificationTest[
  CosmologicalBaikov011`TangentialCheck011[],
  {0, 0},
  TestID -> "I011-tangential-field-preserves-Lambda"
]

VerificationTest[
  CosmologicalBaikov011`NumericalCrossCheck011[],
  {
    {-4/147, 0, 0, 0, 0},
    {-3/280, -311/2940, 0, 0, 0},
    {0, 0, -4/147, 0, 0},
    {-165/45038, -895/472899, 976/472899, -2526/22519, 264/22519},
    {-72/22519, -1366/472899, -184/157633, 954/22519, -3630/22519}
  },
  TestID -> "I011-matrix-generic-rational-point"
]

VerificationTest[
  Together[
    CosmologicalBaikov011`ReductionData011[]["FourthIBPResidual"] /. {
      CosmologicalBaikov`aK -> 5,
      CosmologicalBaikov`bK -> -1,
      CosmologicalBaikov`cK -> 1,
      CosmologicalBaikov`x2sq -> 1,
      CosmologicalBaikov`x3sq -> 2,
      CosmologicalBaikov`eps -> 1/7
    }
  ],
  0,
  TestID -> "I011-unused-fourth-v-minus-one-IBP"
]

VerificationTest[
  FreeQ[
    CosmologicalBaikov011`DifferentialMatrixA011[],
    _CosmologicalBaikov011`J011
  ],
  True,
  TestID -> "I011-matrix-has-no-unreduced-integrals"
]

VerificationTest[
  Simplify[
    CosmologicalBaikov011`CutGramPolynomial011[#] & /@
      CosmologicalBaikov011`LeadingSingularityLocations011[]
  ],
  {0, 0},
  TestID -> "I011-leading-singularity-locations"
]

VerificationTest[
  Module[{block, p, rules, transformed},
    block = CosmologicalBaikov011`DifferentialMatrixA011[][[4 ;; 5, 4 ;; 5]];
    p = CosmologicalBaikov011`LeadingSingularityBasisMatrix011[];
    rules = {
      CosmologicalBaikov`aK -> 5,
      CosmologicalBaikov`bK -> -1,
      CosmologicalBaikov`cK -> 1,
      CosmologicalBaikov`x2sq -> 1,
      CosmologicalBaikov`x3sq -> 2,
      CosmologicalBaikov`eps -> 0
    };
    transformed = Map[Together, p.block.Inverse[p] /. rules, {2}];
    {transformed[[1, 2]], transformed[[2, 1]]}
  ],
  {0, 0},
  TestID -> "I011-leading-singularity-basis-diagonal-at-epsilon-zero"
]
