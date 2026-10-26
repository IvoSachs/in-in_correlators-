packageDirectory = DirectoryName[$InputFileName];
Get[FileNameJoin[{packageDirectory, "CosmologicalBaikovFlatTriangle.wl"}]];

VerificationTest[
  Together[Total[CosmologicalBaikovFlatTriangle`FlatBubbleCoefficientsA[]]],
  0,
  TestID -> "common-flat-bubble-pole-cancels"
]

VerificationTest[
  CosmologicalBaikovFlatTriangle`FlatFiniteBasis[],
  {
    CosmologicalBaikovFlatTriangle`FlatOne,
    CosmologicalBaikovFlatTriangle`FlatLogY,
    CosmologicalBaikovFlatTriangle`FlatLogX,
    CosmologicalBaikovFlatTriangle`FlatTriangle
  },
  TestID -> "finite-flat-basis"
]

VerificationTest[
  CosmologicalBaikovFlatTriangle`FlatCurvatureChecks[],
  ConstantArray[0, {3, 4, 4}],
  TestID -> "finite-flat-connection-is-integrable"
]

VerificationTest[
  Module[{ax},
    ax = CosmologicalBaikovFlatTriangle`FlatConnectionA[];
    {
      ax[[2, 1]],
      ax[[3, 1]],
      ax[[4, 4]]
    }
  ],
  {
    1/CosmologicalBaikovFlatTriangle`sF,
    1/CosmologicalBaikovFlatTriangle`sF -
      1/CosmologicalBaikovFlatTriangle`xF,
    CosmologicalBaikovFlatTriangle`yF/
      (2 CosmologicalBaikovFlatTriangle`lambdaF)
  },
  TestID -> "a-connection-diagonal-and-log-derivatives"
]
