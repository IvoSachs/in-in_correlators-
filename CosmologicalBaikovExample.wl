Get[FileNameJoin[{DirectoryName[$InputFileName], "CosmologicalBaikov.wl"}]];

Print["Master basis:"];
Print[CosmologicalBaikov`MasterBasis[]];

Print["Correct top-sector identity:"];
Print[CosmologicalBaikov`TopSectorIdentity[]];

Print["Constructing exact differential matrix d_a m = A.m ..."];
A = CosmologicalBaikov`DifferentialMatrixA[];
Print["Dimensions: ", Dimensions[A], "; the first three compact rows are:"];
Print[MatrixForm[Take[A, 3]]];
Print["Rows 4 and 5 are exact but intentionally not printed in raw LinearSolve form."];

Print["Generic rational cross-check:"];
Print[MatrixForm[CosmologicalBaikov`NumericalCrossCheck[]]];
