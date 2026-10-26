testFile = FileNameJoin[{DirectoryName[$InputFileName], "CosmologicalBaikovTests.wlt"}];
report = TestReport[testFile];
Print[report];
If[report["TestsFailed"] == 0 && report["TestsSucceeded"] > 0,
  Exit[0],
  Exit[1]
];
