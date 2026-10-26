report = TestReport[FileNameJoin[{
  DirectoryName[$InputFileName],
  "CosmologicalBaikovFiniteNestedTests.wlt"
}]];
Print[report];
If[report["TestsFailed"] > 0, Exit[1], Exit[0]];
