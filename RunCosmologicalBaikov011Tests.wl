testFile = FileNameJoin[{
  DirectoryName[$InputFileName],
  "CosmologicalBaikov011Tests.wlt"
}];
report = TestReport[testFile];
Print[report];
If[Length[report["TestsFailedWrongResultsKeys"]] > 0 ||
   Length[report["TestsFailedWithMessagesKeys"]] > 0,
  Exit[1],
  Exit[0]
];
