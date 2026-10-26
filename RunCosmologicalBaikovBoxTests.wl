testFile = FileNameJoin[{
  DirectoryName[$InputFileName],
  "CosmologicalBaikovBoxTests.wlt"
}];
Print[TestReport[testFile]];
