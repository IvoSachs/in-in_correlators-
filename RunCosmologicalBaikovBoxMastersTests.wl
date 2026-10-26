testFile = FileNameJoin[{
  DirectoryName[$InputFileName],
  "CosmologicalBaikovBoxMastersTests.wlt"
}];
Print[TestReport[testFile]];
