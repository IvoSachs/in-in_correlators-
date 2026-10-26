testFile = FileNameJoin[{
  DirectoryName[$InputFileName],
  "CosmologicalBaikovFlatTriangleTests.wlt"
}];
Print[TestReport[testFile]];
