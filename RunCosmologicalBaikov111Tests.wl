testFile = FileNameJoin[{
  DirectoryName[$InputFileName],
  "CosmologicalBaikov111Tests.wlt"
}];
Print[TestReport[testFile]];
