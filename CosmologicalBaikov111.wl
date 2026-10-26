(* ::Package:: *)

If[!MemberQ[$Packages, "CosmologicalBaikov`"],
  Get[FileNameJoin[{DirectoryName[$InputFileName], "CosmologicalBaikov.wl"}]]
];
If[!MemberQ[$Packages, "CosmologicalBaikov011`"],
  Get[FileNameJoin[{DirectoryName[$InputFileName], "CosmologicalBaikov011.wl"}]]
];

BeginPackage[
  "CosmologicalBaikov111`",
  {"CosmologicalBaikov`", "CosmologicalBaikov011`"}
];

FBI::usage = "FBI[m1,m2,m3,v] is an inert integral head for Ibar^{m1 m2 m3,v 0 0;0}.";
FullMasterBasis::usage = "FullMasterBasis[] returns the eleven-master basis closing the partial_a I111 equation.";
I101ParameterRules::usage = "I101ParameterRules[] gives the q2=p3-p2, q3=p3 map from the I011 system to the I101 system.";
I101BasisTransform::usage = "I101BasisTransform[] maps the transformed I011 basis to the canonical {I000,I100,I001,I101,I101[-D4]} basis.";
DifferentialMatrixA101::usage = "DifferentialMatrixA101[] returns the five-by-five partial_a connection for the I101 pinch.";
DifferentialMatrixAFull::usage = "DifferentialMatrixAFull[] returns the eleven-by-eleven connection in partial_a mFull=AFull.mFull.";
TopSectorReductionRow::usage = "TopSectorReductionRow[] returns the completely reduced coefficient row for partial_a I111.";
ReducedTopSectorDerivative::usage = "ReducedTopSectorDerivative[] returns partial_a I111 as a linear combination of FullMasterBasis[].";
FullReductionData::usage = "FullReductionData[] returns the component connections, embeddings, full connection, and top-sector row.";
CutGramPolynomial101::usage = "CutGramPolynomial101[z] is the I101 maximal-cut Gram polynomial in the original D4 variable.";
LeadingSingularityLocations101::usage = "LeadingSingularityLocations101[] gives the two original-D4 cut locations for I101.";
LeadingSingularityBasisMatrix101::usage = "LeadingSingularityBasisMatrix101[] maps {I101,I101[-D4]} to branch-adapted combinations.";
LeadingSingularityTransformFull::usage = "LeadingSingularityTransformFull[] returns the eleven-dimensional algebraic cut-basis transformation.";
TriangularConnectionEpsilonZero::usage = "TriangularConnectionEpsilonZero[rules] evaluates the full eps=0 connection at exact kinematics and returns it in lower-triangular leading-singularity form. Substitution before the similarity transform avoids symbolic expression swell.";
ExactTriangularTransformFull::usage = "ExactTriangularTransformFull[{rho110,rho101,rho011},var] returns the finite-eps Riccati basis transformation.";
ExactTriangularRiccatiEquationsFull::usage = "ExactTriangularRiccatiEquationsFull[{rho110,rho101,rho011},var,rules] returns the three Riccati equations for exact finite-eps triangularization. Supplying exact kinematic rules while leaving var free keeps the formulas compact.";
ExactTriangularConnectionFull::usage = "ExactTriangularConnectionFull[{rho110,rho101,rho011},var] returns the exactly lower-triangular finite-eps connection after imposing the Riccati equations. Use the option \"KinematicRules\" to substitute fixed parameters before transforming.";
StrictUpperEntries::usage = "StrictUpperEntries[matrix] returns all entries strictly above the diagonal.";
OverlapConsistencyChecks::usage = "OverlapConsistencyChecks[] compares the duplicated vacuum and tadpole derivative rows among the three pinch systems.";
GenericityConditionsFull::usage = "GenericityConditionsFull[] lists the elementary nonvanishing conditions used by the combined reduction.";
NumericalCrossCheckFull::usage = "NumericalCrossCheckFull[] evaluates the exact eleven-by-eleven connection at a generic rational point.";
TildeLambda0Abbreviation::usage = "TildeLambda0Abbreviation is an inert atomic abbreviation for the kinematic polynomial TildeLambda0[].";
TildeLambda0AbbreviationRule::usage = "TildeLambda0AbbreviationRule[] replaces the explicit kinematic polynomial TildeLambda0[] by TildeLambda0Abbreviation.";
TildeLambda0ExpansionRule::usage = "TildeLambda0ExpansionRule[] expands TildeLambda0Abbreviation back to the exact kinematic polynomial.";
ExpandTildeLambda0Abbreviation::usage = "ExpandTildeLambda0Abbreviation[expr] restores the explicit polynomial represented by TildeLambda0Abbreviation.";
CompactifyConnectionEpsilonZero::usage = "CompactifyConnectionEpsilonZero[matrix] applies eps->0, abbreviates TildeLambda0, and performs the Simplify/Factor/FullSimplify rational compaction entrywise. This is the reproducible but potentially slow construction.";
CompactConnectionEpsilonZero::usage = "CompactConnectionEpsilonZero[] returns the cached exact compact eps=0 full connection.";
CompactTriangularConnectionEpsilonZero::usage = "CompactTriangularConnectionEpsilonZero[] returns the cached exact compact eps=0 connection in the full leading-singularity basis.";

Begin["`Private`"];

ClearAll[FBI];
$packageDirectory111 = DirectoryName[$InputFileName];

FullMasterBasis[] := {
  FBI[0, 0, 0, 0],
  FBI[0, 1, 0, 0],
  FBI[1, 0, 0, 0],
  FBI[0, 0, 1, 0],
  FBI[1, 1, 0, 0],
  FBI[1, 1, 0, -1],
  FBI[1, 0, 1, 0],
  FBI[1, 0, 1, -1],
  FBI[0, 1, 1, 0],
  FBI[0, 1, 1, -1],
  FBI[1, 1, 1, 0]
};

(* Put q2=p3-p2 and q3=p3.  Then
     Lambda=(2b-c) q2.q3-b q2^2-(a+b-c) q3^2+lambda0,
   while D2(q)=D1(p), D3(q)=D3(p), and
     D4(q)=D3(p)-D4(p)+b-x3sq. *)
I101ParameterRules[] := {
  aK -> aK + bK - cK,
  cK -> 2 bK - cK,
  x2sq -> x0sq
};

(* If j is the transformed I011 basis and n is the canonical I101 basis,
   n=T101.j.  In particular,
     j5=I100+(b-x3sq) I101-I101[-D4]. *)
I101BasisTransform[] := {
  {1, 0, 0, 0, 0},
  {0, 1, 0, 0, 0},
  {0, 0, 1, 0, 0},
  {0, 0, 0, 1, 0},
  {0, 1, 0, bK - x3sq, -1}
};

Options[DifferentialMatrixA101] = {"SimplifyFunction" -> Identity};
DifferentialMatrixA101[OptionsPattern[]] := Module[
  {f, mapped, transform},
  f = OptionValue["SimplifyFunction"];
  mapped = CosmologicalBaikov011`DifferentialMatrixA011[] /.
    I101ParameterRules[];
  transform = I101BasisTransform[];
  Map[f, transform.mapped.Inverse[transform], {2}]
];

ClearAll[embedRow];
embedRow[row_List, positions_List, dimension_Integer] := Module[{answer},
  answer = ConstantArray[0, dimension];
  answer[[positions]] = row;
  answer
];

ClearAll[transformDisjointBlocks];
transformDisjointBlocks[
    connection_?MatrixQ, blocks_List, transforms_List, var_: None] := Module[
  {result, inverseTransforms},
  result = connection;
  inverseTransforms = Inverse /@ transforms;
  Do[
    result[[blocks[[i]], All]] =
      transforms[[i]].result[[blocks[[i]], All]],
    {i, Length[blocks]}
  ];
  Do[
    result[[All, blocks[[i]]]] =
      result[[All, blocks[[i]]]].inverseTransforms[[i]],
    {i, Length[blocks]}
  ];
  If[var =!= None,
    Do[
      result[[blocks[[i]], blocks[[i]]]] +=
        D[transforms[[i]], var].inverseTransforms[[i]],
      {i, Length[blocks]}
    ]
  ];
  result
];

ClearAll[buildFullReductionData];
buildFullReductionData[] := Module[
  {
    a110, a101, a011, matrix, topRow, tilde,
    indices110, indices101, indices011
  },
  a110 = CosmologicalBaikov`DifferentialMatrixA[];
  a101 = DifferentialMatrixA101[];
  a011 = CosmologicalBaikov011`DifferentialMatrixA011[];

  indices110 = {1, 2, 3, 5, 6};
  indices101 = {1, 3, 4, 7, 8};
  indices011 = {1, 2, 4, 9, 10};

  matrix = ConstantArray[0, {11, 11}];

  (* Vacuum and the three one-denominator sectors. *)
  matrix[[1]] = embedRow[a110[[1]], indices110, 11];
  matrix[[2]] = embedRow[a110[[2]], indices110, 11];
  matrix[[3]] = embedRow[a110[[3]], indices110, 11];
  matrix[[4]] = embedRow[a011[[3]], indices011, 11];

  (* The three two-denominator sectors. *)
  matrix[[5]] = embedRow[a110[[4]], indices110, 11];
  matrix[[6]] = embedRow[a110[[5]], indices110, 11];
  matrix[[7]] = embedRow[a101[[4]], indices101, 11];
  matrix[[8]] = embedRow[a101[[5]], indices101, 11];
  matrix[[9]] = embedRow[a011[[4]], indices011, 11];
  matrix[[10]] = embedRow[a011[[5]], indices011, 11];

  (* Equation (6) of the supplied PDF, with the corrected second line of
     equation (7) used for tilde Lambda0. *)
  tilde = TildeLambda0[];
  topRow = (
    1/2 UnitVector[11, 5]
    + (aK - cK/2) matrix[[5]]
    + (bK - cK/2) matrix[[7]]
    + cK/2 matrix[[9]]
  )/tilde;
  topRow[[11]] = (bK - x3sq)/(2 tilde);
  matrix[[11]] = topRow;

  <|
    "Masters" -> FullMasterBasis[],
    "I110Indices" -> indices110,
    "I101Indices" -> indices101,
    "I011Indices" -> indices011,
    "A110" -> a110,
    "A101" -> a101,
    "A011" -> a011,
    "TopSectorRow" -> topRow,
    "RawMatrix" -> matrix
  |>
];

FullReductionData[] := FullReductionData[] = buildFullReductionData[];

Options[DifferentialMatrixAFull] = {"SimplifyFunction" -> Identity};
DifferentialMatrixAFull[OptionsPattern[]] := Module[{f},
  f = OptionValue["SimplifyFunction"];
  Map[f, FullReductionData[]["RawMatrix"], {2}]
];

TopSectorReductionRow[] := FullReductionData[]["TopSectorRow"];

ReducedTopSectorDerivative[] :=
  TopSectorReductionRow[].FullMasterBasis[];

CutGramPolynomial101[z_] :=
  x0sq x3sq - (z - cK/2 + x3sq)^2;

LeadingSingularityLocations101[] := {
  cK/2 - x3sq + Sqrt[x0sq x3sq],
  cK/2 - x3sq - Sqrt[x0sq x3sq]
};

LeadingSingularityBasisMatrix101[] := Module[
  {z = LeadingSingularityLocations101[]},
  {{-z[[1]], 1}, {-z[[2]], 1}}
];

LeadingSingularityTransformFull[] := Module[{transform},
  transform = IdentityMatrix[11];
  transform[[5 ;; 6, 5 ;; 6]] =
    CosmologicalBaikov`LeadingSingularityBasisMatrix[];
  transform[[7 ;; 8, 7 ;; 8]] =
    LeadingSingularityBasisMatrix101[];
  transform[[9 ;; 10, 9 ;; 10]] =
    CosmologicalBaikov011`LeadingSingularityBasisMatrix011[];
  transform
];

TriangularConnectionEpsilonZero[kinematicRules_List] := Module[
  {transform, connection},
  connection = Map[
    Together,
    DifferentialMatrixAFull[] /. eps -> 0 /. kinematicRules,
    {2}
  ];
  transform = LeadingSingularityTransformFull[] /. kinematicRules;
  Map[Together, transform.connection.Inverse[transform], {2}]
];

ExactTriangularTransformFull[
    rhos : {_Symbol, _Symbol, _Symbol}, var_Symbol] := Module[{transform},
  transform = IdentityMatrix[11];
  transform[[5 ;; 6, 5 ;; 6]] =
    TriangularizingBasisMatrix[rhos[[1]], var];
  transform[[7 ;; 8, 7 ;; 8]] =
    TriangularizingBasisMatrix[rhos[[2]], var];
  transform[[9 ;; 10, 9 ;; 10]] =
    TriangularizingBasisMatrix[rhos[[3]], var];
  transform
];

ExactTriangularRiccatiEquationsFull[
    rhos : {_Symbol, _Symbol, _Symbol}, var_Symbol,
    kinematicRules_List : {}] := Module[
  {connection, blocks},
  connection = DifferentialMatrixAFull[];
  blocks = connection[[#, #]] & /@ {
    Range[5, 6], Range[7, 8], Range[9, 10]
  };
  If[kinematicRules =!= {},
    blocks = (Map[Together, # /. kinematicRules, {2}] &) /@ blocks
  ];
  MapThread[
    TriangularizingRiccatiEquation[#1, #2, var] &,
    {blocks, rhos}
  ]
];

Options[ExactTriangularConnectionFull] = {
  "SimplifyFunction" -> Identity,
  "KinematicRules" -> {}
};
ExactTriangularConnectionFull[
    rhos : {_Symbol, _Symbol, _Symbol}, var_Symbol,
    OptionsPattern[]] := Module[
  {
    f, kinematicRules, connection, transform, transforms, blocks,
    equations, derivativeRules, result
  },
  f = OptionValue["SimplifyFunction"];
  kinematicRules = OptionValue["KinematicRules"];
  connection = DifferentialMatrixAFull[];
  If[kinematicRules =!= {},
    connection = Map[Together, connection /. kinematicRules, {2}]
  ];
  transform = ExactTriangularTransformFull[rhos, var];
  If[kinematicRules =!= {}, transform = transform /. kinematicRules];
  blocks = {Range[5, 6], Range[7, 8], Range[9, 10]};
  transforms = transform[[#, #]] & /@ blocks;
  equations = ExactTriangularRiccatiEquationsFull[
    rhos, var, kinematicRules
  ];
  derivativeRules = MapThread[
    First@First@Solve[#1, #2'[var]] &,
    {equations, rhos}
  ];
  result = transformDisjointBlocks[
      connection, blocks, transforms, var
    ] /. derivativeRules;
  Map[f, result, {2}]
];

StrictUpperEntries[matrix_?MatrixQ] := Flatten@Table[
  matrix[[i, j]],
  {i, 1, Length[matrix] - 1},
  {j, i + 1, Length[matrix]}
];

OverlapConsistencyChecks[] := Module[
  {a110, a101, a011, i110, i101, i011},
  a110 = FullReductionData[]["A110"];
  a101 = FullReductionData[]["A101"];
  a011 = FullReductionData[]["A011"];
  i110 = FullReductionData[]["I110Indices"];
  i101 = FullReductionData[]["I101Indices"];
  i011 = FullReductionData[]["I011Indices"];
  {
    embedRow[a110[[1]], i110, 11] - embedRow[a011[[1]], i011, 11],
    embedRow[a110[[1]], i110, 11] - embedRow[a101[[1]], i101, 11],
    embedRow[a110[[2]], i110, 11] - embedRow[a011[[2]], i011, 11],
    embedRow[a110[[3]], i110, 11] - embedRow[a101[[2]], i101, 11],
    embedRow[a011[[3]], i011, 11] - embedRow[a101[[3]], i101, 11]
  }
];

GenericityConditionsFull[] := And[
  GenericityConditions[],
  GenericityConditions011[],
  (GenericityConditions011[] /. I101ParameterRules[]),
  TildeLambda0[] != 0
];

NumericalCrossCheckFull[] := Together[
  DifferentialMatrixAFull[] /. {
    aK -> 5, bK -> -1, cK -> 1,
    x0sq -> 4, x2sq -> 1, x3sq -> 2,
    eps -> 1/7
  }
];

TildeLambda0AbbreviationRule[] :=
  TildeLambda0[] -> TildeLambda0Abbreviation;

TildeLambda0ExpansionRule[] :=
  TildeLambda0Abbreviation -> TildeLambda0[];

ExpandTildeLambda0Abbreviation[expr_] :=
  expr /. TildeLambda0ExpansionRule[];

CompactifyConnectionEpsilonZero[matrix_?MatrixQ] := Module[
  {abbreviated, numerators, denominators},
  abbreviated = Simplify[
    (matrix /. eps -> 0) /. TildeLambda0AbbreviationRule[]
  ];
  numerators = Factor[Numerator[abbreviated]];
  denominators = Factor[Denominator[abbreviated]];
  FullSimplify[numerators/denominators]
];

CompactConnectionEpsilonZero[] :=
  CompactConnectionEpsilonZero[] = Module[{dataFile, stored},
    dataFile = FileNameJoin[{
      $packageDirectory111,
      "CosmologicalBaikov111ConnectionEpsilonZeroCompact.expr"
    }];
    If[
      FileExistsQ[dataFile],
      stored = Get[dataFile];
      stored /. Global`TL0 -> TildeLambda0Abbreviation,
      CompactifyConnectionEpsilonZero[DifferentialMatrixAFull[]]
    ]
  ];

CompactTriangularConnectionEpsilonZero[] :=
  CompactTriangularConnectionEpsilonZero[] = Module[
    {dataFile, stored, transform, raw},
    dataFile = FileNameJoin[{
      $packageDirectory111,
      "CosmologicalBaikov111TriangularConnectionEpsilonZeroCompact.expr"
    }];
    If[
      FileExistsQ[dataFile],
      stored = Get[dataFile];
      stored /. Global`TL0 -> TildeLambda0Abbreviation,
      transform = LeadingSingularityTransformFull[];
      raw = Simplify[
        transform.CompactConnectionEpsilonZero[].Inverse[transform]
      ];
      FullSimplify[Factor[Numerator[raw]]/Factor[Denominator[raw]]]
    ]
  ];

End[];
EndPackage[];
