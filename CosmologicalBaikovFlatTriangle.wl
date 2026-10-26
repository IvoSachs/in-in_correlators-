(* ::Package:: *)

BeginPackage["CosmologicalBaikovFlatTriangle`"];

xF::usage = "xF denotes Q2^2.";
yF::usage = "yF denotes Q3^2.";
zF::usage = "zF denotes Q2.Q3.";
sF::usage = "sF denotes (Q2-Q3)^2=xF+yF-2 zF.";
lambdaF::usage = "lambdaF denotes the flat Gram determinant zF^2-xF yF.";

FlatOne::usage = "FlatOne is the constant component of the finite flat-space basis.";
FlatLogY::usage = "FlatLogY denotes log(sF/yF), equivalently the finite bubble difference T101-T011.";
FlatLogX::usage = "FlatLogX denotes log(sF/xF), equivalently the finite bubble difference T110-T011.";
FlatTriangle::usage = "FlatTriangle denotes the finite off-shell massless triangle T111.";

FlatFiniteBasis::usage = "FlatFiniteBasis[] returns {FlatOne,FlatLogY,FlatLogX,FlatTriangle}.";
FlatBubbleCoefficientsA::usage = "FlatBubbleCoefficientsA[] returns the coefficients of {T011,T101,T110} in the a derivative of T111 before taking finite differences.";
FlatSourceCoefficientsA::usage = "FlatSourceCoefficientsA[] returns the coefficients of {FlatLogY,FlatLogX} in the finite a-derivative source.";
FlatConnectionX::usage = "FlatConnectionX[] is the connection for partial_x in the finite flat basis.";
FlatConnectionY::usage = "FlatConnectionY[] is the connection for partial_y in the finite flat basis.";
FlatConnectionC::usage = "FlatConnectionC[] is the connection for partial_c=(1/2) partial_z in the finite flat basis.";
FlatConnectionA::usage = "FlatConnectionA[] is FlatConnectionX[], since partial_a xF=1 at fixed p2,p3,b,c.";
FlatCurvatureChecks::usage = "FlatCurvatureChecks[] returns the three zero-curvature matrices for the x,y,c connection.";
FlatSourceIntegrandA::usage = "FlatSourceIntegrandA[] returns the finite bubble-source integrand in equation (31), including its factor 1/2.";
ExtendedBlockPattern::usage = "ExtendedBlockPattern[] displays the forced lower-block-triangular structure of the nested connection in the order {rational,LogY,LogX,triangle}.";

Begin["`Private`"];

sF := xF + yF - 2 zF;
lambdaF := zF^2 - xF yF;

FlatFiniteBasis[] := {FlatOne, FlatLogY, FlatLogX, FlatTriangle};

FlatBubbleCoefficientsA[] := {
  -zF/xF,
  zF/xF + (yF - zF)/sF,
  (zF - yF)/sF
};

FlatSourceCoefficientsA[] := {
  zF/xF + (yF - zF)/sF,
  (zF - yF)/sF
};

FlatConnectionX[] := Module[{cy, cx},
  {cy, cx} = FlatSourceCoefficientsA[];
  {
    {0, 0, 0, 0},
    {1/sF, 0, 0, 0},
    {1/sF - 1/xF, 0, 0, 0},
    {0, cy/(2 lambdaF), cx/(2 lambdaF), yF/(2 lambdaF)}
  }
];

FlatConnectionY[] := Module[{cy, cx},
  cy = (zF - xF)/sF;
  cx = zF/yF + (xF - zF)/sF;
  {
    {0, 0, 0, 0},
    {1/sF - 1/yF, 0, 0, 0},
    {1/sF, 0, 0, 0},
    {0, cy/(2 lambdaF), cx/(2 lambdaF), xF/(2 lambdaF)}
  }
];

FlatConnectionC[] := Module[{cy, cx},
  cy = (zF - yF)/sF;
  cx = (zF - xF)/sF;
  {
    {0, 0, 0, 0},
    {-1/sF, 0, 0, 0},
    {-1/sF, 0, 0, 0},
    {0, cy/(2 lambdaF), cx/(2 lambdaF), -zF/(2 lambdaF)}
  }
];

FlatConnectionA[] := FlatConnectionX[];

FlatCurvatureChecks[] := Module[{ax, ay, ac, xy, xc, yc},
  ax = FlatConnectionX[];
  ay = FlatConnectionY[];
  ac = FlatConnectionC[];
  xy = D[ay, xF] - D[ax, yF] + ay.ax - ax.ay;
  xc = D[ac, xF] - 1/2 D[ax, zF] + ac.ax - ax.ac;
  yc = D[ac, yF] - 1/2 D[ay, zF] + ac.ay - ay.ac;
  Map[Together, {xy, xc, yc}, {3}]
];

FlatSourceIntegrandA[] := Module[{coefficients},
  coefficients = FlatSourceCoefficientsA[];
  1/2 (
    coefficients[[1]] FlatLogY
    + coefficients[[2]] FlatLogX
  )
];

ExtendedBlockPattern[] := {
  {"A-rational", 0, 0, 0},
  {"C-LogY-to-rational", "A-LogY", 0, 0},
  {"C-LogX-to-rational", 0, "A-LogX", 0},
  {0, "C-triangle-to-LogY", "C-triangle-to-LogX", "A-triangle"}
};

End[];
EndPackage[];
