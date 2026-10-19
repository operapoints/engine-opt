(* ::Package:: *)

(* LORd, part 1 of 2. Run this file first in a fresh kernel.
   This file clears Global` so reloading the pair always starts cleanly.

   Solid circular shafts: outer radii rC,rS,rT; lengths LC,LS,LT;
   Young's moduli EC,ES,ET; center-shaft material density rhoS.
   Overhang distributed kinetic energies are omitted, so their densities
   are not needed. mC,mT,JdC,JpC,JdT,JpT denote the separate end ROTORS.

   Main outputs:
     Ktt, Ktu, Kuu          geometry-expanded 4 x 4 stiffness blocks
     KttCompact, ...        same blocks with IC,IS,IT abbreviations
     geometryAssignments   ordered reusable scalar definitions for C++
     geometryRules         fully expanded versions of those definitions
     h,Q0,Q1,Q2            shared interpolation/integration products
     solveS[parameterRules] numerical S, without any symbolic inversion

   Ordering for each plane: u={u0,u1,u2,u3}, theta={th0,th1,th2,th3}.
   theta_z=dy/dx; theta_y=dz/dx (theta_y is sign-reversed RH rotation).
   Ktt S = -Ktu. Row/column 1 corresponds to station x0, not x1.
*)

ClearAll["Global`*"];

lengths = {LC, LS, LT};
radii = {rC, rS, rT};
geometryInputParameters = {EC, ES, ET, rC, rS, rT, LC, LS, LT, rhoS};
geometryAssumptions = And @@ Thread[geometryInputParameters > 0];

(* Ordered definitions: evaluate in this order if emitting scalar C++.
   IC,IS,IT are AREA moments, with units length^4.
   jdS,jpS are LOCAL MASS inertia per unit length.
   JdSlocal is their integrated local diametral inertia.
   JdS is the WHOLE center shaft's diametral inertia about its midpoint,
   consistently with the meaning of a rigid rotor's JdC/JdT.
   No subtraction of whole-body inertia is used in the beam energies. *)
geometryAssignments = {
  AC -> Pi rC^2, AS -> Pi rS^2, AT -> Pi rT^2,
  IC -> Pi rC^4/4, IS -> Pi rS^4/4, IT -> Pi rT^4/4,
  muS -> rhoS AS,
  mS -> muS LS,
  jdS -> rhoS IS,
  jpS -> 2 jdS,
  JdSlocal -> jdS LS,
  JpS -> jpS LS,
  JdS -> JdSlocal + mS LS^2/12
};
expandGeometry[expr_] := expr //. geometryAssignments;
geometryRules = (#[[1]] -> expandGeometry[#[[2]]]) & /@
  geometryAssignments;
(* Symbols on the left remain UNASSIGNED. This intentionally preserves
   abbreviations in compact expressions. Use /. geometryRules to expand. *)

zeroQ[expr_] := AllTrue[Flatten[{expr}],
  (TrueQ[Cancel[Together[#]] == 0] &)];

(* Derive unloaded cubic Hermite interpolation, using local coordinate s. *)
p = a0 + a1 s + a2 s^2 + a3 s^3;
coeffRules = First@Solve[{
  (p /. s -> 0) == uL,
  (D[p, s] /. s -> 0) == tL,
  (p /. s -> ell) == uR,
  (D[p, s] /. s -> ell) == tR
}, {a0, a1, a2, a3}];
h = Factor /@ D[p /. coeffRules, {{uL, tL, uR, tR}}];
dh = D[h, s];
ddh = D[h, {s, 2}];
(* For element d={uL,thetaL,uR,thetaR}, u(s)=h.d.
   Integrals of displacement, slope, and curvature products respectively. *)
Q0 = Map[Factor, Integrate[Outer[Times, h, h],
  {s, 0, ell}, Assumptions -> ell > 0], {2}];
Q1 = Map[Factor, Integrate[Outer[Times, dh, dh],
  {s, 0, ell}, Assumptions -> ell > 0], {2}];
Q2 = Map[Factor, Integrate[Outer[Times, ddh, ddh],
  {s, 0, ell}, Assumptions -> ell > 0], {2}];

(* Assemble the full one-plane stiffness before slope condensation.
   Full ordering is {u0,u1,u2,u3,th0,th1,th2,th3}. *)
U = {u0, u1, u2, u3};
Theta = {th0, th1, th2, th3};
fullPlaneCoordinates = Join[U, Theta];
rigiditiesCompact = {EC IC, ES IS, ET IT};
eye8 = IdentityMatrix[8];
elementSelectors = Table[
  eye8[[{e, 4 + e, e + 1, 5 + e}]], {e, 1, 3}];
elementStiffnessCompact = Table[
  rigiditiesCompact[[e]] (Q2 /. ell -> lengths[[e]]), {e, 1, 3}];
KplaneCompact = Total[Table[
  Transpose[elementSelectors[[e]]] . elementStiffnessCompact[[e]] .
    elementSelectors[[e]], {e, 1, 3}]];
VplaneCompact = fullPlaneCoordinates . KplaneCompact . fullPlaneCoordinates/2;

(* Extract the requested blocks by differentiating the energy. *)
KttCompact = D[VplaneCompact, {Theta, 2}];
KtuCompact = D[D[VplaneCompact, {Theta}], {U}];
KuuCompact = D[VplaneCompact, {U, 2}];
KutCompact = Transpose[KtuCompact];
Ktt = KttCompact /. geometryRules;
Ktu = KtuCompact /. geometryRules;
Kuu = KuuCompact /. geometryRules;

(* The only solve here is numerical and explicitly invoked by the user.
   Required rules: EC,ES,ET,rC,rS,rT,LC,LS,LT; rhoS is not used by S.
   No LinearSolve is evaluated with symbolic input. *)
solveS::nnum = "Supply numerical values for all stiffness inputs.";
solveS[parameterRules_List] := Module[{a, b},
  a = N[Ktt /. parameterRules];
  b = N[Ktu /. parameterRules];
  If[!MatrixQ[a, NumericQ] || !MatrixQ[b, NumericQ],
    Message[solveS::nnum]; Return[$Failed]];
  -LinearSolve[a, b]
];

offBandEntries = Flatten[Table[
  If[Abs[i - j] > 1, KttCompact[[i, j]], 0], {i, 1, 4}, {j, 1, 4}]];
stiffnessChecks = <|
  "Hermite endpoint conditions" -> zeroQ[{
    h /. s -> 0, dh /. s -> 0,
    h /. s -> ell, dh /. s -> ell} - IdentityMatrix[4]],
  "Energy Hessian matches assembled stiffness" ->
    zeroQ[D[VplaneCompact, {fullPlaneCoordinates, 2}] - KplaneCompact],
  "Ktt is symmetric" -> zeroQ[KttCompact - Transpose[KttCompact]],
  "Ktt is tridiagonal" -> zeroQ[offBandEntries],
  "Uniform translation has zero bending" -> zeroQ[
    KplaneCompact . Join[ConstantArray[1, 4], ConstantArray[0, 4]]],
  "Straight-line displacement has zero bending" -> zeroQ[
    KplaneCompact . Join[{0, LC, LC + LS, LC + LS + LT},
      ConstantArray[1, 4]]],
  "Whole-shaft inertia has expected cylinder value" -> zeroQ[
    expandGeometry[JdS - mS (LS^2/12 + rS^2/4)]]
|>;

lordGeometryReady = True;
stiffnessChecks

(* Inspect/export:
   Ktt // MatrixForm
   Ktu // MatrixForm
   geometryAssignments // TableForm
   JdS /. geometryRules
   Snum = solveS[parameterRules];
*)


Ktu


(* LORd, part 2 of 2. Run 01_shaft_stiffness.wl FIRST in the same kernel.
   This file REUSES its geometry definitions and Q0,Q1,Q2 integrals.
   It never computes a symbolic inverse or substitutes a symbolic S solve.

   S is a 4 x 4 matrix of independent SCALAR symbols s00,...,s33.
   ThetaZ=S.Y and ThetaY=S.Z, with the same S in both bending planes.
   S and all material/geometry parameters are constant in time.
   At runtime supply Snum obtained from Ktt Snum = -Ktu.

   Outputs:
     LCompact,H0Compact,MCompact,GCompact,KCompact: reuse IC,...,jpS
     L,H0,M,G,K: the same expressions expanded to geometry input symbols
     M4Compact,K4Compact,W4Compact: three symmetric 4 x 4 shared products
     coefficientAssignments,H0Coefficients: compact C++ generation route
     lagrangianChecks,elResidual: exact symbolic verification

   Distributed kinetic energy: center shaft ONLY. End rotor kinetic
   energies and all three elements' strain energies are retained.
   Gyro sign follows the supplied convention: positive Omega is CW when
   viewed from the compressor; theta_y is the reversed RH y angle.
*)

If[!TrueQ[lordGeometryReady],
  Print["Run 01_shaft_stiffness.wl first in the same kernel."]; Abort[]];

Y = {y0, y1, y2, y3}; Z = {z0, z1, z2, z3};
Ydot = {vy0, vy1, vy2, vy3}; Zdot = {vz0, vz1, vz2, vz3};
Yddot = {ay0, ay1, ay2, ay3}; Zddot = {az0, az1, az2, az3};
q = Join[Y, Z]; qdot = Join[Ydot, Zdot]; qddot = Join[Yddot, Zddot];
X = Join[q, qdot];
S = Table[Symbol["s" <> ToString[i] <> ToString[j]],
  {i, 0, 3}, {j, 0, 3}];
sVariables = Flatten[S];
ThetaZ = S . Y; ThetaY = S . Z;
ThetaZdot = S . Ydot; ThetaYdot = S . Zdot;
totalDt[expr_] := D[expr, {q}] . qdot + D[expr, {qdot}] . qddot;

(* P maps one-plane displacements to {u0..u3,theta0..theta3}.
   Be maps them to the element's {uL,thetaL,uR,thetaR}.
   The row for the outer slope remains S[[1]] or S[[4]]; outer slopes
   have not been independently approximated or discarded. *)
P = Join[IdentityMatrix[4], S];
Be = Table[elementSelectors[[e]] . P, {e, 1, 3}];
Q0S = Q0 /. ell -> LS;
Q1S = Q1 /. ell -> LS;
centerTranslationGram = Transpose[Be[[2]]] . Q0S . Be[[2]];
centerSlopeGram = Transpose[Be[[2]]] . Q1S . Be[[2]];

(* Element bending matrices in the retained displacement basis. Do not
   replace this congruence with a Schur-complement expression while S is
   symbolic: the stationarity relation is only imposed numerically later. *)
elementK4Compact = Table[
  Transpose[Be[[e]]] . elementStiffnessCompact[[e]] . Be[[e]], {e, 1, 3}];
K4Compact = Total[elementK4Compact];
{VCCompact, VSCompact, VTCompact} = Table[
  (Y . elementK4Compact[[e]] . Y + Z . elementK4Compact[[e]] . Z)/2,
  {e, 1, 3}];
By = DiagonalMatrix[{0, KbyC, KbyT, 0}];
Bz = DiagonalMatrix[{0, KbzC, KbzT, 0}];
Vb = (Y . By . Y + Z . Bz . Z)/2;
VCompact = Vb + VCCompact + VSCompact + VTCompact;

(* Center kinetic energies. muS,jdS,jpS come from file 1's definitions.
   Local inertia is evaluated directly from geometry; no JdS subtraction. *)
TStranslationalCompact = muS/2 (
  Ydot . centerTranslationGram . Ydot + Zdot . centerTranslationGram . Zdot);
TSrotaryCompact = jdS/2 (
  Ydot . centerSlopeGram . Ydot + Zdot . centerSlopeGram . Zdot);
TSgyroCompact = jpS Omega/2 (
  Y . centerSlopeGram . Zdot - Z . centerSlopeGram . Ydot);
TSCompact = TStranslationalCompact + TSrotaryCompact + TSgyroCompact;

(* Lumped rotor tilt energies use their station slopes. *)
compressorSlopeGram = Outer[Times, S[[1]], S[[1]]];
turbineSlopeGram = Outer[Times, S[[4]], S[[4]]];
rotorTranslationMass = DiagonalMatrix[{mC, 0, 0, mT}];
TCCompact = mC (vy0^2 + vz0^2)/2 +
  JdC/2 (Ydot . compressorSlopeGram . Ydot + Zdot . compressorSlopeGram . Zdot) +
  JpC Omega/2 (Y . compressorSlopeGram . Zdot - Z . compressorSlopeGram . Ydot);
TTCompact = mT (vy3^2 + vz3^2)/2 +
  JdT/2 (Ydot . turbineSlopeGram . Ydot + Zdot . turbineSlopeGram . Zdot) +
  JpT Omega/2 (Y . turbineSlopeGram . Zdot - Z . turbineSlopeGram . Ydot);
TCompact = TSCompact + TCCompact + TTCompact;
LCompact = TCompact - VCompact;

(* Shared 4 x 4 products. W4 is the polar-inertia coupling BEFORE Omega. *)
M4Compact = rotorTranslationMass +
  muS centerTranslationGram + jdS centerSlopeGram +
  JdC compressorSlopeGram + JdT turbineSlopeGram;
W4Compact = jpS centerSlopeGram +
  JpC compressorSlopeGram + JpT turbineSlopeGram;
zero4 = ConstantArray[0, {4, 4}];
MCompact = ArrayFlatten[{{M4Compact, zero4}, {zero4, M4Compact}}];
KCompact = ArrayFlatten[{{K4Compact + By, zero4}, {zero4, K4Compact + Bz}}];
GCompact = Omega ArrayFlatten[{{zero4, -W4Compact}, {W4Compact, zero4}}];

(* Actual Hessian derivation: q and qdot are independent scalar symbols.
   L is homogeneous quadratic, so evaluation at the origin is unnecessary. *)
H0Compact = D[LCompact, {X, 2}];
eomCompact = MCompact . qddot + GCompact . qdot + KCompact . q;
momentum = D[LCompact, {qdot}];
elResidual = (Cancel[Together[#]] &) /@
  (totalDt[momentum] - D[LCompact, {q}] - eomCompact);

(* Geometry-expanded outputs requested for direct code generation.
   These still contain symbolic entries of S, never its rational solve. *)
{L, H0, M, G, K} =
  {LCompact, H0Compact, MCompact, GCompact, KCompact} /. geometryRules;
{VC, VS, VT, TS, TC, TT} =
  {VCCompact, VSCompact, VTCompact, TSCompact, TCCompact, TTCompact} /.
    geometryRules;

(* Optional compact code-generation interface: evaluate the ordered
   geometryAssignments from file 1, then these 30 scalar assignments,
   then the very small H0Coefficients expressions. Each independent
   symmetric coefficient is computed once and reused in all H blocks.
   s00 means row 0, col 0 of S, i.e. the coefficient of u0 in theta0. *)
symmetricSymbols[prefix_String] := Table[
  Symbol[prefix <> ToString[Min[i, j]] <> ToString[Max[i, j]]],
  {i, 0, 3}, {j, 0, 3}];
m4Symbols = symmetricSymbols["mc"];
k4Symbols = symmetricSymbols["kc"];
w4Symbols = symmetricSymbols["wc"];
upperPositions = Flatten[Table[{i, j}, {i, 1, 4}, {j, i, 4}], 1];
coefficientAssignments = Join[
  Thread[Extract[m4Symbols, upperPositions] -> Extract[M4Compact, upperPositions]],
  Thread[Extract[k4Symbols, upperPositions] -> Extract[K4Compact, upperPositions]],
  Thread[Extract[w4Symbols, upperPositions] -> Extract[W4Compact, upperPositions]]
];
Mcoeff = ArrayFlatten[{{m4Symbols, zero4}, {zero4, m4Symbols}}];
Kcoeff = ArrayFlatten[{{k4Symbols + By, zero4}, {zero4, k4Symbols + Bz}}];
Gcoeff = Omega ArrayFlatten[{{zero4, -w4Symbols}, {w4Symbols, zero4}}];
H0Coefficients = ArrayFlatten[{{-Kcoeff, -Gcoeff/2}, {Gcoeff/2, Mcoeff}}];

lagrangianChecks = <|
  "Coordinates and velocities are independent scalar symbols" ->
    AllTrue[Join[X, qddot, sVariables], (MatchQ[#, _Symbol] &)],
  "Velocity variables independent of coordinates" -> zeroQ[D[qdot, {q}]],
  "No derivative functions in model expressions" ->
    (Cases[{LCompact, H0Compact}, _Derivative | _Function,
      Infinity, Heads -> True] === {}),
  "Reduced bending equals full congruence" ->
    zeroQ[K4Compact - Transpose[P] . KplaneCompact . P],
  "Mass symmetry" -> zeroQ[MCompact - Transpose[MCompact]],
  "Stiffness symmetry" -> zeroQ[KCompact - Transpose[KCompact]],
  "Gyroscopic skew symmetry" -> zeroQ[GCompact + Transpose[GCompact]],
  "Hessian agrees with M G K blocks" -> zeroQ[
    H0Compact - ArrayFlatten[{{-KCompact, -GCompact/2}, {GCompact/2, MCompact}}]],
  "Hessian reconstructs L" -> zeroQ[LCompact - X . H0Compact . X/2],
  "Euler-Lagrange reconstruction" -> zeroQ[elResidual],
  "Shared coefficients reconstruct Hessian" -> zeroQ[
    (H0Coefficients /. coefficientAssignments) - H0Compact]
|>;

(* Numerical use after providing the physical parameter rules:
   Snum = solveS[parameterRules];
   sRules = Thread[Flatten[S] -> Flatten[Snum]];
   Hnum = N[H0 /. sRules /. parameterRules];
   Mnum = N[M /. sRules /. parameterRules];

   A fully specified, illustrative test point (SI units):
   parameterRules = {
     EC->200.*10^9, ES->200.*10^9, ET->200.*10^9,
     rC->0.008, rS->0.006, rT->0.008,
     LC->0.04, LS->0.15, LT->0.04, rhoS->7800.,
     mC->0.2, mT->0.15, JdC->0.0001, JdT->0.00008,
     JpC->0.0002, JpT->0.00016,
     KbyC->10.^7, KbzC->10.^7, KbyT->10.^7, KbzT->10.^7,
     Omega->1000.
   };
   Snum = solveS[parameterRules];
   Norm[(Ktt /. parameterRules).Snum + (Ktu /. parameterRules)]
   sRules = Thread[Flatten[S] -> Flatten[Snum]];
   Hnum = N[H0 /. sRules /. parameterRules];

   Inspect/export:
   H0 // MatrixForm
   H0Coefficients // MatrixForm
   coefficientAssignments // TableForm
   ByteCount[H0]
   Variables[H0]
*)

lagrangianChecks


H0 = FullSimplify[H0];


SetDirectory["/home/opera/Documents/Projects/engine-opt/loworder_rd/wolfram_generations_fresh_out_the_oven"]


(* Load AFTER 01_shaft_stiffness.wl and 02_lagrangian_in_S.wl.

   Export an ordinary C++ member function with named const double groups:

     report = exportH0Cpp["LORd_model_H0_generated.cpp"];

   Defaults: class LORd_model, params type LORd_model_params_t,
   header LORd_model.hpp, member H0, method make_H0.
   Requires Eigen and C++11 or newer. Positive lengths/radii/moduli and
   physically valid inertias are assumed, as in the symbolic model.

   Pipeline: parameter copies -> useful geometry definitions -> CSE of
   Ktt/Ktu -> Eigen LLT solve for S -> CSE of 30 shared coefficients ->
   CSE of independent Hessian entries -> stores + symmetric copies.

   This implements deterministic exact-subexpression CSE using documented
   expression operations. No undocumented Experimental` result parsing,
   Compile/LibraryLink runtime, or manually invented algebra is required.
   CSE is algebraically checked before any output files are written.
   It is syntactic CSE, not a claim of globally minimal operation count.

   A standalone generated *_test.cpp is included by default. It contains
   its own minimal class/params definitions and reference Hessians from
   Mathematica, so no project header is needed to compile that test.
   Compile, for example:
     c++ -std=c++11 -O3 -DNDEBUG -I/path/to/eigen \
         LORd_model_H0_generated_test.cpp -o H0_test
     ./H0_test

   Geometry field names are rC,rS,rT and rhoS. IC,IS,IT and center inertias
   are derived locally; the OLD parameter struct must be updated to expose
   the geometry fields. The return report includes the needed declarations.
*)

ClearAll[lorCppSymbols, lorCppExpr, lorCppPower, lorCppName,
  lorArithmeticCSE, lorEmitConsts, lorEmitSet, lorCppTestCode,
  exportH0Cpp];

lorCppSymbols[expr_] := DeleteCases[
  DeleteDuplicates[Cases[expr, _Symbol, Infinity]], Pi];

lorCppName[s_Symbol] := Module[{name = SymbolName[s]},
  If[!StringMatchQ[name, RegularExpression["[A-Za-z_][A-Za-z0-9_]*"]],
    Throw[Failure["InvalidCppIdentifier", <|"Symbol" -> HoldForm[s]|>],
      "LORdCpp"]];
  name
];

(* A scalar printer for this model's arithmetic subset. All numeric
   literals are doubles, so exact Mathematica rationals cannot silently
   become integer division in C++. Small integer powers use products. *)
lorCppExpr[Pi] := "3.141592653589793238462643383279502884";
lorCppExpr[x_Integer] := ToString[x, InputForm] <> ".0";
lorCppExpr[x_Rational] := "(" <> lorCppExpr[Numerator[x]] <> " / " <>
  lorCppExpr[Denominator[x]] <> ")";
lorCppExpr[x_Real] := StringReplace[
  ToString[N[x, MachinePrecision], InputForm],
  {"*^" -> "e", RegularExpression["`[0-9.]*"] -> ""}];
lorCppExpr[x_Symbol] := lorCppName[x];
lorCppExpr[x_Plus] := "(" <>
  StringRiffle[lorCppExpr /@ (List @@ x), " + "] <> ")";
lorCppExpr[x_Times] := "(" <>
  StringRiffle[lorCppExpr /@ (List @@ x), " * "] <> ")";
lorCppPower[b_, 0] := "1.0";
lorCppPower[b_, 1] := lorCppExpr[b];
lorCppPower[b_, n_Integer /; n < 0] :=
  "(1.0 / " <> lorCppPower[b, -n] <> ")";
lorCppPower[b_, n_Integer /; 2 <= n <= 8] := "(" <>
  StringRiffle[ConstantArray["(" <> lorCppExpr[b] <> ")", n], " * "] <> ")";
lorCppPower[b_, n_] := "std::pow(" <> lorCppExpr[b] <> ", " <>
  lorCppExpr[n] <> ")";
lorCppExpr[Power[b_, n_]] := lorCppPower[b, n];
lorCppExpr[expr_] := Throw[Failure["UnsupportedCppExpression",
  <|"Expression" -> HoldForm[expr]|>], "LORdCpp"];

(* Bottom-up extraction gives topologically ordered const declarations.
   Repeated children are extracted before repeated parents. All outputs
   in a stage are optimized together, not one matrix entry at a time. *)
lorArithmeticCSE[expressions_List, prefix_String] := Module[
  {original, out, defs = {}, found, reusable, chosen, tmp, counter = 0,
   restored, verified},
  original = FactorTerms /@ expressions;
  out = original;
  While[True,
    found = Tally[Cases[out, node : (_Plus | _Times | _Power) :> node,
      Infinity]];
    reusable = Select[found,
      (Last[#] >= 2 && LeafCount[First[#]] >= 3) &];
    If[reusable === {}, Break[]];
    chosen = First[First[SortBy[reusable,
      {LeafCount[First[#]] &, ToString[First[#], InputForm] &}]]];
    tmp = Symbol["LORdCppTemps`" <> prefix <> ToString[counter++]];
    AppendTo[defs, tmp -> chosen];
    out = out /. (chosen -> tmp);
  ];
  restored = out //. Reverse[defs];
  verified = And @@ MapThread[zeroQ[#1 - #2] &, {original, restored}];
  If[!TrueQ[verified], Throw[Failure["CSEVerificationFailed",
    <|"Stage" -> prefix|>], "LORdCpp"]];
  <|"Temporaries" -> defs, "Outputs" -> out,
    "Verified" -> verified,
    "InputLeafCount" -> LeafCount[original],
    "OutputLeafCount" -> LeafCount[{Last /@ defs, out}]|>
];

lorEmitConsts[rules_List] := (
  "    const double " <> lorCppName[First[#]] <> " = " <>
    lorCppExpr[Last[#]] <> ";" &) /@ rules;
lorEmitSet[target_String, expr_] :=
  "    " <> target <> " = " <> lorCppExpr[expr] <> ";";

Options[exportH0Cpp] = {
  "ClassName" -> "LORd_model",
  "ParamsType" -> "LORd_model_params_t",
  "Header" -> "LORd_model.hpp",
  "Method" -> "make_H0",
  "Member" -> "H0",
  "WriteTest" -> True
};

exportH0Cpp[path_String, OptionsPattern[]] := Catch[Module[
  {className, paramsType, header, method, member, inputParameters,
   nonzeroPositions, all4, upper16, ktPositions, kuPositions, hPositions,
   preValues, coeffValues, hValues, needed, geometryDefs, used, unknown,
   fields, prePlan, coefficientPlan, hPlan, matrixTarget, lines, body,
   source, outputPath, testPath = None, testSource, smokeCases,
   smokeParams, parameterDeclarations, checks},

  If[!TrueQ[lordGeometryReady] || !ListQ[coefficientAssignments],
    Throw[Failure["MissingModel", <|
      "Message" -> "Load files 01 and 02 before exporting."|>], "LORdCpp"]];
  className = OptionValue["ClassName"];
  paramsType = OptionValue["ParamsType"];
  header = OptionValue["Header"];
  method = OptionValue["Method"];
  member = OptionValue["Member"];
  If[!TrueQ[And @@ (StringQ /@ {className, paramsType, header, method, member})],
    Throw[Failure["InvalidOptions", <||>], "LORdCpp"]];

  ClearAll["LORdCppTemps`*"];
  all4 = Tuples[Range[4], 2];
  nonzeroPositions[matrix_, positions_] := Select[positions,
    (!TrueQ[Extract[matrix, #] == 0] &)];
  ktPositions = nonzeroPositions[KttCompact, all4];
  kuPositions = nonzeroPositions[KtuCompact, all4];
  upper16 = Flatten[Table[{i, j}, {i, 1, 16}, {j, i, 16}], 1];
  hPositions = nonzeroPositions[H0Coefficients, upper16];
  preValues = Join[Extract[KttCompact, ktPositions],
    Extract[KtuCompact, kuPositions]];
  coeffValues = Last /@ coefficientAssignments;
  hValues = Extract[H0Coefficients, hPositions];

  (* Retain only geometry aliases used by the actual outputs, including
     dependencies. This removes reporting-only quantities such as JdS. *)
  needed = FixedPoint[Function[vars, Union[vars,
    lorCppSymbols[Last /@ Select[geometryAssignments,
      (MemberQ[vars, First[#]] &)]]]],
    lorCppSymbols[{preValues, coeffValues, hValues}]];
  geometryDefs = Select[geometryAssignments, (MemberQ[needed, First[#]] &)];
  inputParameters = {EC, ES, ET, rC, rS, rT, LC, LS, LT, rhoS,
    mC, mT, JdC, JdT, JpC, JpT,
    KbyC, KbyT, KbzC, KbzT, Omega};
  used = lorCppSymbols[{Last /@ geometryDefs, preValues, coeffValues, hValues}];
  unknown = Complement[used, Join[inputParameters, First /@ geometryDefs,
    sVariables, First /@ coefficientAssignments]];
  If[unknown =!= {}, Throw[Failure["UnboundModelSymbols",
    <|"Symbols" -> unknown|>], "LORdCpp"]];
  fields = Select[inputParameters, (MemberQ[used, #] &)];

  prePlan = lorArithmeticCSE[preValues, "lorPre"];
  coefficientPlan = lorArithmeticCSE[coeffValues, "lorCoeff"];
  hPlan = lorArithmeticCSE[hValues, "lorStore"];
  matrixTarget[name_String, ij_List] := name <> "(" <>
    ToString[ij[[1]] - 1] <> ", " <> ToString[ij[[2]] - 1] <> ")";

  lines = Join[
    {"    auto& H0 = this->" <> member <> ";"},
    ("    const double " <> lorCppName[#] <> " = params." <>
      lorCppName[#] <> ";" &) /@ fields,
    {"", "    // Geometry, in dependency order."},
    lorEmitConsts[geometryDefs],
    {"", "    // Slope stiffness blocks; no symbolic inversion.",
     "    Eigen::Matrix4d Ktt = Eigen::Matrix4d::Zero();",
     "    Eigen::Matrix4d Ktu = Eigen::Matrix4d::Zero();"},
    lorEmitConsts[prePlan["Temporaries"]],
    MapThread[lorEmitSet, {
      Join[matrixTarget["Ktt", #] & /@ ktPositions,
        matrixTarget["Ktu", #] & /@ kuPositions], prePlan["Outputs"]}],
    {"    Eigen::LLT<Eigen::Matrix4d> slopeFactor(Ktt);",
     "    assert(slopeFactor.info() == Eigen::Success);",
     "    const Eigen::Matrix4d S = -slopeFactor.solve(Ktu);", ""},
    Flatten[Table[
      "    const double " <> lorCppName[S[[i, j]]] <> " = " <>
        matrixTarget["S", {i, j}] <> ";", {i, 1, 4}, {j, 1, 4}]],
    {"", "    // CSE across all 30 independent matrix coefficients."},
    lorEmitConsts[coefficientPlan["Temporaries"]],
    lorEmitConsts[Thread[(First /@ coefficientAssignments) ->
      coefficientPlan["Outputs"]]],
    {"", "    // Zero all entries, then fill the upper triangle and mirror.",
     "    H0.setZero();"},
    lorEmitConsts[hPlan["Temporaries"]],
    Flatten[MapThread[Function[{ij, value},
      If[ij[[1]] == ij[[2]],
        {lorEmitSet[matrixTarget["H0", ij], value]},
        {lorEmitSet[matrixTarget["H0", ij], value],
         "    " <> matrixTarget["H0", Reverse[ij]] <> " = " <>
           matrixTarget["H0", ij] <> ";"}]],
      {hPositions, hPlan["Outputs"]}]]
  ];
  body = "void " <> className <> "::" <> method <> "(const " <>
    paramsType <> " params) {\n" <> StringRiffle[lines, "\n"] <> "\n}\n";
  source = "// Generated by 03_export_H0_cpp.wl; edit the model, then regenerate.\n" <>
    "#include <Eigen/Cholesky>\n#include <cassert>\n#include <cmath>\n" <>
    If[header === "", "", "#include \"" <> header <> "\"\n"] <> "\n" <> body;
  outputPath = ExpandFileName[path];
  parameterDeclarations = ("    double " <> lorCppName[#] <> ";" &) /@ fields;

  (* References are evaluated directly from H0, independently of the CSE
     declarations. Check both spin directions, rest, and different geometry. *)
  If[TrueQ[OptionValue["WriteTest"]],
    smokeParams = {
      EC -> 200.*10^9, ES -> 190.*10^9, ET -> 180.*10^9,
      rC -> .008, rS -> .006, rT -> .007,
      LC -> .04, LS -> .15, LT -> .05, rhoS -> 7800.,
      mC -> .2, mT -> .15, JdC -> .0001, JdT -> .00008,
      JpC -> .0002, JpT -> .00016,
      KbyC -> 1.*10^7, KbzC -> 1.2*10^7,
      KbyT -> .9*10^7, KbzT -> 1.1*10^7, Omega -> 1000.
    };
    smokeCases = {
      smokeParams,
      Join[{Omega -> 0.}, smokeParams],
      Join[{Omega -> -1700.}, smokeParams],
      Join[{rS -> .009, LS -> .11, LC -> .06, ES -> 120.*10^9}, smokeParams]
    };
    testSource = lorCppTestCode[body, className, paramsType, method, member,
      fields, parameterDeclarations, smokeCases];
    testPath = FileNameJoin[{DirectoryName[outputPath],
      FileBaseName[outputPath] <> "_test.cpp"}];
  ];
  checks = <|"Stiffness CSE" -> prePlan["Verified"],
    "Coefficient CSE" -> coefficientPlan["Verified"],
    "Hessian-store CSE" -> hPlan["Verified"]|>;
  If[Export[outputPath, source, "Text"] === $Failed,
    Throw[Failure["SourceWriteFailed", <|"Path" -> outputPath|>], "LORdCpp"]];
  If[StringQ[testPath] && Export[testPath, testSource, "Text"] === $Failed,
    Throw[Failure["TestWriteFailed", <|"Path" -> testPath|>], "LORdCpp"]];
  <|"Source" -> outputPath, "TestSource" -> testPath,
    "RequiredParameterFields" -> StringRiffle[parameterDeclarations, "\n"],
    "Checks" -> checks,
    "TemporaryCounts" -> <|
      "Stiffness" -> Length[prePlan["Temporaries"]],
      "Coefficients" -> Length[coefficientPlan["Temporaries"]],
      "Stores" -> Length[hPlan["Temporaries"]]|>,
    "SourceBytes" -> StringLength[source]|>
], "LORdCpp"];

lorCppTestCode[body_String, className_String, paramsType_String,
  method_String, member_String, fields_List, declarations_List,
  cases_List] := Module[{testBlocks, rules, sn, hn, vals, block, i},
  testBlocks = Table[
    rules = cases[[i]];
    sn = solveS[rules];
    If[sn === $Failed, Throw[Failure["TestSlopeSolveFailed", <||>], "LORdCpp"]];
    hn = N[H0 /. Thread[Flatten[S] -> Flatten[sn]] /. rules];
    If[!MatrixQ[hn, NumericQ],
      Throw[Failure["TestReferenceNotNumeric", <||>], "LORdCpp"]];
    vals = Flatten[hn];
    block = Join[
      {"    {", "        " <> paramsType <> " params{};"},
      ("        params." <> lorCppName[#] <> " = " <>
        lorCppExpr[# /. rules] <> ";" &) /@ fields,
      {"        model." <> method <> "(params);",
       "        const double reference[256] = {\n            " <>
        StringRiffle[lorCppExpr /@ vals, ", "] <> "\n        };",
       "        for (int r = 0; r < 16; ++r) {",
       "            for (int c = 0; c < 16; ++c) {",
       "                const double expected = reference[16*r+c];",
       "                const double actual = model." <> member <> "(r,c);",
       "                const double error = std::abs(actual - expected);",
       "                const double tolerance = 1e-8 + 1e-10*std::abs(expected);",
       "                if (!std::isfinite(actual) || error > tolerance) {",
       "                    std::cerr << \"Hessian mismatch: case " <> ToString[i] <>
         ", entry (\" << r << \",\" << c << \")\\n\";",
       "                    return 1;",
       "                }", "            }", "        }", "    }"}
    ];
    StringRiffle[block, "\n"], {i, Length[cases]}];
  "// Standalone reference test: no project header required.\n" <>
  "#include <Eigen/Core>\n#include <Eigen/Cholesky>\n" <>
  "#include <cassert>\n#include <cmath>\n#include <iostream>\n\n" <>
  "struct " <> paramsType <> " {\n" <>
    StringRiffle[declarations, "\n"] <> "\n};\n\n" <>
  "struct " <> className <> " {\n" <>
  "    Eigen::Matrix<double,16,16> " <> member <> ";\n" <>
  "    void " <> method <> "(const " <> paramsType <> " params);\n};\n\n" <>
  body <> "\nint main() {\n    " <> className <> " model;\n" <>
  StringRiffle[testBlocks, "\n"] <>
  "\n    std::cout << \"All generated Hessian tests passed.\\n\";\n" <>
  "    return 0;\n}\n"
];


   Get["01_shaft_stiffness.wl"];
   Get["02_lagrangian_in_S.wl"];
   Get["03_export_H0_cpp.wl"];
   report = exportH0Cpp["LORd_model_H0_generated.cpp"];
   report["Checks"]
   report["RequiredParameterFields"]

   Custom project header:
   exportH0Cpp["LORd_model_H0_generated.cpp", "Header" -> "lord_model.h"]



report["Source"]
report["TestSource"]



