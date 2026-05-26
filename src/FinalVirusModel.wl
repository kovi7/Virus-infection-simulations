(* ::Package:: *)

(* ============================================================ *)
(* FinalVirusModel.wl                                           *)
(* Projekt: przestrzenne modelowanie infekcji wirusowej          *)
(* Wolfram Mathematica                                           *)
(* ============================================================ *)
(*
   Cel pliku:
   - Jeden wspolny kod modelu dla wszystkich notebookow.
   - Model prosty: automat komorkowy bez wolnych wirionow.
   - Model zaawansowany: komorki dyskretne + pola ciagle V, C, A.
   - Kazda symulacja z klatkami startuje z wirusem w centrum siatki.
     Jesli InitialInfected > 1, jedna komorka jest w centrum, reszta losowo.

   Uwaga merytoryczna:
   To jest model dydaktyczny fragmentu 2D tkanki, nie model calego organizmu
   ani kalibracja kliniczna konkretnej choroby.
*)

ClearAll["Global`*"];

(* ============================================================ *)
(* 1. Stale stanow komorek                                      *)
(* ============================================================ *)

healthy = 0;          (* zdrowa komorka podatna na infekcje *)
latent = 1;           (* zakazona, faza lag, jeszcze nie produkuje wirionow *)
infectious = 2;       (* zakazona, produkuje wolne wiriony *)
dead = 3;             (* martwa po infekcji *)
protected = 4;        (* komorka chroniona, np. szczepienie *)
mutant = 5;           (* wariant ucieczkowy, slabo neutralizowany *)
deadDrug = 6;         (* smierc toksyczna po leczeniu, opcjonalna *)

stateColorRules = <|
   healthy -> White,
   latent -> Orange,
   infectious -> Red,
   dead -> Black,
   protected -> LightGreen,
   mutant -> Yellow,
   deadDrug -> Purple
|>;

stateLegendRules[] := {
   healthy -> "healthy",
   latent -> "lag",
   infectious -> "infectious",
   dead -> "dead",
   protected -> "protected",
   mutant -> "mutant",
   deadDrug -> "drug-toxic"
};

(* ============================================================ *)
(* 2. Narzedzia numeryczne i pomocnicze                         *)
(* ============================================================ *)

ClearAll[MergeParameters, maskEqual, maskMember, boundedNeighbors, mooreCount,
  laplacian2D, clipNonnegative, hill, normalizePairs, centerPosition,
  initialInfectionPositions, safeRandomSample];

MergeParameters[base_Association, overrides_: <||>] := Join[base, overrides];

maskEqual[m_, value_] := Map[If[# === value, 1.0, 0.0] &, m, {2}];
maskMember[m_, values_List] := Map[If[MemberQ[values, #], 1.0, 0.0] &, m, {2}];

boundedNeighbors[{i_, j_}, size_Integer] := Select[
   ({i, j} + #) & /@ DeleteCases[Tuples[{-1, 0, 1}, 2], {0, 0}],
   1 <= #[[1]] <= size && 1 <= #[[2]] <= size &
];

mooreKernel = {{1, 1, 1}, {1, 0, 1}, {1, 1, 1}};
lapKernel = {{0, 1, 0}, {1, -4, 1}, {0, 1, 0}};

mooreCount[m_] := ListConvolve[mooreKernel, ArrayPad[m, 1, 0]];
laplacian2D[m_] := ListConvolve[lapKernel, ArrayPad[m, 1, "Fixed"]];

clipNonnegative[m_] := Map[If[# < 0, 0.0, N[#]] &, m, {2}];

hill[x_, k_, n_] := If[x <= 0, 0.0, x^n/(k^n + x^n + 10^-12)];

normalizePairs[pairs_] := Module[{ys, ymin, ymax},
   ys = pairs[[All, 2]];
   ymin = Min[ys];
   ymax = Max[ys];
   If[Abs[ymax - ymin] < 10^-12,
      ({#[[1]], 0.0} & /@ pairs),
      ({#[[1]], (#[[2]] - ymin)/(ymax - ymin)} & /@ pairs)
   ]
];

centerPosition[size_Integer] := {Ceiling[size/2], Ceiling[size/2]};

safeRandomSample[list_, n_Integer] := If[Length[list] == 0 || n <= 0, {}, RandomSample[list, Min[n, Length[list]]]];

(*
   Wazne zalozenie projektu:
   Kazda symulacja wizualna startuje od centrum siatki.
   Jesli InitialInfected > 1, jedna komorka jest w centrum, a reszta losowo.
*)
initialInfectionPositions[p_Association] := Module[
   {size, n, placement, center, all, rest, radius, disk},
   size = p["Size"];
   n = p["InitialInfected"];
   placement = Lookup[p, "InitialPlacement", "CenterPlusRandom"];
   center = centerPosition[size];
   all = DeleteCases[Tuples[Range[size], 2], center];
   Switch[placement,
      "Center",
         Take[SortBy[Tuples[Range[size], 2], EuclideanDistance[#, center] &], n],
      "CenterDisk",
         radius = Lookup[p, "InitialRadius", 2];
         disk = DeleteCases[Select[Tuples[Range[size], 2], EuclideanDistance[#, center] <= radius &], center];
         Join[{center}, safeRandomSample[disk, n - 1]],
      "Random",
         (* opcja techniczna, ale w notebookach uzywamy CenterPlusRandom *)
         RandomSample[Tuples[Range[size], 2], n],
      _,
         Join[{center}, safeRandomSample[all, n - 1]]
   ]
];


(* ============================================================ *)
(* 3. Parametry modeli                                          *)
(* ============================================================ *)

ClearAll[DefaultSimpleParameters, DefaultAdvancedParameters, ScenarioParameters,
  SimpleScenarioParameters];

DefaultSimpleParameters[] := <|
   "Size" -> 80,
   "Steps" -> 100,
   "SnapshotEvery" -> 2,
   "InitialInfected" -> 1,
   "InitialPlacement" -> "CenterPlusRandom",
   "InitialRadius" -> 2,
   "BetaContact" -> 0.16,
   "DeathTime" -> 18,
   "DeathProbability" -> 0.0,
   "VaccineCoverage" -> 0.0,
   "MutationProb" -> 0.0,
   "Seed" -> 123
|>;

DefaultAdvancedParameters[] := <|
   "Size" -> 80,
   "Steps" -> 100,
   "SnapshotEvery" -> 2,
   "dt" -> 1.0,
   "dx" -> 1.0,
   "InitialInfected" -> 1,
   "InitialPlacement" -> "CenterPlusRandom",
   "InitialRadius" -> 2,
   "InitialVirus" -> 15.0,
   "InitialRNA" -> 0.10,
   
   (* transmisja *)
   "BetaFree" -> 0.030,
   "BetaContact" -> 0.070,
   "Vhalf" -> 8.0,
   
   (* wiriony V(t,x,y) *)
   "Dv" -> 0.12,
   "VirusClearance" -> 0.045,
   "ReleaseRate" -> 2.8,
   
   (* wewnatrzkomorkowa replikacja RNA *)
   "LagTime" -> 7,
   "RNAReplicationRate" -> 0.36,
   "RNACarryingCapacity" -> 1.0,
   "RNAThreshold" -> 0.45,
   
   (* smierc komorek *)
   "DeathTime" -> 32,
   "DeathProbability" -> 0.006,
   
   (* pole chemotaktyczne C(t,x,y) *)
   "UseChemokine" -> True,
   "Dc" -> 0.18,
   "ChemokineDecay" -> 0.10,
   "ChemokineFromInfected" -> 0.80,
   "ChemokineFromDead" -> 0.25,
   
   (* komorki odpornosciowe *)
   "UseImmuneCells" -> False,
   "ImmuneInitialCount" -> 0,
   "ImmuneActivationTime" -> 30,
   "ImmuneRecruitmentPerStep" -> 0,
   "ImmuneMovement" -> "Chemotaxis",       (* "Random" albo "Chemotaxis" *)
   "ChemotaxisStrength" -> 0.85,
   "ImmuneKillProbability" -> 0.25,
   "ImmuneVirusClearance" -> 0.09,
   
   (* przeciwciala A(t,x,y) *)
   "UseAntibodies" -> False,
   "Da" -> 0.20,
   "AntibodyDecay" -> 0.060,
   "AntibodyFromVirus" -> 0.22,
   "AntibodyFromInfected" -> 0.12,
   "AntibodyHillK" -> 0.020,
   "AntibodyHillN" -> 2,
   "AntibodyVirusClearance" -> 0.70,
   "InitialAntibody" -> 0.0,
   
   (* leczenie *)
   "UseDrug" -> False,
   "DrugStart" -> 99999,
   "DrugDose" -> 0.0,
   "DrugElimination" -> 0.06,
   "DoseInterval" -> 0,
   "DrugToxicity" -> 0.0,
   "DrugProductionEffect" -> 1.0,
   "DrugInfectionEffect" -> 0.75,
   "DrugInfectedKilling" -> 0.0,
   
   (* szczepienie *)
   "VaccineCoverage" -> 0.0,
   "VaccineProtection" -> 0.92,
   
   (* mutacje *)
   "UseMutation" -> False,
   "MutationProb" -> 0.0,
   "MutantAntibodyEscape" -> 0.75,
   "MutantBetaMultiplier" -> 1.15,
   
   "Seed" -> 123
|>;

SimpleScenarioParameters[name_String] := Module[{p = DefaultSimpleParameters[]},
   Switch[name,
      "SimpleBaseline",
         p,
      "SimpleFastProfile",
         Join[p, <|"BetaContact" -> 0.28, "DeathTime" -> 13, "Steps" -> 120, "Seed" -> 11|>],
      "SimpleSlowDamagingProfile",
         Join[p, <|"BetaContact" -> 0.075, "DeathTime" -> 34, "Steps" -> 180, "Seed" -> 12|>],
      _, p
   ]
];

ScenarioParameters[name_String] := Module[{p = DefaultAdvancedParameters[]},
   Switch[name,
      "BaselineAdvanced",
         p,
      
      "SARSCoV2Like",
         Join[p, <|
            "BetaFree" -> 0.040, "BetaContact" -> 0.070, "Dv" -> 0.16,
            "ReleaseRate" -> 3.2, "LagTime" -> 5, "DeathTime" -> 24,
            "InitialInfected" -> 1, "Seed" -> 21
         |>],
      "HCVLike",
         Join[p, <|
            "BetaFree" -> 0.010, "BetaContact" -> 0.115, "Dv" -> 0.07,
            "ReleaseRate" -> 1.6, "LagTime" -> 16, "DeathTime" -> 55,
            "RNAReplicationRate" -> 0.22, "InitialInfected" -> 1, "Seed" -> 22
         |>],
      
      "CellFreeDominant",
         Join[p, <|"BetaFree" -> 0.075, "BetaContact" -> 0.010, "Dv" -> 0.20, "ReleaseRate" -> 3.8, "Seed" -> 31|>],
      "CellToCellDominant",
         Join[p, <|"BetaFree" -> 0.005, "BetaContact" -> 0.130, "Dv" -> 0.05, "ReleaseRate" -> 1.6, "Seed" -> 32|>],
      "BothTransmissionModes",
         Join[p, <|"BetaFree" -> 0.030, "BetaContact" -> 0.080, "Dv" -> 0.12, "ReleaseRate" -> 2.8, "Seed" -> 33|>],
      
      "LowDiffusion",
         Join[p, <|"Dv" -> 0.03, "Seed" -> 41|>],
      "MediumDiffusion",
         Join[p, <|"Dv" -> 0.14, "Seed" -> 42|>],
      "HighDiffusion",
         Join[p, <|"Dv" -> 0.55, "Seed" -> 43|>],
      
      "ShortLag",
         Join[p, <|"LagTime" -> 2, "RNAReplicationRate" -> 0.55, "Seed" -> 51|>],
      "LongLag",
         Join[p, <|"LagTime" -> 18, "RNAReplicationRate" -> 0.18, "Seed" -> 52|>],
      
      "RandomImmune",
         Join[p, <|
            "UseImmuneCells" -> True, "ImmuneInitialCount" -> 45, "ImmuneActivationTime" -> 18,
            "ImmuneMovement" -> "Random", "UseChemokine" -> True, "Seed" -> 61
         |>],
      "ChemotaxisImmune",
         Join[p, <|
            "UseImmuneCells" -> True, "ImmuneInitialCount" -> 45, "ImmuneActivationTime" -> 18,
            "ImmuneMovement" -> "Chemotaxis", "UseChemokine" -> True, "Seed" -> 61
         |>],
      
      "FastImmuneProfile",
         Join[p, <|
            "UseImmuneCells" -> True, "UseAntibodies" -> True,
            "ImmuneActivationTime" -> 10, "ImmuneInitialCount" -> 70,
            "ImmuneKillProbability" -> 0.32, "AntibodyFromVirus" -> 0.33,
            "AntibodyVirusClearance" -> 0.95, "Seed" -> 71
         |>],
      "SlowImmuneProfile",
         Join[p, <|
            "UseImmuneCells" -> True, "UseAntibodies" -> True,
            "ImmuneActivationTime" -> 38, "ImmuneInitialCount" -> 18,
            "ImmuneKillProbability" -> 0.16, "AntibodyFromVirus" -> 0.12,
            "AntibodyVirusClearance" -> 0.42, "Seed" -> 72
         |>],
      
      "FirstExposure",
         Join[p, <|
            "UseImmuneCells" -> True, "UseAntibodies" -> True,
            "ImmuneActivationTime" -> 32, "ImmuneInitialCount" -> 25,
            "InitialAntibody" -> 0.0, "AntibodyFromVirus" -> 0.14, "Seed" -> 81
         |>],
      "SecondExposureMemory",
         Join[p, <|
            "UseImmuneCells" -> True, "UseAntibodies" -> True,
            "ImmuneActivationTime" -> 8, "ImmuneInitialCount" -> 55,
            "InitialAntibody" -> 0.10, "AntibodyFromVirus" -> 0.32, "Seed" -> 81
         |>],
      
      "NoDrug",
         Join[p, <|
            "BetaContact" -> 0.080,
            "BetaFree" -> 0.055,
            "Dv" -> 0.12,
            "ReleaseRate" -> 3.8,
            "VirusClearance" -> 0.035,
            "LagTime" -> 3,
            "InitialRNA" -> 0.75,
            "RNAThreshold" -> 0.18,
            "RNAReplicationRate" -> 0.45,
            "DeathTime" -> 28,
            "DeathProbability" -> 0.001,
            "UseImmuneCells" -> False,
            "UseAntibodies" -> False,
            "UseChemokine" -> False,
            "UseMutation" -> False,
            "UseDrug" -> False,
            "Seed" -> 91
         |>],
         
      "EarlyDrug",
         Join[p, <|
            "BetaContact" -> 0.080,
            "BetaFree" -> 0.055,
            "Dv" -> 0.12,
            "ReleaseRate" -> 3.8,
            "VirusClearance" -> 0.035,
            "LagTime" -> 3,
            "InitialRNA" -> 0.75,
            "RNAThreshold" -> 0.18,
            "RNAReplicationRate" -> 0.45,
            "DeathTime" -> 28,
            "DeathProbability" -> 0.001,
            "UseImmuneCells" -> False,
            "UseAntibodies" -> False,
            "UseChemokine" -> False,
            "UseMutation" -> False,
            "UseDrug" -> True,
		    "DrugStart" -> 8,
		    "DrugDose" -> 0.96,
		    "DrugElimination" -> 0.005,
		    "DoseInterval" -> 0,
		    "DrugToxicity" -> 0.0003,
		    "DrugProductionEffect" -> 1.0,
		    "DrugInfectionEffect" -> 0.99,
		    "DrugInfectedKilling" -> 0.04,
		    "Seed" -> 91
         |>],
         
      "LateDrug",
         Join[p, <|
            "BetaContact" -> 0.080,
            "BetaFree" -> 0.055,
            "Dv" -> 0.12,
            "ReleaseRate" -> 3.8,
            "VirusClearance" -> 0.035,
            "LagTime" -> 3,
            "InitialRNA" -> 0.75,
            "RNAThreshold" -> 0.18,
            "RNAReplicationRate" -> 0.45,
            "DeathTime" -> 28,
            "DeathProbability" -> 0.001,
            "UseImmuneCells" -> False,
            "UseAntibodies" -> False,
            "UseChemokine" -> False,
            "UseMutation" -> False,
            "UseDrug" -> True,
	        "DrugStart" -> 40,
            "DrugDose" -> 0.98,
	        "DrugElimination" -> 0.005,
	        "DoseInterval" -> 0,
	        "DrugToxicity" -> 0.0003,
	        "DrugProductionEffect" -> 1.0,
	        "DrugInfectionEffect" -> 0.99,
	        "DrugInfectedKilling" -> 0.04,
            "Seed" -> 91
         |>],
         
      "SingleDose",
         Join[p, <|
            "BetaContact" -> 0.080,
            "BetaFree" -> 0.055,
            "Dv" -> 0.12,
            "ReleaseRate" -> 3.8,
            "VirusClearance" -> 0.035,
            "LagTime" -> 3,
            "InitialRNA" -> 0.75,
            "RNAThreshold" -> 0.18,
            "RNAReplicationRate" -> 0.45,
            "DeathTime" -> 28,
            "DeathProbability" -> 0.001,
            "UseImmuneCells" -> False,
            "UseAntibodies" -> False,
            "UseChemokine" -> False,
            "UseMutation" -> False,
            "UseDrug" -> True,
            "DrugStart" -> 10,
	        "DrugDose" -> 0.95,
	        "DrugElimination" -> 0.012,
	        "DoseInterval" -> 0,
	        "DrugToxicity" -> 0.0005,
	        "DrugProductionEffect" -> 1.0,
	        "DrugInfectionEffect" -> 0.90,
	        "DrugInfectedKilling" -> 0.030,
            "Seed" -> 92
         |>],
         
      "CyclicDose",
         Join[p, <|
            "BetaContact" -> 0.080,
            "BetaFree" -> 0.055,
            "Dv" -> 0.12,
            "ReleaseRate" -> 3.8,
            "VirusClearance" -> 0.035,
            "LagTime" -> 3,
            "InitialRNA" -> 0.75,
            "RNAThreshold" -> 0.18,
            "RNAReplicationRate" -> 0.45,
            "DeathTime" -> 28,
            "DeathProbability" -> 0.001,
            "UseImmuneCells" -> False,
            "UseAntibodies" -> False,
            "UseChemokine" -> False,
            "UseMutation" -> False,
            "UseDrug" -> True,
            "DrugStart" -> 10,
	        "DrugDose" -> 0.95,
	        "DrugElimination" -> 0.025,
	        "DoseInterval" -> 25,
	        "DrugToxicity" -> 0.0005,
	        "DrugProductionEffect" -> 1.0,
	        "DrugInfectionEffect" -> 0.95,
	        "DrugInfectedKilling" -> 0.035,
            "Seed" -> 92
         |>],
      
      "NoVaccination",
         Join[p, <|"VaccineCoverage" -> 0.0, "BetaContact" -> 0.25, "BetaFree" -> 0.15,  "Seed" -> 101|>],
      "Vaccination",
         Join[p, <|"VaccineCoverage" -> 0.38, "InitialAntibody" -> 0.05, "UseAntibodies" -> True, "Seed" -> 101|>],
      
      "StableVirusWithVaccine",
		   Join[p, <|
		      "BetaContact" -> 0.35,           
              "BetaFree" -> 0.25,
		      "UseMutation" -> False,
		      "VaccineCoverage" -> 0.50,
		      "VaccineProtection" -> 0.95,
		      "UseAntibodies" -> True,
		      "InitialAntibody" -> 0.16,
		      "AntibodyFromVirus" -> 0.30,
		      "AntibodyVirusClearance" -> 0.90,
		      "Seed" -> 111
		   |>],
		"MutationEscape",
		   Join[p, <|
		      "BetaContact" -> 0.35,           
              "BetaFree" -> 0.25,
		      "UseMutation" -> True,
		      "MutationProb" -> 0.30,
		      "MutantAntibodyEscape" -> 0.99,
		      "MutantBetaMultiplier" -> 5.0,
		      "VaccineCoverage" -> 0.50,
		      "VaccineProtection" -> 0.95,
		      "UseAntibodies" -> True,
		      "InitialAntibody" -> 0.16,
		      "AntibodyFromVirus" -> 0.30,
		      "AntibodyVirusClearance" -> 0.90,
		      "Seed" -> 111
		   |>],
      _, p
   ]
];

(* ============================================================ *)
(* 4. Model prosty                                              *)
(* ============================================================ *)

ClearAll[RunSimpleCellularModel, RunSimpleScenario];

RunSimpleCellularModel[overrides_: <||>] := Module[
   {p, size, steps, snapshotEvery, state, age, positions, pos, i, j, t,
    newState, newAge, infMask, nInf, prob, hist, snapshots, susceptibleQ},
   
   p = MergeParameters[DefaultSimpleParameters[], overrides];
   SeedRandom[p["Seed"]];
   size = p["Size"];
   steps = p["Steps"];
   snapshotEvery = p["SnapshotEvery"];
   
   state = ConstantArray[healthy, {size, size}];
   age = ConstantArray[0, {size, size}];
   
   If[p["VaccineCoverage"] > 0,
      Do[If[RandomReal[] < p["VaccineCoverage"], state[[i, j]] = protected], {i, size}, {j, size}]
   ];
   
   positions = initialInfectionPositions[p];
   Do[
      state[[Sequence @@ pos]] = latent;
      age[[Sequence @@ pos]] = 0,
      {pos, positions}
   ];
   
   hist = {};
   snapshots = {};
   
   Do[
      If[Mod[t - 1, snapshotEvery] == 0,
         AppendTo[snapshots, <|"Time" -> t - 1, "State" -> state, "Virus" -> ConstantArray[0.0, {size, size}], "ImmuneCells" -> {}|>]
      ];
      
      infMask = maskMember[state, {latent, infectious, mutant}];
      nInf = mooreCount[infMask];
      newState = state;
      newAge = age;
      
      Do[
         Switch[state[[i, j]],
            healthy,
               prob = 1 - (1 - p["BetaContact"])^nInf[[i, j]];
               If[RandomReal[] < prob,
                  newState[[i, j]] = latent;
                  newAge[[i, j]] = 0;
               ],
            protected,
               Null,
            latent | infectious | mutant,
               newAge[[i, j]] = age[[i, j]] + 1;
               If[p["MutationProb"] > 0 && state[[i, j]] =!= mutant && RandomReal[] < p["MutationProb"],
                  newState[[i, j]] = mutant
               ];
               If[newAge[[i, j]] >= p["DeathTime"] || RandomReal[] < p["DeathProbability"],
                  newState[[i, j]] = dead
               ],
            _, Null
         ],
         {i, size}, {j, size}
      ];
      
      state = newState;
      age = newAge;
      
      AppendTo[hist, <|
         "Time" -> t,
         "Healthy" -> Count[Flatten[state], healthy],
         "Latent" -> Count[Flatten[state], latent],
         "Infectious" -> Count[Flatten[state], infectious],
         "Mutant" -> Count[Flatten[state], mutant],
         "Dead" -> Count[Flatten[state], dead],
         "Protected" -> Count[Flatten[state], protected],
         "TotalVirus" -> 0.0,
         "Drug" -> 0.0
      |>],
      {t, 1, steps}
   ];
   
   AppendTo[snapshots, <|"Time" -> steps, "State" -> state, "Virus" -> ConstantArray[0.0, {size, size}], "ImmuneCells" -> {}|>];
   
   <|"Model" -> "Simple", "Parameters" -> p, "History" -> hist, "Snapshots" -> snapshots,
     "FinalState" -> state|>
];

RunSimpleScenario[name_String, overrides_: <||>] := RunSimpleCellularModel[Join[SimpleScenarioParameters[name], overrides]];

(* ============================================================ *)
(* 5. Model zaawansowany hybrydowy                              *)
(* ============================================================ *)

ClearAll[updateField, drugEfficacyStep, initializeImmuneCells, immuneField,
  moveImmuneCell, moveImmuneCells, RunAdvancedVirusModel, RunScenario];

updateField[field_, source_, diffusion_, decay_, dt_, dx_] := Module[{nSub, h, x},
   nSub = Max[1, Ceiling[4 diffusion dt/dx^2]];
   h = dt/nSub;
   x = field;
   Do[
      x = clipNonnegative[x + h (diffusion laplacian2D[x]/dx^2 - decay x + source)],
      {nSub}
   ];
   x
];

drugEfficacyStep[current_, t_, p_] := Module[{d = current},
   If[TrueQ[p["UseDrug"]] && t >= p["DrugStart"],
      If[p["DoseInterval"] > 0,
         If[Mod[t - p["DrugStart"], p["DoseInterval"]] == 0, d = Min[1.0, d + p["DrugDose"]]],
         If[t == p["DrugStart"], d = Min[1.0, d + p["DrugDose"]]]
      ]
   ];
   d (1 - p["DrugElimination"])
];

initializeImmuneCells[p_Association] := Module[{size, n},
   size = p["Size"];
   n = If[p["UseImmuneCells"], p["ImmuneInitialCount"], 0];
   safeRandomSample[Tuples[Range[size], 2], n]
];

immuneField[immuneCells_, size_Integer] := Module[{field = ConstantArray[0.0, {size, size}], pos},
   Do[field[[Sequence @@ pos]] += 1.0, {pos, immuneCells}];
   field
];

moveImmuneCell[pos_, chemokine_, p_] := Module[
   {size, neigh, values, best, movement, strength},
   size = p["Size"];
   movement = p["ImmuneMovement"];
   strength = p["ChemotaxisStrength"];
   neigh = boundedNeighbors[pos, size];
   If[Length[neigh] == 0, Return[pos]];
   If[movement === "Chemotaxis" && RandomReal[] < strength,
      values = chemokine[[Sequence @@ #]] & /@ neigh;
      best = Pick[neigh, values, Max[values]];
      RandomChoice[best],
      RandomChoice[neigh]
   ]
];

moveImmuneCells[cells_, chemokine_, p_] := moveImmuneCell[#, chemokine, p] & /@ cells;

RunAdvancedVirusModel[overrides_: <||>] := Module[
   {p, size, steps, snapshotEvery, dt, dx, state, age, rna, virus, chem, antibody,
    positions, pos, immuneCells, drug, t, i, j, newState, newAge, newRNA,
    infMask, infectiousMask, mutantMask, deadMask, localInf, localImmune, sourceVirus,
    sourceChem, sourceAntibody, effectiveClearance, pFree, pContact, pTotal,
    mutEscape, drugProductionFactor, drugInfectionFactor, hist, snapshots, 
    recruit, neigh, killTargets, target},
   
   p = MergeParameters[DefaultAdvancedParameters[], overrides];
   SeedRandom[p["Seed"]];
   size = p["Size"];
   steps = p["Steps"];
   snapshotEvery = p["SnapshotEvery"];
   dt = p["dt"];
   dx = p["dx"];
   
   state = ConstantArray[healthy, {size, size}];
   age = ConstantArray[0, {size, size}];
   rna = ConstantArray[0.0, {size, size}];
   virus = ConstantArray[0.0, {size, size}];
   chem = ConstantArray[0.0, {size, size}];
   antibody = ConstantArray[p["InitialAntibody"], {size, size}];
   drug = 0.0;
   
   If[p["VaccineCoverage"] > 0,
      Do[If[RandomReal[] < p["VaccineCoverage"], state[[i, j]] = protected], {i, size}, {j, size}]
   ];
   
   positions = initialInfectionPositions[p];
   Do[
      state[[Sequence @@ pos]] = latent;
      age[[Sequence @@ pos]] = 0;
      rna[[Sequence @@ pos]] = p["InitialRNA"];
      virus[[Sequence @@ pos]] = p["InitialVirus"],
      {pos, positions}
   ];
   
   immuneCells = initializeImmuneCells[p];
   hist = {};
   snapshots = {};
   
   Do[
      If[Mod[t - 1, snapshotEvery] == 0,
         AppendTo[snapshots, <|"Time" -> t - 1, "State" -> state, "Virus" -> virus,
            "Chemokine" -> chem, "Antibody" -> antibody, "RNA" -> rna,
            "ImmuneCells" -> immuneCells, "Drug" -> drug|>]
      ];
      
      drug = drugEfficacyStep[drug, t, p];
      drugProductionFactor = Clip[1 - p["DrugProductionEffect"] drug, {0, 1}];
      drugInfectionFactor = Clip[1 - p["DrugInfectionEffect"] drug, {0, 1}];
      
      (* Maski stanow *)
      infMask = maskMember[state, {latent, infectious, mutant}];
      infectiousMask = maskMember[state, {infectious, mutant}];
      mutantMask = maskEqual[state, mutant];
      deadMask = maskMember[state, {dead, deadDrug}];
      localInf = mooreCount[infectiousMask];
      localImmune = mooreCount[immuneField[immuneCells, size]];
      
      (* Wewnetrzna replikacja RNA: wzrost logistyczny w zakazonych komorkach *)
      newRNA = rna;
      Do[
         If[MemberQ[{latent, infectious, mutant}, state[[i, j]]],
            newRNA[[i, j]] = rna[[i, j]] + dt p["RNAReplicationRate"] rna[[i, j]] (1 - rna[[i, j]]/p["RNACarryingCapacity"]);
            newRNA[[i, j]] = Min[p["RNACarryingCapacity"], Max[0.0, newRNA[[i, j]]]],
            newRNA[[i, j]] = 0.0
         ],
         {i, size}, {j, size}
      ];
      rna = newRNA;
      
      (* latent -> infectious po fazie lag albo po przekroczeniu progu RNA *)
      newState = state;
      newAge = age;
      Do[
         If[state[[i, j]] === latent,
            newAge[[i, j]] = age[[i, j]] + 1;
            If[newAge[[i, j]] >= p["LagTime"] || rna[[i, j]] >= p["RNAThreshold"],
               newState[[i, j]] = infectious
            ]
         ];
         If[MemberQ[{infectious, mutant}, state[[i, j]]], newAge[[i, j]] = age[[i, j]] + 1],
         {i, size}, {j, size}
      ];
      state = newState;
      age = newAge;
      
      infectiousMask = maskMember[state, {infectious, mutant}];
      mutantMask = maskEqual[state, mutant];
      infMask = maskMember[state, {latent, infectious, mutant}];
      deadMask = maskMember[state, {dead, deadDrug}];
      
      (* Pole chemotaktyczne C: zrodlo z zakazonych i martwych komorek *)
      sourceChem = p["ChemokineFromInfected"] infMask + p["ChemokineFromDead"] deadMask;
      chem = If[TrueQ[p["UseChemokine"]],
         updateField[chem, sourceChem, p["Dc"], p["ChemokineDecay"], dt, dx],
         ConstantArray[0.0, {size, size}]
      ];
      
      (* Ruch i dzialanie komorek odpornosciowych *)
      If[TrueQ[p["UseImmuneCells"]],
         If[t >= p["ImmuneActivationTime"] && p["ImmuneRecruitmentPerStep"] > 0,
            recruit = safeRandomSample[Tuples[Range[size], 2], p["ImmuneRecruitmentPerStep"]];
            immuneCells = Join[immuneCells, recruit]
         ];
         immuneCells = moveImmuneCells[immuneCells, chem, p];
         Do[
            neigh = Join[{pos}, boundedNeighbors[pos, size]];
            killTargets = Select[neigh, MemberQ[{latent, infectious, mutant}, state[[Sequence @@ #]]] &];
            If[Length[killTargets] > 0,
               target = RandomChoice[killTargets];
               If[RandomReal[] < p["ImmuneKillProbability"],
                  state[[Sequence @@ target]] = dead;
                  rna[[Sequence @@ target]] = 0.0
               ]
            ],
            {pos, immuneCells}
         ];
      ];
      localImmune = mooreCount[immuneField[immuneCells, size]];
      
      (* Przeciwciala A: zrodlo przez funkcje Hilla *)
      sourceAntibody = If[TrueQ[p["UseAntibodies"]],
         MapThread[
            p["AntibodyFromVirus"] hill[#1, p["AntibodyHillK"], p["AntibodyHillN"]] +
            p["AntibodyFromInfected"] hill[#2, 1.0, p["AntibodyHillN"]] &,
            {virus, infMask}, 2],
         ConstantArray[0.0, {size, size}]
      ];
      antibody = If[TrueQ[p["UseAntibodies"]],
         updateField[antibody, sourceAntibody, p["Da"], p["AntibodyDecay"], dt, dx],
         ConstantArray[0.0, {size, size}]
      ];
      
      (* Produkcja, dyfuzja i usuwanie wolnych wirionow V *)
      (* Mutant jest czesciowo odporny na przeciwciala. *)
      mutEscape = 1 - p["MutantAntibodyEscape"] mutantMask;
      effectiveClearance = p["VirusClearance"] +
         p["ImmuneVirusClearance"] localImmune +
         p["AntibodyVirusClearance"] antibody mutEscape;
      
      sourceVirus = p["ReleaseRate"] drugProductionFactor infectiousMask rna;
      virus = clipNonnegative[
         virus + dt (p["Dv"] laplacian2D[virus]/dx^2 - effectiveClearance virus + sourceVirus)
      ];
        
      (* Infekcja nowych komorek: cell-free + cell-to-cell *)
      localInf = mooreCount[infectiousMask];
      newState = state;
      newAge = age;
      Do[
         If[state[[i, j]] === healthy || state[[i, j]] === protected,
            pFree = drugInfectionFactor p["BetaFree"] virus[[i, j]]/(p["Vhalf"] + virus[[i, j]] + 10^-12);
            pContact = 1 - (1 - drugInfectionFactor p["BetaContact"])^localInf[[i, j]];
            pTotal = 1 - (1 - pFree) (1 - pContact);
            If[state[[i, j]] === protected, pTotal = (1 - p["VaccineProtection"]) pTotal];
            If[RandomReal[] < pTotal,
               newState[[i, j]] = latent;
               newAge[[i, j]] = 0;
               rna[[i, j]] = p["InitialRNA"];
            ]
         ],
         {i, size}, {j, size}
      ];
      state = newState;
      age = newAge;
      
      (* Mutacje i smierc komorek zakazonych *)
      Do[
			If[
		      MemberQ[{latent, infectious}, state[[i, j]]] &&
		      TrueQ[p["UseMutation"]] &&
		      RandomReal[] < p["MutationProb"],
		      
		      state[[i, j]] = mutant
		   ];
		   
		   (* Bezposrednie dzialanie leku na komorki zakazone.
		      Komorki zabite przez lek oznaczamy jako deadDrug, czyli fioletowe. *)
		   If[
		      MemberQ[{latent, infectious, mutant}, state[[i, j]]] &&
		      TrueQ[p["UseDrug"]] &&
		      p["DrugInfectedKilling"] > 0 &&
		      RandomReal[] < p["DrugInfectedKilling"] drug,
		      
		      state[[i, j]] = deadDrug;
		      rna[[i, j]] = 0.0
		   ];
		   
		   If[
		      MemberQ[{latent, infectious, mutant}, state[[i, j]]] &&
		      (age[[i, j]] >= p["DeathTime"] || RandomReal[] < p["DeathProbability"]),
		      
		      state[[i, j]] = dead;
		      rna[[i, j]] = 0.0
		   ];
		   
		   If[
		      state[[i, j]] === healthy &&
		      p["DrugToxicity"] > 0 &&
		      RandomReal[] < p["DrugToxicity"] drug,
		      
		      state[[i, j]] = deadDrug
		   ],
		   
		   {i, size}, {j, size}
		];
      
      AppendTo[hist, <|
         "Time" -> t,
         "Healthy" -> Count[Flatten[state], healthy],
         "Latent" -> Count[Flatten[state], latent],
         "Infectious" -> Count[Flatten[state], infectious],
         "Mutant" -> Count[Flatten[state], mutant],
         "Dead" -> Count[Flatten[state], dead],
         "Protected" -> Count[Flatten[state], protected],
         "TotalVirus" -> Total[Flatten[virus]],
         "MeanVirus" -> Mean[Flatten[virus]],
         "MeanChemokine" -> Mean[Flatten[chem]],
         "MeanAntibody" -> Mean[Flatten[antibody]],
         "MeanRNA" -> Mean[Flatten[rna]],
         "Drug" -> drug,
         "ImmuneCount" -> Length[immuneCells]
      |>],
      {t, 1, steps}
   ];
   
   AppendTo[snapshots, <|"Time" -> steps, "State" -> state, "Virus" -> virus,
      "Chemokine" -> chem, "Antibody" -> antibody, "RNA" -> rna,
      "ImmuneCells" -> immuneCells, "Drug" -> drug|>];
   
   <|"Model" -> "AdvancedHybrid", "Parameters" -> p, "History" -> hist,
     "Snapshots" -> snapshots, "FinalState" -> state, "FinalVirus" -> virus,
     "FinalChemokine" -> chem, "FinalAntibody" -> antibody|>
];

RunScenario[name_String, overrides_: <||>] := RunAdvancedVirusModel[Join[ScenarioParameters[name], overrides]];

(* ============================================================ *)
(* 6. Wizualizacje                                              *)
(* ============================================================ *)

ClearAll[virusPointsFromField, GridStatePlot, AnimateResult, CompareAnimations,
  CompareInfectedPlot, CompareVirusPlot, CompactSummaryPlot, ParameterTable,
  ViroscapeStyleFrame, ViroscapeAnimation, HillFunctionPlot, HillFunctionManipulate, ExperimentHillImpact];

Options[GridStatePlot] = {
   "ShowVirus" -> True,
   "ShowImmune" -> True,
   "ShowChemokineHint" -> False,
   "MaxVirions" -> 300,
   "Title" -> Automatic,
   ImageSize -> 420
};

virusPointsFromField[v_, maxVirions_Integer, size_Integer] := Module[
   {vmax, threshold, candidates, sample},
   vmax = Max[Flatten[N[v]]];
   If[vmax <= 0, Return[{}]];
   threshold = 0.12 vmax;
   candidates = Select[Tuples[Range[size], 2], v[[#[[1]], #[[2]]]] > threshold &];
   sample = safeRandomSample[candidates, maxVirions];
   ({#[[2]] - 0.5, size - #[[1]] + 0.5} + RandomReal[{-0.33, 0.33}, 2]) & /@ sample
];

GridStatePlot[result_Association, frameIndex_: -1, opts : OptionsPattern[]] := Module[
   {snap, state, virus, immuneCells, size, title, imageSize, rects, colorRules,
    vPoints, immunePoints, time, drug, subtitle, virusLayer, immuneLayer},
   
   snap = If[frameIndex === -1, Last[result["Snapshots"]], result["Snapshots"][[frameIndex]]];
   state = snap["State"];
   virus = Lookup[snap, "Virus", ConstantArray[0.0, Dimensions[state]]];
   immuneCells = Lookup[snap, "ImmuneCells", {}];
   size = Length[state];
   time = snap["Time"];
   drug = Lookup[snap, "Drug", 0.0];
   imageSize = OptionValue[ImageSize];
   title = OptionValue["Title"];
   If[title === Automatic, title = "t = " <> ToString[time]];
   subtitle = If[drug > 0, "   lek = " <> ToString[NumberForm[drug, {3, 2}]], ""];
   
   (* Jesli lek jest obecny, zdrowe tlo robi sie lekko niebieskie.
      To jest wizualizacja terapii, nie zmiana biologii modelu. *)
   colorRules = If[
      drug > 0,
      Join[
         stateColorRules,
         <|
            healthy -> Blend[{White, RGBColor[0.72, 0.92, 1.0]}, Min[1.0, 1.6 drug]],
            protected -> Blend[{LightGreen, RGBColor[0.55, 0.88, 0.90]}, Min[1.0, 1.6 drug]]
         |>
      ],
      stateColorRules
   ];
   
   rects = Flatten[
      Table[
         {EdgeForm[Directive[GrayLevel[0.82], Thin]], colorRules[state[[i, j]]],
          Rectangle[{j - 1, size - i}, {j, size - i + 1}]},
         {i, size}, {j, size}], 1];
   
   vPoints = If[TrueQ[OptionValue["ShowVirus"]], virusPointsFromField[virus, OptionValue["MaxVirions"], size], {}];
   immunePoints = If[TrueQ[OptionValue["ShowImmune"]], ({#[[2]] - 0.5, size - #[[1]] + 0.5} & /@ immuneCells), {}];
   virusLayer = If[Length[vPoints] > 0, {Purple, PointSize[0.0065], Point[vPoints]}, {}];
   immuneLayer = If[Length[immunePoints] > 0, {Blue, PointSize[0.014], Point[immunePoints]}, {}];
   
   Graphics[
      {rects, virusLayer, immuneLayer},
      Frame -> True,
      FrameTicks -> None,
      PlotRange -> {{0, size}, {0, size}},
      ImageSize -> imageSize,
      PlotLabel -> Style[title <> subtitle, 13, Bold]
   ]
];

AnimateResult[result_Association, opts : OptionsPattern[GridStatePlot]] := ListAnimate[
   Table[GridStatePlot[result, k, opts], {k, 1, Length[result["Snapshots"]]}],
   AnimationRunning -> False,
   DefaultDuration -> 14
];

CompareAnimations[result1_Association, result2_Association, title1_String, title2_String, opts : OptionsPattern[GridStatePlot]] := Module[{n},
   n = Min[Length[result1["Snapshots"]], Length[result2["Snapshots"]]];
   ListAnimate[
      Table[
         Grid[{{Style[title1, 14, Bold], Style[title2, 14, Bold]},
               {GridStatePlot[result1, k, opts, ImageSize -> 360], GridStatePlot[result2, k, opts, ImageSize -> 360]}},
            Spacings -> {3, 1}],
         {k, n}],
      AnimationRunning -> False,
      DefaultDuration -> 14
   ]
];

(* ------------------------------------------------------------ *)
(* Bezpieczne funkcje do czytania historii symulacji              *)
(* ------------------------------------------------------------ *)

ClearAll[toHistoryAssociation, hvalue, cleanSeries, infectedSeries, virusSeries,
  healthySeries, deadSeries, protectedSeries, antibodySeries, drugSeries,
  PlotCounts, CompareInfectedPlot, CompareVirusPlot, CompactSummaryPlot];

toHistoryAssociation[row_] := If[AssociationQ[row], row, Association[row]];

hvalue[row_, key_String, default_: 0] := Module[{a},
   a = toHistoryAssociation[row];
   Lookup[a, key, default]
];

cleanSeries[data_] := Select[data, MatchQ[#, {_?NumericQ, _?NumericQ}] &];

infectedSeries[result_Association] := cleanSeries[
   ({hvalue[#, "Time"], hvalue[#, "Latent"] + hvalue[#, "Infectious"] + hvalue[#, "Mutant"]} & /@ result["History"])
];

virusSeries[result_Association] := cleanSeries[
   ({hvalue[#, "Time"], hvalue[#, "TotalVirus"]} & /@ result["History"])
];

healthySeries[result_Association] := cleanSeries[
   ({hvalue[#, "Time"], hvalue[#, "Healthy"]} & /@ result["History"])
];

deadSeries[result_Association] := cleanSeries[
   ({hvalue[#, "Time"], hvalue[#, "Dead"]} & /@ result["History"])
];

protectedSeries[result_Association] := cleanSeries[
   ({hvalue[#, "Time"], hvalue[#, "Protected"]} & /@ result["History"])
];

antibodySeries[result_Association] := cleanSeries[
   ({hvalue[#, "Time"], hvalue[#, "MeanAntibody"]} & /@ result["History"])
];

drugSeries[result_Association] := cleanSeries[
   ({hvalue[#, "Time"], hvalue[#, "Drug"]} & /@ result["History"])
];

PlotCounts[result_Association] := Module[{series, labels},
   series = {
      healthySeries[result],
      infectedSeries[result],
      deadSeries[result],
      protectedSeries[result]
   };
   labels = {"zdrowe", "zakazone", "martwe", "odporne/szczepione"};
   ListLinePlot[
      series,
      PlotLegends -> labels,
      Frame -> True,
      FrameLabel -> {"czas", "liczba komorek"},
      PlotLabel -> "Dynamika liczby komorek",
      PlotRange -> All,
      ImageSize -> Large,
      PlotTheme -> "Detailed"
   ]
];

CompareInfectedPlot[results_Association] := Module[{names, lines},
   names = Keys[results];
   lines = infectedSeries /@ Values[results];
   ListLinePlot[
      lines,
      PlotLegends -> names,
      Frame -> True,
      FrameLabel -> {"czas", "liczba zakazonych komorek"},
      PlotLabel -> "Porownanie liczby zakazonych komorek",
      PlotRange -> All,
      ImageSize -> Large,
      PlotTheme -> "Detailed"
   ]
];

CompareVirusPlot[results_Association] := Module[{names, lines},
   names = Keys[results];
   lines = virusSeries /@ Values[results];
   ListLinePlot[
      lines,
      PlotLegends -> names,
      Frame -> True,
      FrameLabel -> {"czas", "suma V(t,x,y)"},
      PlotLabel -> "Calkowita ilosc wolnych wirionow w modelu",
      PlotRange -> All,
      ImageSize -> Large,
      PlotTheme -> "Detailed"
   ]
];

CompactSummaryPlot[result_Association] := Module[
   {infected, deadLine, virusLine, antibodyLine, drugLine},
   infected = infectedSeries[result];
   deadLine = deadSeries[result];
   virusLine = virusSeries[result];
   antibodyLine = antibodySeries[result];
   drugLine = drugSeries[result];
   ListLinePlot[
      {normalizePairs[infected], normalizePairs[deadLine], normalizePairs[virusLine],
       normalizePairs[antibodyLine], drugLine},
      PlotLegends -> {"zakazone", "martwe", "wolne wiriony V", "przeciwciala A", "lek"},
      Frame -> True,
      FrameLabel -> {"czas", "wartosc znormalizowana"},
      PlotLabel -> "Skrocony wykres dynamiki modelu",
      PlotRange -> All,
      ImageSize -> Large,
      PlotTheme -> "Detailed"
   ]
];

ViroscapeStyleFrame[result_Association, frameIndex_: -1, maxVirions_: 450] := GridStatePlot[
   result, frameIndex,
   "ShowVirus" -> True,
   "ShowImmune" -> True,
   "MaxVirions" -> maxVirions,
   "Title" -> "Widok komorki + wolne wiriony, inspiracja Figure 2a",
   ImageSize -> 520
];

ViroscapeAnimation[result_Association] := AnimateResult[
   result,
   "ShowVirus" -> True,
   "ShowImmune" -> True,
   "MaxVirions" -> 450,
   ImageSize -> 520
];

HillFunctionPlot[k_: 0.05, n_: 2] := Plot[
   x^n/(k^n + x^n), {x, 0, 0.25},
   Frame -> True,
   FrameLabel -> {"lokalne stezenie wirusa V", "znormalizowana odpowiedz"},
   PlotLabel -> "Funkcja Hilla: nieliniowa i nasycajaca odpowiedz immunologiczna",
   PlotRange -> {0, 1.05},
   ImageSize -> Large,
   PlotTheme -> "Detailed"
];

HillFunctionManipulate[] := Manipulate[
   Plot[hill[x, k, n], {x, 0, 0.25},
      PlotRange -> {{0, 0.25}, {0, 1.05}},
      Frame -> True,
      FrameLabel -> {"Lokalne st\:0119\:017cenie wirusa V", "Znormalizowana odpowied\:017a immunologiczna"},
      PlotLabel -> Style["Interaktywna charakterystyka funkcji Hilla w modelu", 12, Bold],
      PlotTheme -> "Detailed",
      ImageSize -> Large,
      Epilog -> {
         {Dashed, Gray, Line[{{k, 0}, {k, 0.5}}], Line[{{0, 0.5}, {k, 0.5}}]},
         {Red, PointSize[Large], Point[{k, 0.5}]},
         Text[Style[Row[{"Punkt zwrotny (k) = ", k}], 10, Darker[Gray], Background -> White], {k + 0.03, 0.45}]
      }
   ],
   {{k, 0.02, "Pr\[OAcute]g aktywacji (k)"}, 0.005, 0.1, 0.005, Appearance -> "Labeled"},
   {{n, 2, "Stromo\:015b\[CAcute] odpowiedzi (n)"}, 1, 5, 0.5, Appearance -> "Labeled"},
   TrackedSymbols :> {k, n}
];

ParameterTable[] := Grid[
   Prepend[
      {
       {"BetaFree", "zakazenie przez wolne wiriony V"},
       {"BetaContact", "zakazenie przez kontakt z zakazonymi sasiadami"},
       {"Dv", "wspolczynnik dyfuzji wolnych wirionow"},
       {"LagTime", "opoznienie przed produkcja wirionow"},
       {"RNAReplicationRate", "tempo wzrostu wewnatrzkomorkowego RNA"},
       {"ChemokineFromInfected", "sila sygnalu chemotaktycznego od zakazonych"},
       {"AntibodyHillK", "prog polowy funkcji Hilla"},
       {"AntibodyHillN", "stromosc funkcji Hilla"},
       {"DrugDose", "skutecznosc dawki leku"},
       {"VaccineCoverage", "udzial chronionych komorek"},
       {"MutationProb", "prawdopodobienstwo mutacji ucieczkowej"}
      },
      {Style["Parametr", Bold], Style["Znaczenie", Bold]}
   ],
   Frame -> All,
   Alignment -> Left,
   Background -> {None, {LightYellow, White}}
];

(* ============================================================ *)
(* 7. Gotowe funkcje eksperymentow                              *)
(* ============================================================ *)

ClearAll[ExperimentSimpleVsAdvanced, ExperimentVirusTypes, ExperimentTransmissionModes,
  ExperimentDiffusion, ExperimentLagPhase, ExperimentChemotaxis,
  ExperimentImmuneProfiles, ExperimentMemory, ExperimentDrugVaccination,
  ExperimentDosing, ExperimentMutation, ExperimentHill];

ExperimentSimpleVsAdvanced[] := Module[{simple, adv, wspolnaBeta},
	wspolnaBeta = 0.10;
   simple = RunSimpleScenario["SimpleBaseline", <|"BetaContact" -> wspolnaBeta, "InitialInfected" -> 1, "Seed" -> 201|>];
   adv = RunScenario["BaselineAdvanced", <|"BetaContact" -> wspolnaBeta,"BetaFree" -> 0.025,"Dv" -> 0.12,"InitialInfected" -> 1, "UseImmuneCells" -> False, "UseAntibodies" -> False, "Seed" -> 201|>];
   <|
      "Simple" -> simple,
      "Advanced" -> adv,
      "Animation" -> CompareAnimations[simple, adv, "Model prosty", "Model zaawansowany", "ShowVirus" -> True, "ShowImmune" -> False],
      "InfectedPlot" -> CompareInfectedPlot[<|"prosty" -> simple, "zaawansowany" -> adv|>],
      "VirusPlot" -> CompareVirusPlot[<|"zaawansowany" -> adv|>]
   |>
];

ExperimentVirusTypes[] := Module[{sars, hcv, fast, slow},
   sars = RunScenario["SARSCoV2Like"];
   hcv = RunScenario["HCVLike"];
   covid = RunSimpleScenario["SimpleFastProfile"];
   hanta = RunSimpleScenario["SimpleSlowDamagingProfile"];
   <|
      "SARSCov2Like" -> sars,
      "HCVLike" -> hcv,
      "AnimationSARSvsHCV" -> CompareAnimations[sars, hcv, "SARS-Cov-2-like", "HCV-like", "ShowVirus" -> True, "ShowImmune" -> False],
      "PlotSARSvsHCV" -> CompareInfectedPlot[<|"SARS-Cov-2-like" -> sars, "HCV-like" -> hcv|>],
      "SarsCov2Like" -> covid,
      "HantaLike" -> hanta,
      "AnimationPhenotypes" -> CompareAnimations[covid, hanta, "SARS-Cov-2-like", "Hantavirus-like", "ShowVirus" -> False, "ShowImmune" -> False],
      "PlotPhenotypes" -> CompareInfectedPlot[<|"SARS-Cov-2-like" -> covid, "Hantavirus-like" -> hanta|>]
   |>
];

ExperimentTransmissionModes[] := Module[{free, contact, both},
   free = RunScenario["CellFreeDominant"];
   contact = RunScenario["CellToCellDominant"];
   both = RunScenario["BothTransmissionModes"];
   <|
      "CellFree" -> free,
      "CellToCell" -> contact,
      "Both" -> both,
      "Animation" -> CompareAnimations[free, contact, "dominuje cell-free", "dominuje cell-to-cell", "ShowVirus" -> True, "ShowImmune" -> False],
      "InfectedPlot" -> CompareInfectedPlot[<|"cell-free" -> free, "cell-to-cell" -> contact, "oba mechanizmy" -> both|>],
      "VirusPlot" -> CompareVirusPlot[<|"cell-free" -> free, "cell-to-cell" -> contact, "oba" -> both|>]
   |>
];

ExperimentDiffusion[] := Module[{low, med, high},
   low = RunScenario["LowDiffusion"];
   med = RunScenario["MediumDiffusion"];
   high = RunScenario["HighDiffusion"];
   <|
      "Low" -> low,
      "Medium" -> med,
      "High" -> high,
      "AnimationLowHigh" -> CompareAnimations[low, high, "niska dyfuzja", "wysoka dyfuzja", "ShowVirus" -> True, "ShowImmune" -> False],
      "VirusPlot" -> CompareVirusPlot[<|"niska Dv" -> low, "srednia Dv" -> med, "wysoka Dv" -> high|>]
   |>
];

ExperimentLagPhase[] := Module[{short, long},
   short = RunScenario["ShortLag"];
   long = RunScenario["LongLag"];
   <|
      "ShortLag" -> short,
      "LongLag" -> long,
      "Animation" -> CompareAnimations[short, long, "krotka faza lag", "dluga faza lag", "ShowVirus" -> True, "ShowImmune" -> False],
      "InfectedPlot" -> CompareInfectedPlot[<|"krotki lag" -> short, "dlugi lag" -> long|>],
      "VirusPlot" -> CompareVirusPlot[<|"krotki lag" -> short, "dlugi lag" -> long|>]
   |>
];

ExperimentChemotaxis[] := Module[{random, chemo},
   random = RunScenario["RandomImmune"];
   chemo = RunScenario["ChemotaxisImmune"];
   <|
      "Random" -> random,
      "Chemotaxis" -> chemo,
      "Animation" -> CompareAnimations[random, chemo, "ruch losowy", "chemotaksja po polu C", "ShowVirus" -> True, "ShowImmune" -> True],
      "Plot" -> CompareInfectedPlot[<|"ruch losowy" -> random, "chemotaksja" -> chemo|>]
   |>
];

ExperimentImmuneProfiles[] := Module[{fast, slow},
   fast = RunScenario["FastImmuneProfile"];
   slow = RunScenario["SlowImmuneProfile"];
   <|
      "Fast" -> fast,
      "Slow" -> slow,
      "Animation" -> CompareAnimations[fast, slow, "mlodsza osoba", "starsza osoba", "ShowVirus" -> True, "ShowImmune" -> True],
      "Plot" -> CompareInfectedPlot[<|"szybka/silna odpornosc" -> fast, "opozniona/slabsza" -> slow|>]
   |>
];

ExperimentMemory[] := Module[{first, second},
   first = RunScenario["FirstExposure"];
   second = RunScenario["SecondExposureMemory"];
   <|
      "First" -> first,
      "Second" -> second,
      "Animation" -> CompareAnimations[first, second, "pierwszy kontakt", "wtorny kontakt, pamiec", "ShowVirus" -> True, "ShowImmune" -> True],
      "Plot" -> CompareInfectedPlot[<|"pierwszy kontakt" -> first, "wtorny kontakt" -> second|>]
   |>
];

ExperimentDrugVaccination[] := Module[{noDrug, early, noVax, vax},
   noDrug = RunScenario["NoDrug"];
   early = RunScenario["EarlyDrug"];
   noVax = RunScenario["NoVaccination"];
   vax = RunScenario["Vaccination"];
   <|
      "NoDrug" -> noDrug,
      "EarlyDrug" -> early,
      "AnimationDrug" -> CompareAnimations[noDrug, early, "bez leku", "wczesny lek", "ShowVirus" -> True, "ShowImmune" -> False],
      "PlotDrug" -> CompareInfectedPlot[<|"bez leku" -> noDrug, "wczesny lek" -> early|>],
      "NoVaccination" -> noVax,
      "Vaccination" -> vax,
      "AnimationVaccination" -> CompareAnimations[noVax, vax, "brak szczepienia", "wyszczepienie", "ShowVirus" -> True, "ShowImmune" -> False],
      "PlotVaccination" -> CompareInfectedPlot[<|"brak szczepienia" -> noVax, "wyszczepienie" -> vax|>]
   |>
];

ExperimentDosing[] := Module[{early, late, single, cyclic},
   early = RunScenario["EarlyDrug"];
   late = RunScenario["LateDrug"];
   single = RunScenario["SingleDose"];
   cyclic = RunScenario["CyclicDose"];
   <|
      "Early" -> early,
      "Late" -> late,
      "AnimationEarlyLate" -> CompareAnimations[early, late, "wczesna interwencja", "spozniona interwencja", "ShowVirus" -> True, "ShowImmune" -> False],
      "PlotEarlyLate" -> CompareInfectedPlot[<|"wczesna" -> early, "spozniona" -> late|>],
      "Single" -> single,
      "Cyclic" -> cyclic,
      "AnimationSingleCyclic" -> CompareAnimations[single, cyclic, "dawka jednorazowa", "dawkowanie cykliczne", "ShowVirus" -> True, "ShowImmune" -> False],
      "PlotSingleCyclic" -> CompareInfectedPlot[<|"jednorazowa" -> single, "cykliczna" -> cyclic|>]
   |>
];

ExperimentMutation[] := Module[{stable, mut},
   stable = RunScenario["StableVirusWithVaccine"];
   mut = RunScenario["MutationEscape"];
   <|
      "Stable" -> stable,
      "Mutation" -> mut,
      "Animation" -> CompareAnimations[stable, mut, "wirus stabilny", "mutacja ucieczkowa", "ShowVirus" -> True, "ShowImmune" -> False],
      "Plot" -> CompareInfectedPlot[<|"stabilny" -> stable, "mutujacy" -> mut|>]
   |>
];

(* ============================================================ *)
(* 8. Therapeutic window heatmap                                *)
(* ============================================================ *)

ClearAll[TestTherapyAdvanced, ExperimentTherapyHeatmap];

TestTherapyAdvanced[dTime_, eff_] := Module[
   {numRepeats, results, meanHealthy, res},
   numRepeats = 10;
   results = Table[
      res = RunAdvancedVirusModel[<|
            "BetaContact" -> 0.1,
            "BetaFree" -> 0.01,
            "DeathTime" -> 8,
            "LagTime" -> 3,
            "UseDrug" -> True,
            "DrugStart" -> dTime,
            "DrugDose" -> eff,
            "DrugElimination" -> 0.01,
            "DrugInfectionEffect" -> eff,
            "UseImmuneCells" -> False,
            "UseAntibodies" -> False,
            "UseChemokine" -> False,
            "UseMutation" -> False,
            "Steps" -> 100,
            "Size" -> 60,
            "Seed" -> RandomInteger[{1, 10000}]
          |>];
      N[res["History"][[-1]]["Healthy"]],
      {numRepeats}
   ];
   meanHealthy = Mean[results];
   Return[N[(meanHealthy/3600)*100]]
];

ExperimentTherapyHeatmap[] := Module[
   {data, plot},
   data = Flatten[
      ParallelTable[  
         {time, efficacy, TestTherapyAdvanced[time, efficacy]},
         {time, Range[5, 45, 10]},
         {efficacy, Range[0.2, 1.0, 0.2]}
      ], 1
   ];
   plot = ListDensityPlot[data, 
      ColorFunction -> "ThermometerColors",
      PlotLegends -> BarLegend[Automatic, LegendLabel -> "% Zdrowej Tkanki"],
      FrameLabel -> {"Czas podania leku (krok)", "Skuteczno\:015b\[CAcute] leku"},
      PlotLabel -> "Statystyczne okno terapeutyczne (\:015arednia z 10 pr\[OAcute]b)",
      InterpolationOrder -> 1, 
      PlotTheme -> "Detailed", 
      PlotRange -> All,
      PlotRangePadding -> None, 
      Method -> {"InterpolationPoints" -> 100},
      ImageSize -> Large
   ];
   <|"Data" -> data, "Plot" -> plot|>
];




(* ============================================================ *)
(* Hill Function Impact *)
(* ============================================================ *)
ExperimentHill[] := <|"Plot" -> HillFunctionPlot[]|>;

ExperimentHillImpact[] := Module[{baseNoCells, baseWithCells, lowK, highK, lowKPhase, highKPhase, phasePlot, anim},
   
   baseNoCells = Join[ScenarioParameters["FastImmuneProfile"], <|"UseImmuneCells" -> False|>];
   
   lowK = RunAdvancedVirusModel[Join[baseNoCells, <|"AntibodyHillK" -> 2.0, "Seed" -> 100|>]];
   highK = RunAdvancedVirusModel[Join[baseNoCells, <|"AntibodyHillK" -> 45.0, "Seed" -> 100|>]];
   
   baseWithCells = Join[ScenarioParameters["FastImmuneProfile"], <|"UseImmuneCells" -> True|>];
   
   lowKPhase = RunAdvancedVirusModel[Join[baseWithCells, <|"AntibodyHillK" -> 2.0, "Seed" -> 100|>]];
   highKPhase = RunAdvancedVirusModel[Join[baseWithCells, <|"AntibodyHillK" -> 45.0, "Seed" -> 100|>]];

   (* Wykres fazowy zasilany osobnymi danymi (lowKPhase i highKPhase) *)
   phasePlot = ListLinePlot[{
      Table[{krok["MeanVirus"], krok["MeanAntibody"]}, {krok, lowKPhase["History"]}],
      Table[{krok["MeanVirus"], krok["MeanAntibody"]}, {krok, highKPhase["History"]}]
   },
      PlotLegends -> {"Niski prog (k=2.0)", "Wysoki prog (k=45.0)"},
      Frame -> True,
      FrameLabel -> {"Srednia ilosc wirusa w tkance (V)", "Sredni poziom przeciwcial (A)"},
      PlotLabel -> "Dynamika odpornosci - dane z symulacji",
      PlotTheme -> "Detailed",
      PlotStyle -> Thickness[0.003],
      ImageSize -> 600
   ];

   anim = CompareAnimations[
      lowK, 
      highK, 
      "Szybka reakcja (Niski prog k = 2.0)", 
      "Spozniona reakcja (Wysoki prog k = 45.0)", 
      "ShowVirus" -> True, 
      "ShowImmune" -> True
   ];
   
   <|
      "LowKData" -> lowK,
      "HighKData" -> highK,
      "InfectedPlot" -> CompareInfectedPlot[<|
         "Szybka odpowiedz (k = 2.0)" -> lowK,
         "Spozniona odpowiedz (k = 45.0)" -> highK
      |>],
      "VirusPlot" -> CompareVirusPlot[<|
         "Szybka odpowiedz (k = 2.0)" -> lowK,
         "Spozniona odpowiedz (k = 45.0)" -> highK
      |>],
      "PhasePlot" -> phasePlot,
      "Animation" -> anim
   |>
];


(* ============================================================ *)
(*Bariera*)
(* ============================================================ *)
ClearAll[bbbZones, inBarrier, inBrain, RunBBBModel,
  GridBBBPlot, AnimateBBBResult, CompareBBBAnimations,
  CompareBBBInfectedPlot, ExperimentBBB];

(* Wspolne agresywne parametry wirusa \[LongDash] infekcja musi dotrzec do bariery *)
bbbBase[] := Join[DefaultAdvancedParameters[], <|
   "BloodFraction" -> 0.20, "BarrierWidth" -> 2,
   "BetaFree" -> 0.055, "BetaContact" -> 0.110,
   "ReleaseRate" -> 3.5, "Dv" -> 0.18,
   "LagTime" -> 3, "InitialRNA" -> 0.75,
   "RNAThreshold" -> 0.20, "RNAReplicationRate" -> 0.50,
   "DeathTime" -> 35, "DeathProbability" -> 0.001,
   "BrainBetaMultiplier" -> 1.6, "BrainDeathTime" -> 55,
   "BrainRNARate" -> 0.28,
   "BBBPermeability" -> 0.0, "BBBCellCrossing" -> 0.0,
   "Steps" -> 120, "Seed" -> 501
|>];

bbbZones[p_] := Module[{jB = Floor[p["BloodFraction"] p["Size"]]},
   <|"jBlood" -> jB, "jBarrier" -> jB + p["BarrierWidth"],
     "jBrain" -> jB + p["BarrierWidth"] + 1|>];

inBarrier[j_, z_] := j > z["jBlood"] && j <= z["jBarrier"];
inBrain[j_, z_]   := j > z["jBarrier"];

RunBBBModel[overrides_: <||>] := Module[
   {p, size, steps, se, dt, dx, state, age, rna, virus, chem, antibody,
    drug, zones, t, i, j, newState, newAge, newRNA,
    infMask, infectiousMask, deadMask, localInf, localImmune,
    effectiveClearance, pFree, pContact, pTotal,
    dpf, dif, hist, snapshots, beta, rnaRate},

   p    = Join[bbbBase[], overrides];
   SeedRandom[p["Seed"]];
   size = p["Size"]; steps = p["Steps"]; se = p["SnapshotEvery"];
   dt   = p["dt"];   dx    = p["dx"];    zones = bbbZones[p];

   state    = ConstantArray[healthy, {size, size}];
   age      = ConstantArray[0,       {size, size}];
   rna      = ConstantArray[0.0,     {size, size}];
   virus    = ConstantArray[0.0,     {size, size}];
   chem     = ConstantArray[0.0,     {size, size}];
   antibody = ConstantArray[0.0,     {size, size}];
   drug     = 0.0;

   (* Bariera = protected *)
   Do[If[inBarrier[j,zones], state[[i,j]] = protected],
      {i,size},{j,zones["jBlood"]+1,zones["jBarrier"]}];

   (* Ognisko: lewa 1/3 strefy krwi, srodek pionowy *)
   Module[{jS = Max[1, Floor[zones["jBlood"]*0.3]], iS = Ceiling[size/2]},
      state[[iS,jS]] = latent; rna[[iS,jS]] = p["InitialRNA"];
      virus[[iS,jS]] = p["InitialVirus"]];

   hist = {}; snapshots = {};

   Do[
      If[Mod[t-1,se]==0,
         AppendTo[snapshots, <|"Time"->t-1,"State"->state,"Virus"->virus,
            "ImmuneCells"->{},"Drug"->drug,"Zones"->zones|>]];

      drug = drugEfficacyStep[drug, t, p];
      dpf  = Clip[1 - p["DrugProductionEffect"] drug, {0,1}];
      dif  = Clip[1 - p["DrugInfectionEffect"]  drug, {0,1}];

      infMask        = maskMember[state, {latent,infectious,mutant}];
      infectiousMask = maskMember[state, {infectious,mutant}];
      deadMask       = maskMember[state, {dead,deadDrug}];
      localInf       = mooreCount[infectiousMask];
      localImmune    = ConstantArray[0.0, {size,size}];

      (* RNA *)
      newRNA = rna;
      Do[If[MemberQ[{latent,infectious,mutant}, state[[i,j]]],
            rnaRate = If[inBrain[j,zones], p["BrainRNARate"], p["RNAReplicationRate"]];
            newRNA[[i,j]] = Clip[rna[[i,j]] + dt rnaRate rna[[i,j]]
               (1-rna[[i,j]]/p["RNACarryingCapacity"]), {0.,p["RNACarryingCapacity"]}],
            newRNA[[i,j]] = 0.],
         {i,size},{j,size}];
      rna = newRNA;

      (* latent -> infectious *)
      newState = state; newAge = age;
      Do[If[state[[i,j]]===latent,
            newAge[[i,j]] = age[[i,j]]+1;
            If[newAge[[i,j]]>=p["LagTime"] || rna[[i,j]]>=p["RNAThreshold"],
               newState[[i,j]] = infectious]];
         If[MemberQ[{infectious,mutant},state[[i,j]]], newAge[[i,j]]=age[[i,j]]+1],
         {i,size},{j,size}];
      state = newState; age = newAge;
      infectiousMask = maskMember[state,{infectious,mutant}];
      infMask        = maskMember[state,{latent,infectious,mutant}];

      (* Wiriony: dyfuzja + filtr bariery *)
      effectiveClearance = p["VirusClearance"] + p["AntibodyVirusClearance"] antibody;
      virus = clipNonnegative[virus + dt(
         p["Dv"] laplacian2D[virus]/dx^2
         - effectiveClearance virus
         + p["ReleaseRate"] dpf infectiousMask rna)];
      Do[virus[[All,j]] *= p["BBBPermeability"],
         {j, zones["jBlood"]+1, zones["jBarrier"]}];

      (* Nowe infekcje *)
      newState = state; newAge = age;
      Do[If[(state[[i,j]]===healthy) && !inBarrier[j,zones],
            beta     = If[inBrain[j,zones], p["BrainBetaMultiplier"]*p["BetaContact"], p["BetaContact"]];
            pFree    = dif p["BetaFree"] virus[[i,j]]/(p["Vhalf"]+virus[[i,j]]+10^-12);
            pContact = 1-(1-dif beta)^localInf[[i,j]];
            pTotal   = 1-(1-pFree)(1-pContact);
            If[RandomReal[]<pTotal,
               newState[[i,j]]=latent; newAge[[i,j]]=0; rna[[i,j]]=p["InitialRNA"]]],
         {i,size},{j,size}];
      state = newState; age = newAge;

      (* Smierc + lek *)
      Do[If[MemberQ[{latent,infectious,mutant},state[[i,j]]],
            If[TrueQ[p["UseDrug"]] && p["DrugInfectedKilling"]>0 &&
               RandomReal[]<p["DrugInfectedKilling"]*drug,
               state[[i,j]]=deadDrug; rna[[i,j]]=0.];
            Module[{dt2=If[inBrain[j,zones],p["BrainDeathTime"],p["DeathTime"]]},
               If[age[[i,j]]>=dt2 || RandomReal[]<p["DeathProbability"],
                  state[[i,j]]=dead; rna[[i,j]]=0.]]],
         {i,size},{j,size}];

      AppendTo[hist, <|"Time"->t,
         "Healthy"       -> Count[Flatten@state, healthy],
         "Infectious"    -> Count[Flatten@state, infectious|latent|mutant],
         "Dead"          -> Count[Flatten@state, dead|deadDrug],
         "TotalVirus"    -> Total@Flatten@virus,
         "Drug"          -> drug,
         "InfectedBrain" -> Count[Flatten@Table[state[[i,j]],
            {i,size},{j,zones["jBrain"],size}], latent|infectious|mutant]|>],
      {t,1,steps}];

   AppendTo[snapshots,<|"Time"->steps,"State"->state,"Virus"->virus,
      "ImmuneCells"->{},"Drug"->drug,"Zones"->zones|>];
   <|"History"->hist,"Snapshots"->snapshots,"Parameters"->p,"Zones"->zones|>
];

(* ---- Wizualizacja ---- *)

GridBBBPlot[result_, k_:-1, imageSize_:480] := Module[
   {snap,state,virus,zones,size,drug,title,
    rects,cr,jB,jR,vPts,barrierOverlay,brainFrame,labels},

   snap  = If[k===-1, Last@result["Snapshots"], result["Snapshots"][[k]]];
   state = snap["State"]; virus = snap["Virus"];
   zones = snap["Zones"]; size  = Length@state;
   drug  = snap["Drug"];  title = "t = "<>ToString@snap["Time"];
   jB = zones["jBlood"]; jR = zones["jBarrier"];

   cr = If[drug>0,
      Join[stateColorRules,<|healthy->Blend[{White,RGBColor[.72,.92,1.]},Min[1.,1.6drug]]|>],
      stateColorRules];

   rects = Flatten[Table[
      {EdgeForm[Directive[GrayLevel[.82],Thin]], cr[state[[i,j]]],
       Rectangle[{j-1,size-i},{j,size-i+1}]},
      {i,size},{j,size}],1];

   vPts = virusPointsFromField[virus, 400, size];

   barrierOverlay = {FaceForm[Directive[RGBColor[1.,.75,.15],Opacity[.6]]],
      EdgeForm[Directive[RGBColor[.85,.55,0.],Opacity[.9],Thick]],
      Rectangle[{jB,0},{jR,size}]};

   brainFrame = {EdgeForm[Directive[RGBColor[.2,.45,.8],Thick,Opacity[.7]]],
      FaceForm[None], Rectangle[{jR,0},{size,size}]};

   labels = {
      Style[Text["Krew", {jB/2,          size+1.3}], 11,Bold,RGBColor[.6,0,0]],
      Style[Text["BBB",  {(jB+jR)/2,     size+1.3}],  9,Bold,RGBColor[.75,.4,0]],
      Style[Text["Mozg", {(jR+size)/2,   size+1.3}], 11,Bold,RGBColor[.15,.35,.7]]};

   Graphics[{rects, barrierOverlay, brainFrame,
      If[Length@vPts>0,{Purple,PointSize[.006],Point@vPts},{}], labels},
      Frame->True, FrameTicks->None,
      PlotRange->{{0,size},{0,size+2}},
      ImageSize->imageSize, PlotLabel->Style[title,13,Bold]]
];

AnimateBBBResult[result_, imageSize_:520] :=
   ListAnimate[Table[GridBBBPlot[result,k,imageSize],
      {k,1,Length@result["Snapshots"]}],
      AnimationRunning->False, DefaultDuration->14];

CompareBBBAnimations[r1_,r2_,t1_,t2_,imageSize_:400] :=
   ListAnimate[Table[Grid[{
      {Style[t1,14,Bold],      Style[t2,14,Bold]},
      {GridBBBPlot[r1,k,imageSize], GridBBBPlot[r2,k,imageSize]}},
      Spacings->{3,1}],
      {k,Min[Length@r1["Snapshots"],Length@r2["Snapshots"]]}],
      AnimationRunning->False, DefaultDuration->14];

CompareBBBInfectedPlot[r1_,r2_,r3_,labels_] :=
   ListLinePlot[
      Flatten[{
         ({#["Time"],#["InfectedBrain"]}&/@#["History"])&/@{r1,r2,r3},
         ({#["Time"],#["Infectious"]}&/@#["History"])&/@{r1,r2,r3}
      }, 1],
      PlotLegends -> Join[#<>" [mozg]"&/@labels, labels],
      Frame->True,
      FrameLabel->{"czas","liczba zakazonych komorek"},
      PlotLabel->"Neuroinfekcja vs. calkowita infekcja",
      PlotStyle->{Dashed,Dashed,Dashed,Automatic,Automatic,Automatic},
      PlotRange->All, ImageSize->Large, PlotTheme->"Detailed"];

ExperimentBBB[] := Module[{intact,broken,drug},
   intact = RunBBBModel[<|"BBBPermeability"->0.0|>];
   broken = RunBBBModel[<|"BBBPermeability"->0.55|>];
   drug   = RunBBBModel[<|"BBBPermeability"->0.55,
      "UseDrug"->True,"DrugStart"->8,"DrugDose"->0.90,
      "DrugElimination"->0.007,"DoseInterval"->0,
      "DrugInfectionEffect"->0.80,"DrugInfectedKilling"->0.03,
      "DrugToxicity"->0.0002|>];
   <|"Intact"->intact,"Broken"->broken,"Drug"->drug|>
];


Print["FinalVirusModel.wl loaded. Use RunSimpleScenario, RunScenario or Experiment... functions."];

