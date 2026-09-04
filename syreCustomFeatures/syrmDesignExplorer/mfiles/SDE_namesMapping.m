% Copyright 2026
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [map2gui,gui2map] = SDE_namesMapping()


names = [
  "A"             , "Linear current density - maximum (A/m)";
  "AHWC"          , "Linear current density - ASC from base speed (A/m)";
  "Abar"          , "Flux barriers cross-section (mm2)";
  "Ach"           , "Linear current density - characteristic (A/m)";
  "Ad"            , "Linear current density - d-axis (A/m)";
  "Aq"            , "Linear current density - q-axis (A/m)";
  % "Aqirr"         , "(wip) q-axis demagnetization limit current density (A/m)"
  "Ar"            , "Shaft radius (mm)";
  "ArLim"         , "Maximum shaft radius (mm)";
  "Aslots"        , "Slots cross-section area (mm2)";
  "Bmin0"         , "Demagnetization @ max current - minimum flux density (T)";
  "BminHWC"       , "Demagnetization @ base speed - minimum flux density (T)";
  "BminUGO"       , "Demagnetization @ UGO speed - minimum flux densisy (T)";
  "BrAvg"         , "Average PM flux density (T)";
  "BrMax"         , "Maximum PM flux density (T)";
  "BrMin"         , "Minimum PM flux density (T)";
  "F0_Ns"         , "Normalized maximum flux linkage (Vs/turns)";
  "J"             , "Current density (Arms/mm2)";
  "Lbase"         , "Inductance - base unsaturated (H)";
  "Lcqpu"         , "Inductance - circulating component (pu)";
  "Lfqpu"         , "Inductance - flow-through component (pu)";
  "Lmd"           , "Inductance - magnetizing d-axis (H)";
  "Lrpu"          , "Inductance - ribs effect (pu)";
  "Lspu"          , "Inductance - slot leakage (pu)";
  "MaxDef"        , "Structural - maximum airgap deformation (mm)";
  "MaxStress"     , "Structural - maximum stress (MPa)";
  "NsI0"          , "Normalized maximum current (A turns)";
  "NsIHWC"        , "Normalized ASC current from base speed (A turns)";
  "NsIUGO"        , "Normalized ASC current from UGO speed (A turns)";
  "NsIch"         , "Normalized characteristic current (A turns)";
  "PF"            , "Power Factor";
  "Pbase"         , "Base power (W)";
  "PrcRadStress"  , "Structural - 99th percentile radial ribs stress (MPa)";
  "PrcTanStress"  , "Structural - 99th percentile tangential ribs stress (MPa)";
  "RadRibStress"  , "Structural - maximum radial ribs stress (MPa)";
  "Rs"            , "Phase resistance (Ohm)";
  "T"             , "Torque (Nm)";
  "TanRibStress"  , "Structural - maximum tangential ribs stress (MPa)";
  "agclear"       , "Structural - airgap at maximum speed (mm)";
  "bRaw"          , "FEAfix - b points";
  "bb"            , "b";
  "dPM0"          , "Demagnetization @ max current - demagnetized PM volume (pu)";
  "dPMHWC"        , "Demagnetization @ base speed - demagnetized PM volume (pu)";
  "dPMUGO"        , "Demagnetization @ UGO speed - demagnetized PM volume (pu)";
  "dTempCu"       , "Steady-state temperature difference between copper and outer diameter (C)";
  "dg"            , "FEAfix - gamma correction factor (deg)";
  "dx"            , "geometry - dx parameter";
  "fM"            , "Flux linkage - PM component (Vs)";
  "fUGOpu"        , "UGO-base speed ratio";
  "fd"            , "Flux linkage - d-axis component (Vs)";
  "fq"            , "Flux linkage - q-axis component (Vs)";
  "gamma"         , "Current angle (deg)";
  "hc"            , "geometry - hc parameter";
  "hcMin"         , "minimum barrier thickness (mm)";
  "hc_pu"         , "geometry - hc_pu parameter";
  "i0"            , "Maximum current (A)";
  "iAmp"          , "Current amplitude (A)";
  "iHWC"          , "ASC current from base speed (A)";
  "iUGO"          , "ASC current from UGO speed (A)";
  "ich"           , "Characteristic current (A)";
  "id"            , "Current - d-axis component (A)";
  "iq"            , "Current - q-axis component (A)";
  "k0"            , "FEAfix - field-weakening correction factor (pu)";
  "kBmin0"        , "FEAfix - demagnetization at max current - Bmin correctio factor (pu)";
  "kBminHWC"      , "FEAfix - demagnetization at base speed - Bmin correctio factor (pu)";
  "kBminUGO"      , "FEAfix - demagnetization at UGO speed - Bmin correctio factor (pu)";
  "kMaxDef"       , "FEAfix - maximum deformation correction factor (pu)";
  "kMaxStress"    , "FEAfix - maximum stress correction factor (pu)";
  "kPrcRadStress" , "FEAfix - 99th percentile rad stress correction factor (pu)";
  "kPrcTanStress" , "FEAfix - 99th percentile tan stress correction factor (pu)";
  "kRadRibStress" , "FEAfix - maximum rad stress correction factor (pu)";
  "kTanRibStress" , "FEAfix - maximum tan stress correction factor (pu)";
  "kUGO"          , "PM by maximum flux linkage ratio (pu)";
  "kagclear"      , "FEAfix - airgap clearence correction factor (pu)";
  "kd"            , "FEAfix - d-axis inductance correction factor (pu)";
  "kdPM0"         , "FEAfix - demagnetization at max current - demag volume correction factor (pu)";
  "kdPMHWC"       , "FEAfix - demagnetization at base speed - demag volume correction factor (pu)";
  "kdPMUGO"       , "FEAfix - demagnetization at UGO speed - demag volume correction factor (pu)";
  "kdq"           , "wip - 1-Lmq/Lmd";
  "kg"            , "FEAfix - current angle correction factor (pu)";
  "kiHWC"         , "FEAfix - ASC current from base speed correction factor (pu)";
  "kiUGO"         , "FEAfix - ASC current from UGO speed correction factor (pu)";
  "kich"          , "FEAfix - characteristic current correction factor (pu)";
  "kj"            , "Thermal loading (W/m2)";
  "km"            , "FEAfix - PM flux linkage correction factor (pu)";
  "kmPM"          , "FEAfix - PM mass correction factor (pu)";
  "kq"            , "FEAfix - q-axis inductance correction factor (pu)";
  "ksat"          , "Saturation factor (pu)";
  "la"            , "Total barrier thickness (mm)";
  "lend"          , "End-winding length (mm)";
  "lt"            , "Tooth length (mm)";
  "ly"            , "Yoke length (mm)";
  "mCu"           , "Mass - stator winding (kg)";
  "mFeR"          , "Mass - rotor iron (kg)";
  "mFeS"          , "Mass - stator iron (kg)";
  "mPM"           , "Mass - rotor PM (kg)";
  "nbase"         , "Base speed (rpm)";
  "ws"            , "Slot width (mm)";
  "wt"            , "Tooth width (mm)";
  "xRaw"          , "FEAfix - x points";
  "xx"            , "x";
  "ir"            , "Rotor current (A)";
  "M"             , "Inductance - rotor-stator (H)";
  "wp"            , "Rotor pole width (mm)";
  "dalpha_pu"     , "Pole width (pu)";
  "dalpha"        , "Pole width (deg)";
  "lyr"           , "Rotor yoke length (mm)";
  "hpb"           , "Rotor pole body height (mm)";
  "hph"           , "Rotor pole head height (mm)";
  "wp"            , "Rotor pole width (mm)";
  "wb"            , "Rotor coil width (mm)";
  "hb"            , "Rotor coil height (mm)";
  "r_fillet"      , "Rotor pole head fillet (mm)";
  "r_bfillet"     , "Rotor pole yoke fillet (mm)";
  "th_Head_deg"   , "Rotor pole head angle (deg)";
  "th_Yoke_deg"   , "Rotor pole yoke angle (deg)";
  "pont0"         , "Minimum mechanical tolerance (mm)";
  "g"             , "Airgap length (mm)";
  "Acoilf"        , "Rotor coil cross-section (mm2)";
  "Ld"            , "Inductance - d-axis total (H)";
  "Lq"            , "Inductance - q-axis total (H)";
  "NrIr"          , "Normalized rotor current (A turns)";
  "mCoRo"         , "Mass - rotor winding (kg)";
  "Rr_Nr2"        , "Normalized rotor resistance (Ohm/turns^2)";
  "tempCuAvg10s"  , "Temperatures - Copper average 10s (°C)";
  "tempCuAvg30s"  , "Temperatures - Copper average 30s (°C)";
  "tempCuMax10s"  , "Temperatures - Copper hotspot 10s (°C)";
  "tempCuMax30s"  , "Temperatures - Copper hotspot 30s (°C)";
  "tempPMAvg10s"  , "Temperatures - PM average 10s (°C)";
  "tempPMAvg30s"  , "Temperatures - PM average 30s (°C)";
  "tempPMMax10s"  , "Temperatures - PM hotspot 10s (°C)";
  "tempPMMax30s"  , "Temperatures - PM hotspot 30s (°C)";
];


map2gui = containers.Map(names(:,1),names(:,2));
gui2map = containers.Map(names(:,2),names(:,1));