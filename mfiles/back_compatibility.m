% Copyright 2018
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

function [dataSet,geo,per,mat] = back_compatibility(dataSet,geo,per,Dflag)
% 
% [dataSet,geo,per,mat] = back_compatibility(dataSet,geo,per,Dflag)
% 
% Update the project to newest SyR-e version
% INPUT : dataSet
%         geo
%         per
%         Dflag =1-->disp the modification / =0-->don't plot anything
% OUTPUT: dataSet
%         geo (geometry)
%         per (performance)
%         mat (material)

flag = 0;

if nargin()<4
    Dflag=1;
end

%% geo and per
% new fields of geo and per
if ~isfield(geo,'mesh_K')
    if Dflag
        disp('rev420 - new names in geo, mesh_K, mesh_K_MOOA');
    end
    flag=1;
end

if isfield(geo,'delta_sim_MOOA')
    if Dflag
        disp('rev420 - nsim and delta_sim fields moved from geo to per');
    end
    flag=1;
end

if ~isfield(geo,'win')
    if Dflag
        disp('rev420 - new names in geo, mesh_K, mesh_K_MOOA');
    end
    flag=1;
    geo.win.avv = geo.avv;
    geo.win.avv_flag = geo.avv_flag;
    geo.win.Ns = geo.Ns;
    geo.win.kcu = geo.kcu;
    geo.win.Nbob = geo.Nbob;
    geo.win.n3phase = geo.n3phase;
    geo.win.slot_layer_pos = geo.slot_layer_pos;
    geo.win.kracc = geo.kracc;
    fields = {'avv','avv_flag','Ns','kcu','Nbob','n3phase','slot_layer_pos'};
    geo = rmfield(geo,fields);
    
end


%% dataSet

%% from version 1.1 to version 1.4 (november 2015)
% disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')
% disp('update to version 1.4')
if ~isfield(dataSet,'TargetCopperTemp')
    dataSet.TargetCopperTemp = dataSet.CopperTemp;
    dataSet=rmfield(dataSet,'CopperTemp');
    if Dflag
        disp('v1.4 - renamed CopperTemp')
    end
    flag = 1;
end

if ~isfield(dataSet,'HousingTemp')
    dataSet.HousingTemp = 50;
    if Dflag
        disp('v1.4 - added housing temperature')
    end
    flag = 1;
end

if ~isfield(dataSet,'EstimatedCopperTemp')
    dataSet.EstimatedCopperTemp = per.tempcu;
    if Dflag
        disp('v1.4 - added estimated copper temperature')
    end
    flag = 1;
end

if ~isfield(dataSet,'MagLoadingYoke')
    dataSet.MagLoadingYoke = 0.5;
    if Dflag
        disp('v1.4 - added magnetic loading of the yoke')
    end
    flag = 1;
end

if ~isfield(dataSet,'MagLoadingTooth')
    dataSet.MagLoadingTooth = 1;
    if Dflag
        disp('v1.4 - added magnetic loading of the tooth')
    end
    flag = 1;
end

if ~isfield(dataSet,'ThicknessOfPM')
    dataSet.ThicknessOfPM = geo.g*6;
    if Dflag
        disp('v1.4 - added thickness of PM (SPM)')
    end
    flag = 1;
end

if ~isfield(dataSet,'AngleSpanOfPM')
    dataSet.AngleSpanOfPM = 150;
    if Dflag
        disp('v1.4 - added angle span of PM (SPM)')
    end
    flag = 1;
end

if ~isfield(dataSet,'DepthOfBarrier')
    dataSet.DepthOfBarrier = ones(1,geo.nlay);
    if Dflag
        disp('v1.4 - added depth of barrier')
    end
    flag = 1;
end

if ~isfield(dataSet,'StatorSlotOpenBou')
    dataSet.StatorSlotOpenBou = [0.2 0.3];
    dataSet.StatorSlotOpenBouCheck = 0;
    if Dflag
        disp('v1.4 - added boundaries for acs optimization')
    end
    flag = 1;
    
end

if ~isfield(dataSet,'ToothTangDepthBou')
    dataSet.ToothTangDepthBou = geo.g*[1 3];
    dataSet.ToothTangDepthBouCheck = 0;
    if Dflag
        disp('v1.4 - added boundaries for ttd optimization')
    end
    flag = 1;
end

%% check the dimension of rotor parameters

if length(dataSet.HCpu)~=dataSet.NumOfLayers && ~strcmp(dataSet.TypeOfRotor,'SPM')
    dataSet.HCpu = dataSet.HCpu(1:dataSet.NumOfLayers);
    if Dflag
        disp('v1.4 - correct dataSet.HCpu')
    end
    flag = 1;
end

if length(dataSet.ALPHApu)~=dataSet.NumOfLayers && ~strcmp(dataSet.TypeOfRotor,'SPM')
    dataSet.ALPHApu = dataSet.ALPHApu(1:dataSet.NumOfLayers);
    if Dflag
        disp('v1.4 - correct dataSet.ALPHApu')
    end
    flag = 1;
end

if length(dataSet.DepthOfBarrier)~=dataSet.NumOfLayers && ~strcmp(dataSet.TypeOfRotor,'SPM')
    dataSet.DepthOfBarrier = geo.dx;
    dataSet.DepthOfBarrier = dataSet.DepthOfBarrier(1:dataSet.NumOfLayers);
    if Dflag
        disp('v1.4 - correct dataSet.DepthOfBarrier')
    end
    flag = 1;
end

% from version 1.4 to rev 260
% disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')
% disp('update to rev 260')

% if ~isfield(dataSet,'LossEvaluationCheck')
%     dataSet.LossEvaluationCheck = 0;
%     dataSet.HysteresisLossFactor = 0;
%     dataSet.HysteresisFrequencyFactor = 0;
%     dataSet.HysteresisFluxDenFactor = 0;
%     dataSet.EddyCurLossFactor = 0;
%     dataSet.EddyCurLossFactorEdit = 0;
%     dataSet.EvalSpeed = 0;
%     dataSet.IronMassDen = 0;
%     if Dflag
%         disp('rev260 - added parameters for iron loss evaluation')
%     end
%     flag = 1;
% end

if ~isfield(dataSet,'EvalSpeed')
    dataSet.EvalSpeed = 1500;
    if Dflag
        disp('rev260 - added rotor speed')
    end
    flag = 1;
end

if ~isfield(dataSet,'TorqueOptCheck')
    dataSet.TorqueOptCheck = 1;
    dataSet.TorRipOptCheck = 1;
    if Dflag
        disp('rev260 - added check for optimization objective')
    end
    flag = 1;
end

if ~isfield(dataSet,'Qs')
    Q = 6*geo.p*geo.q;
    t2 = gcd(round(dataSet.NumOfSlots*6*dataSet.NumOfPolePairs),2*dataSet.NumOfPolePairs);
    dataSet.Qs = Q/t2;
    clear t2;
    if Dflag
        disp('rev260 - added number of simulated stator slot')
    end
    flag = 1;
end

if ~isfield(dataSet,'xRange')
    dataSet.xRange = [0.5 0.8];
    dataSet.bRange = [0.3 0.7];
    if Dflag
        disp('rev260 - added (x,b) range for syrmDesign')
    end
    flag = 1;
end

% if ~isfield(dataSet,'BarFillFac')
%     dataSet.BarFillFac = 0;
% %     if Dflag
% %         disp('rev260 - added barrier filling factor')
% %     end
% %     flag = 1;
% end

% from rev 260 to 261
% disp('-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-')
% disp('update to 261')

if ~isfield(dataSet,'MassCuOptCheck')
    dataSet.MassCuOptCheck = 0;
    dataSet.MaxCuMass = 0;
    if Dflag
        disp('rev261 - added copper mass in optimization objective')
    end
    flag = 1;
end

if ~isfield(dataSet,'SlotLayerPosCheck')
    dataSet.SlotLayerPosCheck = 0;
    if Dflag
        disp('rev261 - added slot layer position in GUI')
    end
    flag = 1;
end

if ~isfield(dataSet,'RadRibCheck')
    dataSet.RadRibCheck = 0;
    dataSet.RadRibEdit = zeros(1,geo.nlay);
    if Dflag
        disp('rev261 - added radial ribs in GUI')
    end
    flag = 1;
end

%% Bfe and kt
if ~isfield(dataSet,'Bfe')
    dataSet.Bfe=1.5;
    dataSet.kt=1;
    if Dflag
        disp('rev289 - added Bfe and kt')
    end
    flag = 1;
end

% %% PMs area in Seg and ISeg for bonded-to-real PMs (Walter Ventura Master Thesis)
% if ~isfield(dataSet,'Areavert')
%     dataSet.Areavert  = zeros(1,4);
%     dataSet.Areaob    = zeros(1,4);
%     dataSet.Areatot   = zeros(1,4);
%     dataSet.Areavert0 = zeros(1,4);
%     dataSet.Areaob0   = zeros(1,4);
%     dataSet.dob       = ones(1,4);
%     dataSet.dvert     = ones(1,4);
%     %if Dflag
%     %    disp('rev309 - added PMs area evaluation for Seg and ISeg')
%     %end
%     %flag = 1;
% end

% Multi 3-phase option (Simone Adamo)
if ~isfield(dataSet,'Num3PhaseCircuit') %AS
    dataSet.Num3PhaseCircuit=1;
    dataSet.WinFlag=[1 1 1];
    geo.win.n3phase=1;  % add for i0 computation at line 443
    if Dflag
        disp('rev319 - added multi 3-phase windings')
    end
    flag = 1;
end

% Added parameters in dataSet for Vtype rotor geometry (Marco Gallo Master Thesis)
% if ~isfield(dataSet,'SlopeBarrier')
% 
%     dataSet.SlopeBarrier=60;
%     dataSet.SlopeBarrBou=[10 89];
%     dataSet.SlopeBarrBouCheck=0;
%     
%     if Dflag
%         disp('rev325 - added Vtype rotor geometry')
%     end
%     flag = 1;
% end

% Added parameters in dataSet for optimization Vtype rotor (Marco Gallo Master Thesis)
if ~isfield(dataSet,'MaxPMMass')
    dataSet.MaxPMMass= 1.58; %Massa volume magnete macchina di riferimento (mot.232 Pareto OUT_20180308T163427)
    dataSet.MassPMOptCheck=0;
    if Dflag
        disp('rev325 - added PM Mass in optimization objective')
    end
    flag = 1;
end

% Split radial ribs
if ~isfield(dataSet,'RadRibSplit')
    dataSet.RadRibSplit=0;
    if Dflag
        disp('rev326 - added splitted radial ribs for Seg geometry')
    end
    flag=1;
end

% FEAfix syrmDesign
if ~isfield(dataSet,'FEAfixN')
    dataSet.FEAfixN=1;
    if Dflag
        disp('rev326 - added FEAfix for syrmDesign tool')
    end
    flag=1;
end

% Flux Barrier Shift 
if ~isfield(dataSet,'thetaFBS')
    dataSet.thetaFBS=0;
    dataSet.ThetaFBSBouCheck=0;
    dataSet.ThetaFBSBou=[0 360/(6*dataSet.NumOfPolePairs*dataSet.NumOfSlots)];
    if Dflag
        disp('rev326 - added Flux Barrier Shift for Circular and Seg geometries')
    end
    flag=1;
end

% % BrDesign
% if ~isfield(dataSet,'BrDesign')
%     matTmp=material_properties_layer('Bonded-Magnet');
%     dataSet.BrDesign=matTmp.Br;
% end

%% Thermal Loading kj
if ~isfield(dataSet,'ThermalLoadKj')
    kj=dataSet.AdmiJouleLosses/(2*pi*dataSet.StatorOuterRadius/1000*dataSet.StackLength/1000);
    dataSet.ThermalLoadKj=round(kj,0);
    if Dflag
        disp('rev334 - added thermal loading factor kj')
    end
    flag=1;
end

% Parallel slot
if ~isfield(dataSet,'ParallelSlotCheck')
    dataSet.ParallelSlotCheck=0;
    if Dflag
        disp('rev337 - added parallel slot')
    end
    flag=1;
end

% Post-processing temperature
if ~isfield(dataSet,'tempPP')
    dataSet.tempPP = 20;
    if Dflag
        disp('rev337 - added PMs temperature for post-processing')
    end
    flag=1;
end

% Tangential ribs edit
if ~isfield(dataSet,'TanRibEdit')
    dataSet.TanRibEdit = dataSet.MinMechTol*ones(1,geo.nlay);
    if Dflag
        disp('rev343 - added variables tangential ribs')
    end
    flag=1;
end

% Sliding Gap
if ~isfield(dataSet,'slidingGap')
    dataSet.slidingGap = 1;
    if Dflag
        disp('rev347 - FEMM sliding gap feature added')
    end
    flag=1;
end

% Simulated current in post-processing tab
if ~isfield(dataSet,'SimulatedCurrent')
    dataSet.SimulatedCurrent = dataSet.CurrLoPP*calc_io(geo,per);
    if Dflag
        disp('rev355 - added simulated current')
    end
    flag=1;
end

% Evaluation type in post processing
if ~isfield(dataSet,'EvalType')
    dataSet.EvalType = 'singt';
    if dataSet.GammaPP==1000
        dataSet.GammaPP=45;
    end
    if Dflag
        disp('rev367 - added evaluation type in post-processing tab')
    end
    flag=1;
end

%% New PM definition (Seg and Circ)
if ~isfield(dataSet,'PMdim')
%     dataSet = rmfield(dataSet,'BrDesign');
    dataSet.PMtemp = 20;
    if strcmp(dataSet.TypeOfRotor,'Seg')
        if isfield(dataSet,'Areavert')
            dataSet.PMdim      = [dataSet.Areaob;dataSet.Areavert];
            dataSet.PMdim      = dataSet.PMdim(:,1:dataSet.NumOfLayers)./[geo.hc;geo.hc];
            dataSet.PMdim(2,1) = 0;
            dataSet = rmfield(dataSet,'dob');
            dataSet = rmfield(dataSet,'dvert');
            dataSet = rmfield(dataSet,'Areavert');
            dataSet = rmfield(dataSet,'Areaob');
            dataSet = rmfield(dataSet,'Areatot');
            dataSet = rmfield(dataSet,'Areavert0');
            dataSet = rmfield(dataSet,'Areaob0');
        else
            dataSet.PMdim      = -ones(2,dataSet.NumOfLayers);
            dataSet.PMdim(2,1) = 0;
        end
    elseif strcmp(dataSet.TypeOfRotor,'Circular')
        if isfield(dataSet,'BarFillFac')
            dataSet.PMdim = -dataSet.BarFillFac*[ones(1,dataSet.NumOfLayers);zeros(1,dataSet.NumOfLayers)];
        else
            dataSet.PMdim = -ones(1,dataSet.NumOfLayers);zeros(1,dataSet.NumOfLayers);
        end
    else
        dataSet.PMdim = zeros(2,dataSet.NumOfLayers);
    end
    if Dflag
        disp('rev371 - updated PM definition for Seg and Circ')
    end
    flag=1;
end

%% PM design
if ~isfield(dataSet,'CurrPM')
    dataSet.CurrPM = 1;
    dataSet.PMdesign.geo.rotor  = [];
    dataSet.PMdesign.geo.stator = [];
    if Dflag
        disp('rev373 - added target characteristic current and PM design structure')
    end
    flag=1;
end

%% PM optimization
if ~isfield(dataSet,'PMdimBou')
    dataSet.PMdimBou = [0 1];
    dataSet.PMdimBouCheck = 0;
    if Dflag
        disp('rev390 - added PM dimensions optimization')
    end
    flag=1;
end

% Rated current and Rs in dataSet
if ~isfield(dataSet,'RatedCurrent')
    dataSet.RatedCurrent = NaN;
    dataSet.Rs = NaN;
    if Dflag
        disp('rev393 - added rated current and phase resistance to dataSet')
    end
    flag=1;
end

% PM shape factor (beta)
if ~isfield(dataSet,'betaPMshape')
    if strcmp(dataSet.TypeOfRotor,'Vtype')
        dataSet.betaPMshape = (1-dataSet.SlopeBarrier*pi/180*2/pi)*dataSet.NumOfPolePairs/(dataSet.NumOfPolePairs-1);
        dataSet.PMdim       = -ones(1,dataSet.NumOfLayers);
    else
        if isfield(dataSet,'BarFillFac')
            dataSet.betaPMshape = dataSet.BarFillFac;
        else
            dataSet.betaPMshape = 1;
        end
    end
    if isfield(dataSet,'BarFillFac')
        dataSet=rmfield(dataSet,'BarFillFac');
    end
    if isfield(dataSet,'Barfillfac')
        dataSet=rmfield(dataSet,'Barfillfac');
    end
    dataSet.BetaPMshapeBouCheck = dataSet.SlopeBarrBouCheck;
    dataSet.BetaPMshapeBou       = dataSet.SlopeBarrBou;
    dataSet = rmfield(dataSet,'SlopeBarrBou');
    dataSet = rmfield(dataSet,'SlopeBarrBouCheck');
    
    
    if Dflag
        disp('rev398 - removed BarFillFac and V-type angle and added betaPMshape factor')
    end
    flag=1;
end

if ~isfield(dataSet,'Lend')
    if Dflag
        disp('rev421 - end turn inductance (Lend) added');
    end
    flag=1;
    dataSet.Lend = calc_Lend(geo);
%     dataSet.Lend = 0;
end

% error in dataSet.currentpathname
if strfind(dataSet.currentpathname,'\\')
    dataSet.currentpathname = strrep(dataSet.currentpathname,'\\','\');
    if Dflag
        disp('rev420 - dataSet.currentpathname corrected from \\')
    end
    flag=1;
end

% slot model for skin effect
if ~isfield(dataSet,'SlotConductorType')
    dataSet.SlotConductorType       = 'Round';
    dataSet.SlotConductorInsulation = 0.1;
    dataSet.SlotConductorShape      = 1;
    dataSet.SlotConductorRadius     = dataSet.MinMechTol;
    dataSet.SlotConductorWidth      = 2*dataSet.MinMechTol;
    dataSet.SlotConductorHeight     = 2*dataSet.MinMechTol;
    dataSet.SlotConductorNumber     = geo.win.Nbob*2;
    dataSet.SlotConductorFrequency  = [1 10 50 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500];
    if Dflag
        disp('rev432 - added slot model on GUI')
    end
    flag=1;
end

% Motor-CAD tab
if ~isfield(dataSet,'HousingType')
    dataSet.InletTemperature = 40;
    dataSet.HousingType = 'Axial fins (Servo)';
    dataSet.TransientPeriod = 60;
    dataSet.TransientTimeStep = 50;
    if Dflag
        disp('rev445 - added Motor-CAD interface')
    end
    flag=1;
    
end

% Current density value
if ~isfield(dataSet,'CurrentDensity')
    dataSet.CurrentDensity = [];
    if Dflag
        disp('rev446 - added current density on GUI')
    end
    flag=1;
end

% Rotor fillet value
if ~isfield(dataSet,'RotorFillet')
    dataSet.RotorFillet = nan(1,dataSet.NumOfLayers);
    if Dflag
        disp('rev465 - added rotor fillet on GUI')
    end
    flag=1;
end

% PM clearance
if ~isfield(dataSet,'PMclear')
    dataSet.PMclear = zeros(size(dataSet.PMdim));
    if Dflag
        disp('rev473 - added PM clearance')
    end
    flag=1;
end

% Temperature for slot model
if ~isfield(dataSet,'SlotConductorTemperature')
    dataSet.SlotConductorTemperature = dataSet.EstimatedCopperTemp;
    if Dflag
        disp('rev486 - added temperature for slot model evaluation')
    end
    flag=1;
end

% Mass and inertia
if ~isfield(dataSet,'MassWinding')
    dataSet.MassWinding    = 0;
    dataSet.MassMagnet     = 0;
    dataSet.MassStatorIron = 0;
    dataSet.MassRotorIron  = 0;
    dataSet.RotorInertia   = 0;
    if Dflag
        disp('rev487 - added mass of active part and rotor inertia')
    end
    flag=1;
end

% Radial Ribs parameters 
if ~isfield(dataSet,'pontRangEdit')
    dataSet.pontRangEdit    = zeros(1,dataSet.NumOfLayers);
    dataSet.pontRoffsetEdit = zeros(1,dataSet.NumOfLayers);
    dataSet.RotorFilletRadEdit = dataSet.MinMechTol*ones(1,dataSet.NumOfLayers);
    if Dflag
        disp('2020 12 09 - Added improved radial ribs parameters')
    end
    flag=1;
end


% SEG updates 
if ~isfield(dataSet,'RadShiftInner')
    dataSet.RadShiftInner = zeros(1,dataSet.NumOfLayers);
    dataSet.NarrowFactor = ones(1,dataSet.NumOfLayers);
    dataSet.RadRibSplit = dataSet.RadRibSplit*ones(1,dataSet.NumOfLayers);
    if Dflag
        disp('2021 02 23 - Added improved Seg barriers parameters')
    end
    flag=1;
end

% Fillet update 
if ~isfield(dataSet,'RotorFilletIn')
    dataSet.RotorFilletIn = dataSet.MinMechTol*ones(1,dataSet.NumOfLayers);
    dataSet.RotorFilletOut = dataSet.MinMechTol*ones(1,dataSet.NumOfLayers);
    dataSet = rmfield(dataSet,'RotorFilletRadEdit');
    if Dflag
        disp('2021 02 23 - Added improved rotor fillet parameters')
    end
    flag=1;
end

% Barrier Shrink 
if ~isfield(dataSet,'CentralShrink')
    dataSet.CentralShrink = zeros(1,dataSet.NumOfLayers);

    if Dflag
        disp('2021 03 22 - Added Central Barriers Shrink')
    end
    flag=1;
end


% Motor-CAD tab Update
if ~isfield(dataSet,'th_eval_type')
    dataSet.FlowRate = 6;
    dataSet.Fluid = 'W/G 50/50';
    dataSet.th_eval_type = 'Steady State';
    dataSet.MachineTemperature = 80;
    dataSet.TransientPeriod = 60;
    if Dflag
        disp('2021 03 30 - Updated Motor-CAD Tab')
    end
    flag=1; 
end

% Optimization - New Bounds and Obj   
if ~isfield(dataSet,'TanRibBou')
    dataSet.RadRibBouCheck        = 0;
    dataSet.RadRibBou             = [0 0];
    dataSet.TanRibBouCheck        = 0;
    dataSet.TanRibBou             = [0 0];
    dataSet.CentralShrinkBouCheck = 0;
    dataSet.CentralShrinkBou      = [0 0];
    dataSet.RadShiftInnerBouCheck = 0;
    dataSet.RadShiftInnerBou      = [0 0];
    dataSet.PowerFactorOptCheck   = 0;
    dataSet.MinExpPowerFactor     = 0;
    dataSet.NoLoadFluxOptCheck    = 0;
    dataSet.MaxExpNoLoadFlux      = 0; 
    dataSet.MechStressOptCheck    = 0;
    dataSet.MaxExpMechStress      = 0;
    if Dflag
        disp('2021 04 20 - Added new bounds and objectives in the Optimization tab')
    end
    flag=1;
end

% Custom geometry
if ~isfield(dataSet,'custom')
    dataSet.pShape.rotor  = [];
    dataSet.pShape.stator = [];
    dataSet.pShape.magnet = [];
    dataSet.pShape.slot   = [];
    dataSet.pShape.flag   = 0;
    dataSet.custom        = 0;

    if Dflag
        disp('2021 04 28 - Added Custom Geometry')
    end
    flag=1;
end

% Slot model - Bottom conductor gap 
if ~isfield(dataSet,'SlotConductorBottomGap')
    dataSet.SlotConductorBottomGap  = 0;
    geo.win.gapBotCond = 0;
    
    if Dflag
        disp('2021 05 03 - Added bottom conductor gap')
    end
    flag=1;
end

% Double tangential rotor fillet  
if ~isfield(geo,'RotorFilletTan1')
    dataSet.RotorFilletTan1  = nan(1,dataSet.NumOfLayers);
    dataSet.RotorFilletTan2  = nan(1,dataSet.NumOfLayers);
    dataSet.TanRibCheck = 0; 
    %geo = rmfield(geo,'RotorFillet');
    if Dflag
        disp('2021 05 05 - Added double tangential rotor fillet')
    end
    flag=1;
end

% Scaling
if ~isfield(dataSet,'ScaleCheck')
%     dataSet.ScaleFactor  = 1;
%     dataSet.ScaleFactorAxial  = 1;
    dataSet.ScaleCheck = 0;
  
    if Dflag
        disp('2021 05 12 - Added Scale Machine')
    end
    flag=1;
end

% Slot width display
if ~isfield(dataSet,'SlotWidth')
    dataSet.SlotWidth = 0;
    if Dflag
        disp('2021 05 14 - Added slot width display')
    end
    flag=1;
end

% End windigs
if ~isfield(dataSet,'EndWindingsLength')
    [~,~,geo] = calc_io(geo,per);
    dataSet.EndWindingsLength = geo.lend;
    if Dflag
        disp('2021 06 14 - Added End-windings length')
    end
    flag=1;
end

if ~isfield(dataSet,'MapQuadrants')
    if strcmp(dataSet.FluxBarrierMaterial,'Air')
        dataSet.MapQuadrants = 1;
    else
        dataSet.MapQuadrants = 2;
    end
    if Dflag
        disp('2021 07 14 - Added quadrants selection for flux maps')
    end
    flag=1;
end

% MotorCAD Thermal Limits
if ~isfield(dataSet,'AmbTemp')
    dataSet.TempCuLimit = 180;
    dataSet.InitTemp    = 45;
    dataSet.TransTime   = 30;
    dataSet.SimIth0     = NaN;
    dataSet.SimIthpk    = NaN;
    dataSet.AmbTemp     = 50;
    
    if Dflag
        disp('2021 07 25 - Added Thermal Limits Calculation')
    end
    flag=1;
end

% Display the number of parallel
if ~isfield(dataSet,'SlotConductorParallel')
    dataSet.SlotConductorParallel = round(geo.p*geo.q*geo.win.nCond/geo.win.Ns,2);    
    if Dflag
        disp('2021 08 30 - Added the number of parallel display')
    end
    flag=1;
end

%% remove old fields of dataSet
flagClear = 0;
if isfield(dataSet,'BarFillFac')
    dataSet = rmfield(dataSet,'BarFillFac');
    flagClear = 1;
end
if isfield(dataSet,'randFactor')
    dataSet = rmfield(dataSet,'randFactor');
    flagClear = 1;
end

if isfield(dataSet,'BrDesign')
    dataSet = rmfield(dataSet,'BrDesign');
    flagClear = 1;
end

if isfield(dataSet,'dob')
    dataSet = rmfield(dataSet,'dob');
    dataSet = rmfield(dataSet,'dvert');
    dataSet = rmfield(dataSet,'Areavert');
    dataSet = rmfield(dataSet,'Areaob');
    dataSet = rmfield(dataSet,'Areatot');
    dataSet = rmfield(dataSet,'Areavert0');
    dataSet = rmfield(dataSet,'Areaob0');
    flagClear = 1;
end

if isfield(dataSet,'IronMassDen')
    dataSet = rmfield(dataSet,'HysteresisFluxDenFactor');
    dataSet = rmfield(dataSet,'HysteresisFrequencyFactor');
    dataSet = rmfield(dataSet,'HysteresisLossFactor');
    dataSet = rmfield(dataSet,'EddyCurLossFactor');
    dataSet = rmfield(dataSet,'EddyCurLossFactorEdit');
    dataSet = rmfield(dataSet,'IronMassDen');
    dataSet = rmfield(dataSet,'LossEvaluationCheck');
    flagClear = 1;
end

if isfield(dataSet,'SlopeBarrier')
    dataSet = rmfield(dataSet,'SlopeBarrier');
    flagClear = 1;
end

if isfield(dataSet,'XFEMMOpt')
    dataSet = rmfield(dataSet,'XFEMMOpt');
    flagClear = 1;
end

if isfield(dataSet,'XFEMMPPMot')
    dataSet = rmfield(dataSet,'XFEMMPPMot');
    flagClear = 1;
end

if isfield(dataSet,'MachineTemperature')
    dataSet = rmfield(dataSet,'MachineTemperature');
    flagClear = 1;
end

if isfield(dataSet,'MaxExpMechStress')
    dataSet = rmfield(dataSet,'MaxExpMechStress');
    flagClear = 1;
end

if isfield(dataSet,'RotorFillet')
    dataSet = rmfield(dataSet,'RotorFillet');
    flagClear = 1;
end

if flagClear
    disp('Removed old dataSet fields')
end

if ~isfield(dataSet,'RMVTmp')
    dataSet.RMVTmp = 'ON';
end


% rewriting geo, per and mat (and check if mat exist)
% [bounds, ~, geo, per, mat] = data0(dataSet);
[bounds, ~, geo, per, mat] = data0(dataSet);
geo.RQ=buildDefaultRQ(bounds);

% check if RQ and RQnames are correct
if (length(dataSet.RQ)~=length(dataSet.RQnames) || length(geo.RQ)~=length(dataSet.RQ))
    dataSet.RQ=buildDefaultRQ(bounds);
    dataSet.RQnames=geo.RQnames;
    if Dflag
        disp('rev274 - correct RQ and RQnames')
    end
%     flag=1;
end

% message in command window if some data are added
if flag && Dflag
    msg = 'This project was created with an older version of SyR-e: proceed to SAVE MACHINE to update to latest version';
    title = 'WARNING';
%     f = warndlg(msg,title,'modal');
end

end



