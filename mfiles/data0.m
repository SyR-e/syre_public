% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [bounds, objs, geo, per, mat] = data0(dataIn)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data0.m:
% historically: manual input data used as default by the graphic input GUI
% after GUI introduction: translates dataSet (dataIn) into geo, per, mat, etc ..

% per: performance
% geo: geometry
% bounds: bounds of optimization inputs RQ

%% READ INPUTS FROM THE GUI

% main performance target
per.Loss = dataIn.AdmiJouleLosses;             % admitted Joule loss [W]
per.tempcu = dataIn.TargetCopperTemp;           % Target Copper Temperature [C]
%     per.Vdc = dataIn.DCVoltage;              % dc-link voltage [V]
per.overload = dataIn.CurrOverLoad;           % current overload factor used for optimization (1 means Joule loss = per.Loss)
per.BrPP = dataIn.BrPP;                       % Br used for postprocessing [T]
per.tempPP = dataIn.tempPP;                   % PMs temperature in postprocessing [°C]
per.temphous = dataIn.HousingTemp;            % Housing Temperature [C]
per.tempcuest = dataIn.EstimatedCopperTemp;   % Estimated Copper Temperatue [C]

per.kj = dataIn.ThermalLoadKj;
per.i0 = dataIn.RatedCurrent;
per.Rs = dataIn.Rs;
per.Lend = dataIn.Lend;

% MOOA goals penalization tresholds
per.min_exp_torque = dataIn.MinExpTorque;      % minimum expected torque [Nm]
per.max_exp_ripple = dataIn.MaxRippleTorque;    % maximum expected torque ripple in pu during opimization
per.max_Cu_mass    = dataIn.MaxCuMass;         % maximum expected copper mass [kg]
per.max_PM_mass    = dataIn.MaxPMMass;        % massima massa magneti totale [kg]
per.min_pf         = dataIn.MinExpPowerFactor;
per.max_fdq0       = dataIn.MaxExpNoLoadFlux;
% per.max_mechstress = dataIn.MaxExpMechStress;
per.EvalSpeed      = dataIn.EvalSpeed;

% number of simulated rotor positions
% MOOA means during the optimization
per.nsim_MOOA = dataIn.SimPoMOOA;         % simulated positions (6-1)
% geo.randFactor = dataIn.randFactor;       % Noise factor for position number reduction
per.delta_sim_MOOA = dataIn.RotPoMOOA;    % rotor position span [elt degrees]
% evalx means the re-evaluation stage, with a finer resolution
% per.nsim_singt = dataIn.SimPoFine;        % simulated positions (16-1)
% per.delta_sim_singt = dataIn.RotPoFine;   % rotor position span [elt degrees]
per.nsim_singt = dataIn.NumOfRotPosPP;        % simulated positions (16-1)
per.delta_sim_singt = dataIn.AngularSpanPP;   % rotor position span [elt degrees]

geo.BLKLABELSmaterials = {
    'Air';
    'Air';
    dataIn.SlotMaterial;
    dataIn.StatorMaterial;
    dataIn.RotorMaterial;
    dataIn.FluxBarrierMaterial;
    dataIn.ShaftMaterial;
    dataIn.RotorCondMaterial};

geo.pont0 = dataIn.MinMechTol;  % thickness of the structural bridges at the airgap [mm]

% Geometry
geo.RotType = dataIn.TypeOfRotor;
% 'Circular' is the Circular barrier type of rotor, for any number of barriers
% 'ISeg' draws a rotor with the external I-shaped barrier and other
%        Segmented (U-shaped) barriers
% 'Seg' draws all segmented barriers
% 'Fluid' draws barriers shaped according to fluid mechanics
% 'SPM'     : Surface mounted permanent magnet motor
% 'Vtype'   : Vtype barriers
% geo.RemoveTMPfile = dataIn.RMVTmp;      % 'ON' of 'OFF' for remuving the motor folders in tmp

% geo.RaccBarrier='OFF';
%geo.DTrasl=0;

geo.p  = dataIn.NumOfPolePairs;   % pole pairs
geo.R  = dataIn.StatorOuterRadius;% stator outer radius [mm]
geo.r  = dataIn.AirGapRadius;     % airgap radius [mm]
geo.Ar = dataIn.ShaftRadius;      % shaft radius [mm]
geo.g  = dataIn.AirGapThickness;  % airgap [mm]
geo.l  = dataIn.StackLength;      % stack length [mm]

% stator
geo.q   = dataIn.NumOfSlots;      % stator slots per pole per phase
geo.lt  = dataIn.ToothLength;    % tooth length [mm]
geo.acs = dataIn.StatorSlotOpen; % stator slot opening [p.u.]
geo.wt  = dataIn.ToothWidth;       % tooth width [mm]
geo.st  = dataIn.SlotWidth;

if dataIn.ParallelSlotCheck
    geo.parallel_slot = 1;
else
    geo.parallel_slot = 0;
end

geo.betaPMshape = dataIn.betaPMshape;               % barrier filling factor

geo.ttd = dataIn.ToothTangDepth;   % tooth tang depth [mm]
geo.tta = dataIn.ToothTangAngle;   % tooth tang angle (mech degree)
geo.SFR = dataIn.FilletCorner;     % fillet at the back corner of the slot [mm]

% rotor

if strcmp(geo.RotType,'SPM')
    geo.nlay = 1;
else
    geo.nlay  = dataIn.NumOfLayers;    % number of layers
end

geo.racc_pont = 1 * geo.pont0;                  % radius of the fillet at the sides of inner bridges (if any) [mm]
% geo.ang_pont0 = geo.pont0 / geo.r * 180/pi;    % span of the arc corresponding to pont0 at the airgap radius [deg]
geo.hfe_min   = 2*geo.pont0;                      % min tickness of each steel flux guide

% winding description
geo.win.kcu = dataIn.SlotFillFactor;                % slot filling factor (net copper/slot area)
geo.win.avv = dataIn.WinMatr;
geo.win.avv_flag   = dataIn.WinFlag; %AS
geo.win.n3phase    = dataIn.Num3PhaseCircuit; %AS stator 3-phase circuits number
geo.win.kracc      = dataIn.PitchShortFac;       % pitch shortening factor (for end connections length estimation)
geo.win.Ns         = dataIn.TurnsInSeries;          % turns in series per phase (entire motor, one way, all poles in series)
geo.win.Nbob       = geo.win.Ns/geo.p/(geo.q)/size(geo.win.avv,1);  % conductors in slot per layer
if dataIn.SlotLayerPosCheck      % stator slot layer position, two solution are possible, 'over_under', and 'side_by_side'
    geo.win.slot_layer_pos = 'side_by_side';
else
    geo.win.slot_layer_pos = 'over_under';
end
geo.win.Lend = dataIn.Lend; % end turn inductance

% slot model
geo.win.condType   = dataIn.SlotConductorType;
geo.win.condIns    = dataIn.SlotConductorInsulation;
geo.win.condShape  = dataIn.SlotConductorShape;
geo.win.rCond      = dataIn.SlotConductorRadius;
geo.win.wCond      = dataIn.SlotConductorWidth;
geo.win.hCond      = dataIn.SlotConductorHeight;
geo.win.nCond      = dataIn.SlotConductorNumber;
geo.win.pCond      = round(geo.p*geo.q*geo.win.nCond/geo.win.Ns,2);
geo.win.gapBotCond = dataIn.SlotConductorBottomGap;
per.slotModelFreq  = dataIn.SlotConductorFrequency;
per.slotModelTemp  = dataIn.SlotConductorTemperature;

geo.Qs         = dataIn.Qs;                     % number of stator slots in the FEMM simulation

geo.nmax = dataIn.OverSpeed; % overspeed [rpm]

geo.lm = dataIn.ThicknessOfPM;
geo.phi = dataIn.AngleSpanOfPM;

geo.lend = dataIn.EndWindingsLength;

% direction of magnetization in PMs of SPM with multiple segments of PM
geo.PMdir = 'p';    % parallel direction
geo.PMdir = 'r';    % radial direction

% Mesh ratio (all_motor/air-gap)
geo.mesh_K_MOOA = dataIn.Mesh_MOOA;    % optimization
geo.mesh_K = dataIn.Mesh;         % post-processing and manual design

% Rotor
geo.x0 = geo.r/cos(pi/2/geo.p);
geo.dalpha_pu = dataIn.ALPHApu;
geo.dalpha = geo.dalpha_pu*(90/geo.p);   % [mec degrees]
geo.hc_pu = dataIn.HCpu;
geo.dx = dataIn.DepthOfBarrier;
geo.dxIB = dataIn.RadShiftInner;
geo.kOB = dataIn.NarrowFactor;
geo.hcShrink = dataIn.CentralShrink;

% PM sizing
geo.PMdim     = dataIn.PMdim;
geo.flagPMFBS = 0;  % PMFBS==0 --> PMs are not deformed during FBS
% PMFBS==1 --> PMs are deformend during FBS according to the barrier area
geo.PMclear = dataIn.PMclear;
% tangential and radial ribs
geo.pontT = dataIn.TanRibEdit;
geo.pontR = dataIn.RadRibEdit;
geo.radial_ribs_eval = dataIn.RadRibCheck;
geo.radial_ribs_split = dataIn.RadRibSplit;
geo.RotorFilletTan1 = dataIn.RotorFilletTan1;
geo.RotorFilletTan2 = dataIn.RotorFilletTan2;
geo.pontRoffset = dataIn.pontRoffsetEdit;
geo.pontRang = dataIn.pontRangEdit;
geo.RotorFillet1 = dataIn.RotorFilletIn;
geo.RotorFillet2 = dataIn.RotorFilletOut;

% Flux Barrier Shift (mechanical radians)
geo.th_FBS = dataIn.thetaFBS*pi/180;  % th_FBS is the shift angle in mechanical radians

% % Sliding airgap definition
geo.slidingGap = dataIn.slidingGap;

%% bounds: limits of the search space


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTINGS OF MODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RQ and bounds
% RQ: vector of the n optimization inputs
% bounds: n x 2 vector containing the boundaries of each input and the flag
% if enabled (third column used just in the creation, then filtered and
% deleted)

rr = 1;
for ii=1:geo.nlay
    RQnames{rr} = ['dalpha_pu(' int2str(ii) ')'];
    if ii==1
        bounds(rr,:) = [dataIn.Alpha1Bou dataIn.Dalpha1BouCheck];
    else
        bounds(rr,:) = [dataIn.DeltaAlphaBou dataIn.DalphaBouCheck];
    end
    rr=rr+1;
end
for ii=1:geo.nlay
    RQnames{rr} = ['hc_pu(' int2str(ii) ')'];
    if strcmp(geo.RotType,'SPM')
        bounds(rr,:) = [geo.lm*dataIn.hcBou dataIn.hcBouCheck];
    else
        bounds(rr,:) = [dataIn.hcBou dataIn.hcBouCheck];
    end
    rr=rr+1;
end
for ii=1:geo.nlay
    RQnames{rr} = ['dx(' int2str(ii) ')'];
    bounds(rr,:) = [dataIn.DfeBou dataIn.DxBouCheck];
    if ~(strcmp(geo.RotType,'Fluid') || strcmp(geo.RotType,'Seg') || strcmp(geo.RotType,'Circular') || strcmp(geo.RotType,'Vtype'))
        bounds(rr,3) = 0; % dx not included for SPM and ISeg
    end
    rr=rr+1;
end
for ii=1:geo.nlay
    RQnames{rr} = ['Br(' int2str(ii) ')'];
    bounds(rr,:) = [dataIn.BrBou dataIn.BrBouCheck];
    rr=rr+1;
end

RQnames{rr} = 'g';          % airgap
bounds(rr,:) = [dataIn.GapBou dataIn.GapBouCheck];
rr = rr+1;
RQnames{rr} = 'r';          % rotor radius
bounds(rr,:) = [dataIn.GapRadiusBou dataIn.AirgapRadiusBouCheck];
rr = rr+1;
RQnames{rr} = 'wt';         % tooth width
bounds(rr,:) = [dataIn.ToothWiBou dataIn.ToothWidthBouCheck];
rr = rr+1;
RQnames{rr} = 'lt';         % tooth length
bounds(rr,:) = [dataIn.ToothLeBou dataIn.ToothLengthBouCheck];
rr = rr+1;
RQnames{rr} = 'acs';        % stator slot opening [p.u.]
bounds(rr,:) = [dataIn.StatorSlotOpenBou dataIn.StatorSlotOpenBouCheck];
rr = rr+1;
RQnames{rr} = 'ttd';        % tooth tang depth [mm]
bounds(rr,:) = [dataIn.ToothTangDepthBou dataIn.ToothTangDepthBouCheck];
rr = rr+1;
RQnames{rr} = 'th_FBS';     % flux barrier shift [mech deg]
bounds(rr,:) = [dataIn.ThetaFBSBou dataIn.ThetaFBSBouCheck];
rr = rr+1;

for ii=1:geo.nlay
    RQnames{rr} = ['betaPMshape(' int2str(ii) ')'];
    bounds(rr,:) = [dataIn.BetaPMshapeBou dataIn.BetaPMshapeBouCheck];
    rr=rr+1;
end

for ii=1:length(geo.PMdim(:))
    RQnames{rr} = ['PMdim(' int2str((ii)) ')'];
    bounds(rr,:) = [dataIn.PMdimBou dataIn.PMdimBouCheck];
    if (strcmp(geo.RotType,'Circular')||strcmp(geo.RotType,'Vtype'))
        if rem(ii,2)==0
            bounds(rr,3) = 0; % just central PM for Circular and Vtype
        end
    end
    rr=rr+1;
end

for ii=1:geo.nlay
    RQnames{rr} = ['pontR(' int2str(ii) ')'];
    bounds(rr,:) = [dataIn.RadRibBou dataIn.RadRibBouCheck];
    rr=rr+1;
end

for ii=1:geo.nlay
    RQnames{rr} = ['pontT(' int2str(ii) ')'];
    bounds(rr,:) = [dataIn.TanRibBou dataIn.TanRibBouCheck];
    rr=rr+1;
end

for ii=1:geo.nlay
    RQnames{rr} = ['hcShrink(' int2str(ii) ')'];
    bounds(rr,:) = [dataIn.CentralShrinkBou dataIn.CentralShrinkBouCheck];
    rr=rr+1;
end

for ii=1:geo.nlay
    RQnames{rr} = ['dxIB(' int2str(ii) ')'];
    bounds(rr,:) = [dataIn.RadShiftInnerBou dataIn.RadShiftInnerBouCheck];
    rr=rr+1;
end

RQnames{rr} = 'gamma';
bounds(rr,:) = [dataIn.PhaseAngleCurrBou dataIn.GammaBouCheck];


filt_bounds = (bounds(:,3)==1);
if geo.nlay == 1
    filt_bounds(2) = [];
    bounds(2,:) = [];
end
bounds = bounds(filt_bounds,1:2);
RQnames = RQnames(filt_bounds);
geo.RQnames = RQnames;

%% OBJECTIVES
objs = [
    per.min_exp_torque  dataIn.TorqueOptCheck           0
    per.max_exp_ripple  dataIn.TorRipOptCheck           0
    per.max_Cu_mass     dataIn.MassCuOptCheck           0
    per.max_PM_mass     dataIn.MassPMOptCheck           0
    per.min_pf          dataIn.PowerFactorOptCheck      0.1
    per.max_fdq0        dataIn.NoLoadFluxOptCheck       0
    ];
per.MechStressOptCheck =  0;

filt_objs = (objs(:,2)==1);
objs = objs(objs(:,2)==1,:);
per.objs=objs;

% end

% names of the MODE objectives
OBJnames{1} = 'Torque';
OBJnames{2} = 'TorRip';
OBJnames{3} = 'MassCu';
OBJnames{4} = 'MassPM';
OBJnames{5} = 'PF';
OBJnames{6} = 'Fdq0';
% OBJnames{7} = 'MechStress';

% eliminate unnecessary OBJnames
OBJnames = OBJnames(filt_objs);
geo.OBJnames = OBJnames;

%% Machine periodicity
Q = round(6*geo.q*geo.p*geo.win.n3phase);          % number of slots
geo.t = gcd(Q,geo.p);    % number of periods
t2=gcd(Q,2*geo.p);
QsCalc=Q/t2;
psCalc=2*geo.p/t2;

if isfield(geo,'Qs')  % Qs set in the GUI
    geo.ps = round(psCalc*geo.Qs/QsCalc);
    % Check for periodicity: if ps is odd, the simulated portion of motor
    % is anti-periodic; if ps is even, the simulated portion is periodic.
    if rem(geo.ps,2)==0
        geo.periodicity = 4;   % periodic configuration
    else
        geo.periodicity = 5;   % anti-periodic configuration
    end
    geo.tempWinTable = geo.win.avv;
else
    geo.Qs = QsCalc;
    geo.ps = psCalc;
end

mat = assign_mat_prop(dataIn);

