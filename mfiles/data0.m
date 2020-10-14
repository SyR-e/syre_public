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
per.max_Cu_mass = dataIn.MaxCuMass;         % maximum expected copper mass [kg]
per.max_PM_mass= dataIn.MaxPMMass;        % massima massa magneti totale [kg]

per.EvalSpeed = dataIn.EvalSpeed;

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

geo.RaccBarrier='OFF';
% geo.DTrasl=0;

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

    if dataIn.ParallelSlotCheck
        geo.parallel_slot = 1;
    else
        geo.parallel_slot = 0;
    end
    
geo.betaPMshape=dataIn.betaPMshape;               % barrier filling factor

geo.ttd = dataIn.ToothTangDepth;   % tooth tang depth [mm]
geo.tta = dataIn.ToothTangAngle;   % tooth tang angle (mech degree)
geo.SFR = dataIn.FilletCorner;     % fillet at the back corner of the slot [mm]

% rotor
% Parameter for Vtype rotor geometry: slope of semi-barrier - rev.Gallo
% if strcmp(geo.RotType,'Vtype')
%     geo.VanglePM=dataIn.SlopeBarrier*pi/180;
% else
%     geo.VanglePM=NaN;
% end

if strcmp(geo.RotType,'SPM')
    geo.nlay = 1;
    %dataIn.DepthOfBarrier = 1;        %OCT
else
    geo.nlay  = dataIn.NumOfLayers;    % number of layers
end

geo.racc_pont = 1 * geo.pont0;                  % radius of the fillet at the sides of inner bridges (if any) [mm]
geo.ang_pont0 = geo.pont0 / geo.r * 180/pi;    % span of the arc corresponding to pont0 at the airgap radius [deg]
geo.hfe_min   = 2*geo.pont0;                      % min tickness of each steel flux guide

% winding description
geo.win.kcu = dataIn.SlotFillFactor;                % slot filling factor (net copper/slot area)
geo.win.avv = dataIn.WinMatr;
geo.win.avv_flag   = dataIn.WinFlag; %AS
geo.win.n3phase = dataIn.Num3PhaseCircuit; %AS stator 3-phase circuits number
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
geo.win.condType  = dataIn.SlotConductorType;
geo.win.condIns   = dataIn.SlotConductorInsulation;
geo.win.condShape = dataIn.SlotConductorShape;
geo.win.rCond     = dataIn.SlotConductorRadius;
geo.win.wCond     = dataIn.SlotConductorWidth;
geo.win.hCond     = dataIn.SlotConductorHeight;
geo.win.nCond     = dataIn.SlotConductorNumber;
per.slotModelFreq = dataIn.SlotConductorFrequency;
per.slotModelTemp = dataIn.SlotConductorTemperature;

geo.Qs         = dataIn.Qs;                     % number of stator slots in the FEMM simulation

geo.nmax = dataIn.OverSpeed; % overspeed [rpm]

geo.lm = dataIn.ThicknessOfPM;
geo.phi = dataIn.AngleSpanOfPM;
% direction of magnetization in PMs of SPM with multiple segments of PM
geo.PMdir = 'p';    % parallel direction
geo.PMdir = 'r';    % radial direction

% Mesh ratio (all_motor/air-gap)
geo.mesh_K_MOOA = dataIn.Mesh_MOOA;    % optimization
geo.mesh_K = dataIn.Mesh;         % post-processing and manual design

geo.x0 = geo.r/cos(pi/2/geo.p);
geo.dalpha_pu = dataIn.ALPHApu;
% rev.Gallo
geo.dalpha = geo.dalpha_pu*(90/geo.p);   % [mec degrees]
geo.hc_pu = dataIn.HCpu;
geo.dx = dataIn.DepthOfBarrier;

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
geo.RotorFillet=dataIn.RotorFillet;


% Flux Barrier Shift (mechanical radians)
geo.th_FBS=dataIn.thetaFBS*pi/180;  % th_FBS is the shift angle in mechanical radians

% % Sliding airgap definition
geo.slidingGap=dataIn.slidingGap;

%% bounds: limits of the search space
% dalpha1 [p.u.]
bounds_dalpha_1 = [dataIn.Alpha1Bou(1) dataIn.Alpha1Bou(2)];   % first angle [deg]
% dalpha(2:nlay) [p.u.]
if geo.nlay == 1
    bounds_dalpha = [dataIn.DeltaAlphaBou(1) dataIn.DeltaAlphaBou(2)]; % other angles [p.u.]
else
    bounds_dalpha = ones(geo.nlay-1,1) * [dataIn.DeltaAlphaBou(1) dataIn.DeltaAlphaBou(2)]; % other angles [p.u.]
end

% Slope Barrier bounds Vtype rotor geometry - rev.Gallo
% if strcmp(geo.RotType, 'Vtype')
%     bounds_angleDEG = [dataIn.SlopeBarrBou(1) dataIn.SlopeBarrBou(2)];   %bound slope barrier [deg]
% else
%     bounds_angleDEG = [0 0];
% end
if strcmp(geo.RotType, 'SPM')
    % thickness of PM
    bounds_hc = geo.lm * [dataIn.hcBou(1) dataIn.hcBou(2)];
%     geo.hc_pu = mean(bounds_hc);
else
    % barrier ticknesses [p.u.]
    bounds_hc = ones(geo.nlay,1) * [dataIn.hcBou(1)  dataIn.hcBou(2)];
end
% barriers radial offset [p.u.]
bounds_dx = ones(geo.nlay,1) * [dataIn.DfeBou(1) dataIn.DfeBou(2)];
% remanence of the PMs
bounds_Br = ones(geo.nlay,1) * [dataIn.BrBou(1) dataIn.BrBou(2)];
% airgap
bounds_g = [dataIn.GapBou(1)  dataIn.GapBou(2)];
% Rotor radius
bounds_xr = [dataIn.GapRadiusBou(1) dataIn.GapRadiusBou(2)];
% Tooth width [mm]
bounds_wt = [dataIn.ToothWiBou(1) dataIn.ToothWiBou(2)];
% Tooth length [mm]
bounds_lt = [dataIn.ToothLeBou(1) dataIn.ToothLeBou(2)];
% stator slot opening [p.u.]
bounds_acs = [dataIn.StatorSlotOpenBou(1) dataIn.StatorSlotOpenBou(2)];
% tooth tang depth [mm]
bounds_ttd = [dataIn.ToothTangDepthBou(1) dataIn.ToothTangDepthBou(2)];
% phase angle of the current vector
bounds_gamma = [dataIn.PhaseAngleCurrBou(1) dataIn.PhaseAngleCurrBou(2)];
% FBS
bounds_thFBS = [dataIn.ThetaFBSBou(1) dataIn.ThetaFBSBou(2)];
% PM shape
bounds_beta = ones(geo.nlay,1) * [dataIn.BetaPMshapeBou(1) dataIn.BetaPMshapeBou(2)];
% PMdim
if (strcmp(dataIn.TypeOfRotor,'Seg')||strcmp(dataIn.TypeOfRotor,'ISeg'))
    indexPM = [1:1:geo.nlay;geo.nlay+(1:1:geo.nlay)];
elseif (strcmp(dataIn.TypeOfRotor,'Circular')||strcmp(dataIn.TypeOfRotor,'Vtype'))
    indexPM = 1:1:geo.nlay;geo.nlay;
else
    indexPM = zeros(1,geo.nlay);
end

bounds_PMdim = indexPM(:)./indexPM(:) * [dataIn.PMdimBou(1) dataIn.PMdimBou(2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTINGS OF MODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RQ and bounds
% RQ: vector of the n optimization inputs
% bounds: n x 2 vector containing the boundaries of each input
yes_vector = ones(geo.nlay,1);
no_vector  = zeros(geo.nlay,1);
% Added case Vtype - rev.Gallo
if (strcmp(geo.RotType,'Fluid') || strcmp(geo.RotType,'Seg') || strcmp(geo.RotType,'Circular') || strcmp(geo.RotType,'Vtype'))
    flag_dx = yes_vector*dataIn.DxBouCheck;
else
    flag_dx = no_vector*dataIn.DxBouCheck;
end
% Re-define bounds case nlay=1 - rev.Gallo
if geo.nlay == 1
    bounds = [
        bounds_dalpha_1 dataIn.Dalpha1BouCheck
        bounds_dalpha   yes_vector(1)*dataIn.DalphaBouCheck
        bounds_hc       dataIn.hcBouCheck
        bounds_dx       flag_dx
        bounds_Br       dataIn.BrBouCheck
        bounds_beta     dataIn.BetaPMshapeBouCheck
        bounds_g        dataIn.GapBouCheck
        bounds_xr       dataIn.AirgapRadiusBouCheck
        bounds_wt       dataIn.ToothWidthBouCheck
        bounds_lt       dataIn.ToothLengthBouCheck
        bounds_acs      dataIn.StatorSlotOpenBouCheck
        bounds_ttd      dataIn.ToothTangDepthBouCheck
        bounds_thFBS    dataIn.ThetaFBSBouCheck
        bounds_PMdim    dataIn.PMdimBouCheck*(indexPM(:)./indexPM(:))
        bounds_gamma    dataIn.GammaBouCheck];
    
else
    bounds = [
        bounds_dalpha_1 dataIn.Dalpha1BouCheck
        bounds_dalpha   yes_vector(1:end-1)*dataIn.DalphaBouCheck
        bounds_hc       yes_vector*dataIn.hcBouCheck
        bounds_dx       flag_dx
        bounds_Br       yes_vector*dataIn.BrBouCheck
        bounds_beta     yes_vector*dataIn.BetaPMshapeBouCheck
        bounds_g        dataIn.GapBouCheck
        bounds_xr       dataIn.AirgapRadiusBouCheck
        bounds_wt       dataIn.ToothWidthBouCheck
        bounds_lt       dataIn.ToothLengthBouCheck
        bounds_acs      dataIn.StatorSlotOpenBouCheck
        bounds_ttd      dataIn.ToothTangDepthBouCheck
        bounds_thFBS    dataIn.ThetaFBSBouCheck
        bounds_PMdim    dataIn.PMdimBouCheck*(indexPM(:)./indexPM(:))
        bounds_gamma    dataIn.GammaBouCheck];
    
end

filt_bounds = (bounds(:,3)==1);
if geo.nlay == 1
    filt_bounds(2) = [];
    bounds(2,:) = [];
end
bounds = bounds(filt_bounds,1:2);

%% OBJECTIVES
objs = [per.min_exp_torque  dataIn.TorqueOptCheck
    per.max_exp_ripple      dataIn.TorRipOptCheck
    per.max_Cu_mass         dataIn.MassCuOptCheck
    per.max_PM_mass         dataIn.MassPMOptCheck];

filt_objs = (objs(:,2)==1);
objs = objs(objs(:,2)==1,:);
per.objs=objs;

% end

% names of the variables in RQ
RQnames{1} = 'dalpha';
if geo.nlay > 1
    for k = 2:geo.nlay
        RQnames{k} = 'dalpha';
    end
    kend = k;
    for k = kend+1:kend+geo.nlay
        RQnames{k} = 'hc';
    end
    kend = k;
    for k = kend+1:kend+geo.nlay
        RQnames{k} = 'dx';
    end
    kend = k;
    for k = kend+1:kend+geo.nlay
        RQnames{k} = 'Br';
    end
    kend = k;
    for k = kend+1:kend+geo.nlay
        RQnames{k} = 'betaPMshape';
    end
else
    RQnames{2} = 'hc';
    RQnames{3} = 'dx';
    RQnames{4} = 'Br';
    RQnames{5} = 'betaPMshape';
    k = 5;
end
RQnames{k+1} = 'g';         % airgap
RQnames{k+2} = 'r';         % rotor radius
RQnames{k+3} = 'wt';        % tooth width
RQnames{k+4} = 'lt';        % tooth length
RQnames{k+5} = 'acs';       % stator slot opening [p.u.]
RQnames{k+6} = 'ttd';       % tooth tang depth [mm]
%RQnames{k+7} = 'VanglePM';  % slope barrier [deg]
RQnames{k+7} = 'th_FBS';    % flux barrier shift [mech deg]

k=k+8;

for ii=1:numel(indexPM)
    RQnames{k} = ['PMdim(' int2str(indexPM(ii)) ')'];
    k=k+1;
end

RQnames{k} = 'gamma';     % idq current phase angle

% eliminate unnecessary RQnames
RQnames = RQnames(filt_bounds);
geo.RQnames = RQnames;

% names of the MODE objectives
OBJnames{1} = 'Torque';
OBJnames{2} = 'TorRip';
OBJnames{3} = 'MassCu';
OBJnames{4} = 'MassPM';


% eliminate unnecessary OBJnames
OBJnames = OBJnames(filt_objs);
geo.OBJnames = OBJnames;

% Machine periodicity
Q = 6*geo.q*geo.p*geo.win.n3phase;          % number of slots
geo.t = gcd(Q,geo.p);    % number of periods
t2=gcd(Q,2*geo.p);
QsCalc=Q/t2;
psCalc=2*geo.p/t2;
% if rem(psCalc,2)==0 %AS
%     %periodic machine
% %     geo.tempWinTable = geo.defaultavv;
%     geo.tempWinTable = geo.win.avv;
%     geo.periodicity = 4;
% else
%     %anti-periodic machine
% %     geo.tempWinTable = [geo.defaultavv -geo.defaultavv];
%     geo.tempWinTable = [geo.win.avv -geo.win.avv];
%     geo.periodicity = 5;
% end
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

