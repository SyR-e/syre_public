%% manual_dataSet.m
% Used for build dataSet and allow the manual input of SyR-e

%% Main data
dataSet.NumOfPolePairs    = 2;              % number of pole pairs
dataSet.NumOfSlots        = 2;              % number of slot per pole per phase
dataSet.AirGapThickness   = 0.5;            % airgap thichness [mm]
dataSet.StatorOuterRadius = 67.5;           % stator outer radius [mm]
dataSet.AirGapRadius      = 40.55;          % rotor outer radius [mm]
dataSet.ShaftRadius       = 17;             % shaft radius [mm]
dataSet.StackLength       = 101;            % stack length [mm]
dataSet.TypeOfRotor       = 'Circular';     % type of rotor: Circular, Seg, ISeg, Fluid, SPM

% parameter for parametric design
dataSet.xRange       = [0.5 0.8];           % x=rotor/stator radius
dataSet.bRange       = [0.3 0.7];           % b=iron/copper
dataSet.CurrOverLoad = 2;                   % current overload [pu]
dataSet.Bfe          = 1.5;                 % Steel loading [T]
dataSet.kt           = 1;                   % tooth factor

%% Geometrical data

dataSet.MagLoadingYoke  = 0.5;              % magnetic loading in the yoke
dataSet.MagLoadingTooth = 1;                % magnetic loading tooth

% stator
dataSet.ToothLength    = 16.81;             % tooth length [mm]
dataSet.ToothWidth     = 5.31;              % tooth width [mm]
dataSet.StatorSlotOpen = 0.25;              % stator slot open [pu]
dataSet.ToothTangDepth = 0.75;              % tooth shoe depth [mm]
dataSet.ToothTangAngle = 25;                % tooth shoe angle [°]
dataSet.FilletCorner   = 1;                 % radius of the fillet at the back iron corner of the slot [mm]

% rotor
dataSet.NumOfLayers    = 3;                 % number of layer
dataSet.ALPHApu        = [0.4 0.25 0.25];   % angle of the flux barrier [pu]
dataSet.HCpu           = [0.5 0.5 0.5];     % width of flux barrier [pu]
dataSet.DepthOfBarrier = [0 0 0];           % barrier translation on q-axis [pu]
dataSet.RadRibCheck    = 0;                 % check for manual radial ribs input
dataSet.RadRibEdit     = [0 0 0];           % manual input of radial ribs [mm]
%dataSet.SlopeBarrier   = 60;               %angle in degree - rev.Gallo

% SPM
if strcmp(dataSet.TypeOfRotor,'SPM')
    dataSet.DepthOfBarrier = 1;             % number of segment
end
dataSet.ThicknessOfPM  = 5;                 % thickness of PM [mm]
dataSet.AngleSpanOfPM  = 150;               % angle span of PM [electrical degree]

%% Other option data
dataSet.AdmiJouleLosses     = 1000;         % Permitted Joule losses [W]
dataSet.TargetCopperTemp    = 130;          % Target copper temperature [°C]
dataSet.HousingTemp         = 70;           % Housing temperature [°C]
dataSet.OverSpeed           = 4000;         % maximum speed [rpm]
dataSet.Mesh                = 2;            % mesh parameter
dataSet.Mesh_MOOA           = 10;           % mesh parameter
dataSet.MinMechTol          = 0.4;          % minimum mechanical dimension [mm]
dataSet.betaPMshape         = 1;            % barrier filling factor [pu] (only for circular geometry)
dataSet.Br                  = 0.4;          % remanence of PM [T] (used only for Bonded-Magnet)
dataSet.EstimatedCopperTemp = 130;          % estimated copper temperature
dataSet.Rs                  = NaN;          % phase resistance (evaluated from scripts)


%Seg and ISeg PMs
dataSet.Areaob    = zeros(1,4);
dataSet.Areavert  = zeros(1,4);
dataSet.Areatot   = zeros(1,4);
dataSet.dob       = ones(1,4);
dataSet.dvert     = ones(1,4);
dataSet.Areaob0   = zeros(1,4);
dataSet.Areavert0 = zeros(1,4);

%% Windings data
dataSet.SlotFillFactor    = 0.407;          % slot filling factor [pu]
dataSet.TurnsInSeries     = 122;            % turns in series per phase
dataSet.PitchShortFac     = 1;              % shortening factor
dataSet.SlotLayerPosCheck = 1;              % if 1, in the fractional slot, the slot layers are side-by-side
dataSet.Qs                = 6;              % number of simulated slots
dataSet.WinMatr = [1 1 -3 -3 2 2
                   1 1 -3 -3 2 2];          % winding matrix

%% Materials data
dataSet.SlotMaterial        = 'Copper';         % slot material
dataSet.StatorMaterial      = 'M600-50A';       % stator core material
dataSet.RotorMaterial       = 'M600-50A';       % rotor core material
dataSet.FluxBarrierMaterial = 'Bonded-Magnet';  % flux barrier material
dataSet.ShaftMaterial       = 'M600-50A';       % shaft material
dataSet.RotorCondMaterial   = 'Copper';         % rotor conductor material

%% Optimization data

% optimization parameters
dataSet.MaxGen    = 3;                 % number of generation
dataSet.XPop      = 4;                 % number of individuals for each generation
dataSet.SimPoMOOA = 5;                  % number of simulated position during evolution
dataSet.RotPoMOOA = 30;                 % rotor angular excursion during evolution [elt degree]
dataSet.SimPoFine = 31;                 % number of simulated position during Pareto front re-evaluation
dataSet.RotPoFine = 60;                 % rotor angular excursion during Pareto front re-evaluaton [elt degree]

dataSet.XFEMMOpt = 'N';                 % use XFEMM for optimization
dataSet.RMVTmp   = 'ON';                % remove temp files (only for XFEMM)

dataSet.randFactor = 0;                 % noise factor for position number reduction

% variables check
dataSet.Dalpha1BouCheck        = 1;     % first (outer) barrier angle
dataSet.DalphaBouCheck         = 1;     % other barrier angle
dataSet.hcBouCheck             = 1;     % flux barrier width
dataSet.DxBouCheck             = 1;     % flux barrier translation
dataSet.GapBouCheck            = 0;     % airgap thickness
dataSet.BrBouCheck             = 0;     % remanence of PM
dataSet.AirgapRadiusBouCheck   = 0;     % rotor radius
dataSet.ToothWidthBouCheck     = 0;     % tooth width
dataSet.ToothLengthBouCheck    = 0;     % tooth length
dataSet.StatorSlotOpenBouCheck = 0;     % stator slot open
dataSet.ToothTangDepthBouCheck = 0;     % tooth shoe depth
dataSet.GammaBouCheck          = 0;     % current phase angle

% boundary
dataSet.Alpha1Bou         = [0.25 0.5];     % first (outer) barrier angle
dataSet.DeltaAlphaBou     = [0.17 0.5];     % other barrier angle
dataSet.hcBou             = [0.2 1];        % flux barrier width
dataSet.DfeBou            = [-0.75 0.75];   % flux barrier translation
dataSet.GapBou            = [0.4 0.8];      % airgap thickness
dataSet.BrBou             = [0.3 0.4];      % remanence of PM
dataSet.GapRadiusBou      = [52 78];        % rotor radius
dataSet.ToothWiBou        = [3.8 6.3];      % tooth width
dataSet.ToothLeBou        = [15 22.2];      % tooth length
dataSet.StatorSlotOpenBou = [0.2 0.3];      % stator slot open
dataSet.ToothTangDepthBou = [0.8 1.2];      % tooth shoe depth
dataSet.PhaseAngleCurrBou = [40 75];        % current phase angle

% objectives check
dataSet.TorqueOptCheck = 1;     % maximize torque
dataSet.TorRipOptCheck = 1;     % minimize peak-to-peak torque ripple
dataSet.MassCuOptCheck = 0;     % minimize the copper mass

% objective values
dataSet.MinExpTorque    = 10;   % minimum expected torque
dataSet.MaxRippleTorque = 2;    % max admitted torque ripple
dataSet.MaxCuMass       = 10;   % max admitted copper mass

%% Post-processing data

dataSet.XFEMMPPMot = 'N';       % post-process with XFEMM
dataSet.AngularSpanPP = 60;     % rotor angular excursion [elt degree]
dataSet.GammaPP       = 45;     % phase current angle [elt degree]
dataSet.CurrLoPP      = 1;      % current load [pu]
dataSet.NumOfRotPosPP = 7;      % number of rotor position
dataSet.BrPP          = 0.4;    % remanence of PM [T]
dataSet.NumGrid       = 10;     % number of points for 1 side of the map

% iron loss

dataSet.LossEvaluationCheck = 0;    % evaluate also the iron losses
dataSet.EvalSpeed           = 1000; % evaluation speed [rpm]






