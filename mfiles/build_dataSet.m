function [dataSet] = build_dataSet(geo,per)
%
% Create dataSet for the project saved with older version of SyR-e


%% Geometry (from geo)

dataSet.NumOfPolePairs = geo.p;                             % number of pole pairs
dataSet.AirGapThickness = geo.g;                            % airgap thickness [mm]
dataSet.StatorOuterRadius = geo.R;                          % stator outer radius [mm]
dataSet.AirGapRadius = geo.r;                               % machine airgap radius [mm]
dataSet.ShaftRadius = geo.Ar;                               % shaft radius [mm]
dataSet.StackLength = geo.l;                                % stack length [mm]
dataSet.TypeOfRotor = geo.RotType;                          % type of rotor
dataSet.NumOfSlots = geo.q;                                 % number of slots per pole per phase
dataSet.ToothLength = round(geo.lt,2);                      % tooth length [mm]
dataSet.StatorSlotOpen = round(geo.acs,2);                  % stator slot open in [p.u.]
dataSet.ToothWidth = round(geo.wt,2);                       % tooth width [mm]
dataSet.ToothTangDepth = geo.ttd;                           % tooth tang depth [mm]
dataSet.ToothTangAngle = geo.tta;                           % tooth tang angle (mech degree)
dataSet.FilletCorner = geo.SFR;                             % fillet at the back corner of the slot [mm]
dataSet.SlotMaterial = geo.BLKLABELSmaterials{3};           % slot material
dataSet.StatorMaterial = geo.BLKLABELSmaterials{4};         % stator material
dataSet.RotorMaterial = geo.BLKLABELSmaterials{5};          % rotor material
dataSet.FluxBarrierMaterial = geo.BLKLABELSmaterials{6};    % flux barrier material
dataSet.ShaftMaterial = geo.BLKLABELSmaterials{7};          % shaft material
dataSet.RotorCondMaterial = geo.BLKLABELSmaterials{8};      % rotor conductor material
dataSet.SlotFillFactor = geo.win.kcu;                           % slot fill factor
dataSet.PitchShortFac = round(geo.win.kracc,2);                 % pitch short factor
dataSet.TurnsInSeries= geo.win.Ns;                              % turns in series
dataSet.NumOfLayers = geo.nlay;                             % Number of Layers
dataSet.OverSpeed = geo.nmax;                               % max mechanical speed [rpm]
dataSet.Br = geo.Br;                                        % Br
dataSet.Hc = geo.Hc;                                        % Hc
dataSet.MinMechTol = geo.pont0;                             % minimum mechanical tolerance [mm]
dataSet.ALPHApu = geo.dalpha_pu;                            % angular position of rotor barriers [pu]
dataSet.HCpu = round(geo.hc_pu,2);                          % dimension of the rotor barriers [pu]
dataSet.DepthOfBarrier = geo.dx;                            % depth of rotor barrier or number of PM segments in SPM
if isfield(geo,'lm')
    dataSet.ThicknessOfPM = geo.lm;                         % length of PM in SPM machines [mm]
else
    dataSet.ThicknessOfPM = 6*geo.g;
end

if isfield(geo,'phi')
    dataSet.AngleSpanOfPM = geo.phi;                        % angle span of PM in SPM machines [elt degrees]
else
    dataSet.AngleSpanOfPM = 150;
end

dataSet.WinMatr = geo.win.avv;                                  % winding matrix

%% Performance (form per)

dataSet.AdmiJouleLosses = per.Loss;                         % max admitted losses [W]
dataSet.TargetCopperTemp = per.tempcu;                      % target copper temperature [°C]
if isfield(per,'temphous')
    dataSet.HousingTemp = per.temphous;                     % housing temperature [°C]
else
    dataSet.HousingTemp = 50;
    disp('Added parameter : Housing Temperature in other option tab')
end

if isfield(per,'tempcuest')
    dataSet.EstimatedCopperTemp = per.tempcuest;            % estimate copper temperature [°C]
else
    dataSet.EstimatedCopperTemp = per.tempcu;
end

dataSet.CurrOverLoad = per.overload;                        % overload current [pu]

%% Objective
dataSet.MinExpTorque = per.min_exp_torque;                  % minimum expected torque in the optimization [Nm]
dataSet.MaxRippleTorque = per.max_exp_ripple;               % maximum expected ripple torque [pu]
dataSet.TorqueOptCheck = 1;
dataSet.TorRipOptCheck = 1;


%% Mesh setting

if isfield(geo,'mesh_K_MOOA')
    dataSet.Mesh = geo.mesh_K;                                  % for post-processing and manual design
    dataSet.Mesh_MOOA = geo.mesh_K_MOOA;                        % for optimization
end
if isfield(geo,'K_mesh_MOOA')
    dataSet.Mesh = geo.K_mesh;                                  % for post-processing and manual design (before rev 419)
    dataSet.Mesh_MOOA = geo.K_mesh_MOOA;                        % for optimization
end

%% MOOA setting

dataSet.SimPoMOOA = 20;                     % default number of rotor position
dataSet.RotPoMOOA = 60;                     % default rotor position span [elt degrees]
dataSet.SimPoFine = 20;                     % default simulated position in the re-evaluation of Pareto front
dataSet.RotPoFine = 60;                     % default rotor position span in the re-evaluation of Pareto front [elt degrees]

% if isfield(geo,'randFactor')
%     dataSet.randFactor = geo.randFactor;                    % noise factor for position number reduction
% else
%     dataSet.randFactor = 0;
% end

dataSet.RQnames = [];                              % RQ names
dataSet.RQ = [];

dataSet.MaxGen = 60;                                        % max numer of generation
dataSet.XPop = 60;                                          % number of population
dataSet.XFEMMOpt ='N';                                      % use XFEMM for optimization


%% Boundary setting

% boundary value
dataSet.Alpha1Bou = [0.25 0.5];                             % alpha of the first barrier
dataSet.DeltaAlphaBou = [0.17 0.5];                         % alpha of the other barriers
dataSet.hcBou = [0.2 1];                                    % width of barriers
dataSet.DfeBou = [-0.75 0.75];                              % dx
dataSet.BrBou = [0.3 0.4];                                  % remeance of the PM
dataSet.GapBou = round(geo.g*[0.8 1.2],2);                  % airgap
dataSet.GapRadiusBou = round(geo.r*[0.8 1.2],2);            % rotor radius
dataSet.ToothWiBou = round(geo.wt*[0.8 1.2],2);             % tooth width;
dataSet.StatorSlotOpenBou = [0.2 0.3];                      % stator slot open (geo.acs)
dataSet.ToothTangDepthBou = [0.8 1.2];                      % tooth tang depth (geo.ttd)
dataSet.ToothLeBou = round(geo.lt*[0.8 1.2],2);             % tooth length
dataSet.PhaseAngleCurrBou = [40 75];                        % gamma

% check
%% Added bounds SlopeBarrBouCheck - rev.Gallo
if isfield(geo,'RQnames')
    dataSet.Dalpha1BouCheck = 0;
    dataSet.DalphaBouCheck = 0;
    dataSet.hcBouCheck = 0;
    dataSet.DxBouCheck = 0;
    dataSet.BrBouCheck = 0;
    dataSet.GapBouCheck = 0;
    dataSet.AirgapRadiusBouCheck = 0;
    dataSet.ToothWidthBouCheck = 0;
    dataSet.StatorSlotOpenBouCheck = 0;
    dataSet.ToothTangDepthBouCheck = 0;
    dataSet.ToothLengthBouCheck = 0;
    dataSet.GammaBouCheck = 0;
    dataSet.SlopeBarrBouCheck=0;
    for ii=1:length(geo.RQnames)
        if isequal(geo.RQnames{ii},'dalpha')
            dataSet.Dalpha1BouCheck = 1;
            dataSet.DalphaBouCheck = 1;
        elseif isequal(geo.RQnames{ii},'hc')
            dataSet.hcBouCheck = 1;
        elseif isequal(geo.RQnames{ii},'dx')
            dataSet.DxBouCheck = 1;
        elseif isequal(geo.RQnames{ii},'Br')
            dataSet.BrBouCheck = 1;
        elseif isequal(geo.RQnames{ii},'g')
            dataSet.GapBouCheck = 1;
        elseif isequal(geo.RQnames{ii},'r')
            dataSet.AirgapRadiusBouCheck = 1;
        elseif isequal(geo.RQnames{ii},'wt')
            dataSet.ToothWidthBouCheck = 1;
        elseif isequal(geo.RQnames{ii},'lt')
            dataSet.ToothLengthBouCheck = 1;
        elseif isequal(geo.RQnames{ii},'acs')
            dataSet.StatorSlotOpenBouCheck = 1;
        elseif isequal(geo.RQnames{ii},'ttd')
            dataSet.ToothTangDepthBouCheck = 1;
        elseif isequal(geo.RQnames{ii},'gamma')
            dataSet.GammaBouCheck = 1;
        elseif isequal(geo.RQnames{ii},'angleDEG')
            dataSet.SlopeBarrBouCheck = 1;
        end
    end
else
    dataSet.Dalpha1BouCheck = 0;
    dataSet.DalphaBouCheck = 0;
    dataSet.hcBouCheck = 0;
    dataSet.DxBouCheck = 0;
    dataSet.BrBouCheck = 0;
    dataSet.GapBouCheck = 0;
    dataSet.AirgapRadiusBouCheck = 0;
    dataSet.ToothWidthBouCheck = 0;
    dataSet.StatorSlotOpenBouCheck = 0;
    dataSet.ToothTangDepthBouCheck = 0;
    dataSet.ToothLengthBouCheck = 0;
    dataSet.GammaBouCheck = 0;
    dataSet.SlopeBarrBouCheck=0;
    disp('Setted all boudary condition in optimization tab')
end



%% Post-processing setting

dataSet.CurrLoPP = per.overload;                            % current overload [pu]
if strcmp(geo.RotType,'SPM')
    dataSet.GammaPP = 90;                                   % gamma [elt degrees]
else
    dataSet.GammaPP = 45;
end

dataSet.BrPP = geo.Br;                                      % PM remeance [T]
dataSet.NumGrid = 10;                                       % number of current level in map evaluation
dataSet.NumOfRotPosPP = 10;                                 % number of rotor position
dataSet.AngularSpanPP = 60;                                 % angular span on rotor [mech degrees]
dataSet.XFEMMPPMot = 'N';                                   % post-process with XFEMM

disp('Setted all parameters in post-processing tab')


%% Additional parameters

dataSet.MagLoadingYoke = 0.5;
dataSet.MagLoadingTooth = 1;
dataSet.RMVTmp = geo.RemoveTMPfile;







