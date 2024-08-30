warning off;
InverterModel = 'Switching';    
%InverterModel = 'Average';     

load motorModel;

Slx_name = [motorModel.data.motorName '_ctrl_INST'];

motorModelType = motorModel.SyreDrive.modelSetup.motorModelType;
inverterModelType = motorModel.SyreDrive.modelSetup.modelType;

% S-fun Parameters
Tstep = 2e-6;
Ts    = 1/motorModel.SyreDrive.Converter.fPWM;


%% ----------------Machine and Converter  Parameters------------------------%
VDC     = motorModel.data.Vdc;
Rs      = motorModel.data.Rs;
p       = motorModel.data.p;
J       = motorModel.data.J;
Rfe     = 1e5;
accel   = 10000; % rpm/s
Bm      = 0;                     % damping constant (Nm/(rad/sec))
Tf      = 0;                 % friction loss (Nm)
Tv      = 0;%P0/(nmax*pi/30)^3;     % ventilation loss coefficient (W/(rad/s)^3 = Nm/(rad/s)^2)

MTPA    = motorModel.controlTrajectories.MTPA;
i0      = motorModel.data.i0;

switch motorModel.SyreDrive.modelSetup.Ctrl_type
    case 'Current control'
        idRef = motorModel.WaveformSetup.CurrAmpl*cosd(motorModel.WaveformSetup.CurrAngle);
        iqRef = motorModel.WaveformSetup.CurrAmpl*sind(motorModel.WaveformSetup.CurrAngle);
        nRef  = motorModel.WaveformSetup.EvalSpeed;
    otherwise
        idRef = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.id,i0);
        iqRef = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.iq,i0);
        nRef  = motorModel.data.n0;
end
Ld_inic = interp2(motorModel.IncInductanceMap_dq.Id,motorModel.IncInductanceMap_dq.Iq,motorModel.IncInductanceMap_dq.Ldd,idRef,iqRef);
Lq_inic = interp2(motorModel.IncInductanceMap_dq.Id,motorModel.IncInductanceMap_dq.Iq,motorModel.IncInductanceMap_dq.Lqq,idRef,iqRef);

T0(~isnan(motorModel.data.T0)) = motorModel.data.T0;
T0(isempty(T0)) = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.T,i0);
% n0(~isnan(motorModel.data.n0)) = motorModel.data.n0;
% n0(isempty(n0)) = 1000;
if nRef==0
    nRef = 1000;
end

clear MTPA

% Converter
V0 = motorModel.SyreDrive.Converter.V0;             % power semiconductors ON treshold [V]
Rd = motorModel.SyreDrive.Converter.Rd;             % power semiconductors incremental resistance [Ohm]
dT = motorModel.SyreDrive.Converter.dT * 1e-6;      % dead time [s]
      
%% ----------------PM Initialization and Convention------------------------%
Fd0 = interp2(motorModel.FluxMap_dq.Id,motorModel.FluxMap_dq.Iq,motorModel.FluxMap_dq.Fd,0,0);
Fq0 = interp2(motorModel.FluxMap_dq.Id,motorModel.FluxMap_dq.Iq,motorModel.FluxMap_dq.Fq,0,0);
th0 = 0; %initial mechanical angle
F0_alpha = Fd0*cos(th0) - Fq0*sin(th0);
F0_beta  = Fd0*sin(th0) + Fq0*cos(th0);

if strcmp(motorModel.data.axisType,'SR') && strcmp(motorModel.data.motorType,'SR')
    Quad_Maps = 0; %SyR Convention - 1st quadrant maps
elseif strcmp(motorModel.data.axisType,'SR') && strcmp(motorModel.data.motorType,'PM')
    Quad_Maps = 1; %PM-SyR - 1st and 4th quadrant maps
elseif strcmp(motorModel.data.axisType,'PM') && strcmp(motorModel.data.motorType,'PM')
    Quad_Maps = 2; %IPM - 1st and 2st quadrant maps
end


%% ---------------Magnet Flux Estimation----------------------------%

Idd=motorModel.FluxMap_dq.Id;
Iqq=motorModel.FluxMap_dq.Iq;
Fdd=motorModel.FluxMap_dq.Fd;
Fqq=motorModel.FluxMap_dq.Fq;

if(Quad_Maps == 0)      %Syr
    Fm = 0;
elseif (Quad_Maps == 1) %PM-Syr
    Fm = abs(interp2(Idd,Iqq,Fqq,0,0));
elseif (Quad_Maps == 2) %IPM
    Fm = interp2(Idd,Iqq,Fdd,0,0);
end

clear Idd Iqq Fdd Fqq

switch(Quad_Maps)
    case 0
        InitIntg_d = 0.0;
        InitIntg_q = 0.0;
    case 1
        InitIntg_d = Fm*sin(th0);
        InitIntg_q = -Fm*cos(th0);
    case 2
        InitIntg_d = Fm*cos(th0);
        InitIntg_q = Fm*sin(th0);
end


%% ----------------dqt Inverse Flux Maps ------------------------%

% dq or dqt model
switch motorModel.SyreDrive.modelSetup.FMapsModel
    case 'dq Model'
        FMapsModel = 1;
        fD_pu_norm = motorModel.FluxMapInv_dq.fD_pu_norm;
        fQ_pu_norm = motorModel.FluxMapInv_dq.fQ_pu_norm;
        iD_pu_norm = motorModel.FluxMapInv_dq.iD_pu_norm;
        iQ_pu_norm = motorModel.FluxMapInv_dq.iQ_pu_norm;

        switch(Quad_Maps)
            case {0,2}
                fD_vct_ref = motorModel.FluxMapInv_dq.fD_vct_ref;
                fQ_vct_max = motorModel.FluxMapInv_dq.fQ_vct_max;
            case 1
                fQ_vct_ref = motorModel.FluxMapInv_dq.fQ_vct_ref;
                fD_vct_max = motorModel.FluxMapInv_dq.fD_vct_max;
        end
        
    case 'dqt Model'
        FMapsModel = -1;
        
        Fd_max = max(motorModel.FluxMapInv_dqt.dataF.Fd,[],'all');
        Fq_max = max(motorModel.FluxMapInv_dqt.dataF.Fq,[],'all');
        Fd_min = min(motorModel.FluxMapInv_dqt.dataF.Fd,[],'all');
        Fq_min = min(motorModel.FluxMapInv_dqt.dataF.Fq,[],'all');
        th_min = min(motorModel.FluxMapInv_dqt.dataF.th,[],'all');
        th_max = max(motorModel.FluxMapInv_dqt.dataF.th,[],'all');
        th_dqt = motorModel.FluxMapInv_dqt.dataF.th;

        Fd_v=linspace(Fd_min,Fd_max,256);
        Fq_v=linspace(Fq_min,Fq_max,256);
        th_v=linspace(min(th_dqt,[],'all'),max(th_dqt,[],'all'),256);
        [Fd_dqt,Fq_dqt,th_dqt]=meshgrid(Fd_v,Fq_v,th_v);

        Id_dqt=interpn(motorModel.FluxMapInv_dqt.dataF.Fd,motorModel.FluxMapInv_dqt.dataF.Fq,motorModel.FluxMapInv_dqt.dataF.th,motorModel.FluxMapInv_dqt.dataF.Id,Fd_dqt,Fq_dqt,th_dqt,'cubic');
        Iq_dqt=interpn(motorModel.FluxMapInv_dqt.dataF.Fd,motorModel.FluxMapInv_dqt.dataF.Fq,motorModel.FluxMapInv_dqt.dataF.th,motorModel.FluxMapInv_dqt.dataF.Iq,Fd_dqt,Fq_dqt,th_dqt,'cubic');
        T_dqt=interpn(motorModel.FluxMapInv_dqt.dataF.Fd,motorModel.FluxMapInv_dqt.dataF.Fq,motorModel.FluxMapInv_dqt.dataF.th,motorModel.FluxMapInv_dqt.dataF.T,Fd_dqt,Fq_dqt,th_dqt,'cubic');

end

%% ----------------SimScape FEM-based PMSM----------------

dqtMap = motorModel.FluxMap_dqt;

fInt.Fd = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th/motorModel.data.p,dqtMap.data.Fd,'spline','none');
fInt.Fq = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th/motorModel.data.p,dqtMap.data.Fq,'spline','none');
fInt.T  = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th/motorModel.data.p,dqtMap.data.T,'spline','none');

IdMax = max(dqtMap.data.Id(:));
IdMin = min(dqtMap.data.Id(:));
IqMax = max(dqtMap.data.Iq(:));
IqMin = min(dqtMap.data.Iq(:));
thMax = max(dqtMap.data.th(:));
thMin = min(dqtMap.data.th(:));
IdStp = abs(dqtMap.Id(2)-dqtMap.Id(1));
IqStp = abs(dqtMap.Iq(2)-dqtMap.Iq(1));
thStp = 360-thMax;

IdMax = 0.95*max([IdMax -IdMin]);
IqMax = 0.95*max([IqMax -IqMin]);


eMotor.Id = linspace(-IdMax,IdMax,31);
eMotor.Iq = linspace(-IqMax,IqMax,31);
eMotor.th = thMin:thStp/motorModel.data.p:(120/motorModel.data.p)-thStp;

[data.Id,data.Iq,data.th] = ndgrid(eMotor.Id,eMotor.Iq,eMotor.th);
data.Fd = fInt.Fd(data.Id,data.Iq,data.th);
data.Fq = fInt.Fq(data.Id,data.Iq,data.th);
data.T  = fInt.T(data.Id,data.Iq,data.th);

if strcmp(motorModel.data.motorType,'SR')
    % symmetr on both axis
    data.Fd(data.Id<IdMin) = -fInt.Fd(-data.Id(data.Id<IdMin),+data.Iq(data.Id<IdMin),+data.th(data.Id<IdMin));
    data.Fq(data.Id<IdMin) = +fInt.Fq(-data.Id(data.Id<IdMin),+data.Iq(data.Id<IdMin),+data.th(data.Id<IdMin));
    data.T(data.Id<IdMin)  = -fInt.T(-data.Id(data.Id<IdMin),+data.Iq(data.Id<IdMin),+data.th(data.Id<IdMin));
    data.Fd(data.Iq<IqMin) = +fInt.Fd(+data.Id(data.Iq<IqMin),-data.Iq(data.Id<IqMin),+data.th(data.Iq<IqMin));
    data.Fq(data.Iq<IqMin) = -fInt.Fq(+data.Id(data.Iq<IqMin),-data.Iq(data.Id<IqMin),+data.th(data.Iq<IqMin));
    data.T(data.Iq<IqMin)  = -fInt.T(+data.Id(data.Iq<IqMin),-data.Iq(data.Iq<IqMin),+data.th(data.Iq<IqMin));
else
    if strcmp(motorModel.data.axisType,'SR')
        % symmetry just along d axis
        data.Fd(data.Id<IdMin) = -fInt.Fd(-data.Id(data.Id<IdMin),data.Iq(data.Id<IdMin),data.th(data.Id<IdMin));
        data.Fq(data.Id<IdMin) = +fInt.Fq(-data.Id(data.Id<IdMin),data.Iq(data.Id<IdMin),data.th(data.Id<IdMin));
        data.T(data.Id<IdMin)  = -fInt.T(-data.Id(data.Id<IdMin),data.Iq(data.Id<IdMin),data.th(data.Id<IdMin));
        data.Fd(data.Iq<IdMin) = NaN;
        data.Fq(data.Iq<IdMin) = NaN;
        data.T(data.Iq<IdMin)  = NaN;
    else
        % symmetry just along q axis
        data.Fd(data.Id<IdMin) = NaN;
        data.Fq(data.Id<IdMin) = NaN;
        data.T(data.Id<IdMin)  = NaN;
        data.Fd(data.Iq<IqMin) = +fInt.Fd(+data.Id(data.Iq<IqMin),-data.Iq(data.Iq<IqMin),+data.th(data.Iq<IqMin));
        data.Fq(data.Iq<IqMin) = -fInt.Fq(+data.Id(data.Iq<IqMin),-data.Iq(data.Iq<IqMin),+data.th(data.Iq<IqMin));
        data.T(data.Iq<IqMin)  = -fInt.T(+data.Id(data.Iq<IqMin),-data.Iq(data.Iq<IqMin),+data.th(data.Iq<IqMin));
    end
end

eMotor.Fd = data.Fd;
eMotor.Fq = data.Fq;
eMotor.T  = data.T;

eMotor.th(end+1) = 120/motorModel.data.p;
eMotor.Fd(:,:,end+1) = eMotor.Fd(:,:,1);
eMotor.Fq(:,:,end+1) = eMotor.Fq(:,:,1);
eMotor.T(:,:,end+1)  = eMotor.T(:,:,1);

eMotor.p  = motorModel.data.p;
eMotor.Rs = motorModel.data.Rs;
eMotor.L0 = 1e-12;
eMotor.J  = motorModel.data.J;

% Iron loss
switch motorModel.SyreDrive.modelSetup.IronLoss
    case 'No'
        set_param([Slx_name '/Motor model/SimScape FEM-based PMSM/SimScape PMSM'],'loss_param','ee.enum.fem_motor.ironLossesExtended.none')
    case 'Yes'
        set_param([Slx_name '/Motor model/SimScape FEM-based PMSM/SimScape PMSM'],'loss_param','ee.enum.fem_motor.ironLossesExtended.tabulated3D')
        ironLoss.n = linspace(0,motorModel.data.nmax,31);
        [ironLossData.Id,ironLossData.Iq,ironLossData.n] = ndgrid(eMotor.Id,eMotor.Iq,ironLoss.n);
        if strcmp(motorModel.data.motorType,'SR')
            % symmetr on both axis
            fdfqTmp.Id = abs(ironLossData.Id(:,:,1));
            fdfqTmp.Iq = abs(ironLossData.Iq(:,:,1));
        else
            if strcmp(motorModel.data.axisType,'SR')
                % symmetry just along d axis
                fdfqTmp.Id = abs(ironLossData.Id(:,:,1));
                fdfqTmp.Iq = ironLossData.Iq(:,:,1);
            else
                % symmetry just along q axis
                fdfqTmp.Id = ironLossData.Id(:,:,1);
                fdfqTmp.Iq = abs(ironLossData.Iq(:,:,1));
            end
        end
        fdfqTmp.Fd = nan(size(fdfqTmp.Id));
        fdfqTmp.Fq = nan(size(fdfqTmp.Id));
        ironLossData.Pfes = nan(size(ironLossData.Id));
        ironLossData.Pfer = nan(size(ironLossData.Id));
        for ii=1:length(ironLoss.n)
            f = ironLoss.n(ii)/60*motorModel.data.p;
            [~,Pfesh,Pfesc,Pferh,Pferc,~] = calcIronLoss(motorModel.IronPMLossMap_dq,fdfqTmp,f);
            ironLossData.Pfes(:,:,ii) = Pfesh+Pfesc;
            ironLossData.Pfer(:,:,ii) = Pferh+Pferc;
        end
        eMotor.ironLoss.Pfes = ironLossData.Pfes;
        eMotor.ironLoss.Pfer = ironLossData.Pfer;
        eMotor.ironLoss.n = ironLoss.n;
end




%% -----------------Iron Loss Model-------------------------------%%

switch motorModel.SyreDrive.modelSetup.IronLoss
    case 'No'
        IronLoss = 0;
    case 'Yes'
        IronLoss = 1;
        Pfe_h = motorModel.IronPMLossMap_dq.Pfes_h + motorModel.IronPMLossMap_dq.Pfer_h;
        Pfe_c = motorModel.IronPMLossMap_dq.Pfes_c + motorModel.IronPMLossMap_dq.Pfer_c;
        Ppm   = motorModel.IronPMLossMap_dq.Ppm;
        n0    = motorModel.IronPMLossMap_dq.n0;
        expH  = motorModel.IronPMLossMap_dq.expH;
        expC  = motorModel.IronPMLossMap_dq.expC;
        expPM = motorModel.IronPMLossMap_dq.expPM;

        Id_fe = motorModel.IronPMLossMap_dq.Id;
        Iq_fe = motorModel.IronPMLossMap_dq.Iq;
end  

%% -----------------AC Loss Model-------------------------------%%

switch motorModel.SyreDrive.modelSetup.WindingLossAC
    case 'No'
        ACloss.enable = 0;
        ACloss.n = linspace(-motorModel.data.nmax*1.2,motorModel.data.nmax*1.2,3);
        ACloss.f = abs(ACloss.n)/60*motorModel.data.p;
        [ACloss.Rac,ACloss.kAC] = calcRsTempFreq(motorModel.data.Rs,motorModel.data.tempCu,motorModel.data.l,motorModel.data.lend,[],'0',motorModel.data.tempCu,ACloss.f);
        ACloss.Rac = ACloss.Rac*ones(size(ACloss.f));
        ACloss.kAC = ones(size(ACloss.f));
    case 'Yes'
        ACloss.enable = 1;
        ACloss.n = linspace(-motorModel.data.nmax*1.2,motorModel.data.nmax*1.2,501);
        ACloss.f = abs(ACloss.n)/60*motorModel.data.p;
        [ACloss.Rac,ACloss.kAC] = calcRsTempFreq(motorModel.data.Rs,motorModel.data.tempCu,motorModel.data.l,motorModel.data.lend,motorModel.acLossFactor,'LUT',motorModel.data.tempCu,ACloss.f);
end

%% --------------------User Settings-------------------------%

% Ctrl settings
switch motorModel.SyreDrive.modelSetup.Ctrl_type
    case 'Current control'
        Ctrl_type = 0;
    case 'Torque control'
        Ctrl_type = 2;
    case 'Speed control'
        Ctrl_type = 3;
end

% Ctrl Strategy 
switch motorModel.SyreDrive.modelSetup.Ctrl_strategy
    case 'FOC'
        Ctrl_strategy = 0;
    case 'DFVC'
        Ctrl_strategy = 1;
end


% Sensorless on or off
switch motorModel.SyreDrive.SS_on
    case 'Off'
        SS_on = 0;
    case 'On'
        SS_on = 1;
end

% Injected waveform
switch motorModel.SyreDrive.SS_settings.inj_waveform
    case 'Sinusoidal'
        inj_waveform = 0;
    case 'Squarewave'
        inj_waveform = 1;
end 

% Demodulation technique
switch motorModel.SyreDrive.SS_settings.dem
    case 'Current'
        dem = 0;
    case 'Flux'
        dem = 1;
end

% High speed position error estimation technique
switch motorModel.SyreDrive.SS_settings.HS_ctrl
    case 'Active Flux'
        HS_ctrl = 0;
    case 'APP'
        HS_ctrl = 1;
end


 %% coordinate transformations
%3-ph to 2-ph
Clarke = 2/3*[1 -0.5 -0.5;0 sqrt(3)/2 -sqrt(3)/2];
%2-ph to 3-ph
Clarke_inv = [1 0;-0.5 sqrt(3)/2;-0.5 -sqrt(3)/2];
%% Inverter model



set_param([Slx_name '/Inverter Model/Converter (Three-Phase)'],'diode_Vf','0.8');
set_param([Slx_name '/Inverter Model/Converter (Three-Phase)'],'diode_Ron','1e-4');
set_param([Slx_name '/Inverter Model/Converter (Three-Phase)'],'diode_Goff','1e-5');

switch(inverterModelType)
    
    case 'Istantaneous'
        InvModel = 1;
        set_param([Slx_name '/Inverter Model/Solver Configuration'],'UseLocalSolver','off');
        set_param([Slx_name '/Inverter Model/Converter (Three-Phase)'],'device_type','ee.enum.converters.switchingdevice.ideal');
        set_param([Slx_name '/Inverter Model/Converter (Three-Phase)'],'Ron','Rd');
        set_param([Slx_name '/Inverter Model/Converter (Three-Phase)'],'Goff','1e-4');
        set_param([Slx_name '/Inverter Model/Converter (Three-Phase)'],'Vth','0.001');
        
    case 'Average'
        InvModel = 0;
        set_param([Slx_name '/Inverter Model/Solver Configuration'],'UseLocalSolver','on');
        set_param([Slx_name '/Inverter Model/Solver Configuration'],'LocalSolverSampleTime','Ts/10');
        set_param([Slx_name '/Inverter Model/Converter (Three-Phase)'],'device_type','ee.enum.converters.switchingdevice.averaged');
        set_param([Slx_name '/Inverter Model/Converter (Three-Phase)'],'Ron','Rd');
end
