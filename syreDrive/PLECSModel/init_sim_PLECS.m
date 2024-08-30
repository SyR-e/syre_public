
clear
warning off;

load motorModel


% S-fun Parameters
Tstep = 2e-6;
Ts    = 100e-6;


%% ----------------Machine and Converter  Parameters------------------------%
VDC     = motorModel.data.Vdc;
Rs      = motorModel.data.Rs;
p       = motorModel.data.p;
J       = motorModel.data.J;
Rfe     = 1e5;
accel   = 5000; % rpm/s
Bm      = 0.;                     % damping constant (Nm/(rad/sec))
Tf      = 0;                      % friction loss (Nm)
Tv      = 0;%P0/(nmax*pi/30)^3;   % ventilation loss coefficient (W/(rad/s)^3 = Nm/(rad/s)^2)

MTPA    = motorModel.controlTrajectories.MTPA;
i0      = motorModel.data.i0;

id_MTPA = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.id,i0);
iq_MTPA = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.iq,i0);

Ld_inic = interp2(motorModel.IncInductanceMap_dq.Id,motorModel.IncInductanceMap_dq.Iq,motorModel.IncInductanceMap_dq.Ldd,id_MTPA,iq_MTPA);
Lq_inic = interp2(motorModel.IncInductanceMap_dq.Id,motorModel.IncInductanceMap_dq.Iq,motorModel.IncInductanceMap_dq.Lqq,id_MTPA,iq_MTPA);
T0(~isnan(motorModel.data.T0)) = motorModel.data.T0;
T0(isempty(T0)) = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.T,i0);
n0(~isnan(motorModel.data.n0)) = motorModel.data.n0;
n0(isempty(n0)) = 1000;
clear MTPA

% Converter
V0 = motorModel.SyreDrive.Converter.V0;             % power semiconductors ON treshold [V]
Rd = motorModel.SyreDrive.Converter.Rd;             % power semiconductors incremental resistance [Ohm]
dT = motorModel.SyreDrive.Converter.dT * 1e-6;      % dead time [s]

switch motorModel.SyreDrive.modelType
    case 'Istantaneous'
        InverterModel = 1;
    case 'Average'
        InverterModel = 0;
end

      
%% ----------------PM Initialization and Convention------------------------%
Fd0 = interp2(motorModel.FluxMap_dq.Id,motorModel.FluxMap_dq.Iq,motorModel.FluxMap_dq.Fd,0,0);
Fq0 = interp2(motorModel.FluxMap_dq.Id,motorModel.FluxMap_dq.Iq,motorModel.FluxMap_dq.Fq,0,0);
th0 = 0; %initial electrical angle
F0_alpha = Fd0*cos(th0) - Fq0*sin(th0);
F0_beta  = Fd0*sin(th0) + Fq0*cos(th0);

if motorModel.data.axisType == 'SR' & motorModel.data.motorType == 'SR'
    Quad_Maps = 0; %SyR Convention - 1st quadrant maps
elseif motorModel.data.axisType == 'SR' & motorModel.data.motorType == 'PM'
    Quad_Maps = 1; %PM-SyR - 1st and 4th quadrant maps
elseif motorModel.data.axisType == 'PM' & motorModel.data.motorType == 'PM'
    Quad_Maps = 2; %IPM - 1st and 2st quadrant maps
end
%% ---------------Magnets Estimation----------------------------%

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

clear FluxMap_dq
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


%% ---------------- Inverse Flux Maps ------------------------%
 
switch motorModel.SyreDrive.FMapsModel
    case 'dq Model'
        FMapsModel = 1;
    case 'dqt Model'
        FMapsModel = -1;
end       
% 'dq Model'
        
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
        
%  'dqt Model'
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


%% -----------------Iron Loss Model-------------------------------%%

switch motorModel.SyreDrive.IronLoss
    case 'No'
        IronLoss = 0;        
    case 'Yes'
        IronLoss = 1;
        
end  

Pfe_h = motorModel.IronPMLossMap_dq.Pfes_h + motorModel.IronPMLossMap_dq.Pfer_h;
Pfe_c = motorModel.IronPMLossMap_dq.Pfes_c + motorModel.IronPMLossMap_dq.Pfer_c;
Ppm   = motorModel.IronPMLossMap_dq.Ppm;
n0    = motorModel.IronPMLossMap_dq.n0;
expH  = motorModel.IronPMLossMap_dq.expH;
expC  = motorModel.IronPMLossMap_dq.expC;
expPM = motorModel.IronPMLossMap_dq.expPM;

Id_fe = motorModel.IronPMLossMap_dq.Id;
Iq_fe = motorModel.IronPMLossMap_dq.Iq;

id_fe=Id_fe(1,:);
iq_fe=Iq_fe(:,1);
save('SimMatFiles\IronLosses','id_fe','iq_fe','Pfe_c','Ppm','Pfe_h','expH','expC','expPM','IronLoss');


%% --------------------User Settings-------------------------%

% Ctrl settings
switch motorModel.SyreDrive.Ctrl_type
    case 'Current control'
        Ctrl_type = 0;
    case 'Torque control'
        Ctrl_type = 2;
    case 'Speed control'
        Ctrl_type = 3;
end

% Ctrl Strategy 
switch motorModel.SyreDrive.Ctrl_strategy
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



%% Save MatFiles

save('SimMatFiles\Parameters','Tstep','Ts','VDC','Rs','p','J','Rfe','accel','Bm','Tf','Tv','Rd','dT','T0','n0','Fm','th0','Clarke','Clarke_inv');
switch(Quad_Maps)
    case {0,2}
        save('SimMatFiles\InvFluxMapdq','fD_pu_norm','fQ_pu_norm','iD_pu_norm','iQ_pu_norm','fQ_pu_norm','fD_vct_ref','fQ_vct_max');
    case 1
        save('SimMatFiles\InvFluxMapdq','fD_pu_norm','fQ_pu_norm','iD_pu_norm','iQ_pu_norm','fQ_pu_norm','fQ_vct_ref','fD_vct_max');
end
save('SimMatFiles\InvFluxMapdqt','Fd_max','Fq_max','Fd_min','Fq_min','th_min','th_max','th_dqt','Fd_v','Fq_v','th_v','Id_dqt','Iq_dqt','T_dqt');
save('SimMatFiles\UserSettings','Ctrl_type','Ctrl_strategy','SS_on','inj_waveform','dem','HS_ctrl','FMapsModel','InverterModel','Quad_Maps','id_MTPA','iq_MTPA','Ld_inic','Lq_inic','InitIntg_d','InitIntg_q','IronLoss');


 

