warning off;
InverterModel = 'Switching';    
%InverterModel = 'Average';     

load motorModel;

% S-fun Parameters
Tstep = 2e-6;
Ts    = 1/motorModel.SyreDrive.Converter.fPWM;
n_set = motorModel.data.n3phase;
n_phase = n_set*3;

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

id_MTPA = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.id,i0);
iq_MTPA = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.iq,i0);
Ld_inic = interp2(motorModel.IncInductanceMap_dq.Id,motorModel.IncInductanceMap_dq.Iq,motorModel.IncInductanceMap_dq.Ldd,id_MTPA,iq_MTPA);
Lq_inic = interp2(motorModel.IncInductanceMap_dq.Id,motorModel.IncInductanceMap_dq.Iq,motorModel.IncInductanceMap_dq.Lqq,id_MTPA,iq_MTPA);
L_sigma = Lq_inic*0.1;

T0(~isnan(motorModel.data.T0)) = motorModel.data.T0;
T0(isempty(T0)) = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.T,i0);
n0(~isnan(motorModel.data.n0)) = motorModel.data.n0;
n0(isempty(n0)) = 1000;
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
switch motorModel.SyreDrive.FMapsModel
    case 'dq Model'
        FMapsModel = 1;
         %------------------dq Inverse Flux Maps------------------------%  
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


%% -----------------Iron Loss Model-------------------------------%%

switch motorModel.SyreDrive.IronLoss
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


 %% coordinate transformations
T = motorModel.data.T_VSD;
T_inv = inv(T);

%% Inverter model

Slx_name = [motorModel.data.motorName '_ctrl_INST'];

switch(InverterModel)
    
    case 'Switching'
        InvModel = 1;
        for i=1:n_set
            Inverter_Name = ['Inverter Model ' num2str(i)];
            set_param([Slx_name '/' Inverter_Name '/Solver Configuration'],'UseLocalSolver','off');
            set_param([Slx_name '/' Inverter_Name '/Converter (Three-Phase)'],'device_type','ee.enum.converters.switchingdevice.ideal');
            set_param([Slx_name '/' Inverter_Name '/Converter (Three-Phase)'],'Ron','Rd');
            set_param([Slx_name '/' Inverter_Name '/Converter (Three-Phase)'],'Goff','1e-4');
            set_param([Slx_name '/' Inverter_Name '/Converter (Three-Phase)'],'Vth','0.001');            
            set_param([Slx_name '/' Inverter_Name '/Converter (Three-Phase)'],'diode_Vf','0.8');
            set_param([Slx_name '/' Inverter_Name '/Converter (Three-Phase)'],'diode_Ron','1e-4');
            set_param([Slx_name '/' Inverter_Name '/Converter (Three-Phase)'],'diode_Goff','1e-5');
        
        end

    case 'Average'
        InvModel = 0;
        for i=1:n_set
            Inverter_Name = ['Inverter Model ' num2str(i)];
            set_param([Slx_name '/' Inverter_Name '/Solver Configuration'],'UseLocalSolver','on');
            set_param([Slx_name '/' Inverter_Name '/Solver Configuration'],'LocalSolverSampleTime','Ts/10');
            set_param([Slx_name '/' Inverter_Name '/Converter (Three-Phase)'],'device_type','ee.enum.converters.switchingdevice.averaged');
            set_param([Slx_name '/' Inverter_Name '/Converter (Three-Phase)'],'Ron','Rd'); 
            set_param([Slx_name '/' Inverter_Name '/Converter (Three-Phase)'],'diode_Vf','0.8');
            set_param([Slx_name '/' Inverter_Name '/Converter (Three-Phase)'],'diode_Ron','1e-4');
            set_param([Slx_name '/' Inverter_Name '/Converter (Three-Phase)'],'diode_Goff','1e-5');
        end
        
end
