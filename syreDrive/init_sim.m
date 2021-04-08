
clear
load motorModel.mat

warning off;

% S-fun Parameters
Tstep = 2e-6;
Ts    = 100e-6;

% Machine  Parameters
VDC     = motorModel.data.Vdc;
Rs      = motorModel.data.Rs;
P       = motorModel.data.p;
B       = 0;
J       = motorModel.data.J;
accel   = 5000; % rpm/s

MTPA    = motorModel.AOA.MTPA;
i0      = motorModel.data.i0;
id_MTPA = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.id,i0);
iq_MTPA = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.iq,i0);
Ld_inic = interp2(motorModel.Inductance.Id,motorModel.Inductance.Iq,motorModel.Inductance.Ldd,id_MTPA,iq_MTPA);
Lq_inic = interp2(motorModel.Inductance.Id,motorModel.Inductance.Iq,motorModel.Inductance.Lqq,id_MTPA,iq_MTPA);
T0(~isnan(motorModel.data.T0)) = motorModel.data.T0;
T0(isempty(T0)) = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.T,i0);
clear MTPA

% Converter
V0 = motorModel.SyreDrive.Converter.V0;             % power semiconductors ON treshold [V]
Rd = motorModel.SyreDrive.Converter.Rd;             % power semiconductors incremental resistance [Ohm]
dT = motorModel.SyreDrive.Converter.dT * 1e-6;      % dead time [s]
      
%----------------dqt Inverse Flux Maps ------------------------%

Fd_max = max(motorModel.dqtMapF.dataF.Fd,[],'all');
Fq_max = max(motorModel.dqtMapF.dataF.Fq,[],'all');
Fd_min = min(motorModel.dqtMapF.dataF.Fd,[],'all');
Fq_min = min(motorModel.dqtMapF.dataF.Fq,[],'all');
th_min = min(motorModel.dqtMapF.dataF.th,[],'all');
th_max = max(motorModel.dqtMapF.dataF.th,[],'all');
th_dqt = motorModel.dqtMapF.dataF.th;

Fd_v=linspace(Fd_min,Fd_max,256);
Fq_v=linspace(Fq_min,Fq_max,256);
th_v=linspace(min(th_dqt,[],'all'),max(th_dqt,[],'all'),256);
[Fd_dqt,Fq_dqt,th_dqt]=meshgrid(Fd_v,Fq_v,th_v);

Id_dqt=interpn(motorModel.dqtMapF.dataF.Fd,motorModel.dqtMapF.dataF.Fq,motorModel.dqtMapF.dataF.th,motorModel.dqtMapF.dataF.Id,Fd_dqt,Fq_dqt,th_dqt,'cubic');
Iq_dqt=interpn(motorModel.dqtMapF.dataF.Fd,motorModel.dqtMapF.dataF.Fq,motorModel.dqtMapF.dataF.th,motorModel.dqtMapF.dataF.Iq,Fd_dqt,Fq_dqt,th_dqt,'cubic');
T_dqt=interpn(motorModel.dqtMapF.dataF.Fd,motorModel.dqtMapF.dataF.Fq,motorModel.dqtMapF.dataF.th,motorModel.dqtMapF.dataF.T,Fd_dqt,Fq_dqt,th_dqt,'cubic');

%------------------dq Inverse Flux Maps------------------------%

Fd     = motorModel.idiq.Fd;
Fq     = motorModel.idiq.Fq;
Id     = motorModel.idiq.Id;
Iq     = motorModel.idiq.Iq;
T      = motorModel.idiq.T;

%--------------------User Settings-------------------------%

% Ctrl settings
switch motorModel.SyreDrive.Ctrl_type
    case 'Current control'
        Ctrl_type = 0;
    case 'Torque control'
        Ctrl_type = 2;
    case 'Speed control'
        Ctrl_type = 3;
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

% dq or dqt model
switch motorModel.SyreDrive.FMapsModel
    case 'dq Model'
        FMapsModel = 1;
        
    case 'dqt Model'
        FMapsModel = -1;
end

 