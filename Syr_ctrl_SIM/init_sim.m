load motorModel.mat

% S-fun Parameters
Tstep = 2e-6;
Ts    = 200e-6;

% Machine  Parameters
VDC     = motorModel.data.Vdc;
Rs      = motorModel.data.Rs;
P       = motorModel.data.p;
MTPA    = motorModel.AOA.MTPA;
i0      = motorModel.dataSet.RatedCurrent;
id_MTPA = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.id,i0);
iq_MTPA = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.iq,i0);
Ld_inic = interp2(motorModel.Inductance.Id,motorModel.Inductance.Iq,motorModel.Inductance.Ldd,id_MTPA,iq_MTPA);
Lq_inic = interp2(motorModel.Inductance.Id,motorModel.Inductance.Iq,motorModel.Inductance.Lqq,id_MTPA,iq_MTPA);
B       = 0;
J       = motorModel.data.J;
accel   = 1000; % rpm/s

% Converter
V0 = motorModel.SyreDrive.Converter.V0;             % power semiconductors ON treshold [V]
Rd = motorModel.SyreDrive.Converter.Rd;             % power semiconductors incremental resistance [Ohm]
dT = motorModel.SyreDrive.Converter.dT;             % dead time [s]
      
motorModel = MMM_dqtMapF_rescale(motorModel);

Fd_dqt = motorModel.dqtMapF.dataF.Fd;
Fq_dqt = motorModel.dqtMapF.dataF.Fq;
Id_dqt = motorModel.dqtMapF.dataF.Id;
Iq_dqt = motorModel.dqtMapF.dataF.Iq;
T_dqt  = motorModel.dqtMapF.dataF.T;
th_dqt = motorModel.dqtMapF.dataF.th;
        
Fd_max = max(Fd_dqt,[],'all');
Fq_max = max(Fq_dqt,[],'all');
Fd_min = min(Fd_dqt,[],'all');
Fq_min = min(Fq_dqt,[],'all');

Fd     = motorModel.idiq.Fd;
Fq     = motorModel.idiq.Fq;
Id     = motorModel.idiq.Id;
Iq     = motorModel.idiq.Iq;
T      = motorModel.idiq.T;
        
% Ctrl settings
switch motorModel.SyreDrive.Ctrl_type
    case 'Current control'
        Ctrl_type = 0;
    case 'Torque control'
        Ctrl_type = 2;
    case 'Speed control'
        Ctrl_type = 3;
end

% dq or dqt model
switch motorModel.SyreDrive.FMapsModel
    case 'dq Model'
        FMapsModel = 1;
        
    case 'dqt Model'
        FMapsModel = -1;
end

 