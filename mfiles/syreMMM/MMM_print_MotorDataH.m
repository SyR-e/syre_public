% Copyright 2020
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

function MMM_print_MotorDataH(motorModel)
    MotorDataH_path = [motorModel.data.pathname motorModel.data.motorName '_ctrl\User_functions\Inc\MotorData.h'];
    
    MTPA = motorModel.AOA.MTPA;
    MTPV = motorModel.AOA.MTPV;
    MTPA.T(isnan(MTPA.T)) = 0;
    
    i0       = motorModel.data.i0;
    Imax_mot = motorModel.data.Imax;
    nmax_mot = motorModel.data.nmax;
    T_rated(~isnan(motorModel.data.T0)) = motorModel.data.T0;
    T_rated(isempty(T_rated)) = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.T,i0);
    Tmax_mot = 2 * T_rated;
    RS       = motorModel.data.Rs;
    PP       = motorModel.data.p;
    J        = motorModel.data.J;
    
    deadtime = motorModel.SyreDrive.Converter.dT;
    
    id_MTPA  = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.id,i0);
    iq_MTPA  = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.iq,i0);
    Ld_inic  = interp2(motorModel.Inductance.Id,motorModel.Inductance.Iq,motorModel.Inductance.Ldd,id_MTPA,iq_MTPA);
    Lq_inic  = interp2(motorModel.Inductance.Id,motorModel.Inductance.Iq,motorModel.Inductance.Lqq,id_MTPA,iq_MTPA);
    ld_inic  = interp2(motorModel.Inductance.Id,motorModel.Inductance.Iq,motorModel.Inductance.Ldd,id_MTPA,iq_MTPA);
    lq_inic  = interp2(motorModel.Inductance.Id,motorModel.Inductance.Iq,motorModel.Inductance.Lqq,id_MTPA,iq_MTPA);
    ldq_inic = interp2(motorModel.Inductance.Id,motorModel.Inductance.Iq,motorModel.Inductance.Ldq,id_MTPA,iq_MTPA);
 
    fid = fopen(MotorDataH_path,'w');
        fprintf(fid,'#define I_rated            %4.2f\n',i0);
        fprintf(fid,'#define T_rated            %4.2f\n',T_rated);
        fprintf(fid,'#define Tmax_mot           %4.2f\n',Tmax_mot);
        fprintf(fid,'#define Imax_mot           %4.2f\n',Imax_mot);
        fprintf(fid,'#define nmax_mot           %d\n',nmax_mot);
        fprintf(fid,'#define RS                 %4.4f\n',RS);
        fprintf(fid,'#define PP                 %d // pole pairs\n',PP);
        fprintf(fid,'#define ONE_P              1/PP\n');
        fprintf(fid,'#define J                  %4.4f\n',J);
        fprintf(fid,'#define ENCODER_RESOLUTION 2048 //1024 //512\n');
        fprintf(fid,' \n');

        fprintf(fid,'#define deadtime %4.2fe-6\n',deadtime);
        
        fprintf(fid,' \n');
        fprintf(fid,'#define Ld_inic %4.3f\n',Ld_inic);
        fprintf(fid,'#define Lq_inic %4.3f\n',Lq_inic);
        fprintf(fid,'#define ld_inic %4.3f\n',ld_inic);
        fprintf(fid,'#define lq_inic %4.3f\n',lq_inic);
        fprintf(fid,'#define ldq_inic %4.3f\n',ldq_inic);
        fprintf(fid,' \n');
    fclose(fid);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % LUT
    m = 1;  % # of lines
    n = 20; % # of columns (table size is 1 x n)

    Tmax = MTPA.T(end);
    step = Tmax/n;
    T_set = 0:step:Tmax;

    id_set    = interp1(MTPA.T,MTPA.id,T_set);
    id_set(1) = 0;
    iq_set    = interp1(MTPA.T,MTPA.iq,T_set);
    iq_set(1) = 0.2*i0;
    fd_set    = interp1(MTPA.T,MTPA.fd,T_set);
    fd_set(1) = 0;
    fq_set    = interp1(MTPA.T,MTPA.fq,T_set);
    fq_set(1) = 0;
    f_set     = abs(fd_set+j*fq_set);
    f_set(1)  = 0;

    % print txt file (MTPA)
    fid = fopen(MotorDataH_path,'a');
        % fprintf(fid,'//SIGLA MOTORE: %s\n',motor_name);
        fprintf(fid,['//' date '\n']);
        fprintf(fid,'float TMIN    = 0;\n');
        fprintf(fid,'float TMAX    = %4.3f; //Nm\n',Tmax);
        fprintf(fid,'float DT      = %4.4f; //Nm\n',step);
        fprintf(fid,'float INV_DT  = %4.4f; //Nm^-1\n',1/step);

        StampaVarg(fid,id_set,m,n+1,'ID_REF','//MTPA - id','%6.3f')
        StampaVarg(fid,iq_set,m,n+1,'IQ_REF','//MTPA - iq','%6.3f')
        StampaVarg(fid,fd_set,m,n+1,'FD_REF','//MTPA - fd','%6.3f')
        StampaVarg(fid,fq_set,m,n+1,'FQ_REF','//MTPA - fq','%6.3f')
        StampaVarg(fid,f_set,m,n+1,'F_REF','//MTPA - flux amplitude','%6.3f')
        fprintf(fid,' \n');

        FMTPV = linspace(min(abs(MTPV.fd+1i*MTPV.fq)),max(abs(MTPV.fd+1i*MTPV.fq)),n+1);
        TMTPV = interp1(abs(MTPV.fd+1i*MTPV.fq),MTPV.T,FMTPV);


        fprintf(fid,'float FMIN    = %4.3f; //Vs\n',min(FMTPV));
        fprintf(fid,'float FMAX    = %4.3f; //Vs\n',max(FMTPV));
        step = (max(FMTPV)-min(FMTPV))/n;
        fprintf(fid,'float DF      = %4.4f; //Vs\n',step);
        fprintf(fid,'float INV_DF  = %4.4f; //Vs^-1\n',1/step);
        StampaVarg(fid,TMTPV,m,n+1,'T_MTPV','//MTPV - max torque vs Vs','%6.3f')
        fprintf(fid,' \n');
    fclose(fid);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Id = motorModel.fdfq.Id;
    Iq = motorModel.fdfq.Iq;
    Fd = motorModel.fdfq.Fd;
    Fq = motorModel.fdfq.Fq;

    % LUTs for flux observer
    id_tab_min = min(min(Id));
    id_tab_max = max(max(Id));
    iq_tab_min = min(min(Iq));
    iq_tab_max = max(max(Iq));

    % LUT dimension
    m = 51;    % rows
    n = 51;   % columns

    % Fd table: current steps
    Didd = (id_tab_max-id_tab_min)/(n-1);
    Diqd = (iq_tab_max-iq_tab_min)/(m-1);
    % Fq table: current steps
    Diqq = (iq_tab_max-iq_tab_min)/(n-1);
    Didq = (id_tab_max-id_tab_min)/(m-1);

    [idd,iqd]=meshgrid(linspace(id_tab_min,id_tab_max,n),linspace(iq_tab_min,iq_tab_max,m));
    [idq,iqq]=meshgrid(linspace(id_tab_min,id_tab_max,m),linspace(iq_tab_min,iq_tab_max,n));

    fd=interp2(Id,Iq,Fd,idd,iqd);
    fq=interp2(Id,Iq,Fq,idq,iqq);
    fq=fq';

    % print to FluxTables.txt
    fid = fopen(MotorDataH_path,'a');
        fprintf(fid,['//' date '\n']);
        fprintf(fid,['float  ID_TAB_MIN = ' num2str(id_tab_min) ' ;\r\n']);
        fprintf(fid,['float  IQ_TAB_MIN = ' num2str(iq_tab_min) ' ;\r\n']);
        fprintf(fid,['float  ID_TAB_MAX = ' num2str(id_tab_max) ' ;\r\n']);
        fprintf(fid,['float  IQ_TAB_MAX = ' num2str(iq_tab_max) ' ;\r\n']);
        fprintf(fid,['float  DIDD       = ' num2str(Didd,4) ' ;\r\n']);
        fprintf(fid,['float  DIQD       = ' num2str(Diqd,4) ' ;\r\n']);
        fprintf(fid,['float  DIQQ       = ' num2str(Diqq,4) ' ;\r\n']);
        fprintf(fid,['float  DIDQ       = ' num2str(Didq,4) ' ;\r\n']);
        fprintf(fid,['float  INV_DIDD   = ' num2str(1/Didd,4) ' ;\r\n']);
        fprintf(fid,['float  INV_DIQD   = ' num2str(1/Diqd,4) ' ;\r\n']);
        fprintf(fid,['float  INV_DIQQ   = ' num2str(1/Diqq,4) ' ;\r\n']);
        fprintf(fid,['float  INV_DIDQ   = ' num2str(1/Didq,4) ' ;\r\n']);
        fprintf(fid,['float  n_size     = ' num2str(m) ' ;\r\n']);
        % fd(:,1)=0*fd(:,1);  % Fd @ id=0
        StampaVarg(fid,fd',m,n,'FD_LUT','//Fluxd(iq,id)','%6.4f')
        StampaVarg(fid,fq',m,n,'FQ_LUT','//Fluxq(id,iq)','%6.4f')
    fclose(fid);
end