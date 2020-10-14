% Copyright 2020
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MMM_print_AOA(motorModel)

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'C code LUTs\'];

if ~exist([pathname resFolder],'dir')
    mkdir([pathname resFolder]);
end

MTPA = motorModel.AOA.MTPA;
MTPV = motorModel.AOA.MTPV;

% LUT
m = 1;  % # of lines
n = 20; % # of columns (table size is 1 x n)

Tmax = MTPA.T(end);
step = Tmax/n;
T_set = 0:step:Tmax;

id_set    = interp1(MTPA.T,MTPA.id,T_set);
id_set(1) = 0;
iq_set    = interp1(MTPA.T,MTPA.iq,T_set);
iq_set(1) = 0;
fd_set    = interp1(MTPA.T,MTPA.fd,T_set);
fq_set    = interp1(MTPA.T,MTPA.fq,T_set);
f_set     = abs(fd_set+j*fq_set);


% print txt file (MTPA)
fid = fopen([pathname resFolder 'tablesMTPA.txt'],'w');
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

if not(isempty(MTPV.iq))
    
    FMTPV = linspace(min(abs(MTPV.fd+j*MTPV.fq)),max(abs(MTPV.fd+j*MTPV.fq)),n+1);
    TMTPV = interp1(abs(MTPV.fd+j*MTPV.fq),MTPV.T,FMTPV);
    
    fprintf(fid,' \n');
    fprintf(fid,'float FMIN    = %4.3f; //Vs\n',min(FMTPV));
    fprintf(fid,'float FMAX    = %4.3f; //Vs\n',max(FMTPV));
    step = (max(FMTPV)-min(FMTPV))/n;
    fprintf(fid,'float DF      = %4.4f; //Vs\n',step);
    fprintf(fid,'float INV_DF  = %4.4f; //Vs^-1\n',1/step);
    StampaVarg(fid,TMTPV,m,n+1,'T_MTPV','//MTPV - max torque vs Vs','%6.3f')
    
else
    FMTPV = [];
    TMTPV = [];
end
fclose(fid);

hfig(1) = figure();
figSetting();
hax(1) = axes('OuterPosition',[0 0 1 1]);
xlabel('$T_{set}$ [$Nm$]')
ylabel('$i_{dq}$ [$A$]')
set(hfig(1),'FileName',[pathname resFolder 'MTPA current LUTs.fig'])
plot(MTPA.T,MTPA.id,'-b','DisplayName','$i_d$')
plot(MTPA.T,MTPA.iq,'-r','DisplayName','$i_q$')
plot(T_set,id_set,'kx','DisplayName','LUT')
plot(T_set,iq_set,'kx','DisplayName','LUT')
hleg = legend('show','Location','northwest','Orientation','horizontal');
set(hleg,'NumColumns',2);

hfig(2) = figure();
figSetting();
hax(2) = axes('OuterPosition',[0 0 1 1]);
xlabel('$T_{set}$ [$Nm$]')
ylabel('$\lambda_{dq}$ [$A$]')
set(hfig(2),'FileName',[pathname resFolder 'MTPA flux linkages LUTs.fig'])
plot(MTPA.T,MTPA.fd,'-b','DisplayName','$\lambda_d$')
plot(MTPA.T,MTPA.fq,'-r','DisplayName','$\lambda_q$')
plot(MTPA.T,abs(MTPA.fd+j*MTPA.fq),'-','Color',[0 0.8 0],'DisplayName','$|\lambda|$ [$Vs$]')
plot(T_set,fd_set,'kx','DisplayName','LUT')
plot(T_set,fq_set,'kx','DisplayName','LUT')
plot(T_set,f_set,'kx','DisplayName','LUT')
hleg = legend('show','Location','northwest','Orientation','horizontal');
set(hleg,'NumColumns',3);

hfig(3) = figure();
figSetting();
hax(3) = axes('OuterPosition',[0 0 1 1]);
xlabel('$\lambda_{set}$ [$Nm$]')
ylabel('$T_{set}$ [$A$]')
set(hfig(3),'FileName',[pathname resFolder 'MTPV LUTs.fig'])
plot(abs(MTPV.fd+j*MTPV.fq),MTPV.T,'-b','DisplayName','curve')
plot(FMTPV,TMTPV,'kx','DisplayName','LUT');
hleg = legend('show','Location','northwest');

for ii=1:length(hfig)
    savePrintFigure(hfig(ii))
end



