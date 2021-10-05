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

function MMM_print_fdfq(motorModel)

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'C code LUTs\'];

if ~exist([pathname resFolder],'dir')
    mkdir([pathname resFolder]);
end

fdfq = motorModel.FluxMap_dq;
Id = fdfq.Id;
Iq = fdfq.Iq;
Fd = fdfq.Fd;
Fq = fdfq.Fq;

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
fid = fopen([pathname resFolder 'FluxTables.txt'],'w');
fprintf(fid,['//' date '\n']);
fprintf(fid,['float  ID_TAB_MIN = ' num2str(id_tab_min) ' ;\r\n']);
fprintf(fid,['float  IQ_TAB_MIN = ' num2str(iq_tab_min) ' ;\r\n']);
fprintf(fid,['float  ID_TAB_MAX = ' num2str(id_tab_max) ' ;\r\n']);
fprintf(fid,['float  IQ_TAB_MAX = ' num2str(iq_tab_max) ' ;\r\n']);
fprintf(fid,['float  DIDD = ' num2str(Didd,4) ' ;\r\n']);
fprintf(fid,['float  DIQD = ' num2str(Diqd,4) ' ;\r\n']);
fprintf(fid,['float  DIQQ = ' num2str(Diqq,4) ' ;\r\n']);
fprintf(fid,['float  DIDQ = ' num2str(Didq,4) ' ;\r\n']);
fprintf(fid,['float  INV_DIDD = ' num2str(1/Didd,4) ' ;\r\n']);
fprintf(fid,['float  INV_DIQD = ' num2str(1/Diqd,4) ' ;\r\n']);
fprintf(fid,['float  INV_DIQQ = ' num2str(1/Diqq,4) ' ;\r\n']);
fprintf(fid,['float  INV_DIDQ = ' num2str(1/Didq,4) ' ;\r\n']);
% fd(:,1)=0*fd(:,1);  % Fd @ id=0
StampaVarg(fid,fd',m,n,'FD_LUT','//Fluxd(iq,id)','%6.4f')
StampaVarg(fid,fq',m,n,'FQ_LUT','//Fluxq(id,iq)','%6.4f')
fclose(fid);

hfig = figure();
figSetting()
xlabel('$i_{dq}$ [$A$]')
ylabel('$\lambda_{dq}$ [$Vs$]')
title('Magnetic Model LUT')
set(hfig,'FileName',[pathname resFolder 'Magnetic Model LUTs.fig'])
plot(idd',fd')
plot(iqq ,fq')
plot(idd',fd','kx')
plot(iqq ,fq','kx')

savePrintFigure(hfig)



