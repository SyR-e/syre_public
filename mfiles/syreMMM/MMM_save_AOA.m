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

function MMM_save_AOA(motorModel)

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'AOA - ' int2str(motorModel.data.tempPM) 'deg\'];

if ~exist([pathname resFolder],'dir')
    mkdir([pathname resFolder]);
end

AOA  = motorModel.AOA;
MTPA = motorModel.AOA.MTPA;
MTPV = motorModel.AOA.MTPV;
save([pathname resFolder 'AOAcurves.mat'],'AOA')
save([pathname resFolder 'MTPAcurve.mat'],'MTPA')
save([pathname resFolder 'MTPVcurve.mat'],'MTPV')

id_KtMax    = MTPA.id;
iq_KtMax    = MTPA.iq;
T_KtMax     = MTPA.T;
fd_KtMax    = MTPA.fd;
fq_KtMax    = MTPA.fq;
dTpp_KtMax  = MTPA.dTpp;
I_KtMax     = abs(id_KtMax+j*iq_KtMax);
gamma_KtMax = angle(id_KtMax+j*iq_KtMax)*180/pi;
F_KtMax     = abs(fd_KtMax+j*fq_KtMax);
delta_KtMax = angle(fd_KtMax+j*fq_KtMax)*180/pi;

id_KvMax    = MTPV.id;
iq_KvMax    = MTPV.iq;
T_KvMax     = MTPV.T;
fd_KvMax    = MTPV.fd;
fq_KvMax    = MTPV.fq;
F_KvMax     = abs(fd_KvMax+j*fq_KvMax);
delta_KvMax = angle(fd_KvMax+j*fq_KvMax)*180/pi;

save([pathname resFolder 'ktMax_idiq.mat'],'id_KtMax','iq_KtMax','I_KtMax','T_KtMax','dTpp_KtMax','gamma_KtMax');
save([pathname resFolder 'ktMax_FdFq.mat'],'fd_KtMax','fq_KtMax','F_KtMax','delta_KtMax');
save([pathname resFolder 'kvMax_idiq.mat'],'id_KvMax','iq_KvMax','F_KvMax','T_KvMax');
save([pathname resFolder 'kvMax_FdFq.mat'],'fd_KvMax','fq_KvMax','F_KvMax','delta_KvMax');











