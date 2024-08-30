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

function MMM_save_ironLoss(motorModel)

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'Iron Loss Model - ' int2str(motorModel.data.tempPM) 'deg\'];

if ~exist([pathname resFolder],'dir')
    mkdir([pathname resFolder])
end

fdfq     = motorModel.FluxMap_dq;
ironLoss = motorModel.IronPMLossMap_dq;

Id            = ironLoss.Id;
Iq            = ironLoss.Iq;
Fd            = interp2(fdfq.Id,fdfq.Iq,fdfq.Fd,Id,Iq);
Fq            = interp2(fdfq.Id,fdfq.Iq,fdfq.Fq,Id,Iq);
T             = interp2(fdfq.Id,fdfq.Iq,fdfq.T,Id,Iq);
dT            = interp2(fdfq.Id,fdfq.Iq,fdfq.dT,Id,Iq);
dTpp          = interp2(fdfq.Id,fdfq.Iq,fdfq.dTpp,Id,Iq);
Pfes_h        = ironLoss.Pfes_h;
Pfes_c        = ironLoss.Pfes_c;
Pfer_h        = ironLoss.Pfer_h;
Pfer_c        = ironLoss.Pfer_c;
Ppm           = ironLoss.Ppm;
Pfe           = Pfes_h+Pfes_c+Pfer_h+Pfer_c;
velDim        = ironLoss.n0;
factors.expC  = ironLoss.expC;
factors.expH  = ironLoss.expH;
factors.expPM = ironLoss.expPM;
factors.f0    = ironLoss.f0;

per.tempPP = motorModel.data.tempPM;
dataSet.axisType = motorModel.data.axisType;

save([pathname resFolder 'fdfq_idiq_n256_ironLoss.mat'],...
    'Id','Iq','Fd','Fq','T','dT','dTpp',...
    'Pfes_h','Pfes_c','Pfer_h','Pfer_c','Ppm','Pfe','velDim','factors','per','dataSet');



