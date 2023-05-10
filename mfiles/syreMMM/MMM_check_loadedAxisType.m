% Copyright 2022
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

function [motorModel] = MMM_check_loadedAxisType(motorModel,dataSet,fdfq,dqtMap,ironLoss)

if ~isfield(dataSet,'axisType')
   if strcmp(dataSet.TypeOfRotor,'SPM')||strcmp(dataSet.TypeOfRotor,'Vtype')
       dataSet.axisType = 'PM';
   else
       dataSet.axisType = 'SR';
   end
end

if ~strcmp(dataSet.axisType,motorModel.data.axisType)
%    if ~isempty(fdfq)
%        tmp.FluxMap_dq = fdfq;
%    else
%        tmp.Id = nan(2);
%        tmp.Iq = nan(2);
%        tmp.Fd = nan(2);
%        tmp.Fq = nan(2);
%    end
   tmp.FluxMap_dq       = fdfq;
   tmp.FluxMap_dqt      = dqtMap;
   tmp.IronPMLossMap_dq = ironLoss;
   if strcmp(dataSet.axisType,'SR')
       tmp = MMM_sr2pm(tmp);
       disp('Axis changed from SR to PM')
   else
       tmp = MMM_pm2sr(tmp);
       disp('Axis changed from PM to SR')
   end
   if ~isempty(fdfq)
       motorModel.FluxMap_dq = tmp.FluxMap_dq;
   end
   if ~isempty(dqtMap)
       motorModel.FluxMap_dqt = tmp.FluxMap_dqt;
   end
   if ~isempty(ironLoss)
       motorModel.IronPMLossMap_dq = tmp.IronPMLossMap_dq;
   end
else
    if ~isempty(fdfq)
       motorModel.FluxMap_dq = fdfq;
   end
   if ~isempty(dqtMap)
       motorModel.FluxMap_dqt = dqtMap;
   end
   if ~isempty(ironLoss)
       motorModel.IronPMLossMap_dq = ironLoss;
   end
end