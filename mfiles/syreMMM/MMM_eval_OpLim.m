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

function [Plim,resPath] = MMM_eval_OpLim(motorModel,saveFlag)

nCurr = motorModel.data.nCurr;
Imax  = motorModel.data.Imax;
i0    = motorModel.data.i0;
Vdc   = motorModel.data.Vdc;

% if length(nCurr)==1
%     if nCurr==1
%         Ivect = i0;
%     else
%         Ivect = (1/nCurr:1/nCurr:1)*Imax;
%     end
% else
%     Ivect = nCurr*i0;
% end

Ivect = nCurr*i0;

Plim = cell(size(Ivect));

for ii=1:length(Plim)
    Plim{ii} = OpLimEval(motorModel,Ivect(ii),Vdc);
end

[hfig,resPath] = OpLimPlot(Plim,Ivect,motorModel);

%% Save figures
if nargin()==1
    answer = 'No';
    answer = questdlg('Save figures?','Save','Yes','No',answer);
else
    if saveFlag
        answer = 'Yes';
    else
        answer = 'No';
    end
end

if strcmp(answer,'Yes')
    if ~exist(resPath,'dir')
        mkdir(resPath);
    end
    
    for ii=1:length(hfig)
        savePrintFigure(hfig(ii));
    end
    
    save([resPath 'OpLimResults.mat'],'Plim');
    
end



