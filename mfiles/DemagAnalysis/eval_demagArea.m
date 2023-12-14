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

function [DemagArea,resFolder] = eval_demagArea(dataIn)

load([dataIn.currentpathname dataIn.currentfilename])

RatedCurrent = dataIn.RatedCurrent;
CurrLoPP = dataIn.CurrLoPP(1);
% SimulatedCurrent = dataIn.SimulatedCurrent;
GammaPP  = dataIn.GammaPP;
BrPP = dataIn.BrPP;
% NumOfRotPosPP = dataIn.NumOfRotPosPP;
% AngularSpanPP = dataIn.AngularSpanPP;
% NumGrid = dataIn.NumGrid;
tempPP = dataIn.tempPP(1);

per.EvalSpeed = dataIn.EvalSpeed;

clc;

if ~isfield(geo,'axisType')
    if strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype')
        geo.axisType = 'PM';
    else
        geo.axisType = 'SR';
    end
end

if ~strcmp(geo.axisType,dataIn.axisType)
    %geo.axisType = dataIn.axisType;
    if strcmp(dataIn.axisType,'PM')
        geo.th0 = geo.th0 - 90;
    else
        geo.th0 = geo.th0 + 90;
    end
end


per.overload = CurrLoPP;
per.i0       = RatedCurrent;
per.BrPP     = BrPP;
per.tempPP   = tempPP;
per.gamma    = GammaPP;


Br = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,per.tempPP);
per.BrPP = Br;
mat.LayerMag.Hc = Br/(4*pi*1e-7*mat.LayerMag.mu);

pathname = dataIn.currentpathname;
filename = dataIn.currentfilename;

outFolder = [filename(1:end-4) '_results\FEA results\'];
if ~exist([pathname outFolder],'dir')
    mkdir([pathname outFolder]);
end

resFolder = ['demagArea_' int2str(per.tempPP) 'degC_' int2str(per.i0*per.overload) 'Amp_' int2str(per.gamma) 'degGamma' dataIn.axisType '\'];
mkdir([pathname outFolder resFolder]);
resFolder = [pathname outFolder resFolder];

eval_type    = 'demagArea';

% if strcmp(dataIn.axisType,'PM')
%     per.gamma = 180;
% else
%     per.gamma = 90;
% end

[~,~,~,out,~] = FEMMfitness([],geo,per,mat,eval_type,[pathname filename(1:end-4) '.fem']);



DemagArea.current     = per.i0*per.overload;
DemagArea.temperature = per.tempPP;
DemagArea.dPM         = out.SOL.dPM;
DemagArea.Bmin        = out.SOL.Bmin;
DemagArea.SOL         = out.SOL;
DemagArea.xyC         = out.SOL.xyDemagTmpC;
DemagArea.xyV         = out.SOL.xyDemagTmpV;
DemagArea.xyB         = out.SOL.xyDemagTmpB;


save([resFolder 'demagAreaResults.mat'],'DemagArea','geo','per','mat','dataSet');

figure()
figSetting(20,10)
for ii=1:2
    hax(ii) = axes(...
        'OuterPosition',[0+0.5*(ii-1) 0 0.5 1],...
        'XLim',[-1 geo.r+1],...
        'YLim',[-1 geo.r+1],...
        'DataAspectRatio',[1 1 1],...
        'XTick',[],...
        'YTick',[]);
    GUI_Plot_Machine(hax(ii),geo.rotor);
    hchild = get(hax(ii),'Children');
    for jj=1:length(hchild)
        set(hchild(jj),'Color','k','LineWidth',1.5);
    end
    
    title([num2str(DemagArea.dPM(ii)*100) '\% of the PMs demagnetized'])
    
    if ~isempty(DemagArea.xyV{ii})
        vTmp = [DemagArea.xyV{ii};DemagArea.xyV{ii}(1,:)];
        cTmp = DemagArea.xyC{ii};
    else
        vTmp = [0];
        cTmp = [0];
    end
    
    for jj=1:length(vTmp(1,:))
        fill(hax(ii),real(vTmp(:,jj)),imag(vTmp(:,jj)),'r');
        %         plot(hax(ii),real(cTmp(jj)),imag(cTmp(jj)),'r.');
    end
end

saveas(gcf,[resFolder 'Demagnetized_PM_Area.fig']);





