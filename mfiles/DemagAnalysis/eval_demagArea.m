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

function [DemagArea] = eval_demagArea(dataIn)

load([dataIn.currentpathname dataIn.currentfilename])

RatedCurrent = dataIn.RatedCurrent;
CurrLoPP = dataIn.CurrLoPP(1);
% SimulatedCurrent = dataIn.SimulatedCurrent;
% GammaPP  = dataIn.GammaPP;
BrPP = dataIn.BrPP;
% NumOfRotPosPP = dataIn.NumOfRotPosPP;
% AngularSpanPP = dataIn.AngularSpanPP;
% NumGrid = dataIn.NumGrid;

per.EvalSpeed = dataIn.EvalSpeed;

clc;

per.overload = CurrLoPP;
per.i0       = RatedCurrent;
per.BrPP     = BrPP;
per.tempPP   = dataIn.tempPP(1);


Br = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,per.tempPP);
per.BrPP = Br;
mat.LayerMag.Hc = Br/(4*pi*1e-7*mat.LayerMag.mu);

pathname = dataIn.currentpathname;
filename = dataIn.currentfilename;

outFolder = [filename(1:end-4) '_results\FEA results\'];
if ~exist([pathname outFolder],'dir')
    mkdir([pathname outFolder]);
end

resFolder = ['demagArea_' int2str(per.tempPP) 'deg_' int2str(per.i0*per.overload) 'Amp\'];
mkdir([pathname outFolder resFolder]);
resFolder = [pathname outFolder resFolder];

eval_type    = 'demagArea';

if (strcmp(geo.RotType,'SPM')||strcmp(geo.RotType,'Vtype'))
    per.gamma = 180;
else
    per.gamma = 90;
end

nsim = 2;

% Bd = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Bd,tempVect(tt));

% per.BrPP = Br;


openfemm(1)
opendocument([pathname filename(1:end-4) '.fem'])
mi_saveas([resFolder filename(1:end-4) '.fem']);
mi_close;

dPM  = zeros(1,nsim);
Bmin = zeros(1,nsim);
xyC  = cell(1,nsim);
xyV  = cell(1,nsim);
xyB  = cell(1,nsim);
SOL  = cell(1,nsim);


for pp=1:nsim
    [~,tmpFolder]=createTempDir();
    
    copyfile([resFolder filename(1:end-4) '.fem'],[tmpFolder filename(1:end-4) '.fem']);
    
    per.nsim_singt      = 1;
    per.delta_sim_singt = (0.5*360/(6*geo.q*geo.win.n3phase))*(pp-1);
    
    SOL{pp} = simulate_xdeg(geo,per,mat,eval_type,tmpFolder,[filename(1:end-4) '.fem']);
    
    
    dPM(pp)  = SOL{pp}.dPM;
    Bmin(pp) = SOL{pp}.Bmin;
    xyC{pp}  = SOL{pp}.xyDemagTmpC;
    xyV{pp}  = SOL{pp}.xyDemagTmpV;
    xyB{pp}  = SOL{pp}.xyDemagTmpB;
    
    
end

DemagArea.current     = per.i0*per.overload;
DemagArea.temperature = per.tempPP;
DemagArea.dPM         = dPM;
DemagArea.Bmin        = Bmin;
DemagArea.SOL         = SOL;
DemagArea.xyC         = xyC;
DemagArea.xyV         = xyV;
DemagArea.xyB         = xyB;


save([resFolder 'demagAreaResults.mat'],'DemagArea','geo','per','mat','dataSet');

figure()
figSetting(20,10)
for ii=1:nsim
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
    
    title([num2str(dPM(ii)*100) '\% of the PMs demagnetized'])
    
    if ~isempty(xyV{ii})
        vTmp = [xyV{ii};xyV{ii}(1,:)];
        cTmp = xyC{ii};
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





