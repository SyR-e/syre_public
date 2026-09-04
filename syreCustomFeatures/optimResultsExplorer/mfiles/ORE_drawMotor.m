% Copyright 2026
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

function [dataSet] = ORE_drawMotor(optRes,motorID,flagSave)

if nargin()==1
    prompt = {'Enter motorID:'};
    name = 'Motor ID';
    default = {num2str(round(optRes.motNum(1)))};
    answer = inputdlg(prompt,name,1,default);
    motorID = eval(answer{1});
    flagSave = 1;
elseif nargin()==2
    flagSave = 1;
end



geo = optRes.geo;
% per = optRes.per;
mat = optRes.mat;

dataSet = optRes.dataSet;

fem = dimMesh(geo,'singt');

% motorID = 1256;
% warning('motorID set by script!!!')

index = 1:1:numel(optRes.motNum);
index = index(optRes.motNum==motorID);

[geo,gamma,mat] = interpretRQ(optRes.varParetoFront(index,:),geo,mat);
geo.x0 = geo.r/cos(pi/2/geo.p);
% geo.x0 = (geo.r-geo.hs)/cos(pi/2/geo.p);

[rotor,~,geo,mat] = ROTmatr(geo,fem,mat); % rotor and BLKLABELSrot describe the rotor
geo.rotor = rotor;

[geo,stator,~] = STATmatr(geo,fem); % statore and BLKLABELSstat describe the stator
geo.stator=stator;

dataSet.AirGapThickness = round(geo.g.*100)./100;
dataSet.AirGapRadius    = round(geo.r.*100)./100;
dataSet.ToothLength     = round(geo.lt.*100)./100;
dataSet.StatorSlotOpen  = round(geo.acs.*100)./100;
dataSet.ToothWidth      = round(geo.wt.*100)./100;
dataSet.ToothTangDepth  = round(geo.ttd.*100)./100;
dataSet.Br              = round(mat.LayerMag.Br(1).*100)./100;
dataSet.ALPHApu         = round(geo.dalpha_pu.*100)./100;
dataSet.HCpu            = round(geo.hc_pu.*100)./100;
dataSet.DepthOfBarrier  = round(geo.dx.*100)./100;
dataSet.betaPMshape     = round(geo.betaPMshape.*100)./100;
dataSet.thetaFBS        = round(geo.th_FBS*180/pi*100)/100;
dataSet.TanRibEdit      = round(geo.pontT.*100)./100;
dataSet.RadRibEdit      = round(geo.pontR.*100)./100;
dataSet.RotorFilletIn   = round(geo.RotorFillet1.*100)./100;
dataSet.RotorFilletOut  = round(geo.RotorFillet2.*100)./100;
dataSet.RotorFilletTan1 = round(geo.RotorFilletTan1.*100)./100;
dataSet.RotorFilletTan2 = round(geo.RotorFilletTan2.*100)./100;
dataSet.PMdim           = round(geo.PMdim.*100)./100;
dataSet.DepthOfBarrier  = round(geo.dx.*100)./100;
dataSet.CentralShrink   = round(geo.hcShrink.*100)./100;
dataSet.RadShiftInner   = round(geo.dxIB.*100)./100;
dataSet.betaPMshape     = round(geo.betaPMshape.*100)./100;
dataSet.gammaPP         = round(gamma.*100)./100;

dataSet.currentfilename = ['opt' int2str(motorID) '.mat'];
dataSet.currentpathname = checkPathSyntax([cd '\']);


figure();
figSetting();
title(['motorID: ' int2str(motorID)]);
tmp.dataSet = dataSet;
tmp.AxisGeometry = gca;
tmp = GUI_APP_DrawMachine(tmp);

if flagSave
    button = questdlg('Open the motor in the main SyR-e GUI?','Select','Yes','No','Yes');
    if strcmp(button,'Yes')
        hApp = findall(0,'Name','GUI_Syre');
        if isempty(hApp)
            GUI_Syre(dataSet);
        else
            hApp.RunningAppInstance.dataSet = dataSet;
            hApp.RunningAppInstance.GUI_update;
            figure(hApp)
            clear hApp
        end
    end
end

if nargout()==0
    clear dataSet
end

