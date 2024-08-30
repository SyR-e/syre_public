% Copyright 2021
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

function [dataSet] = SDE_saveMotor(map,saveFlag)

dataSet = map.dataSet;
if numel(map.xx)==1
    map.xSelect = map.xx;
    map.bSelect = map.bb;
    
    dataSet.AirGapRadius    = map.xx*dataSet.StatorOuterRadius;
    dataSet.ShaftRadius     = map.Ar;
    dataSet.ToothLength     = map.lt;
    dataSet.ToothWidth      = map.wt;
    dataSet.GammaPP         = map.gamma;
    dataSet.ThermalLoadKj   = map.kj;
    dataSet.CurrentDensity  = map.J;
    dataSet.AdmiJouleLosses = dataSet.ThermalLoadKj*(2*pi*dataSet.StatorOuterRadius/1000*dataSet.StackLength/1000);
    dataSet.TurnsInSeries   = map.geo.win.Ns;

    dataSet.PMdim = -dataSet.PMdimPU./dataSet.PMdimPU;
    dataSet.PMdim(isnan(dataSet.PMdim)) = 0;
    dataSet.PMdim = dataSet.kPM*dataSet.PMdim;

    dataSet.HCpu = map.hc_pu{1,1};
    dataSet.DepthOfBarrier = map.dx{1,1};
else
    x = map.xSelect;
    b = map.bSelect;
    if dataSet.FEAfixN==0
        dataSet.FEAfixN = 1;
    end

    dataSet.AirGapRadius    = x*dataSet.StatorOuterRadius;
    dataSet.ShaftRadius     = interp2(map.xx,map.bb,map.Ar,x,b);
    dataSet.ToothLength     = interp2(map.xx,map.bb,map.lt,x,b);
    dataSet.ToothWidth      = interp2(map.xx,map.bb,map.wt,x,b);
    dataSet.GammaPP         = interp2(map.xx,map.bb,map.gamma,x,b);
    dataSet.ThermalLoadKj   = interp2(map.xx,map.bb,map.kj,x,b);
    dataSet.CurrentDensity  = interp2(map.xx,map.bb,map.J,x,b);
    dataSet.AdmiJouleLosses = dataSet.ThermalLoadKj*(2*pi*dataSet.StatorOuterRadius/1000*dataSet.StackLength/1000);
    if isfield(map,'Ns')
        dataSet.TurnsInSeries = interp2(map.xx,map.bb,map.Ns,x,b);
    else
        dataSet.TurnsInSeries   = map.geo.win.Ns;
    end

    dataSet.PMdim = -dataSet.PMdimPU./dataSet.PMdimPU;
    dataSet.PMdim(isnan(dataSet.PMdim)) = 0;
    dataSet.PMdim = dataSet.kPM*dataSet.PMdim;

    hc_pu = zeros(1,dataSet.NumOfLayers);
    dx    = zeros(1,dataSet.NumOfLayers);
    [m,n]=size(map.xx);
    hcTmp=zeros(m,n);
    dxTmp=zeros(m,n);
    for ii=1:dataSet.NumOfLayers
        for mm=1:m
            for nn=1:n
                hcTmp(mm,nn)=map.hc_pu{mm,nn}(ii);
                dxTmp(mm,nn)=map.dx{mm,nn}(ii);
            end
        end
        hc_pu(ii)=interp2(map.xx,map.bb,hcTmp,x,b);
        dx(ii)=interp2(map.xx,map.bb,dxTmp,x,b);
    end

    dataSet.HCpu = round(hc_pu*100)/100;
    dataSet.DepthOfBarrier = round(dx*100)/100;
end

%Disable Optimization Check
dataSet.PMdimBouCheck = 0;

figure();
figSetting();
title(['$x=' num2str(map.xSelect,2) '$ / $b=' num2str(map.bSelect,2) '$'])


[dataSet,~,~,~] = back_compatibility(dataSet,[],[],0);
tmp.dataSet = dataSet;
tmp.AxisGeometry = gca;

tmp = GUI_APP_DrawMachine(tmp);

if saveFlag
    dataSet = DrawAndSaveMachine(tmp.dataSet);
else
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





