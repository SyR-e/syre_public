% Copyright 2023
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

function [forceOut] = elab_singt_toothForce(geo,per,out,pathname)
% 
% if nargin==0
%     [filename,pathname] = uigetfile([cd '\*.mat'],'Load results from SyR-e');
% end
% 
% load([pathname filename]);

% spatial symmetry (pole-->full rotor simulation)
nRep = 2*geo.p/geo.ps;
nSlot = (6*geo.p*geo.q*geo.win.n3phase);
rF = (geo.r+geo.g);
l  = geo.l;
thRot = out.SOL.th-out.SOL.th(1);
xdeg = per.delta_sim_singt;

alphaSlot = 360/nSlot;
alphaTooth = alphaSlot*(1-geo.acs);

thTooth = ((0:1:nSlot-1)*alphaSlot)';

if nRep~=1
    thGapStep = out.SOL.Ft(2,1)-out.SOL.Ft(1,1);
    thGap = (thGapStep/2:thGapStep:360-thGapStep/2)';
    pGapR = repmat(out.SOL.Fr(:,2:end),[nRep,1]);
    pGapT = repmat(out.SOL.Ft(:,2:end),[nRep,1]);
else
    thGapStep = out.SOL.Ft(2,1)-out.SOL.Ft(1,1);
    thGap = out.SOL.Ft(:,1);
    pGapR = out.SOL.Fr(:,2:end);
    pGapT = out.SOL.Ft(:,2:end);
end

% time symmetry (half period --> full period)
if xdeg==360
    disp('Simulation done on 360 elt deg rotation. OK')
elseif xdeg==180
    disp('Simulation done on 180 elt deg rotation. Symmetry applied')
    thRot = [thRot thRot+180];
    pGapR = [pGapR pGapR];
    pGapT = [pGapT pGapT];
else
    disp(['Simulation done on ' int2str(xdeg) ' elt deg rotation. Symmetry not possible']);
end

% tooth force computation

FrTooth = zeros(length(thTooth),size(pGapR,2));
FtTooth = zeros(length(thTooth),size(pGapR,2));

for ii=1:length(thTooth)
    thMin = thTooth(ii)-alphaTooth/2;
    % if thMin<0
    %     thMin = thMin+360;
    % end
    thMax = thTooth(ii)+alphaTooth/2;
    % if thMax>360
    %     thMax=thMax-360;
    % end
    thTmp = repmat(thGap,[1,size(pGapR,2)]);
    thTmp = [thTmp-360;thTmp;thTmp+360];
    pRtmp = [pGapR;pGapR;pGapR];
    pRtmp(thTmp<thMin) = NaN;
    pRtmp(thTmp>thMax) = NaN;
    FrTooth(ii,:) = mean(pRtmp,'omitmissing')*(alphaTooth*pi/180*rF/1000*l/1000); % [Nm]
    pTtmp = [pGapT;pGapT;pGapT];
    pTtmp(thTmp<thMin) = NaN;
    pTtmp(thTmp>thMax) = NaN;
    FtTooth(ii,:) = mean(pTtmp,'omitmissing')*(alphaTooth*pi/180*rF/1000*l/1000); % [Nm]
end


forceOut.thGap   = thGap;
forceOut.pGapR   = pGapR;
forceOut.pGapT   = pGapT;
forceOut.thTooth = thTooth;
forceOut.FrTooth = FrTooth;
forceOut.FtTooth = FtTooth;
forceOut.thRot   = thRot;

save([pathname 'toothForceOut.mat'],'geo','per','out','forceOut');

% figures
thGapPlot = repmat(thGap,1,size(pGapR,2));
thRotPlot = repmat(thRot,size(pGapR,1),1);

hfig(1) = figure();
figSetting()
%view(3)
xlabel('$\theta_{gap}$ (mech deg)')
ylabel('$\theta_{rot}$ (elt deg)')
title('$p_r$ (Pa)')
%scatter3(thGapPlot(:),thRotPlot(:),pGapR(:),[],pGapR(:),'filled','MarkerEdgeColor','none');
contourf(thGapPlot,thRotPlot,pGapR,'LineStyle','none');
colorbar()
set(gcf,'FileName',[pathname 'radial_pressure.fig'])

hfig(2) = figure();
figSetting()
%view(3)
xlabel('$\theta_{gap}$ (mech deg)')
ylabel('$\theta_{rot}$ (elt deg)')
title('$p_t$ (Pa)')
%scatter3(thGapPlot(:),thRotPlot(:),pGapT(:),[],pGapT(:),'filled','MarkerEdgeColor','none');
contourf(thGapPlot,thRotPlot,pGapT,'LineStyle','none');
colorbar()
set(gcf,'FileName',[pathname 'tangential_pressure.fig'])

nSlotPlot = repmat(1:1:nSlot,1,size(FrTooth,2));
thRotPlot = repmat(thRot,size(FrTooth,1),1);

hfig(3) = figure();
figSetting()
view(3)
xlabel('slot number')
ylabel('$\theta_{rot}$ (elt deg)')
zlabel('$F_r$ (Nm)')
scatter3(nSlotPlot(:),thRotPlot(:),FrTooth(:),[],FrTooth(:),'filled','MarkerEdgeColor','none');
set(gcf,'FileName',[pathname 'radial_tooth_force.fig'])

hfig(4) = figure();
figSetting()
view(3)
xlabel('slot number')
ylabel('$\theta_{rot}$ (elt deg)')
zlabel('$F_t$ (Nm)')
scatter3(nSlotPlot(:),thRotPlot(:),FtTooth(:),[],FtTooth(:),'filled','MarkerEdgeColor','none');
set(gcf,'FileName',[pathname 'tangential_tooth_force.fig'])

for ii=1:length(hfig)
    savePrintFigure(hfig(ii));
end

if nargout()==0
    clear forceOut;
end

