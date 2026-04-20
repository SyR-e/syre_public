% Copyright 2025
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

function plot_singt_toothForce(geo,per,out,forceOut,pathname)
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
if geo.acs~=1
    alphaTooth = alphaSlot*(1-geo.acs);
else
    alphaTooth = alphaSlot*0.5;
end

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
scatter3(nSlotPlot(:),thRotPlot(:),forceOut.FrTooth(:),[],forceOut.FrTooth(:),'filled','MarkerEdgeColor','none');
set(gcf,'FileName',[pathname 'radial_tooth_force.fig'])

hfig(4) = figure();
figSetting()
view(3)
xlabel('slot number')
ylabel('$\theta_{rot}$ (elt deg)')
zlabel('$F_t$ (Nm)')
scatter3(nSlotPlot(:),thRotPlot(:),forceOut.FtTooth(:),[],forceOut.FtTooth(:),'filled','MarkerEdgeColor','none');
set(gcf,'FileName',[pathname 'tangential_tooth_force.fig'])

for ii=1:length(hfig)
    savePrintFigure(hfig(ii));
end

if nargout()==0
    clear forceOut;
end

