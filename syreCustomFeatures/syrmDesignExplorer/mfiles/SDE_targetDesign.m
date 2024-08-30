% Copyright 2022
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

function [hfig] = SDE_targetDesign(map,setup)

if isfield(map,'xSelect')
    map = rmfield(map,'xSelect');
    map = rmfield(map,'bSelect');
    map = rmfield(map,'dataAvailable');
    map = rmfield(map,'dataSelect');
end


if nargin()==1
    if exist('calc_hairpinTurns.m','file')
        Ns = calc_hairpinTurns(map.dataSet);
        Ns = unique(Ns);
        Ns = Ns';
    else
        Ns = map.dataSet.TurnsInSeries;
    end

    prompt  = {'Peak phase current (Apk)','DC link voltage (V)','Target torque (Nm)','Target power factor','Feasible number of turns','Target power (W)'};
    answers = {num2str(map.dataSet.RatedCurrent),int2str(map.Vdc),'100','NaN',mat2str(Ns),'100000'};
    answers = inputdlg(prompt,'Design Targets',1,answers);

    setup.Imax = eval(answers{1});
    setup.Vdc  = eval(answers{2});
    setup.T    = eval(answers{3});
    setup.PF   = eval(answers{4});
    setup.Ns   = eval(answers{5});
    setup.P    = eval(answers{6});
end

setup.w    = setup.P/setup.T;
setup.n    = setup.w*30/pi;
setup.Fmax = setup.Vdc/(sqrt(3)*map.dataSet.NumOfPolePairs*setup.w);
if isnan(setup.PF)
    setup.PF   = setup.P/(sqrt(3)/2*setup.Vdc*setup.Imax*map.dataSet.Num3PhaseCircuit);
end
nPoints = 1001;

xx = linspace(min(map.xx(:)),max(map.xx(:)),nPoints);
bb = linspace(min(map.bb(:)),max(map.bb(:)),nPoints);
[xx,bb] = meshgrid(xx,bb);

feasible = zeros(size(xx));

hfig(1) = figure();
figSetting(14,14);
hax(1) = axes(...
    'OuterPosition',[0 0 1 1],...
    'XLim',[min(xx(:)) max(xx(:))],...
    'YLim',[min(bb(:)) max(bb(:))]);
xlabel('$x$')
ylabel('$b$')
title('$(x,b)$ Design Plane')
hleg(1) = legend('show','Location','southwest');
set(hfig(1),'UserData',map);

colors = get(hax,'ColorOrder');
colors = [colors;colors;colors;colors];
colors = [colors;colors;colors;colors];

contour(map.xx,map.bb,map.T,'-','EdgeColor',colors(2,:),'LineWidth',1,'ShowText','on','DisplayName','$T$ [Nm]');
contour(map.xx,map.bb,map.PF,'-','EdgeColor',colors(1,:),'LineWidth',1,'ShowText','on','DisplayName','$cos \varphi$');
contour(map.xx,map.bb,map.T,'-','EdgeColor',colors(2,:),'LineWidth',2,'ShowText','on','DisplayName',['$T=' num2str(setup.T) '$ Nm'],'LevelList',setup.T);
contour(map.xx,map.bb,map.PF,'-','EdgeColor',colors(1,:),'LineWidth',2,'ShowText','on','DisplayName',['$cos \varphi=' num2str(setup.PF) '$'],'LevelList',setup.PF);
contour(map.xx,map.bb,map.NsI0./setup.Imax,'-','EdgeColor',colors(6,:),'LineWidth',1,'ShowText','on','DisplayName','$N_s\cdot I_{max}$','LevelList',setup.Ns,'LabelFormat','Ns=%d');
contour(map.xx,map.bb,setup.Fmax./map.F0_Ns,'-','EdgeColor',colors(7,:),'LineWidth',1,'ShowText','on','DisplayName','$\lambda_{max}/N_s$','LevelList',setup.Ns,'LabelFormat','Ns=%d');

T = interp2(map.xx,map.bb,map.T,xx,bb);
feasible(T<setup.T)=NaN;
PF = interp2(map.xx,map.bb,map.PF,xx,bb);
feasible(PF<setup.PF)=NaN;

surf(xx,bb,feasible,'EdgeColor','none','FaceColor','g','FaceAlpha',0.5,'DisplayName','Feasible Designs');
plot(map.xRaw,map.bRaw,'o','Color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'DisplayName','FEAfix');

hfig(2) = figure();
figSetting(14,14);
hax(2) = axes(...
    'OuterPosition',[0 0 1 1],...
    'XLim',[min(xx(:)) max(xx(:))],...
    'YLim',[min(bb(:)) max(bb(:))]);
xlabel('$x$')
ylabel('$b$')
title(['Feasibility area @ $I_{max}=' num2str(setup.Imax) '$ Apk / $V_{DC}=' int2str(setup.Vdc) '$ V'])
hleg(2) = legend('show','Location','southwest');
set(hfig(2),'UserData',map);


contour(map.xx,map.bb,map.T,'-','EdgeColor','k','LineWidth',1,'ShowText','on','DisplayName','$T$ [Nm]');
contour(map.xx,map.bb,map.PF,'--','EdgeColor','k','LineWidth',1,'ShowText','on','DisplayName','$cos \varphi$');
contour(map.xx,map.bb,map.T,'-','EdgeColor','k','LineWidth',2,'ShowText','on','DisplayName',['$T=' num2str(setup.T) '$ Nm'],'LevelList',setup.T);
contour(map.xx,map.bb,map.PF,'--','EdgeColor','k','LineWidth',2,'ShowText','on','DisplayName',['$cos \varphi=' num2str(setup.PF) '$'],'LevelList',setup.PF);
contour(map.xx,map.bb,map.NsI0./setup.Imax,':','EdgeColor','k','LineWidth',1,'ShowText','on','DisplayName','$N_s\cdot I_{max}$','LevelList',setup.Ns,'LabelFormat','Ns=%d');
contour(map.xx,map.bb,setup.Fmax./map.F0_Ns,'-.','EdgeColor','k','LineWidth',1,'ShowText','on','DisplayName','$\lambda_{max}/N_s$','LevelList',setup.Ns,'LabelFormat','Ns=%d');

T = interp2(map.xx,map.bb,map.T,xx,bb);
feasible(T<setup.T)=NaN;
PF = interp2(map.xx,map.bb,map.PF,xx,bb);
feasible(PF<setup.PF)=NaN;

surf(xx,bb,feasible,'EdgeColor','none','FaceColor','k','FaceAlpha',0.5,'DisplayName','$(T,cos\varphi)$ feasible designs');
plot(map.xRaw,map.bRaw,'o','Color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'DisplayName','FEAfix');

for ii=1:length(setup.Ns)
    feasible = zeros(size(xx));
    tmp = interp2(map.xx,map.bb,map.NsI0,xx,bb);
    feasible(tmp>=setup.Ns(ii)*setup.Imax) = NaN;
    tmp = interp2(map.xx,map.bb,map.F0_Ns,xx,bb);
    feasible(tmp>=setup.Fmax/setup.Ns(ii)) = NaN;
    if sum(~isnan(feasible(:)))==0
        plotName = ['$N_s=' int2str(setup.Ns(ii)) '$ - unfeasible'];
    else
        plotName = ['$N_s=' int2str(setup.Ns(ii)) '$ - feasible'];
    end
    surf(xx,bb,feasible,'EdgeColor','none','FaceColor',colors(ii,:),'FaceAlpha',0.5,'DisplayName',plotName)
end








