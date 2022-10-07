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

function MMM_plot_skinEffect(motorModel)

skinEffect = motorModel.acLossFactor;

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'Skin Effect Model\'];

f = skinEffect.f;
k = skinEffect.k;
% p = skinEffect.p;

% fPlot = linspace(0,1.2*max(f),1001);
% kPlot = polyval(p,fPlot);

hfig(1) = figure();
figSetting()
hax(1) = axes('OuterPosition',[0 0 1 1],...
    'XLim',[0 max(f(:))],...
    'YLim',[1 max(k(:))]);
xlabel('$f$ [$Hz$]')
ylabel('$k_{AC} = \frac{R_{AC}}{R_{DC}}$')
set(hfig,'FileName',[pathname resFolder 'skinEffectModel.fig'])

plot(f',k','-o','DisplayName','LUT model')
% plot(fPlot,kPlot,'-r','DisplayName','Fit Model')

% hleg = legend('show','Location','northwest');

if isfield(skinEffect,'T')
    T = skinEffect.T;
    tempVect = unique(T);
    hfig(2) = figure();
    figSetting()
    hax(2) = axes('OuterPosition',[0 0 1 1],...
        'XLim',[0 max(f(:))],...
        'YLim',[min(T(:)),max(T(:))],...
        'ZLim',[1 max(k(:))]);
    view(hax(2),3)
    xlabel('$f$ [Hz]')
    ylabel('$\Theta_{Cu}$ [$^\circ$C]')
    zlabel('$k_{AC} = \frac{R_{AC}}{R_{DC}}$')
    set(hfig(2),'FileName',[pathname resFolder 'skinEffectModel2D.fig'])
    
    surf(hax(2),f,T,k,'FaceColor','interp','EdgeColor','k');
    plot3(f(:),T(:),k(:),'bo');

    cla(hax(1));
%     hax(1) = axes('OuterPosition',[0 0 1 1],...
%         'XLim',[0 max(f(:))],...
%         'YLim',[1 max(k(:))]);
    for ii=1:length(tempVect)
        plot(hax(1),f(ii,:),k(ii,:),'-o','DisplayName',['$\Theta_{Cu}=' int2str(tempVect(ii)) '^\circ$C'])
    end
    legend(hax(1),'show','Location','northwest');
    
    hfig(3) = figure();
    figSetting();
    set(hfig(3),'Filename',[pathname resFolder 'phaseACresistance.fig'])
    hax(3) = axes('OuterPosition',[0 0 1 1],...
        'XLim',[0 max(f(:))]);
    xlabel('$f$ [Hz]')
    ylabel('$R_{AC}$ [$\Omega$]')
    Rs0   = motorModel.data.Rs;
    temp0 = motorModel.data.tempCu;
    l     = motorModel.data.l;
    lend  = motorModel.data.lend;
    Rs = calcRsTempFreq(Rs0,temp0,l,lend,skinEffect,'LUT',T,f);
    for ii=1:length(tempVect)
        plot(hax(3),f(ii,:),Rs(ii,:),'-o','DisplayName',['$\Theta_{Cu}=' int2str(tempVect(ii)) '^\circ$C'])
    end
    legend(hax(3),'show','Location','northeastoutside');


end

%% Save figures
answer = 'No';
answer = questdlg('Save figures?','Save','Yes','No',answer);
if strcmp(answer,'Yes')
    if ~exist([pathname resFolder],'dir')
        mkdir([pathname resFolder]);
    end
    
    for ii=1:length(hfig)
        savePrintFigure(hfig(ii));
    end
end

