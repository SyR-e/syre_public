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

function MMM_plot_VSgamma(motorModel)

flag4Q = 1;
if flag4Q
    warning('Debug mode active...4Q plots')
end

fdfq = motorModel.FluxMap_dq;
% mirror the flux maps, if needed
if strcmp(motorModel.data.axisType,'SR')
    if min(fdfq.Id,[],'all')>=0
        fdfq.Id = [-fliplr(fdfq.Id(:,2:end)) fdfq.Id];
        fdfq.Iq = [fdfq.Iq(:,2:end) fdfq.Iq];
        fdfq.Fd = [-fliplr(fdfq.Fd(:,2:end)) fdfq.Fd];
        fdfq.Fq = [fliplr(fdfq.Fq(:,2:end)) fdfq.Fq];
        fdfq.T  = [-fliplr(fdfq.T(:,2:end)) fdfq.T];
        fdfq.dTpp = [fliplr(fdfq.dTpp(:,2:end)) fdfq.dTpp];
    end
    if strcmp(motorModel.data.motorType,'SR')
        if min(fdfq.Iq,[],'all')>=0
            fdfq.Id = [fdfq.Id(2:end,:);fdfq.Id];
            fdfq.Iq = [-flipud(fdfq.Iq(2:end,:));fdfq.Iq];
            fdfq.Fd = [flipud(fdfq.Fd(2:end,:));fdfq.Fd];
            fdfq.Fq = [-flipud(fdfq.Fq(2:end,:));fdfq.Fq];
            fdfq.T  = [-flipud(fdfq.T(2:end,:));fdfq.T];
            fdfq.dTpp  = [flipud(fdfq.dTpp(2:end,:));fdfq.dTpp];
        end
    end
else
    if min(fdfq.Iq,[],'all')>=0
        fdfq.Id = [fdfq.Id(2:end,:);fdfq.Id];
        fdfq.Iq = [-flipud(fdfq.Iq(2:end,:));fdfq.Iq];
        fdfq.Fd = [flipud(fdfq.Fd(2:end,:));fdfq.Fd];
        fdfq.Fq = [-flipud(fdfq.Fq(2:end,:));fdfq.Fq];
        fdfq.T  = [-flipud(fdfq.T(2:end,:));fdfq.T];
        fdfq.dTpp  = [flipud(fdfq.dTpp(2:end,:));fdfq.dTpp];
    end
end

% load data
Id   = fdfq.Id;
Iq   = fdfq.Iq;
Fd   = fdfq.Fd;
Fq   = fdfq.Fq;
T    = fdfq.T;
dTpp = fdfq.dTpp;

% nCurr = motorModel.data.nCurr;
% i0    = motorModel.data.i0;
% Imax  = motorModel.data.Imax;
axisType = motorModel.data.axisType;

Ivect = motorModel.WaveformSetup.CurrAmpl;

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'plotVSgamma - ' int2str(motorModel.data.tempPM) 'deg\'];

IPF = sin(atan2(Iq,Id)-atan2(Fq,Fd));

% IdMax = max(abs(Id),[],'all');
% IqMax = max(abs(Iq),[],'all');

% Imax = min([Imax,IdMax,IqMax]);

if flag4Q
    numPoints = 500;
    gammaVect = linspace(0,360,numPoints);
else
    numPoints = 46;
    
    if strcmp(axisType,'SR')
        gammaVect = linspace(0,90,numPoints);
    else
        gammaVect = linspace(90,180,numPoints);
    end
end


% create figures
figNames{1} = 'torqueVSgamma';
figNames{2} = 'pwrfctVSgamma';
figNames{3} = 'fluxVSgamma';
figNames{4} = 'deltaVSgamma';
figNames{5} = 'torripVSgamma';

for ii=1:length(figNames)
    hfig(ii) = figure();
    figSetting();
    hax(ii) = axes('OuterPosition',[0 0 1 1],...
        'XLim',[min(gammaVect) max(gammaVect)],...
        'XTick',min(gammaVect):15:max(gammaVect));
    xlabel('$\gamma$ [$^\circ$]')
    set(hfig(ii),'FileName',[pathname resFolder figNames{ii} '.fig']);
    set(hfig(ii),'Name',figNames{ii});
    switch ii
        case 1
            ylabel('$T$ [Nm]')
            legend(hax(ii),'show','Location','northwest')
        case 2
            ylabel('$cos \varphi$')
            legend(hax(ii),'show','Location','northwest')
        case 3
            ylabel('$|\lambda|$ [Vs]')
            legend(hax(ii),'show','Location','southwest')
        case 4
            ylabel('$\delta$ [$^\circ$]')
            legend(hax(ii),'show','Location','northwest')
        case 5
            ylabel('$\Delta T_{pp}$ [Nm]')
            legend(hax(ii),'show','Location','northwest')
    end
end

% Ivect = nCurr*i0;

% Ivect = Ivect(Ivect<=Imax);

for ii=1:length(Ivect)
    TVect     = interp2(Id,Iq,T,Ivect(ii)*cosd(gammaVect),Ivect(ii)*sind(gammaVect));
    IPFVect   = interp2(Id,Iq,IPF,Ivect(ii)*cosd(gammaVect),Ivect(ii)*sind(gammaVect));
    FVect     = interp2(Id,Iq,abs(Fd+j*Fq),Ivect(ii)*cosd(gammaVect),Ivect(ii)*sind(gammaVect));
    deltaVect = interp2(Id,Iq,atan2(Fq,Fd)*180/pi,Ivect(ii)*cosd(gammaVect),Ivect(ii)*sind(gammaVect));
    dTppVect  = interp2(Id,Iq,dTpp,Ivect(ii)*cosd(gammaVect),Ivect(ii)*sind(gammaVect));
    
    plot(hax(1),gammaVect,TVect,'-x','DisplayName',['$I=' num2str(round(Ivect(ii),2)) '$ A'])
    plot(hax(2),gammaVect,IPFVect,'-x','DisplayName',['$I=' num2str(round(Ivect(ii),2)) '$ A'])
    plot(hax(3),gammaVect,FVect,'-x','DisplayName',['$I=' num2str(round(Ivect(ii),2)) '$ A'])
    plot(hax(4),gammaVect,deltaVect,'-x','DisplayName',['$I=' num2str(round(Ivect(ii),2)) '$ A'])
    plot(hax(5),gammaVect,dTppVect,'-x','DisplayName',['$I=' num2str(round(Ivect(ii),2)) '$ A'])
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






