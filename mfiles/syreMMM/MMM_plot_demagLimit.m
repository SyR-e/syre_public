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

function MMM_plot_demagLimit(motorModel)

% load data
demagLimit = motorModel.DemagnetizationLimit;
i0         = motorModel.data.i0;
Imax       = motorModel.data.Imax;
test       = demagLimit.test;
nIter      = 1:1:numel(test.Idemag);

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'Demagnetization Limit\'];

if isfield(motorModel,'mat')
    mat = motorModel.mat.LayerMag;
    if isfield(mat,'kgm3')
        demagLimit.Br      = interp1(mat.temp.temp,mat.temp.Br,demagLimit.tempPM);
        demagLimit.Bd      = interp1(mat.temp.temp,mat.temp.Bd,demagLimit.tempPM);
        demagLimit.test.Br = interp1(mat.temp.temp,mat.temp.Br,demagLimit.test.tempPM);
        demagLimit.test.Bd = interp1(mat.temp.temp,mat.temp.Bd,demagLimit.test.tempPM);
    else
        demagLimit.Br      = nan(size(demagLimit.tempPM));
        demagLimit.Bd      = nan(size(demagLimit.tempPM));
        demagLimit.test.Br = nan(size(demagLimit.test.tempPM));
        demagLimit.test.Bd = nan(size(demagLimit.test.tempPM));
    end
else
    demagLimit.Br      = nan(size(demagLimit.tempPM));
    demagLimit.Bd      = nan(size(demagLimit.tempPM));
    demagLimit.test.Br = nan(size(demagLimit.test.tempPM));
    demagLimit.test.Bd = nan(size(demagLimit.test.tempPM));
end

test = demagLimit.test;

% Prepare figures
for ii=1:3
    hfig(ii) = figure();
    figSetting()
    hax(ii,1) = axes('OuterPosition',[0 0 1 1]);
    xlabel('$\Theta_{PM}$ ($^\circ$C)')
    set(hax(ii,1),'XLim',[min(demagLimit.tempPM) max(demagLimit.tempPM)]);
    switch ii
        case 1
            ylabel('$I_{demag}$ (Apk)')
            set(hax(ii,1),'YLim',[0 max(demagLimit.Idemag)]);
            set(hfig(ii),'FileName',[pathname resFolder 'currentVStemperature.fig']);
            set(hfig(ii),'Name','currentVStemperature');
            legend('show','Location','best')
        case 2
            ylabel('$B_{min}$ (T)')
            set(hax(ii,1),'YLim',[min(demagLimit.Bmin) max(demagLimit.Br)]);
            set(hfig(ii),'FileName',[pathname resFolder 'fluxDensityVStemperature.fig']);
            set(hfig(ii),'Name','fluxDensityVStemperature');
            legend('show','Location','best')
        case 3
            ylabel('Volume of PM demagnetized (\%)')
            set(hax(ii,1),'YLim',[0 5]);
            set(hfig(ii),'FileName',[pathname resFolder 'volumeVStemperature.fig']);
            set(hfig(ii),'Name','volumeVStemperature');
            legend('show','Location','best')
    end
end

hfig(4) = figure();
figSetting(20,20)
hax(4,1) = axes('OuterPosition',[0 3/4 1 1/4]);
hax(4,2) = axes('OuterPosition',[0 2/4 1 1/4]);
hax(4,3) = axes('OuterPosition',[0 1/4 1 1/4]);
hax(4,4) = axes('OuterPosition',[0 0/4 1 1/4]);
set(hfig(4),'FileName',[pathname resFolder 'iterations.fig'])
set(hfig(4),'Name','iterations')

yMax = max([test.Idemag,100,test.tempPM]);


for ii=1:4
    set(hax(4,ii),'XLim',[0 length(test.Idemag)+1]);
    xlabel('iteration')
    switch ii
        case 1
            ylabel(hax(4,ii),'$i_{test}$ (A)')
            set(hax(4,ii),'YLim',[0 max(test.Idemag)])
        case 2
            ylabel(hax(4,ii),'PM demag (\%)')
            set(hax(4,ii),'YLim',[0 100])
        case 3
            ylabel(hax(4,ii),'$\Theta_{PM}$ ($^\circ$C)')
            tmpLim = [min(test.tempPM) max(test.tempPM)];
            if tmpLim(1)>0
                tmpLim(1) = 0;
            end
            set(hax(4,ii),'YLim',tmpLim);
            plot(hax(4,ii),[0 length(test.Idemag)+1],[0 0],'-k','LineWidth',1,'HandleVisibility','off')
        case 4
            ylabel(hax(4,ii),'$B_{min}$ (T)')
            tmpLim = [min([test.Bmin test.Br test.Bd]) max([test.Bmin test.Br test.Bd])];
            if tmpLim(1)>0
                tmpLim(1) = 0;
            end
            set(hax(4,ii),'YLim',tmpLim);
            plot(hax(4,ii),[0 length(test.Idemag)+1],[0 0],'-k','LineWidth',1,'HandleVisibility','off')
    end
%     for jj=1:length(xIndex)
%         plot(hax(ii),xIndex(jj)*[1 1],yMax*[-1 1],'--k','LineWidth',0.5,'HandleVisibility','off')
%     end
end

hfig(5) = figure();
figSetting()
hax(5,1) = axes('OuterPosition',[0 0 1 1]);
xlabel('$\Theta_{PM}$ ($^\circ$C)')
ylabel('$I_{demag}$ (Apk)')
zlabel('PM demag (\%)')
% clim([0 100])
zlim([0 100])
colormap turbo
colorbar
view(3)
set(hfig(5),'FileName',[pathname resFolder 'demagMap.fig']);
set(hfig(5),'Name','demagMap');

% Plot figures

plot(hax(1,1),demagLimit.tempPM,demagLimit.Idemag,'-o','DisplayName','Demagnetization Limit');
plot(hax(1,1),demagLimit.tempPM,Imax*ones(size(demagLimit.tempPM)),'-','DisplayName','Maximum current')
plot(hax(1,1),demagLimit.tempPM,i0*ones(size(demagLimit.tempPM)),'-','DisplayName','Rated current')

plot(hax(2,1),demagLimit.tempPM,demagLimit.Bmin,'-o','DisplayName','$B_{min}$')
plot(hax(2,1),demagLimit.tempPM,demagLimit.Bd,'-x','DisplayName','$B_{d}$')
plot(hax(2,1),demagLimit.tempPM,demagLimit.Br,'-d','DisplayName','$B_{r}$')

plot(hax(3,1),demagLimit.tempPM,demagLimit.dPM*100,'-o','DisplayName','dPM')

stem(hax(4,1),nIter,test.Idemag,'bo','LineWidth',1.5,'DisplayName','$I_{test}$');
stem(hax(4,2),nIter,test.dPM*100,'bo','LineWidth',1.5,'DisplayName','$dPM$');
plot(hax(4,3),nIter,test.tempPM,'bo','DisplayName','$\Theta_{PM}$ [$^\circ$C]');
X = [1 nIter nIter(end) 1];
Y = [min([test.Bmin test.Br test.Bd]) test.Bd min([test.Bmin test.Br test.Bd]) min([test.Bmin test.Br test.Bd])];
fill(hax(4,4),X,Y,'r','FaceAlpha',0.25,'EdgeColor','none');
X = [1 nIter fliplr(nIter) 1];
Y = [min([test.Bmin test.Br test.Bd]) test.Br fliplr(test.Bd) min([test.Bmin test.Br test.Bd])];
fill(hax(4,4),X,Y,'g','FaceAlpha',0.25,'EdgeColor','none');
plot(hax(4,4),nIter,test.Bmin,'-bo','DisplayName','$B_{min}$')
plot(hax(4,4),nIter,test.Bd,'-r','DisplayName','$B_d$')
plot(hax(4,4),nIter,test.Br,'-g','DisplayName','$B_r$')

scatter3(hax(5,1),test.tempPM,test.Idemag,test.dPM*100,[],test.dPM*100,'filled','LineWidth',0.5,'markerEdgeColor','k','DisplayName','Test points')
stem3(hax(5,1),test.tempPM,test.Idemag,test.dPM*100,':k','HandleVisibility','off')
plot3(hax(5,1),demagLimit.tempPM,demagLimit.Idemag,demagLimit.dPM,'-k','DisplayName','demagnetization limit @ 1\%')


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


