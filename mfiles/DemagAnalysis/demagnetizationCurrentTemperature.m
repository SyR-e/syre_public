% Copyright 2018
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
%    See the License for the spcific language governing permissions and
%    limitations under the License.

function [DemagOut]=demagnetizationCurrentTemperature(dataSet)
% 
% [DemagOut]=demagnetizationCurrentTemperature(setup)
% 
% Compute the demagnetization current function of the temperature.

clc
if nargin<1
    
    [filename,pathname,~]=uigetfile([cd '\.mat'],'Select a PM machine');
    load([pathname filename]);
    
    prompt={'Maximum current [A]',...
            'Number of current levels',...
            };
    name='Demagnetization input';
    [i0,~]=calc_io(geo,per);
    if ~isfield(mat.LayerMag,'temp')
        error('No thermal properties for the selected PMs. Please select another PMs or update the material properties')
    end
    defaultanswer={ num2str(i0),...
                    num2str(dataSet.NumGrid),...
                    };
    answer=inputdlg(prompt,name,1,defaultanswer);
    
    dataSet.currentfilename = filename;
    dataSet.currentpathname = pathname;
    dataSet.CurrLoPP        = eval(answer{1})/i0;
    dataSet.NumGrid         = eval(answer{2});
    dataSet.tempPP          = mat.LayerMag.temp.temp;
    
end

load([dataSet.currentpathname dataSet.currentfilename],'geo','per','mat');

%[i0,~]=calc_io(geo,per);

i0 = dataSet.SimulatedCurrent(end)/dataSet.CurrLoPP(end);

filename = dataSet.currentfilename;
pathname = dataSet.currentpathname;
Imax     = dataSet.CurrLoPP*i0;
TempVect = dataSet.tempPP;
ntemp    = length(TempVect);
ncur     = dataSet.NumGrid;

resFolder=[filename(1:end-4) '_Demag_' int2str(Imax) 'A\'];
if ~isfolder([pathname resFolder])
    mkdir(pathname,resFolder);
end

source=[pathname filename];
destin=[pathname resFolder filename];
copyfile(source,destin);

source=[pathname filename(1:end-4) '.fem'];
destin=[pathname resFolder filename(1:end-4) '.fem'];
copyfile(source,destin);

input.filename     = dataSet.currentfilename;
input.pathname     = [dataSet.currentpathname resFolder];
% input.Imax         = Imax;
% input.ncur         = ncur;
input.geo          = geo;
input.mat          = mat;
input.per          = per;
%input.per.overload = linspace(0,Imax,ncur)/i0;
input.dataSet      = dataSet;

Ivect=linspace(0,Imax,ncur+1)/i0;
Ivect=Ivect(2:end); % avoid zero current
input.per.overload=Ivect;

for tt=1:ntemp
    disp(['Evaluation temperature ' int2str(tt) ' of ' int2str(ntemp)]);
    inputTemp            = input;
    inputTemp.per.tempPP = TempVect(tt);
    inputTemp.figFlag    = 1;
    OUT(tt)=demagnetizationDetectionPoint(inputTemp);
end
disp('FEMM evaluations done!')
%close all

save([pathname resFolder 'DemagCurrTemp_RawData.mat'],'dataSet','geo','per','mat','OUT')

Iout      = zeros(ncur,ntemp,2); % dimensions=[current,temperature,position]
dPMout    = zeros(ncur,ntemp,2);
TempOut   = zeros(ncur,ntemp,2);
xyDemag.C = cell(ncur,ntemp,2);
xyDemag.V = cell(ncur,ntemp,2);
xyDemag.B = cell(ncur,ntemp,2);

for tt=1:ntemp
    for ii=1:ncur
        for pp=1:2
            Iout(ii,tt,pp)      = OUT(tt).Ivect(ii);
            VolDem              = OUT(tt).VolDem(ii,pp);
            VolTot              = OUT(tt).VolTot(ii,pp);
            dPMout(ii,tt,pp)    = VolDem./VolTot;
            TempOut(ii,tt,pp)   = TempVect(tt);
            xyDemag.C{ii,tt,pp} = OUT(tt).xyDemag.C{ii,pp};
            xyDemag.V{ii,tt,pp} = OUT(tt).xyDemag.V{ii,pp};
            xyDemag.B{ii,tt,pp} = OUT(tt).xyDemag.B{ii,pp};
        end
    end
end

%DemagOut.setup   = setup;
DemagOut.I       = Iout;
DemagOut.Temp    = TempOut;
DemagOut.dPMpos  = dPMout;
DemagOut.dPMtot  = sum(dPMout,3)/2;
DemagOut.xyDemag = xyDemag;

DemagOut.TempDemag = TempVect;
DemagOut.IDemag    = Imax*ones(size(TempVect));

dPMtot = max(dPMout,[],3);
for tt=1:length(DemagOut.TempDemag)
    Ivect = OUT(tt).Ivect;
    index=find(dPMtot(:,tt)<=0,1,'last');
    DemagOut.IDemag(tt)=Ivect(index);
end

save([pathname resFolder 'DemagCurrentTemperature.mat'],'DemagOut','dataSet','geo','per','mat','OUT')

if length(TempVect)>1
    figure()
    figSetting()
    xlabel('$I_{demag}$ [A]')
    ylabel('$\theta_{PM}$ [$^\circ$C]')
    zlabel('$PM_{demag}$ [pu]')
    view(3)
    set(gca,'XLim',[0 Imax],'YLim',[min(TempVect) max(TempVect)],'ZLim',[0 1]);
    set(gca,'PlotBoxAspectRatio',[1 1 0.8]);
    surf(Iout(:,:,1),TempOut(:,:,1),dPMout(:,:,1),'FaceColor',[0 0 1],'EdgeColor',0.5*[0 0 1],'DisplayName','Position 1');
    surf(Iout(:,:,2),TempOut(:,:,2),dPMout(:,:,2),'FaceColor',[1 0 0],'EdgeColor',0.5*[1 0 0],'DisplayName','Position 2');
    saveas(gcf,[pathname resFolder 'DemagSurface.fig']);
    
    figure()
    figSetting()
    xlabel('$\theta_{PM}$ [$^\circ$C]')
    ylabel('$I_{demag}$ [A]')
    title('Demagnetization Limit $I_{demag}=f(\theta_{PM})$')
    set(gca,'XLim',[min(TempVect) max(TempVect)],'YLim',[0 Imax]);
    plot(DemagOut.TempDemag,DemagOut.IDemag,'-ro','DisplayName','Demagnetization Limit')
    saveas(gcf,[pathname resFolder 'DemagLimit.fig']);
end

figure()
figSetting()
xlabel('$I_{demag}$ [A]')
ylabel('$PM_{demag}$ [pu]')
set(gca,'XLim',[0 Imax],'YLim',[0 1]);
for tt=1:ntemp
    plot(Iout(:,tt,1),dPMtot(:,tt),'-','DisplayName',['$\theta_{PM}=' num2str(TempVect(tt)) '^\circ C$'])
end
legend('show','Location','NorthWest');
saveas(gcf,[pathname resFolder 'DemagCurves.fig']);


if sum(sum(sum(dPMout)))>0
    mkdir([pathname resFolder 'critical machines\'])
end

for ii=1:size(dPMout,1)
    for tt=1:size(dPMout,2)
        for pp=1:size(dPMout,3)
            setup.pathname=[pathname resFolder];
            setup.filename='DemagCurrentTemperature.mat';
            setup.iInd=ii;
            setup.tInd=tt;
            setup.pInd=pp;
            if dPMout(ii,tt,pp)>0
                plotDemagResults(setup);
                saveas(gcf,[pathname resFolder 'critical machines\curr' int2str(ii) '_temp' int2str(tt) '_pos' int2str(pp) '.fig']);
            end
        end
    end
end


% figure()
% figSetting()
% xlabel('$\theta_{PM}$ [$^\circ$C]')
% ylabel('$I_{demag}$ [A]')
% title('Demagnetization Limit $I_{demag}=f(\theta_{PM})$')
% set(gca,'XLim',[min(TempVect) max(TempVect)],'YLim',[0 Imax]);
% Tmax=[min(TempOut) min(TempOut) c(1,2:end) max(TempOut) min(TempOut)];
% Imax=[min(Iout)    max(Iout)    c(2,2:end) min(Iout)    min(Iout)];
% fill(Tmax,Imax,'g')
% Tmax=[max(TempOut) c(1,2:end) max(TempOut)];
% Imax=[max(Iout)    c(2,2:end) min(Iout)];
% fill(Tmax,Imax,'r')





