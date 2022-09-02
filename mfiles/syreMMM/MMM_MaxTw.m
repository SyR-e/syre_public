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

function MMM_MaxTw(motorModel,hax)

% load data
Id  = motorModel.FluxMap_dq.Id;
Iq  = motorModel.FluxMap_dq.Iq;
Tem = motorModel.FluxMap_dq.T;

TwData = motorModel.TnSetup;

motorType = motorModel.data.motorType;

nVect = linspace(TwData.nmin,TwData.nmax,TwData.nstep);
TVect = linspace(TwData.Tmin,TwData.Tmax,TwData.Tstep);

% Initialize matrix results
[nmap,Tmap] = meshgrid(nVect,TVect);

TwMap.n     = nmap;             % speed reference [rpm]
TwMap.T     = Tmap;             % torque (mechanical) reference [Nm]
TwMap.Tem   = nan(size(nmap));  % electro-magnetic torque [Nm]
TwMap.Tout  = nan(size(nmap));  % torque produced (mechanical, for operating limits) [Nm]
TwMap.Id    = nan(size(nmap));  % d-axis magnetizing current [A]
TwMap.Iq    = nan(size(nmap));  % q-axis magnetizing current [A]
TwMap.Im    = nan(size(nmap));  % dq magnetizing current [A]
TwMap.Iph   = nan(size(nmap));  % dq phase current [A]
TwMap.Vph   = nan(size(nmap));  % dq phase voltage [A]
TwMap.Fd    = nan(size(nmap));  % d-axis magnetizing flux linkage [Vs]
TwMap.Fq    = nan(size(nmap));  % q-axis magnetizing flux linkage [Vs]
TwMap.Vo    = nan(size(nmap));  % peak line voltage [V]
TwMap.Io    = nan(size(nmap));  % peak phase current [A]
TwMap.PF    = nan(size(nmap));  % power factor
TwMap.P     = nan(size(nmap));  % Output power [W]
TwMap.Ploss = nan(size(nmap));  % Total loss [W]
TwMap.Pjs   = nan(size(nmap));  % Stator Joule loss (total) [W]
TwMap.PjDC  = nan(size(nmap));  % Stator Joule loss (only DC) [W]
TwMap.PjAC  = nan(size(nmap));  % Stator Joule loss (only AC) [W]
TwMap.Pfe   = nan(size(nmap));  % Total iron loss [W]
TwMap.Pfes  = nan(size(nmap));  % Stator iron loss [W]
TwMap.Pfer  = nan(size(nmap));  % Rotor iron loss [W]
TwMap.Ppm   = nan(size(nmap));  % Permanent magnet loss [W]
TwMap.Pjr   = nan(size(nmap));  % Rotor Joule loss [W]
TwMap.Pmech = nan(size(nmap));  % Mechanical loss [W]
TwMap.Eo    = nan(size(nmap));  % Back-emf [V]
TwMap.Ife   = nan(size(nmap));  % Current for iron loss [A]
TwMap.slip  = nan(size(nmap));  % Rotor slip
TwMap.Ir    = nan(size(nmap));  % Rotor current [A]
TwMap.Rs    = nan(size(nmap));  % Stator resistance [Ohm]
TwMap.eff   = nan(size(nmap));  % Efficiency [pu]


if nargin==1
    % no link to an existing axis...create a new axis
    figure()
    figSetting();
    hax = axes('OuterPosition',[0 0 1 1]);
end

cla(hax)

if min(nVect)~=max(nVect)
    set(hax,'XLim',[min(nVect) max(nVect)]);
end

if min(TVect)~=max(TVect)
    set(hax,'YLim',[min(TVect) max(TVect)]);
end

hr=plot(hax,0,0,'r.','LineWidth',1.5);
hg=plot(hax,0,0,'g.','LineWidth',1.5);



disp('Map evaluation in progress...')
fprintf(' %06.2f%%',0)

for ii=1:numel(TwMap.n)
    
    [out] = calcTnPoint(motorModel,Tmap(ii),nmap(ii));
    
    % update figure and matrices
    if ~isnan(out.T)
        xdata = [get(hg,'XData') nmap(ii)];
        ydata = [get(hg,'YData') Tmap(ii)];
        set(hg,'XData',xdata,'YData',ydata);
        drawnow();
        
        TwMap.Tout(ii)  = out.T;
        TwMap.Id(ii)    = out.Id;
        TwMap.Iq(ii)    = out.Iq;
        TwMap.Fd(ii)    = out.Fd;
        TwMap.Fq(ii)    = out.Fq;
        TwMap.Tem(ii)   = out.Tem;
        TwMap.Vo(ii)    = out.Vo;
        TwMap.Io(ii)    = out.Io;
        TwMap.Im(ii)    = out.Im;
        TwMap.Iph(ii)   = out.Iph;
        TwMap.Vph(ii)   = out.Vph;
        TwMap.PF(ii)    = out.PF;
        TwMap.P(ii)     = out.P;
        TwMap.Ploss(ii) = out.Ploss;
        TwMap.Pjs(ii)   = out.Pjs;
        TwMap.PjDC(ii)  = out.PjDC;
        TwMap.PjAC(ii)  = out.PjAC;
        TwMap.Pfe(ii)   = out.Pfe;
        TwMap.Pfes(ii)  = out.Pfes;
        TwMap.Pfer(ii)  = out.Pfer;
        TwMap.Ppm(ii)   = out.Ppm;
        TwMap.Pjr(ii)   = out.Pjr;
        TwMap.Pmech(ii) = out.Pmech;
        TwMap.Eo(ii)    = out.Eo;
        TwMap.Ife(ii)   = out.Ife;
        TwMap.slip(ii)  = out.slip;
        TwMap.Ir(ii)    = out.Ir;
        TwMap.Rs(ii)    = out.Rs;
        TwMap.eff(ii)   = out.eff;
    else
        xdata = [get(hr,'XData') nmap(ii)];
        ydata = [get(hr,'YData') Tmap(ii)];
        set(hr,'XData',xdata,'YData',ydata);
        drawnow();
%         pause(0.01);
    end
    
    fprintf('\b\b\b\b\b\b\b')
    fprintf('%06.2f%%',ii/numel(TwMap.n)*100)

end

disp(' ')
disp('Maps Evaluated');

TwMap.dTpp = interp2(motorModel.FluxMap_dq.Id,motorModel.FluxMap_dq.Iq,motorModel.FluxMap_dq.dTpp,TwMap.Id,TwMap.Iq);

TwMap.T_top_W = max(TwMap.Tout);
TwMap.T_top_W(TwMap.T_top_W<0) = 0;
TwMap.T_bot_W = min(TwMap.Tout);
TwMap.T_bot_W(TwMap.T_bot_W>0) = 0;

%% plot results

pathname = motorModel.data.pathname;
motName = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'TwMap_' datestr(now,30) '\'];

indexFig = 1;

figNames{1}  = 'Torque limits';
figNames{2}  = 'Efficiency map';
figNames{3}  = 'Output power map';
figNames{4}  = 'Total loss map';
figNames{5}  = 'Stator Joule loss map';
figNames{6}  = 'Stator DC Joule loss map';
figNames{7}  = 'Stator AC Joule loss map';
figNames{8}  = 'Iron loss map';
figNames{9}  = 'Stator iron loss map';
figNames{10} = 'Rotor iron loss map';
figNames{11} = 'PM loss map';
figNames{12} = 'Phase current map';
figNames{13} = 'Line voltage map';
figNames{14} = 'Power factor map';
figNames{15} = 'Mechanical loss map';
figNames{16} = 'Iron and mechanical loss map';
figNames{17} = 'Rotor Joule loss map';
figNames{18} = 'Rotor slip map';
figNames{19} = 'Rotor current map';
figNames{20} = 'Torque ripple';

flagPlot = 1:1:numel(figNames);
if strcmp(TwData.IronLossFlag,'No')
    flagPlot(8)  = 0;
    flagPlot(9)  = 0;
    flagPlot(10) = 0;
    flagPlot(11) = 0;
end
if (TwData.MechLoss==0)
    flagPlot(15) = 0;
end
if (strcmp(TwData.IronLossFlag,'No')&&(TwData.MechLoss==0))
    flagPlot(16) = 0;
end
if strcmp(TwData.SkinEffectFlag,'No')
    flagPlot(6)  = 0;
    flagPlot(7)  = 0;
end
if ~strcmp(motorType,'IM')
    flagPlot(17) = 0;
    flagPlot(18) = 0;
    flagPlot(19) = 0;
end

flagPlot = flagPlot(flagPlot~=0);

for ii=1:length(flagPlot)
    hfig(ii) = figure();
    figSetting();
    set(hfig(ii),'FileName',[pathname resFolder figNames{flagPlot(ii)} '.fig']);
    hax(ii) = axes(...
        'XLim',[TwData.nmin TwData.nmax],...
        'YLim',[TwData.Tmin TwData.Tmax]);
    xlabel('$n$ [rpm]')
    ylabel('$T$ [Nm]')
    switch flagPlot(ii)
        case 1
            title('Torque limits')
            plot(unique(TwMap.n),TwMap.T_top_W,'-bx')
            plot(unique(TwMap.n),TwMap.T_bot_W,'-bx')
        case 2
            title('Efficiency map [p.u.]')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.eff,[64:2:86 87:1:100]/100);
            clabel(c,h)
            colorbar
        case 3
            title('Output power [W]')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.P);
            clabel(c,h)
            colorbar
        case 4
            title('Total loss [W]')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.Ploss);
            clabel(c,h)
            colorbar
        case 5
            title('Stator Joule loss [W]')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.Pjs);
            clabel(c,h)
            colorbar
        case 6
            title('Stator Joule DC loss [W]')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.PjDC);
            clabel(c,h)
            colorbar
        case 7
            title('Stator Joule AC loss [W]')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.PjAC);
            clabel(c,h)
            colorbar
        case 8
            title('Iron loss [W]')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.Pfe);
            clabel(c,h)
            colorbar
        case 9
            title('Stator iron loss [W]')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.Pfes);
            clabel(c,h)
            colorbar
        case 10
            title('Rotor iron loss [W]')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.Pfer);
            clabel(c,h)
            colorbar
        case 11
            title('Permanent magnet loss [W]')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.Ppm);
            clabel(c,h)
            colorbar
        case 12
            title('Phase current map [Apk]')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.Io);
            clabel(c,h)
            colorbar
        case 13
            title('Line voltage map [Vpk]')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.Vo);
            clabel(c,h)
            colorbar
        case 14
            title('Power factor map')
            [c,h] = contourf(TwMap.n,TwMap.T,abs(TwMap.PF),[0.4:0.05:1]);
            clabel(c,h)
            colorbar
        case 15
            title('Mechanical loss map')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.Pmech);
            clabel(c,h)
            colorbar
        case 16
            title('Iron + PM + mechanical loss map')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.Pmech+TwMap.Pfes+TwMap.Pfer+TwMap.Ppm);
            clabel(c,h)
            colorbar
        case 17
            title('Rotor Joule loss map')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.Pjr);
            clabel(c,h)
            colorbar
        case 18
            title('Rotor slip map')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.slip,[-1:0.1:1]);
            clabel(c,h)
            colorbar
        case 19
            title('Rotor current map [A]')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.Ir);
            clabel(c,h)
            colorbar
        case 20
            title('Peak-to-peak torque ripple [Nm]')
            [c,h] = contourf(TwMap.n,TwMap.T,TwMap.dTpp);
            clabel(c,h)
            colorbar
    end
    if flagPlot(ii)>1
        plot(unique(TwMap.n),TwMap.T_top_W,'-k')
        plot(unique(TwMap.n),TwMap.T_bot_W,'-k')
    end
end

ii=ii+1;

hfig(ii) = figure();
figSetting();
set(hfig(ii),'FileName',[pathname resFolder 'Control locus.fig']);
hax(ii) = axes(...
    'XLim',[min(Id,[],'all') max(Id,[],'all')],...
    'YLim',[min(Id,[],'all') max(Iq,[],'all')],...
    'PlotBoxAspectRatio',[1 1 1]);
xlabel('$i_d$ [A]')
ylabel('$i_q$ [A]')
if strcmp(TwData.Control,'Maximum efficiency')
    title('Maximum efficiency locus')
elseif strcmp(TwData.Control,'MTPA')
    title('Maximum torque per minimum copper loss locus')
end
[c,h] = contour(Id,Iq,Tem,'k','DisplayName','$T$ [Nm]');
clabel(c,h)
plot(TwMap.Id,TwMap.Iq)


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
    save([pathname resFolder 'TwMap.mat'],'motorModel','TwMap','TwData');
    
end
