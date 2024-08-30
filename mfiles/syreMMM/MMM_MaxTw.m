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

function [TwMap,resFolderOut] = MMM_MaxTw(motorModel,hax,saveFlag)

% load data
Id  = motorModel.FluxMap_dq.Id;
Iq  = motorModel.FluxMap_dq.Iq;
Tem = motorModel.FluxMap_dq.T;

TwData = motorModel.TnSetup;

motorType = motorModel.data.motorType;

nVect = linspace(TwData.nmin,TwData.nmax,TwData.nstep);
TVect = linspace(TwData.Tmin,TwData.Tmax,TwData.Tstep);

kRefineMap = 2;

% Initialize matrix results
[nmap,Tmap] = meshgrid(nVect,TVect);

TwMap.n       = nmap;             % speed reference [rpm]
TwMap.T       = Tmap;             % torque (mechanical) reference [Nm]
TwMap.Tem     = nan(size(nmap));  % electro-magnetic torque [Nm]
TwMap.Tout    = nan(size(nmap));  % torque produced (mechanical, for operating limits) [Nm]
TwMap.Id      = nan(size(nmap));  % d-axis magnetizing current [A]
TwMap.Iq      = nan(size(nmap));  % q-axis magnetizing current [A]
TwMap.Im      = nan(size(nmap));  % dq magnetizing current [A]
TwMap.Iph     = nan(size(nmap));  % dq phase current [A]
TwMap.Vph     = nan(size(nmap));  % dq phase voltage [A]
TwMap.Fd      = nan(size(nmap));  % d-axis magnetizing flux linkage [Vs]
TwMap.Fq      = nan(size(nmap));  % q-axis magnetizing flux linkage [Vs]
TwMap.Vo      = nan(size(nmap));  % peak line voltage [V]
TwMap.Io      = nan(size(nmap));  % peak phase current [A]
TwMap.PF      = nan(size(nmap));  % power factor
TwMap.P       = nan(size(nmap));  % Output power [W]
TwMap.Ploss   = nan(size(nmap));  % Total loss [W]
TwMap.Pjs     = nan(size(nmap));  % Stator Joule loss (total) [W]
TwMap.PjDC    = nan(size(nmap));  % Stator Joule loss (only DC) [W]
TwMap.PjAC    = nan(size(nmap));  % Stator Joule loss (only AC) [W]
TwMap.Pfe     = nan(size(nmap));  % Total iron loss [W]
TwMap.Pfes    = nan(size(nmap));  % Stator iron loss [W]
TwMap.Pfer    = nan(size(nmap));  % Rotor iron loss [W]
TwMap.Ppm     = nan(size(nmap));  % Permanent magnet loss [W]
TwMap.Pjr     = nan(size(nmap));  % Rotor Joule loss [W]
TwMap.Pmech   = nan(size(nmap));  % Mechanical loss [W]
TwMap.Eo      = nan(size(nmap));  % Back-emf [V]
TwMap.Ife     = nan(size(nmap));  % Current for iron loss [A]
TwMap.slip    = nan(size(nmap));  % Rotor slip
TwMap.Ir      = nan(size(nmap));  % Rotor current [A]
TwMap.Rs      = nan(size(nmap));  % Stator resistance [Ohm]
TwMap.dTpp    = nan(size(nmap));  % Peak-to-peak torque ripple [Nm]
TwMap.eff     = nan(size(nmap));  % Efficiency [pu]
TwMap.Idemag  = nan(size(nmap));  % Demagnetization limit [Apk]
TwMap.IHWC    = nan(size(nmap));  % Peak current during ASC (HWC) [Apk]
TwMap.F0      = nan(size(nmap));  % Flusso a vuoto [Vs]
TwMap.VUGO    = nan(size(nmap));  % Tensione a vuoto [Vs]
TwMap.ASCsafe = nan(size(nmap));  % Safe area for ASC
TwMap.UGOsafe = nan(size(nmap));  % Safe area for UGO


if (nargin==1)||isempty(hax)
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

% flux maps refinement
if kRefineMap>1
    nPoints = kRefineMap*size(motorModel.FluxMap_dq.Id,1);
    motorModelFine = motorModel;
    motorModelFine.FluxMap_dq = mapsReInterpolation(motorModel.FluxMap_dq,'Id','Iq',nPoints,'linear');
else
    motorModelFine = motorModel;
end

disp('Map evaluation in progress...')
fprintf(' %06.2f%%',0)

for ii=1:numel(TwMap.n)
    
    [out] = calcTnPoint(motorModelFine,Tmap(ii),nmap(ii));
    
    % update figure and matrices
    if ~isnan(out.T)
        xdata = [get(hg,'XData') nmap(ii)];
        ydata = [get(hg,'YData') Tmap(ii)];
        set(hg,'XData',xdata,'YData',ydata);
        drawnow();
        
        TwMap.Tout(ii)    = out.T;
        TwMap.Id(ii)      = out.Id;
        TwMap.Iq(ii)      = out.Iq;
        TwMap.Fd(ii)      = out.Fd;
        TwMap.Fq(ii)      = out.Fq;
        TwMap.Tem(ii)     = out.Tem;
        TwMap.Vo(ii)      = out.Vo;
        TwMap.Io(ii)      = out.Io;
        TwMap.Im(ii)      = out.Im;
        TwMap.Iph(ii)     = out.Iph;
        TwMap.Vph(ii)     = out.Vph;
        TwMap.PF(ii)      = out.PF;
        TwMap.P(ii)       = out.P;
        TwMap.Ploss(ii)   = out.Ploss;
        TwMap.Pjs(ii)     = out.Pjs;
        TwMap.PjDC(ii)    = out.PjDC;
        TwMap.PjAC(ii)    = out.PjAC;
        TwMap.Pfe(ii)     = out.Pfe;
        TwMap.Pfes(ii)    = out.Pfes;
        TwMap.Pfer(ii)    = out.Pfer;
        TwMap.Ppm(ii)     = out.Ppm;
        TwMap.Pjr(ii)     = out.Pjr;
        TwMap.Pmech(ii)   = out.Pmech;
        TwMap.Eo(ii)      = out.Eo;
        TwMap.Ife(ii)     = out.Ife;
        TwMap.slip(ii)    = out.slip;
        TwMap.Ir(ii)      = out.Ir;
        TwMap.Rs(ii)      = out.Rs;
        TwMap.dTpp(ii)    = out.dTpp;
        TwMap.eff(ii)     = out.eff;
        TwMap.Idemag(ii)  = out.Idemag;
        TwMap.IHWC(ii)    = out.IHWC;
        TwMap.F0(ii)      = out.F0;
        TwMap.VUGO(ii)    = out.VUGO;
        TwMap.ASCsafe(ii) = out.ASCsafe;
        TwMap.UGOsafe(ii) = out.UGOsafe;
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

%TwMap.dTpp = interp2(motorModel.FluxMap_dq.Id,motorModel.FluxMap_dq.Iq,motorModel.FluxMap_dq.dTpp,TwMap.Id,TwMap.Iq);

TwMap.limits.n    = unique(TwMap.n)';
TwMap.limits.Tmax = max(TwMap.Tout);
TwMap.limits.Tmax(TwMap.limits.Tmax<0) = 0;
TwMap.limits.Tmin = min(TwMap.Tout);
TwMap.limits.Tmin(TwMap.limits.Tmin>0) = 0;

%% plot results

pathname = motorModel.data.pathname;
motName = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'TwMap_' datestr(now,30) '\'];

resFolderOut = [pathname resFolder];

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
figNames{21} = 'ASC safe';
figNames{22} = 'UGO safe';
figNames{23} = 'ASC HWC current';
figNames{24} = 'No load voltage';
figNames{25} = 'Turn-off safe state';

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
if (strcmp(TwData.IronLossFlag,'No')&&(sum(TwData.MechLoss==0)))
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
    set(hfig(ii),'Name',figNames{flagPlot(ii)});
    hax(ii) = axes(...
        'XLim',[TwData.nmin TwData.nmax],...
        'YLim',[TwData.Tmin TwData.Tmax]);
    xlabel('$n$ (rpm)')
    ylabel('$T$ (Nm)')
    switch flagPlot(ii)
        case 1
            title('Torque limits')
            plot(TwMap.limits.n,TwMap.limits.Tmax,'-bx')
            plot(TwMap.limits.n,TwMap.limits.Tmin,'-bx')
        case 2
            title('Efficiency map (p.u.)')
            contourf(TwMap.n,TwMap.T,TwMap.eff,[64:2:86 87:1:100]/100,'ShowText','on');
            colorbar
        case 3
            title('Output power (W)')
            contourf(TwMap.n,TwMap.T,TwMap.P,'ShowText','on');
            colorbar
        case 4
            title('Total loss (W)')
            contourf(TwMap.n,TwMap.T,TwMap.Ploss,'ShowText','on');
            colorbar
        case 5
            title('Stator Joule loss (W)')
            contourf(TwMap.n,TwMap.T,TwMap.Pjs,'ShowText','on');
            colorbar
        case 6
            title('Stator Joule DC loss (W)')
            contourf(TwMap.n,TwMap.T,TwMap.PjDC,'ShowText','on');
            colorbar
        case 7
            title('Stator Joule AC loss (W)')
            contourf(TwMap.n,TwMap.T,TwMap.PjAC,'ShowText','on');
            colorbar
        case 8
            title('Iron loss (W)')
            contourf(TwMap.n,TwMap.T,TwMap.Pfe,'ShowText','on');
            colorbar
        case 9
            title('Stator iron loss (W)')
            contourf(TwMap.n,TwMap.T,TwMap.Pfes,'ShowText','on');
            colorbar
        case 10
            title('Rotor iron loss (W)')
            contourf(TwMap.n,TwMap.T,TwMap.Pfer,'ShowText','on');
            colorbar
        case 11
            title('Permanent magnet loss (W)')
            contourf(TwMap.n,TwMap.T,TwMap.Ppm,'ShowText','on');
            colorbar
        case 12
            title('Phase current map (Apk)')
            contourf(TwMap.n,TwMap.T,TwMap.Io,'ShowText','on');
            colorbar
        case 13
            title('Line voltage map (Vpk)')
            contourf(TwMap.n,TwMap.T,TwMap.Vo,'ShowText','on');
            colorbar
        case 14
            title('Power factor map')
            contourf(TwMap.n,TwMap.T,abs(TwMap.PF),[0.4:0.05:1],'ShowText','on');
            colorbar
        case 15
            title('Mechanical loss map (W)')
            contourf(TwMap.n,TwMap.T,TwMap.Pmech,'ShowText','on');
            colorbar
        case 16
            title('Iron + PM + mechanical loss map (W)')
            contourf(TwMap.n,TwMap.T,TwMap.Pmech+TwMap.Pfes+TwMap.Pfer+TwMap.Ppm,'ShowText','on');
            colorbar
        case 17
            title('Rotor Joule loss map (W)')
            contourf(TwMap.n,TwMap.T,TwMap.Pjr,'ShowText','on');
            colorbar
        case 18
            title('Rotor slip map')
            contourf(TwMap.n,TwMap.T,TwMap.slip,[-1:0.1:1],'ShowText','on');
            colorbar
        case 19
            title('Rotor current map (A)')
            contourf(TwMap.n,TwMap.T,TwMap.Ir,'ShowText','on');
            colorbar
        case 20
            title('Peak-to-peak torque ripple (Nm)')
            contourf(TwMap.n,TwMap.T,TwMap.dTpp,'ShowText','on');
            colorbar
        case 21
            title('Active Short Circuit safe state area')
            plot(TwMap.n(TwMap.ASCsafe==1),TwMap.T(TwMap.ASCsafe==1),'g.')
            plot(TwMap.n(TwMap.ASCsafe==0),TwMap.T(TwMap.ASCsafe==0),'r.')
        case 22
            title('Uncontrolled Generator Operation safe state area')
            plot(TwMap.n(TwMap.UGOsafe==1),TwMap.T(TwMap.UGOsafe==1),'g.')
            plot(TwMap.n(TwMap.UGOsafe==0),TwMap.T(TwMap.UGOsafe==0),'r.')
        case 23
            title('Active Short Circuit Hyper-Worst-Case current (Apk)')
            contourf(TwMap.n,TwMap.T,TwMap.IHWC,'ShowText','on');
            colorbar
        case 24
            title('Line no-load voltage in UGO [Vpk]')
            contourf(TwMap.n,TwMap.T,TwMap.VUGO,'ShowText','on');
            colorbar
        case 25
            title('Turn-off safe state')
            plot(TwMap.n(TwMap.UGOsafe==1),TwMap.T(TwMap.UGOsafe==1),'go','DisplayName','UGO safe')
            plot(TwMap.n(TwMap.UGOsafe==0),TwMap.T(TwMap.UGOsafe==0),'ro','DisplayName','UGO unsafe')
            plot(TwMap.n(TwMap.ASCsafe==1),TwMap.T(TwMap.ASCsafe==1),'g.','DisplayName','ASC safe')
            plot(TwMap.n(TwMap.ASCsafe==0),TwMap.T(TwMap.ASCsafe==0),'r.','DisplayName','ASC unsafe')
            legend('show','Location','northeast');
    end
    if flagPlot(ii)>1
        plot(TwMap.limits.n,TwMap.limits.Tmax,'-k','HandleVisibility','off')
        plot(TwMap.limits.n,TwMap.limits.Tmin,'-k','HandleVisibility','off')
    end
end

ii=ii+1;

hfig(ii) = figure();
figSetting();
set(hfig(ii),'FileName',[pathname resFolder 'Control locus.fig']);
set(hfig(ii),'Name','Control locus');
hax(ii) = axes(...
    'XLim',[min(Id,[],'all') max(Id,[],'all')],...
    'YLim',[min(Id,[],'all') max(Iq,[],'all')],...
    'PlotBoxAspectRatio',[1 1 1]);
xlabel('$i_d$ (A)')
ylabel('$i_q$ (A)')
if strcmp(TwData.Control,'Maximum efficiency')
    title('Maximum efficiency locus')
elseif strcmp(TwData.Control,'MTPA')
    title('Maximum torque per minimum copper loss locus')
end
[c,h] = contour(Id,Iq,Tem,'k','DisplayName','$T$ (Nm)');
clabel(c,h)
plot(TwMap.Id,TwMap.Iq)


%% Save figures
if nargin()==2
    answer = 'No';
    answer = questdlg('Save figures?','Save','Yes','No',answer);
else
    if saveFlag
        answer = 'Yes';
    else
        answer = 'No';
    end
end

if strcmp(answer,'Yes')
    if ~exist([pathname resFolder],'dir')
        mkdir([pathname resFolder]);
    end
    
    save([pathname resFolder 'TwMap.mat'],'motorModel','TwMap','TwData');

    for ii=1:length(hfig)
        savePrintFigure(hfig(ii));
    end
end
