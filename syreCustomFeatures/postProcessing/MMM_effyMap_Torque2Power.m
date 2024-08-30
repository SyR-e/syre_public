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

function [PwMap] = MMM_effyMap_Torque2Power(TwMap,pathname)

if nargin()==0
    [filename,pathname] = uigetfile([cd '\*.mat'],'Select the TwMap file');
    load([pathname filename],'TwMap')
end

if isfield(TwMap,'T_top_W')

    limits.n = unique(TwMap.n)';
    limits.Tmax = TwMap.T_top_W;
    limits.Tmin = TwMap.T_bot_W;
    
    TwMap = rmfield(TwMap,'T_top_W');
    TwMap = rmfield(TwMap,'T_bot_W');
else
    limits = TwMap.limits;

    TwMap = rmfield(TwMap,'limits');
end

limits.Tmax = smooth(limits.Tmax)';
limits.Tmin = smooth(limits.Tmin)';

limits.Pmax = limits.Tmax.*limits.n*pi/30;
limits.Pmin = limits.Tmin.*limits.n*pi/30;

fields = fieldnames(TwMap);

nvect = unique(TwMap.n);
Tvect = unique(TwMap.T);

Pvect = linspace(min(TwMap.P(:))*1.1,max(TwMap.P(:))*1.1,numel(Tvect));

[nmap,Pmap] = meshgrid(nvect,Pvect);

x = TwMap.n;
y = TwMap.T;

for ii=1:length(fields)
    if strcmp(fields{ii},'n')
        PwMap.n = nmap;
    elseif strcmp(fields{ii},'P')
        PwMap.P = Pmap;
    else
        z = TwMap.(fields{ii});
        if (min(z(:))==0 && max(z(:))==0)
            PwMap.(fields{ii}) = zeros(size(nmap));
        elseif sum(isnan(z(:)))==numel(z)
            PwMap.(fields{ii}) = zeros(size(nmap));
        else
            xFilt = x(~isnan(z));
            yFilt = y(~isnan(z));
            zFilt = z(~isnan(z));
            fInt = scatteredInterpolant(xFilt(:),yFilt(:),zFilt(:),'linear','linear');
            zTmp = fInt(nmap,Pmap./(nmap*pi/30));
            zTmp(Pmap>limits.Pmax) = NaN;
            zTmp(Pmap<limits.Pmin) = NaN;
            PwMap.(fields{ii}) = zTmp;
        end
    end
end

PwMap.ASCsafe = round(PwMap.ASCsafe);
PwMap.UGOsafe = round(PwMap.UGOsafe);

% PwMap.ASCsafe(PwMap.ASCsafe<=0.5) = 0;
% PwMap.ASCsafe(PwMap.ASCsafe>0.5) = 1;
% PwMap.UGOsafe(PwMap.UGOsafe<=0.5) = 0;
% PwMap.UGOsafe(PwMap.UGOsafe>0.5) = 1;

PwMap.limits = limits;

%% plot figures

resFolder = 'PowerSpeed_plane/';
mkdir([pathname resFolder]);


figNames{1}  = 'Power limits';
figNames{2}  = 'Efficiency map';
figNames{3}  = 'Output torque map';
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

for ii=1:length(figNames)
    hfig(ii) = figure();
    figSetting();
    set(hfig(ii),'FileName',[pathname resFolder figNames{ii} '.fig']);
    set(hfig(ii),'Name',figNames{ii});
    hax(ii) = axes(...
        'XLim',[min(PwMap.n(:)) max(PwMap.n(:))],...
        'YLim',[min(PwMap.P(:)) max(PwMap.P(:))]);
    xlabel('$n$ [rpm]')
    ylabel('$P$ [W]')
    switch ii
        case 1
            title('Power limits')
            plot(PwMap.limits.n,PwMap.limits.Pmax,'-bx')
            plot(PwMap.limits.n,PwMap.limits.Pmin,'-bx')
        case 2
            title('Efficiency map [p.u.]')
            contourf(PwMap.n,PwMap.P,PwMap.eff,[64:2:86 87:1:100]/100,'ShowText','on');
            colorbar
        case 3
            title('Output torque [Nm]')
            contourf(PwMap.n,PwMap.P,PwMap.T,'ShowText','on');
            colorbar
        case 4
            title('Total loss [W]')
            contourf(PwMap.n,PwMap.P,PwMap.Ploss,'ShowText','on');
            colorbar
        case 5
            title('Stator Joule loss [W]')
            contourf(PwMap.n,PwMap.P,PwMap.Pjs,'ShowText','on');
            colorbar
        case 6
            title('Stator Joule DC loss [W]')
            contourf(PwMap.n,PwMap.P,PwMap.PjDC,'ShowText','on');
            colorbar
        case 7
            title('Stator Joule AC loss [W]')
            contourf(PwMap.n,PwMap.P,PwMap.PjAC,'ShowText','on');
            colorbar
        case 8
            title('Iron loss [W]')
            contourf(PwMap.n,PwMap.P,PwMap.Pfe,'ShowText','on');
            colorbar
        case 9
            title('Stator iron loss [W]')
            contourf(PwMap.n,PwMap.P,PwMap.Pfes,'ShowText','on');
            colorbar
        case 10
            title('Rotor iron loss [W]')
            contourf(PwMap.n,PwMap.P,PwMap.Pfer,'ShowText','on');
            colorbar
        case 11
            title('Permanent magnet loss [W]')
            contourf(PwMap.n,PwMap.P,PwMap.Ppm,'ShowText','on');
            colorbar
        case 12
            title('Phase current map [Apk]')
            contourf(PwMap.n,PwMap.P,PwMap.Io,'ShowText','on');
            colorbar
        case 13
            title('Line voltage map [Vpk]')
            contourf(PwMap.n,PwMap.P,PwMap.Vo,'ShowText','on');
            colorbar
        case 14
            title('Power factor map')
            contourf(PwMap.n,PwMap.P,abs(PwMap.PF),[0.4:0.05:1],'ShowText','on');
            colorbar
        case 15
            title('Mechanical loss map')
            contourf(PwMap.n,PwMap.P,PwMap.Pmech,'ShowText','on');
            colorbar
        case 16
            title('Iron + PM + mechanical loss map')
            contourf(PwMap.n,PwMap.P,PwMap.Pmech+PwMap.Pfes+PwMap.Pfer+PwMap.Ppm,'ShowText','on');
            colorbar
        case 17
            title('Rotor Joule loss map')
            contourf(PwMap.n,PwMap.P,PwMap.Pjr,'ShowText','on');
            colorbar
        case 18
            title('Rotor slip map')
            contourf(PwMap.n,PwMap.P,PwMap.slip,[-1:0.1:1],'ShowText','on');
            colorbar
        case 19
            title('Rotor current map [A]')
            contourf(PwMap.n,PwMap.P,PwMap.Ir,'ShowText','on');
            colorbar
        case 20
            title('Peak-to-peak torque ripple [Nm]')
            contourf(PwMap.n,PwMap.P,PwMap.dTpp,'ShowText','on');
            colorbar
        case 21
            title('Active Short Circuit safe state area')
            plot(PwMap.n(PwMap.ASCsafe==1),PwMap.P(PwMap.ASCsafe==1),'g.')
            plot(PwMap.n(PwMap.ASCsafe==0),PwMap.P(PwMap.ASCsafe==0),'r.')
        case 22
            title('Uncontrolled Generator Operation safe state area')
            plot(PwMap.n(PwMap.UGOsafe==1),PwMap.P(PwMap.UGOsafe==1),'g.')
            plot(PwMap.n(PwMap.UGOsafe==0),PwMap.P(PwMap.UGOsafe==0),'r.')
        case 23
            title('Active Short Circuit Hyper-Worst-Case current [Apk]')
            contourf(PwMap.n,PwMap.P,PwMap.IHWC,'ShowText','on');
            colorbar
        case 24
            title('Line no-load voltage in UGO [Vpk]')
            contourf(PwMap.n,PwMap.P,PwMap.VUGO,'ShowText','on');
            colorbar
    end

    if ii>1
        plot(PwMap.limits.n,PwMap.limits.Pmax,'-k')
        plot(PwMap.limits.n,PwMap.limits.Pmin,'-k')
    end
end

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
    
    save([pathname resFolder 'PwMap.mat'],'PwMap');

    for ii=1:length(hfig)
        savePrintFigure(hfig(ii));
    end
end

if nargout()==0
    clear PwMap
end






