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

function [demagLimit,resFolder] = eval_idemag_vol(setup)

clc

if nargin==0
    if exist('LastPath.mat','file')
        load LastPath.mat
    else
        pathname = cd;
    end
    [filename,pathname,~]=uigetfile([pathname '\.mat'],'Select a machine');
    
    load([pathname filename])
    
    % inputs
    prompt={
        'Temperature vector'
        'Figure flag (1=yes/0=no)'};
    name='Demag detection';
    numlines=1;
    answer={mat2str(mat.LayerMag.temp.temp),...
        '1'};
    
    answer=inputdlg(prompt,name,numlines,answer);
    
    setup.filename = filename;
    setup.pathname = pathname;
    setup.tempVect = eval(answer{1});
    setup.figFlag  = eval(answer{2});
else
    
    temp = setup;
    clear setup
    setup.filename = temp.currentfilename;
    setup.pathname = temp.currentpathname;
    setup.tempVect = temp.tempPP;
    setup.figFlag  = 1;
    
    
    load([setup.pathname setup.filename]);

    dataSet.axisType = temp.axisType;
    dataSet.tempPP   = setup.tempVect;
    clear temp

    if ~isfield(dataSet,'axisType')
        if strcmp(dataSet.TypeOfRotor,'SPM') || strcmp(dataSet.TypeOfRotor,'Vtype')
            dataSet.axisType = 'PM';
        else
            dataSet.axisType = 'SR';
        end
    end

    if ~strcmp(geo.axisType,dataSet.axisType)
        %geo.axisType = dataIn.axisType;
        if strcmp(dataSet.axisType,'PM')
            geo.th0 = geo.th0 - 90;
        else
            geo.th0 = geo.th0 + 90;
        end
    end

    


end

pathname = setup.pathname;
filename = setup.filename;
tempVect = setup.tempVect;
figFlag  = setup.figFlag;


per.nsim_singt      = dataSet.NumOfRotPosPP;
per.delta_sim_singt = dataSet.AngularSpanPP;
per.tempPP          = tempVect;

outFolder = [filename(1:end-4) '_results\FEA results\'];
if ~exist([pathname outFolder],'dir')
    mkdir([pathname outFolder]);
end

resFolder = ['demagCurve - ' datestr(now,30) '\'];
mkdir([pathname outFolder],resFolder);
resFolder = [pathname outFolder resFolder];

save([resFolder 'demagCurveSetup.mat'],'setup');

IdemagVect   = zeros(size(tempVect));
dPMdemagVect = zeros(size(tempVect));
BminVect     = zeros(size(tempVect));

disp('Starting FEMM simulations...')
openfemm(1)
opendocument([pathname filename(1:end-4) '.fem'])

mi_saveas([resFolder filename(1:end-4) '.fem']);
mi_close;

i0 = per.i0;

maxIter = 500;
dPMtarget = 0.01;
tolPU = 0.5;
IstepLim = i0/100;

if strcmp(dataSet.axisType,'PM')
    per.gamma = 180;
else
    per.gamma = 90;
end

for tt=1:length(tempVect)
    
    disp(['Temperature ' int2str(tt) ' of ' int2str(length(tempVect)) ' - ' int2str(tempVect(tt)) ' Celsius degree' ])
    
    %nsim = 2;
    Br = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,tempVect(tt));
    Bd = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Bd,tempVect(tt));
    
    %Btol = (Br-Bd)*tolPU;
    
    per.tempPP = tempVect(tt);
    per.BrPP   = Br;
    mat.LayerMag.Hc = Br/(4*pi*1e-7*mat.LayerMag.mu);
    
    done=0;
    
    if figFlag
        hfig=figure();
        figSetting();
        xlabel('$I_{demag} [$A$]$')
        ylabel('PM demag [$\%$]')
        title(['$\theta_{PM}=' int2str(tempVect(tt)) '^\circ$C'])
        legend('show','Location','northeastoutside');
        set(gca,'YLim',[0 100],'YTick',0:10:100)
        hax=gca;
        drawnow
    end
    
    ii=1;
    
    Iiter   = nan(1,maxIter);
    dPMiter = zeros(1,maxIter);
    Biter   = zeros(1,maxIter);
    Titer   = zeros(1,maxIter);
    Istep = i0;

    while ~done
        if (ii==1)
                per.overload=0;
                Istep = i0;
        elseif ii==2
            if tt==1
                per.overload=1;
                Istep = i0;
            else
                per.overload = IdemagVect(tt-1)/i0;
                Istep = per.overload*i0;
                maxIter = 20;
            end
        else
            if max(dPMiter,[],'omitnan')<dPMtarget
                per.overload = max(Iiter,[],'omitnan')/i0*2;
                Istep = per.overload*i0;
            elseif min(dPMiter,[],'omitnan')>dPMtarget
                per.overload = max(Iiter,[],'omitnan')/i0/2;
                Istep = per.overload*i0;
            else
                Istep = Istep/2;
                if dPMiter(ii-1)>dPMtarget
                    per.overload = (Iiter(ii-1)-Istep)/i0;
                else
                    per.overload = (Iiter(ii-1)+Istep)/i0;
                end
            end
            if isnan(per.overload)
                per.overload=0;
            end
            if per.overload<0
                per.overload=-0;
            end
        end
        
        [~,~,~,out,~] = FEMMfitness([],geo,per,mat,'idemag',[pathname filename(1:end-4) '.fem']);
        
        Iiter(ii) = abs(out.id+j*out.iq);
        if out.dPM>=dPMiter(ii)
            dPMiter(ii) = out.dPM;
            Biter(ii)   = out.Bmin;
            Titer(ii)   = out.T;
        end
        
        if figFlag
            plot(hax,Iiter(ii)*[1 1],[0 dPMiter(ii)]*100,'-o','MarkerIndices',2,'DisplayName',['iteration ' int2str(ii)])
            drawnow
        end
        
        disp([' Iteration ' int2str(ii) ' --> ' int2str(dPMiter(ii)*100) '% @ ' int2str(Iiter(ii)) ' A'])
        
        if (dPMiter(ii)<dPMtarget && dPMiter(ii)>(dPMtarget*(1-tolPU))) % PMs demagnetized below the tolerance
            done=1;
        elseif ((Iiter(ii)==0) && dPMiter(ii)>(dPMtarget*(1-tolPU))) % PMs demagnetized at zero current
            done = 1;
            Iiter(ii) = 0;
        elseif ii==maxIter % max number of iteration reached
            done=1;
        elseif ((ii>2)&&(Iiter(ii)==Iiter(ii-1))) % two iterations at the same current level
            done = 1;
            Iiter(ii) = 0;
            dPMiter(ii) = dPMtarget;
        elseif ((ii>2)&&abs(Iiter(ii)-Iiter(ii-1))<IstepLim)
            done = 1;
            Iiter(ii) = NaN;
        else
            done = 0;
            ii = ii+1;
        end
    end
    
    dPMiter = dPMiter(~isnan(Iiter));
    Biter   = Biter(~isnan(Iiter));
    Iiter   = Iiter(~isnan(Iiter));
    if ~isempty(Iiter)
        IdemagVect(tt)   = Iiter(end);
        dPMdemagVect(tt) = dPMiter(end);
        BminVect(tt)     = Biter(end);
    else
        IdemagVect(tt)   = NaN;
        dPMdemagVect(tt) = NaN;
        BminVect(tt)     = NaN;
    end
    eval(['detectionIter.temp_' int2str(tt) '.Iiter=Iiter;'])
    eval(['detectionIter.temp_' int2str(tt) '.dPMiter=dPMiter;'])
    eval(['detectionIter.temp_' int2str(tt) '.Biter=Biter;'])

    if figFlag
        saveas(hfig,[resFolder 'iterations_' int2str(tempVect(tt)) 'C.fig']);
    end
    
end

disp('FEMM simulations done!')

figure()
figSetting();
xlabel('$\theta_{PM}$ [$^\circ$C]')
ylabel('$I_{demag}$ [A]')
plot(tempVect,IdemagVect,'-ro');
title('Demagnetization limit @ 1\%')
saveas(gcf,[resFolder 'Demagnetization limit.fig'])

figure()
figSetting();
xlabel('$\theta_{PM}$ [$^\circ C$]')
ylabel('PM demag [$\%$]')
plot(tempVect,dPMdemagVect,'-ko');
plot(tempVect,dPMtarget*ones(size(tempVect)),'-g','DisplayName','Target')
plot(tempVect,dPMtarget*(1+tolPU)*ones(size(tempVect)),'-r','DisplayName','Target','HandleVisibility','off')
plot(tempVect,dPMtarget*(1-tolPU)*ones(size(tempVect)),'-r','DisplayName','Target','HandleVisibility','off')
title(['Demagnetization limit @ ' int2str(dPMtarget) '\%'])
saveas(gcf,[resFolder 'Demagnetization limit vol.fig'])

figure()
figSetting();
xlabel('$\theta_{PM}$ [$^\circ$C]')
ylabel('$B_{min}$ [T]')
plot(tempVect,BminVect,'-o','DisplayName','$B_{min}$');
plot(tempVect,interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Bd,tempVect),'-','DisplayName','$B_d$')
plot(tempVect,interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,tempVect),'-','DisplayName','$B_r$')
title('Minimum PM flux density')
saveas(gcf,[resFolder 'Demagnetization limit fluxDensity.fig'])

demagLimit.temperature = tempVect;
demagLimit.current     = IdemagVect;
demagLimit.dPM         = dPMdemagVect;
demagLimit.Bmin        = BminVect;
demagLimit.setup       = setup;
demagLimit.iterations  = detectionIter;

save([resFolder 'demagCurve.mat'],'demagLimit','setup','dataSet');


