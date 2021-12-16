% Copyright 2019
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

function [demagLimit] = eval_idemag(setup)

clc

if nargin==0
    load LastPath.mat
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
    clear temp
    
    load([setup.pathname setup.filename]);
    if ~isfield(dataSet,'axisType')
        if strcmp(dataSet.TypeOfRotor,'SPM') || strcmp(dataSet.TypeOfRotor,'Vtype')
            dataSet.axisType = 'PM';
        else
            dataSet.axisType = 'SR';
        end
    end
    
    if ~strcmp(dataSet.axisType,dataIn.axisType)
        if strcmp(dataSet.axisType,'PM')
            geo.th0 = geo.th0 + 90;
        else
            geo.th0 = geo.th0 - 90;
        end
    end
end

pathname = setup.pathname;
filename = setup.filename;
tempVect = setup.tempVect;
figFlag  = setup.figFlag;

dataSet.tempPP      = setup.tempVect;

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

IdemagVect=zeros(size(tempVect));

disp('Starting FEMM simulations...')
openfemm(1)
opendocument([pathname filename(1:end-4) '.fem'])

% the PMs group must be changed, from 201 on, based on the number of PMs
% area (no more true from rev433, when PMs are automatically set with
% different groups)
% tmp=geo.BLKLABELS.rotore.xy;
% index=201;
% for ii=1:length(tmp(:,1))
%     if tmp(ii,3)==6
%         mi_selectlabel(tmp(ii,1),tmp(ii,2));
%         mi_setgroup(index);
%         mi_clearselected;
%         index=index+1;
%     end
% end
% clear index

mi_saveas([resFolder filename(1:end-4) '.fem']);
mi_close;

geo0 = geo;
mat0 = mat;
per0 = per;

% i0=calc_io(geo,per);
i0 = per.i0;

maxIter = 20;
tolPU   = 2/100;

if (strcmp(geo.RotType,'SPM')||strcmp(geo.RotType,'Vtype'))
    per.gamma = 180;
else
    per.gamma = 90;
end

for tt=1:length(tempVect)
    
    disp(['Temperature ' int2str(tt) ' of ' int2str(length(tempVect)) ' - ' int2str(tempVect(tt)) ' Celsius degree' ])
    
    nsim = 2;
    Br = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,tempVect(tt));
    Bd = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Bd,tempVect(tt));
    
    Btol = (Br-Bd)*tolPU;
    
    per.tempPP = tempVect(tt);
    per.BrPP   = Br;
    mat.LayerMag.Hc = Br/(4*pi*1e-7*mat.LayerMag.mu);
    
    done=0;
    
    if figFlag
        hfig=figure();
        figSetting();
        xlabel('$I_{demag} [$A$]$')
        ylabel('$B_r$ [$T$]')
        %plot(0,Br,'go','DisplayName',['$B_r=' num2str(Br) '\,T$'])
        %plot(0,Bd,'ro','DisplayName',['$B_d=' num2str(Bd) '\,T$'])
        title(['$\theta_{PM}=' int2str(tempVect(tt)) '^\circ C$'])
        legend('show','Location','northeastoutside');
        hax=gca;
        drawnow
    end
    
    ii=1;
    
    Iiter = nan(1,maxIter);
    Biter = Br*1.1*ones(1,maxIter);
    
    while ~done
        if (ii==1)
            if (tt==1)
                per.overload=0;
            end
        elseif ii==2
            %             if tt==1
            per.overload=1;
            %             else
            %                 BrOld = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,tempVect(tt-1));
            %                 BrNew = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,tempVect(tt));
            %                 BdOld = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,tempVect(tt-1));
            %                 BdNew = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,tempVect(tt));
            %                 per.overload=IdemagVect(tt-1)/i0*(BrNew/BrOld)*(1-BdNew/BdOld);
            %             end
        else
            per.overload = (Iiter(ii-2)+(Bd-Biter(ii-2))*(Iiter(ii-1)-Iiter(ii-2))/(Biter(ii-1)-Biter(ii-2)))/i0;
            if isnan(per.overload)
                per.overload=0;
            end
            if per.overload<0
                per.overload=-0;
            end
            %             per.overload = Iiter(ii-2)/i0;
        end
        
        for pp=1:nsim
            [~,tmpFolder]=createTempDir();
            
            copyfile([resFolder filename(1:end-4) '.fem'],[tmpFolder filename(1:end-4) '.fem']);
            
            per.nsim_singt      = 1;
            per.delta_sim_singt = (0.5*360/(6*geo.q*geo.win.n3phase))*(pp-1);
            
            SOL = simulate_xdeg(geo,per,mat,'idemag',tmpFolder,[filename(1:end-4) '.fem']);
            
            Iiter(ii)=abs(SOL.id+j*SOL.iq);
            if SOL.Bmin<Biter(ii)
                Biter(ii)=SOL.Bmin;
            end
        end
        
        if figFlag
            plot(hax,Iiter(ii)*[1 1],[Bd Biter(ii)],'-o','MarkerIndices',2,'DisplayName',['iteration ' int2str(ii)])
            drawnow
        end
        
        disp([' Iteration ' int2str(ii) ' --> ' num2str(Biter(ii)) ' T @ ' int2str(Iiter(ii)) ' A - ' int2str((Biter(ii)-Bd)/(Br-Bd)*100) '% @ ' int2str(per.overload) ' p.u.'])
        
        if (((Biter(ii)-Bd)<Btol)&&Biter(ii)>Bd) % PMs not demagnetized and over the tolerance
            done=1;
        elseif ((Iiter(ii)==0)&&Biter(ii)<Bd) % PMs demagnetized at zero current
            done = 1;
            Iiter(ii) = 0;
            %             Biter(ii) = NaN;
        elseif ii==maxIter % max number of iteration reached
            done=1;
        elseif ((ii>2)&&(Iiter(ii)==Iiter(ii-1)))
            done = 1;
            Iiter(ii) = 0;
            Biter(ii) = Bd;
        else
            done = 0;
            ii = ii+1;
        end
        %         if ii==maxIter+1
        %             done=1;
        %             Iiter(ii-1) = NaN;
        %             Biter(ii-1) = NaN;
        %         end
    end
    
    %     Iiter=Iiter(Biter>Bd);
    %     Biter=Biter(Biter>Bd);
    Biter=Biter(~isnan(Iiter));
    Iiter=Iiter(~isnan(Iiter));
    if ~isempty(Iiter)
        IdemagVect(tt)=Iiter(end);
    else
        IdemagVect(tt)=NaN;
    end
    eval(['detectionIter.temp_' int2str(tt) '.Iiter=Iiter;'])
    eval(['detectionIter.temp_' int2str(tt) '.Biter=Biter;'])
    
    if figFlag
        saveas(hfig,[resFolder 'iterations_' int2str(tempVect(tt)) 'C.fig']);
    end
    
end

disp('FEMM simulations done!')

figure()
figSetting();
xlabel('$\theta_{PM}$ [$^\circ C$]')
ylabel('$I_{demag}$ [$A$]')
plot(tempVect,IdemagVect,'-ro');
title('Demagnetization limit')
saveas(gcf,[resFolder 'Demagnetization limit.fig'])

demagLimit.temperature = tempVect;
demagLimit.current     = IdemagVect;
demagLimit.setup       = setup;
demagLimit.iterations  = detectionIter;

save([resFolder 'demagCurve.mat'],'demagLimit','setup','dataSet');


