% Copyright 2014
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

function evalParetoFront(filename,dataSet)

% re-evaluates all the machines of the Pareto front in 15 rotor
% positions (instead of the 5 positions with random offset used by MODE)

% pivot_cost: selects how to sort the machines of the front
% pivot_cost = 1 (default), COST_sorted is ordered according to the values
% of first cost function (normally T)
% pivot_cost = 2, COST_sorted is ordered according to the values of the
% second cost function (normally dT) .. and so on


pivot_cost = 1;

syreRoot = fileparts(which('syre.png')); %OCT
if nargin<1
    [filename, pathname_ini] = uigetfile('results\*.mat', 'Pick a file');
    load([pathname_ini filename]);
    if ~exist('dataSet')
        dataSet=build_dataSet(geo0,per);
        [dataSet,geo0,per,mat] = back_compatibility(dataSet,geo0,per,1);
    end
    save('dataSet','dataSet');
else
    if isoctave()
        pathname_ini=[syreRoot filesep 'results' filesep];
        load([pathname_ini filename]);
    else
        pathname_ini=[syreRoot filesep 'results' filesep];
        filename=[filename,'.mat'];
        load([pathname_ini filename]);
    end
end

if ~exist('mat')
    [dataSet,geo0,per,mat] = back_compatibility(dataSet,geo0,per,1);
    [bounds, objs, geo0, per, mat] = data0(dataSet);
    per.objs=objs;
end

% dir_name = strrep(filename,'end_','');
dir_name = strrep(filename,'.mat','');
pathname_res = [pathname_ini dir_name '\'];
[~,MESSAGE,~] = mkdir(pathname_res);

eval_type = 'singt';

if isempty(MESSAGE)
    runcase = 'Yes';
else
    runcase = questdlg('Overwrite existing files??','Warning','Yes','No','No');
end
if ~strcmp(runcase,'Yes')
    return
end

geo=geo0;       % assign to geo intial geometric data (same in data0)
per_temp=per;   % re-assign because matlabpool don't work...
mat=mat;
COST = [];
[PSet,iA] = unique(OUT.PSet,'rows');
x = PSet;

% RQ variables, set of optimization inputs

% for j = 1:length(geo0.RQnames)
%     nameTemp{j} = upper(eval(['geo0.RQnames{j}']));
%     %     if strcmp(nameTemp{j},'DALPHA')
%     %         k = k +1;
%     %         nameTemp{j} = [nameTemp{j} num2str(k)];
%     %     end
%     %     if strcmp(nameTemp{j},'HC')
%     %         y = y +1;
%     %         nameTemp{j} = [nameTemp{j} num2str(y)];
%     %     end
%     %     if strcmp(nameTemp{j},'DX')
%     %         t = t +1;
%     %         nameTemp{j} = [nameTemp{j} num2str(t)];
%     %     end
%     %     if strcmp(nameTemp{j},'BR')
%     %         u = u +1;
%     %         nameTemp{j} = [nameTemp{j} num2str(u)];
%     %     end
%     %     if strfind(nameTemp{j},'PMDIM')
%     %         u = u +1;
%     %         nameTemp{j} = [nameTemp{j}(1:5) int2str(u)];
%     %     end
%     %     if strcmp(nameTemp{j},'BETAPMSHAPE')
%     %         u = u +1;
%     %         nameTemp{j} = [nameTemp{j} num2str(u)];
%     %     end
%     %  eval([nameTemp{j} ' = zeros(size(x,1),1)']);
% end

% goals
T   = zeros(size(x,1),1);
dT  = zeros(size(x,1),1);
MCu = zeros(size(x,1),1);
MPM = zeros(size(x,1),1);
% per_temp.MechStressOptCheck = 1;

Tini=now();
parfor m = 1:size(x,1)
    
    Tstart=now;
    mot_index = num2str(m);
    if m<10
        mot_index = ['0' mot_index];
    end
    disp(['Evaluating machine ' mot_index ' of ' num2str(size(x,1))])
    
    % FEA evaluates each machine of the front
    
    [cost{m},geometry{m},material{m},output{m},~] = FEMMfitness(x(m,:)',geo,per_temp,mat,eval_type,strrep(filename,'.mat','.fem'));
    
    
    if per_temp.MechStressOptCheck
        warning 'off'
        mkdir([pathname_res '\MechStress']);
        warning 'on'
        figure();
        figSetting();
        set(gca,'DataAspectRatio',[1 1 1]);
        pdeplot(output{m}.structModel,'XYData',output{m}.sVonMises/1e6,'ZData',output{m}.sVonMises/1e6)
        colormap turbo
        view(2)
        title(['Von Mises Stress [MPa] - mot\_' mot_index])
        sVM_max = max(output{m}.sVonMises);
        set(gca,'CLim',[0 sVM_max]/1e6)
        set(gcf,'FileName',[pathname_res '\MechStress\mot_' mot_index '_mech.fig'])
        savePrintFigure(gcf)
        close
    end
    % debug
    %     cost{m} = x(m,end-1:end);
    %     geometry{m} = geo;
    % debug
    
    Tend = now;
    EvalTime = datevec(Tend-Tstart);
    endTime=datenum(EvalTime)*size(x,1)/4+Tini;
    disp(['Evaluated machine ' mot_index ' of ' num2str(size(x,1)) ' in ' int2str(EvalTime(4)) ' h ' num2str(EvalTime(5)) ' min ' num2str(round(EvalTime(6))) 'sec']);
    disp(['End of Pareto evaluation : ' datestr(endTime)]);
end



per_temp.MechStressOptCheck = 0;

disp('Saving Pareto front machines...')

% remove optimization checks from dataSet of the saved machines
dataSet.Dalpha1BouCheck        = 0;
dataSet.DalphaBouCheck         = 0;
dataSet.hcBouCheck             = 0;
dataSet.DxBouCheck             = 0;
dataSet.GapBouCheck            = 0;
dataSet.BrBouCheck             = 0;
dataSet.AirgapRadiusBouCheck   = 0;
dataSet.ToothWidthBouCheck     = 0;
dataSet.ToothLengthBouCheck    = 0;
dataSet.StatorSlotOpenBouCheck = 0;
dataSet.ToothTangDepthBouCheck = 0;
dataSet.BetaPMshapeBouCheck    = 0;
dataSet.ThetaFBSBouCheck       = 0;
dataSet.PMdimBouCheck          = 0;
dataSet.GammaBouCheck          = 0;
dataSet.RadRibBouCheck         = 0;
dataSet.CentralShrinkBouCheck  = 0;
dataSet.RadShiftInnerBouCheck  = 0;
dataSet.MechStressOptCheck     = 0;
dataSet.FilletRad1BouCheck     = 0;
dataSet.FilletRad2BouCheck     = 0;
dataSet.FilletTan1BouCheck     = 0;
dataSet.FilletTan2BouCheck     = 0;

dataSet.RQ                     = [];
dataSet.RQnames                = [];

for m=1:size(x,1)
    clear geo
    mot_index = num2str(m);
    disp(['Saving machine ' mot_index ' of ' num2str(size(x,1))])
    if m<10
        mot_index = ['0' mot_index];
    end
    
    geo=geometry{m};
    mat=material{m};
    out=output{m};
    
    %     openfemm(1)
    [geo,gamma,mat] = interpretRQ(x(m,:),geo,mat);
    [geo,mat] = draw_motor_in_FEMM(geo,mat,pathname_res,['mot_'  mot_index '.fem']);
    
%     if isoctave() %OCT
%         %         (to be debugged
%         %  - mot0.fem does not exist any longer
%         %  - pathname changed
%         current_dir=strcat(syreRoot, '\mot0.fem');
%         movefile (current_dir, pathname);
%         old_name= strcat(pathname, '\mot0.fem');
%         new_name= strcat(pathname, '\mot_', mot_index, '.fem');
%         rename(old_name,new_name);
%         clear old_name new_name current_dir
%     end
    
    % Save data geometry mot
    geo.RQ = x(m,:);
    
%     if isoctave()                       %OCT
%         dataSet.AirGapThickness = geo.g; % airgap thickness
%         dataSet.AirGapThickness = round(dataSet.AirGapThickness .*100) ./100;
%         dataSet.AirGapRadius = geo.r; % machine airgap radius
%         dataSet.AirGapRadius = round(dataSet.AirGapRadius .*100) ./100;
%         dataSet.ToothLength = geo.lt; % tooth length
%         dataSet.ToothLength = round(dataSet.ToothLength .*100) ./100;
%         dataSet.StatorSlotOpen = geo.acs; % stator slot open in [p.u.]
%         dataSet.StatorSlotOpen = round(dataSet.StatorSlotOpen .*100) ./100;
%         dataSet.ToothWidth = geo.wt; % Bgap/Btooth (determines tooth width^-1, yoke width^-1)
%         dataSet.ToothWidth = round(dataSet.ToothWidth .*100) ./100;
%         dataSet.ToothTangDepth = geo.ttd; % tooth tang depth [mm]
%         dataSet.ToothTangDepth = round(dataSet.ToothTangDepth .*100) ./100;
%         dataSet.Br = round(mat.LayerMag.Br .*10000) ./10000; % Br
%         dataSet.Br = round(dataSet.Br .*100) ./100;
%         dataSet.ALPHApu = geo.dalpha_pu;
%         dataSet.ALPHApu = round(dataSet.ALPHApu .*100) ./100;
%         dataSet.PMdim   = round(geo.PMdim.*100)./100;
%         dataSet.betaPMshape = round(geo.betaPMshape.*100)./100;
%     else
%         dataSet.AirGapThickness = round(geo.g,2);
%         dataSet.AirGapRadius = round(geo.r,2);
%         dataSet.ToothLength = round(geo.lt,2);
%         dataSet.StatorSlotOpen = round(geo.acs,2);
%         dataSet.ToothWidth = round(geo.wt,2);
%         dataSet.ToothTangDepth = round(geo.ttd,2);
%         dataSet.Br = round(mat.LayerMag.Br,4); % Br
%         dataSet.Br = round(dataSet.Br,2);
%         dataSet.ALPHApu = round(geo.dalpha_pu,2);
%         dataSet.PMdim   = round(geo.PMdim,2);
%         dataSet.betaPMshape = round(geo.betaPMshape,2);
%     end
%     [dataSet] = SaveInformation(geo,mat,dataSet);

    dataSet.AirGapThickness = round(geo.g.*100)./100;
    dataSet.AirGapRadius    = round(geo.r.*100)./100;
    dataSet.ToothLength     = round(geo.lt.*100)./100;
    dataSet.StatorSlotOpen  = round(geo.acs.*100)./100;
    dataSet.ToothWidth      = round(geo.wt.*100)./100;
    dataSet.ToothTangDepth  = round(geo.ttd.*100)./100;
    dataSet.Br              = round(mat.LayerMag.Br.*100)./100;
    dataSet.ALPHApu         = round(geo.dalpha_pu.*100)./100;
    dataSet.HCpu            = round(geo.hc_pu.*100)./100;
    dataSet.DepthOfBarrier  = round(geo.dx.*100)./100;
    dataSet.betaPMshape     = round(geo.betaPMshape.*100)./100;
    dataSet.thetaFBS        = round(geo.th_FBS*180/pi*100)/100;
    dataSet.TanRibEdit      = round(geo.pontT.*100)./100;
    dataSet.RadRibEdit      = round(geo.pontR.*100)./100;
    dataSet.RotorFilletIn   = round(geo.RotorFillet1.*100)./100;
    dataSet.RotorFilletOut  = round(geo.RotorFillet2.*100)./100;
    dataSet.RotorFilletTan1 = round(geo.RotorFilletTan1.*100)./100;
    dataSet.RotorFilletTan2 = round(geo.RotorFilletTan2.*100)./100;
    dataSet.PMdim           = round(geo.PMdim.*100)./100;
    dataSet.DepthOfBarrier  = round(geo.dx.*100)./100;
    dataSet.CentralShrink   = round(geo.hcShrink.*100)./100;
    dataSet.RadShiftInner   = round(geo.dxIB.*100)./100;
    dataSet.betaPMshape     = round(geo.betaPMshape.*100)./100;
    dataSet.gammaPP         = round(gamma.*100)./100;

    dataSet.RQ              = geo.RQ;

    if strcmp(geo.RotType,'SPM')
        dataSet.ThicknessOfPM = geo.hc_pu;
    else
        dataSet.HCpu = geo.hc_pu;
    end
%     if isoctave()                                           %OCT
%         dataSet.HCpu = round(dataSet.HCpu .*100) ./100;
%         dataSet.DepthOfBarrier = geo.dx;                    % the depth of the barriers radial-wise in per unit
%         dataSet.DepthOfBarrier = round(dataSet.DepthOfBarrier .*100) ./100;
%     else
%         dataSet.HCpu = round(dataSet.HCpu,2);
%         dataSet.DepthOfBarrier = round(geo.dx,2);
%     end
%     dataSet.RQ = geo.RQ;
    
    if isoctave()  %OCT
        name_file = strcat(pathname, '\mot_', mot_index, '.mat');
        save ('-mat7-binary', name_file,'geo','cost','per','dataSet','mat');
        clear name_file
    else
        save([pathname_res '\mot_' mot_index '.mat'],'geo','cost','per','dataSet','mat','out');
    end
    
    %%%%%%%%%
    COST = [COST; cost{m}];
    if m<10
        mot_index = ['0' mot_index];
    end
    
end

copyfile(strcat(pathname_ini, filename),strcat(pathname_res, filename));
n_mot = size(COST,1);
[~,I] = sort(COST(:,pivot_cost));
name_case = strrep(filename,'.mat','');

% Pareto front
if OUT.Param.NOBJ==1
    figSetting
    xlabel('Generation number')
    ylabel(['log(' geo0.OBJnames{1} ')'])
    h=gcf();
    if isoctave()   %OCT
        fig_name=strcat(pathname_res, '\goalVSgenerations - ', name_case);
        hgsave(h,[fig_name]);
    else
        saveas(gcf,[pathname_res '\goalVSgenerations - ' name_case '.fig'])
    end
    
elseif OUT.Param.NOBJ==2
    close
    figure(1), clf
    figSetting
    for ii=1:n_mot
        plot(COST(ii,1),COST(ii,2),'x'),
        text(COST(ii,1)+0.1,COST(ii,2)+0.1,num2str(ii));
    end
    
    xlabel(geo0.OBJnames{1})
    ylabel(geo0.OBJnames{2})
    
    [front,idx] = FastParetoEstimation(x,COST);
    nmot_actual_front=find(idx);
    
    for ii=1:length(nmot_actual_front)
        C1_front(ii)=front(ii,end-1);C2_front(ii)=front(ii,end);
        plot(C1_front(ii),C2_front(ii),'ro','LineWidth',2),
        text(C1_front(ii)+0.1,C2_front(ii)+0.1,num2str(nmot_actual_front(ii)));
    end
    [~,ii_tfs]=sort(C1_front);
    plot(C1_front(ii_tfs),C2_front(ii_tfs),'r','LineWidth',2)
    title('Pareto front')
    h=gcf();
    if isoctave()
        fig_name=strcat(pathname_res, 'Pareto - ', name_case);
        hgsave(h,[fig_name]);
    else
        saveas(gcf,[pathname_res 'Pareto - ' name_case '.fig']);
    end
    
    
elseif OUT.Param.NOBJ==3
    close
    figure(1)
    figSetting
    view(3)
    for ii = 1:n_mot
        plot3(COST(ii,1),COST(ii,2),COST(ii,3),'x');
        text(COST(ii,1)+0.1,COST(ii,2)+0.1,COST(ii,3)+0.1,int2str(ii));
    end
    xlabel(geo0.OBJnames{1});
    ylabel(geo0.OBJnames{2});
    zlabel(geo0.OBJnames{3});
    [front,idx] = FastParetoEstimation(x,COST);
    nmot_actual_front=find(idx);
    
    for ii=1:length(nmot_actual_front)
        plot3(front(ii,end-2),front(ii,end-1),front(ii,end),'ro','LineWidth',2),
        text(front(ii,end-2)+0.1,front(ii,end-1)+0.1,front(ii,end)+0.1,int2str(nmot_actual_front(ii)));
    end
    title('Pareto front')
    h=gcf();
    if isoctave()
        fig_name=strcat(pathname_res, '\Pareto - ', name_case);
        hgsave(h,[fig_name]);
    else
        saveas(gcf,[pathname_res '\Pareto - ' name_case '.fig']);
    end
    
end

% bw = 0.7;
% figure()
% 
% 
% if OUT.Param.NOBJ<4
%     subplot(OUT.Param.NOBJ,1,1)
%     index1 = OUT.Param.NOBJ;
%     index2 = 1;
% elseif OUT.Param.NOBJ == 4
%     subplot(2,2,1)
%     index1 = 2;
%     index2 = 2;
% elseif OUT.Param.NOBJ > 4
%     subplot(3,2,1)
%     index1 = 3;
%     index2 = 2;
% end
% 
% figSetting(11*index2,5*index1)
% sgtitle('Objectives')
% 
% for ii=1:OUT.Param.NOBJ
%     subplot(index1,index2,ii)
%     grid on
%     bar(1:length(I),sign(per.objs(ii,1))*COST(I,ii),bw,'FaceColor',[0 0.7 0])
%     set(gca,'XLim',[0 n_mot+1],'xTick',1:1:length(I),'xTickLabel',I);
%     xlabel('Machine Number');
%     ylabel(geo0.OBJnames{ii})
% end
% 
% if isoctave()
%     fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_Objectives');
%     hgsave(h,[fig_name]);
% else
%     saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_Objectives.fig']);
% end

if OUT.eval_type=='MO_OA'
    if OUT.Param.NOBJ>1
        save([pathname_res '\' filename],'front','idx','-append');
    end
end

for ii=1:OUT.Param.NOBJ
    hfig = figure();
    figSetting();
    xlabel('Machine Number')
    ylabel(geo0.OBJnames{ii})
    title(['Objective: ' geo0.OBJnames{ii}])
    set(gca,'XLim',[0 length(I)+1],'XTick',1:1:length(I),'XTickLabel',I);
    bar(1:1:length(I),abs(COST(I,ii)),'BarWidth',1,'FaceColor',[0 0 1],'EdgeColor',[0 0 0.5]);
    figName = [pathname_res 'BarChart - ' name_case '_sort' int2str(pivot_cost) '_' geo0.OBJnames{ii} '.fig'];
    if isoctave()
        hgsave(hfig,figName);
    else
        saveas(hfig,figName);
    end
end

for ii=1:length(geo0.RQnames)
    figName = geo0.RQnames{ii};
    figName(figName=='(') = '';
    figName(figName==')') = '';
    figName(figName=='_') = '-';
    hfig = figure();
    figSetting();
    xlabel('Machine Number')
    ylabel(figName)
    title(['Variable: ' figName])
    set(gca,'XLim',[0 length(I)+1],'XTick',1:1:length(I),'XTickLabel',I);
    bar(1:1:length(I),x(I,ii),'BarWidth',1,'FaceColor',[0 0 1],'EdgeColor',[0 0 0.5]);
    figName = [pathname_res 'BarChart - ' name_case '_sort' int2str(pivot_cost) '_' figName '.fig'];
    if isoctave()
        hgsave(hfig,figName);
    else
        saveas(hfig,figName);
    end
end



% [uniqueRQ,~,~] = unique(geo0.RQnames,'stable');
% for ii=1:length(uniqueRQ)
%     tmp =  eval(upper(uniqueRQ{ii}));
%     figure()
%     figSetting(12,length(tmp(1,:))*5)
%     for kk=1:length(tmp(1,:))
%         subplot(length(tmp(1,:)),1,kk);
%         set(gca,'xTickLabel',I,'xTick',1:1:length(I));
%         bar(1:length(I),tmp(I,kk),0.6,'b')
%         xlim([0 n_mot+1]);
%         xlabel('Machine Number');
%         if length(tmp(1,:))>1
%             ylabel([uniqueRQ{ii} num2str(kk)]);
%         else
%             ylabel(uniqueRQ{ii});
%         end
%     end
%     h=gcf();
%     if isoctave()
%         fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_',uniqueRQ{ii});
%         hgsave(h,[fig_name]);
%     else
%         saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' uniqueRQ{ii} '.fig']);
%     end
%     %disp([uniqueRQ{ii},' = ',num2str(round(eval(upper(uniqueRQ{ii})),2))])
% end

