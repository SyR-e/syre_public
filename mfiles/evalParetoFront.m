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

% filter duplicate solutions
% if isoctave() %OCT
%     [PSet,iA] = unique_oct(OUT.PSet,'rows');
% else
%     [PSet,iA] = unique(OUT.PSet,'rows');
% end

[PSet,iA] = unique(OUT.PSet,'rows');

x = PSet;

% RQ variables, set of optimization inputs
k = 0; y = 0; t = 0; u = 0;
for j = 1:length(geo0.RQnames)
    nameTemp{j} = upper(eval(['geo0.RQnames{j}']));
    if strcmp(nameTemp{j},'DALPHA')
        k = k +1;
        nameTemp{j} = [nameTemp{j} num2str(k)];
    end
    if strcmp(nameTemp{j},'HC')
        y = y +1;
        nameTemp{j} = [nameTemp{j} num2str(y)];
    end
    if strcmp(nameTemp{j},'DX')
        t = t +1;
        nameTemp{j} = [nameTemp{j} num2str(t)];
    end
    if strcmp(nameTemp{j},'BR')
        u = u +1;
        nameTemp{j} = [nameTemp{j} num2str(u)];
    end
    if strfind(nameTemp{j},'PMDIM')
        u = u +1;
        nameTemp{j} = [nameTemp{j}(1:5) int2str(u)];
    end
    if strcmp(nameTemp{j},'BETAPMSHAPE')
        u = u +1;
        nameTemp{j} = [nameTemp{j} num2str(u)];
    end
    eval([nameTemp{j} ' = zeros(size(x,1),1)']);
end

% goals
T   = zeros(size(x,1),1);
dT  = zeros(size(x,1),1);
MCu = zeros(size(x,1),1);
MPM = zeros(size(x,1),1);
% PFES = T; PFER = T; FD90 = T;

% geo.nsim_singt = 2; % debug
Tini=now();
parfor m = 1:size(x,1)
    
    Tstart=now;
    mot_index = num2str(m);
    disp(['Evaluating machine ' mot_index ' of ' num2str(size(x,1))])
    
    % FEA evaluates each machine of the front
    
    [cost{m},geometry{m},material{m},output{m},~] = FEMMfitness(x(m,:)',geo,per_temp,mat,eval_type,strrep(filename,'.mat','.fem'));
    
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
    [geo,~,mat] = interpretRQ(x(m,:),geo,mat);
    [geo,mat] = draw_motor_in_FEMM(geo,mat,pathname_res,['mot_'  mot_index '.fem']);
%     mi_saveas([pathname '\mot_'  mot_index '.fem']);
    
    % Save png file mot
    mi_zoomnatural;
%     mi_savebitmap([pathname_res '\mot_'  mot_index '.bmp']);
    mi_savemetafile([pathname_res '\mot_'  mot_index '.emf'])
    closefemm
    
    if isoctave() %OCT
%         (to be debugged
%  - mot0.fem does not exist any longer
%  - pathname changed
        current_dir=strcat(syreRoot, '\mot0.fem');
        movefile (current_dir, pathname);
        old_name= strcat(pathname, '\mot0.fem');
        new_name= strcat(pathname, '\mot_', mot_index, '.fem');
        rename(old_name,new_name);
        clear old_name new_name current_dir
    end
    
    % Save data geometry mot
    geo.RQ = x(m,:);
    
    if isoctave()          %OCT
        dataSet.AirGapThickness = geo.g; % airgap thickness
        dataSet.AirGapThickness = round(dataSet.AirGapThickness .*100) ./100;
        dataSet.AirGapRadius = geo.r; % machine airgap radius
        dataSet.AirGapRadius = round(dataSet.AirGapRadius .*100) ./100;
        dataSet.ToothLength = geo.lt; % tooth length
        dataSet.ToothLength = round(dataSet.ToothLength .*100) ./100;
        dataSet.StatorSlotOpen = geo.acs; % stator slot open in [p.u.]
        dataSet.StatorSlotOpen = round(dataSet.StatorSlotOpen .*100) ./100;
        dataSet.ToothWidth = geo.wt; % Bgap/Btooth (determines tooth width^-1, yoke width^-1)
        dataSet.ToothWidth = round(dataSet.ToothWidth .*100) ./100;
        dataSet.ToothTangDepth = geo.ttd; % tooth tang depth [mm]
        dataSet.ToothTangDepth = round(dataSet.ToothTangDepth .*100) ./100;
        dataSet.Br = round(mat.LayerMag.Br .*10000) ./10000; % Br
        dataSet.Br = round(dataSet.Br .*100) ./100;
        dataSet.ALPHApu = geo.dalpha_pu;
        dataSet.ALPHApu = round(dataSet.ALPHApu .*100) ./100;
        dataSet.PMdim   = round(geo.PMdim.*100)./100;
        dataSet.betaPMshape = round(geo.betaPMshape.*100)./100;
    else
        dataSet.AirGapThickness = geo.g; % airgap thickness
        dataSet.AirGapThickness = round(dataSet.AirGapThickness,2);
        dataSet.AirGapRadius = geo.r; % machine airgap radius
        dataSet.AirGapRadius = round(dataSet.AirGapRadius,2);
        dataSet.ToothLength = geo.lt; % tooth length
        dataSet.ToothLength = round(dataSet.ToothLength,2);
        dataSet.StatorSlotOpen = geo.acs; % stator slot open in [p.u.]
        dataSet.StatorSlotOpen = round(dataSet.StatorSlotOpen,2);
        dataSet.ToothWidth = geo.wt; % Bgap/Btooth (determines tooth width^-1, yoke width^-1)
        dataSet.ToothWidth = round(dataSet.ToothWidth,2);
        dataSet.ToothTangDepth = geo.ttd; % tooth tang depth [mm]
        dataSet.ToothTangDepth = round(dataSet.ToothTangDepth,2);
        dataSet.Br = round(mat.LayerMag.Br,4); % Br
        dataSet.Br = round(dataSet.Br,2);
        dataSet.ALPHApu = geo.dalpha_pu;
        dataSet.ALPHApu = round(dataSet.ALPHApu,2);
        dataSet.PMdim   = round(geo.PMdim,2);
        dataSet.betaPMshape = round(geo.betaPMshape,2);
    end
    if strcmp(geo.RotType,'SPM')
        dataSet.ThicknessOfPM = geo.hc_pu;
    else
        dataSet.HCpu = geo.hc_pu;
    end
    if isoctave()                                           %OCT
        dataSet.HCpu = round(dataSet.HCpu .*100) ./100;
        dataSet.DepthOfBarrier = geo.dx;    % the depth of the barriers radial-wise in per unit
        dataSet.DepthOfBarrier = round(dataSet.DepthOfBarrier .*100) ./100;
    else
        dataSet.HCpu = round(dataSet.HCpu,2);
        dataSet.DepthOfBarrier = geo.dx;    % the depth of the barriers radial-wise in per unit
        dataSet.DepthOfBarrier = round(dataSet.DepthOfBarrier,2);
    end
    dataSet.RQ = geo.RQ;
    
    if isoctave()  %OCT
        name_file = strcat(pathname, '\mot_', mot_index, '.mat');
        save ('-mat7-binary', name_file,'geo','cost','per','dataSet','mat');
        clear name_file
    else
        save([pathname_res '\mot_' mot_index '.mat'],'geo','cost','per','dataSet','mat','out');
    end
    
    %% %%%%%%%
    COST=[COST; cost{m}];
    if m<10
        mot_index = ['0' mot_index];
    end
    
    %     geo.Br
    % geometry summary of the re evaluated Pareto front
    % replace fields od geo0 with content of RQ
    %     k = 0; y = 0; t = 0; u = 0;
    for j = 1:length(nameTemp)
        eval([nameTemp{j} '(m) = x(m,j);']);
    end
    
    % Ho inserito anche l'obiettivo MassPM (sono variabili che rimangono all'interno della funzione, non hanno output...
    % ma la aggiungo lo stesso -  rev.Gallo 12/04/2018
    % cost functions (modify in the future)
    temp1 = 1;
    if strcmp(geo.OBJnames{temp1},'Torque')
        T(m,1) = COST(m,temp1); %out.Tn;
        %         cost(temp1) = -out.T;
        temp1 = temp1+1;
    end
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'TorRip')
        %         cost(temp1) = out.dTpp;
        dT(m,1) = COST(m,temp1); %out.ripple_pu;
        temp1 = temp1+1;
    end
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'MassCu')
        %         cost(temp1)=calcMassCu(geo);
        MCu(m,1) = COST(m,temp1); %out.MassCu;
        temp1=temp1+1;
    end
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'MassPM')
        MPM(m,1) = COST(m,temp1); %out.MassCu;
        temp1=temp1+1;
    end
end

%[STATUS,MESSAGE,MESSAGEID] = copyfile([pathname_ini,filename],[pathname,'\']);
% mat_name=strcat(pathname_ini, filename); %OCT
% mat_dest=strcat(pathname,'\');
[STATUS,MESSAGE,MESSAGEID] = copyfile(strcat(pathname_ini, filename),strcat(pathname_res, filename));
% clear mat_name mat_dest
out_all = struct('T', T, 'dT', dT, 'MCu', MCu, 'MPM', MPM); %, 'PFES', PFES, 'PFER', PFER, 'FD90', FD90);

% sort the solutions of the Front according to T or dT
n_mot = size(COST,1);
[Y,I] = sort(COST(:,pivot_cost));

COST_sorted = [COST(I,:) I];
x_sorted = x(I,:);

T_s = out_all.T(I,:);
DT_s = out_all.dT(I,:);
MCu_s = out_all.MCu(I,:);
MPM_s = out_all.MPM(I,:);

out_all_sorted = struct('T', T_s, 'dT', DT_s, 'MCu', MCu_s);

name_output = strcat(filename,'_sort_');
name_output = strrep(name_output,'.mat','');
name_case = strrep(filename,'.mat','');

% Pareto front
if OUT.Param.NOBJ==1
    xlabel('generation')
    ylabel(['log(' geo0.OBJnames{1} ')'])
    h=gcf();
    if isoctave()   %OCT
        fig_name=strcat(pathname_res, '\goalVSgenerations - ', name_case);
        hgsave(h,[fig_name]);
    else
        saveas(gcf,[pathname_res '\goalVSgenerations - ' name_case '.fig'])
    end
    figure()
    hold on
    grid on
    box on
    xlabel('Machine number')
    ylabel(geo0.OBJnames)
    set(gca,'XLim',[0 length(COST)+1],'XTick',1:1:length(COST))
    bar(COST)
    h=gcf();
    if isoctave() %OCT
        fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_goal ');
        hgsave(h,[fig_name]);
    else
        saveas(gcf,[pathname_res 'BarChart - ' name_case '_goal.fig'])
    end
elseif OUT.Param.NOBJ==2
    close
    figure(1), clf, hold on
    for ii=1:n_mot
        plot(COST(ii,1),COST(ii,2),'x'),
        text(COST(ii,1)+0.1,COST(ii,2)+0.1,num2str(ii));
    end
    grid on, hold off
    xlabel(geo0.OBJnames{1}), ylabel(geo0.OBJnames{2})
    
    [front,idx] = FastParetoEstimation(x,COST);
    nmot_actual_front=find(idx); %macchine sul fronte vero
    
    figure(1), hold on
    for ii=1:length(nmot_actual_front) %ii = 1:length(x) PERCHE' USARE LENGTH(X) ABBIAMO LA LUNGHEZZA DI T ???????????
        C1_front(ii)=front(ii,end-1);C2_front(ii)=front(ii,end);
        plot(C1_front(ii),C2_front(ii),'ro','LineWidth',2),
        text(C1_front(ii)+0.1,C2_front(ii)+0.1,num2str(nmot_actual_front(ii)));
    end
    [T_front_sorted,ii_tfs]=sort(C1_front);
    plot(C1_front(ii_tfs),C2_front(ii_tfs),'r','LineWidth',2)
    set(gca,'FontName','Arial','FontSize',12)
    grid on, hold off
    xlabel(geo0.OBJnames{1}),
    ylabel(geo0.OBJnames{2})
    figure_title = [geo.RotType ' nlay = ' num2str(geo.nlay) ' p = ' num2str(geo.p) ' - magnete plasto Br = ' num2str(mat.LayerMag.Br) ' T'];
    title(figure_title)
    h=gcf();
    if isoctave()
        fig_name=strcat(pathname_res, 'Pareto - ', name_case);
        hgsave(h,[fig_name]);
    else
        saveas(gcf,[pathname_res 'Pareto - ' name_case '.fig']);
    end
    
    % Bar charts
    % 1) build a bar chart for each element of nameTemp
    % 2) build a bar chart for each goal
    bw = 0.7;   % bar width
    
    % RQ elements: they are not yet sorted
    k = 0; y = 0; t = 0; u = 0;
    for j = 1:length(nameTemp)
        if strfind(nameTemp{j},'DALPHA')
            k = k +1;
            figure(2)
            set(gca,'xTickLabel',I);subplot(geo.nlay,1,k);
            bar(eval(nameTemp{j}),bw,'r'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Machanical Deg');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_DALPHA');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_DALPHA.fig']);
            end
        elseif strfind(nameTemp{j},'HC')
            y = y +1;
            figure(3)
            set(gca,'xTickLabel',I);subplot(geo.nlay,1,y);
            bar(eval(nameTemp{j}),bw,'g'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number');
            if strcmp(geo.RotType,'SPM')
                ylabel('PM thickness [mm]');
            else
                ylabel('Barrier Width [mm]');
            end
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_HC');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_HC.fig']);
            end
        elseif strfind(nameTemp{j},'DX')
            t = t +1;
            figure(4)
            set(gca,'xTickLabel',I);subplot(geo.nlay,1,t);
            bar(eval(nameTemp{j}),bw,'y'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Barrier Depth [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_DX');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_DX.fig']);
            end
        elseif strfind(nameTemp{j},'BR')
            u = u +1;
            figure(5);
            set(gca,'xTickLabel',I); subplot(geo.nlay,1,u);
            bar(eval(nameTemp{j}),bw,'w'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Br [T]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_BR');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_BR.fig']);
            end
        elseif strcmp(nameTemp{j},'G')
            figure;
            set(gca,'xTickLabel',I);
            bar(eval(nameTemp{j}),bw,'c'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Airgap Thickness [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        elseif strcmp(nameTemp{j},'R')
            figure;
            set(gca,'xTickLabel',I);
            bar(eval(nameTemp{j}),bw,'m'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Airgap Radius [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        elseif strcmp(nameTemp{j},'WT')
            figure;
            set(gca,'xTickLabel',I);
            bar(eval(nameTemp{j}),bw,'k'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Tooth Width [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        elseif strcmp(nameTemp{j},'LT')
            figure;
            set(gca,'xTickLabel',I);subplot()
            bar(eval(nameTemp{j}),bw,'y'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Tooth Length [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        elseif strcmp(nameTemp{j},'ACS')
            figure;
            set(gca,'xTickLabel',I);subplot()
            bar(eval(nameTemp{j}),bw,'r'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Slot Open [p.u.]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        elseif strcmp(nameTemp{j},'TTD')
            figure;
            set(gca,'xTickLabel',I);subplot()
            bar(eval(nameTemp{j}),bw,'g'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Tooth Tang Depth [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        else
            figure;
            set(gca,'xTickLabel',I);
            bar(eval(nameTemp{j}),bw,'b'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Electrical Deg');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        end
    end
    
    % cost functions (already sorted)
    figure, subplot(2,1,1);
    bar(sign(per.objs(1,1))*COST(I,1),bw,'y'), grid,
    xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
    xlabel('Machine Number'), ylabel(geo0.OBJnames{1})
    legend(geo0.OBJnames{1});
    % title('Torque in descending order');
    subplot(2,1,2);
    bar(sign(per.objs(2,1))*COST(I,2),bw,'b'), grid,
    xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
    xlabel('Machine Number'), ylabel(geo0.OBJnames{2});
    legend(geo0.OBJnames{2});
    h=gcf();
    if isoctave()
        fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_goal1&goal2');
        hgsave(h,[fig_name]);
    else
        saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_goal1&goal2.fig']);
    end
elseif OUT.Param.NOBJ==3
    close
    figure(1), clf, hold on, view(3)
    for ii=1:n_mot
        plot3(COST(ii,1),COST(ii,2),COST(ii,3),'x');
        text(COST(ii,1)+0.1,COST(ii,2)+0.1,COST(ii,3)+0.1,int2str(ii));
    end
    grid on
    xlabel(geo0.OBJnames{1});
    ylabel(geo0.OBJnames{2});
    zlabel(geo0.OBJnames{3});
    
    [front,idx] = FastParetoEstimation(x,COST);
    nmot_actual_front=find(idx); %macchine sul fronte vero
    
    figure(1)
    for ii=1:length(nmot_actual_front) %ii = 1:length(x) PERCHE' USARE LENGTH(X) ABBIAMO LA LUNGHEZZA DI T ???????????
        plot3(front(ii,end-2),front(ii,end-1),front(ii,end),'ro','LineWidth',2),
        text(front(ii,end-2)+0.1,front(ii,end-1)+0.1,front(ii,end)+0.1,int2str(nmot_actual_front(ii)));
    end
    figure_title = [geo.RotType ' nlay = ' num2str(geo.nlay) ' p = ' num2str(geo.p) ' - magnete plasto Br = ' num2str(mat.LayerMag.Br) ' T'];
    title(figure_title)
    h=gcf();
    if isoctave()
        fig_name=strcat(pathname_res, '\Pareto - ', name_case);
        hgsave(h,[fig_name]);
    else
        saveas(gcf,[pathname_res '\Pareto - ' name_case '.fig']);
    end
    
    % Bar charts
    % 1) build a bar chart for each element of nameTemp
    % 2) build a bar chart for each goal
    bw = 0.7;   % bar width
    
    % RQ elements: they are not yet sorted
    k = 0; y = 0; t = 0; u = 0;
    for j = 1:length(nameTemp)
        if strfind(nameTemp{j},'DALPHA')
            k = k +1;
            figure(2)
            set(gca,'xTickLabel',I);subplot(geo.nlay,1,k);
            bar(eval(nameTemp{j}),bw,'r'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Machanical Deg');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_DALPHA');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_DALPHA.fig']);
            end
        elseif strfind(nameTemp{j},'HC')
            y = y +1;
            figure(3)
            set(gca,'xTickLabel',I);subplot(geo.nlay,1,y);
            bar(eval(nameTemp{j}),bw,'g'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number');
            if strcmp(geo.RotType,'SPM')
                ylabel('PM thickness [mm]');
            else
                ylabel('Barrier Width [mm]');
            end
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_HC');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_HC.fig']);
            end
        elseif strfind(nameTemp{j},'DX')
            t = t +1;
            figure(4)
            set(gca,'xTickLabel',I);subplot(geo.nlay,1,t);
            bar(eval(nameTemp{j}),bw,'y'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Barrier Depth [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_DX');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_DX.fig']);
            end
        elseif strfind(nameTemp{j},'BR')
            u = u +1;
            figure(5);
            set(gca,'xTickLabel',I); subplot(geo.nlay,1,u);
            bar(eval(nameTemp{j}),bw,'w'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Br [T]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_BR');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_BR.fig']);
            end
        elseif strcmp(nameTemp{j},'G')
            figure;
            set(gca,'xTickLabel',I);
            bar(eval(nameTemp{j}),bw,'c'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Airgap Thickness [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        elseif strcmp(nameTemp{j},'R')
            figure;
            set(gca,'xTickLabel',I);
            bar(eval(nameTemp{j}),bw,'m'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Airgap Radius [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        elseif strcmp(nameTemp{j},'WT')
            figure;
            set(gca,'xTickLabel',I);
            bar(eval(nameTemp{j}),bw,'k'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Tooth Width [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        elseif strcmp(nameTemp{j},'LT')
            figure;
            set(gca,'xTickLabel',I);subplot()
            bar(eval(nameTemp{j}),bw,'y'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Tooth Length [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        elseif strcmp(nameTemp{j},'ACS')
            figure;
            set(gca,'xTickLabel',I);subplot()
            bar(eval(nameTemp{j}),bw,'r'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Slot Open [p.u.]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        elseif strcmp(nameTemp{j},'TTD')
            figure;
            set(gca,'xTickLabel',I);subplot()
            bar(eval(nameTemp{j}),bw,'g'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Tooth Tang Depth [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        else
            figure;
            set(gca,'xTickLabel',I);
            bar(eval(nameTemp{j}),bw,'b'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('PU');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        end
    end
    
    % cost functions (already sorted)
    figure, subplot(3,1,1);
    bar(sign(per.objs(1,1))*COST(I,1),bw,'b'), grid,
    xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
    xlabel('Machine Number'), ylabel(geo0.OBJnames{1})
    legend(geo0.OBJnames{1});
    % title('Torque in descending order');
    subplot(3,1,2);
    bar(sign(per.objs(2,1))*COST(I,2),bw,'b'), grid,
    xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
    xlabel('Machine Number'), ylabel(geo0.OBJnames{2});
    legend(geo0.OBJnames{2});
    subplot(3,1,3);
    bar(sign(per.objs(3,1))*COST(I,3),bw,'b'), grid,
    xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
    xlabel('Machine Number'), ylabel(geo0.OBJnames{3});
    legend(geo0.OBJnames{3});
    h=gcf();
    if isoctave()
        fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_goal1&goal2');
        hgsave(h,[fig_name]);
    else
        saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_goal1&goal2.fig']);
    end
else
    % Bar charts
    % 1) build a bar chart for each element of nameTemp
    % 2) build a bar chart for each goal
    bw = 0.7;   % bar width
    
    % RQ elements: they are not yet sorted
    k = 0; y = 0; t = 0; u = 0;
    for j = 1:length(nameTemp)
        if strfind(nameTemp{j},'DALPHA')
            k = k +1;
            figure(2)
            set(gca,'xTickLabel',I);subplot(geo.nlay,1,k);
            bar(eval(nameTemp{j}),bw,'r'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Machanical Deg');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_DALPHA');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_DALPHA.fig']);
            end
        elseif strfind(nameTemp{j},'HC')
            y = y +1;
            figure(3)
            set(gca,'xTickLabel',I);subplot(geo.nlay,1,y);
            bar(eval(nameTemp{j}),bw,'g'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number');
            if strcmp(geo.RotType,'SPM')
                ylabel('PM thickness [mm]');
            else
                ylabel('Barrier Width [mm]');
            end
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_HC');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_HC.fig']);
            end
        elseif strfind(nameTemp{j},'DX')
            t = t +1;
            figure(4)
            set(gca,'xTickLabel',I);subplot(geo.nlay,1,t);
            bar(eval(nameTemp{j}),bw,'y'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Barrier Depth [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_DX');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_DX.fig']);
            end
        elseif strfind(nameTemp{j},'BR')
            u = u +1;
            figure(5);
            set(gca,'xTickLabel',I); subplot(geo.nlay,1,u);
            bar(eval(nameTemp{j}),bw,'w'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Br [T]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_BR');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_BR.fig']);
            end
        elseif strcmp(nameTemp{j},'G')
            figure;
            set(gca,'xTickLabel',I);
            bar(eval(nameTemp{j}),bw,'c'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Airgap Thickness [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        elseif strcmp(nameTemp{j},'R')
            figure;
            set(gca,'xTickLabel',I);
            bar(eval(nameTemp{j}),bw,'m'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Airgap Radius [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        elseif strcmp(nameTemp{j},'WT')
            figure;
            set(gca,'xTickLabel',I);
            bar(eval(nameTemp{j}),bw,'k'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Tooth Width [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        elseif strcmp(nameTemp{j},'LT')
            figure;
            set(gca,'xTickLabel',I);subplot()
            bar(eval(nameTemp{j}),bw,'y'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Tooth Length [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        elseif strcmp(nameTemp{j},'ACS')
            figure;
            set(gca,'xTickLabel',I);subplot()
            bar(eval(nameTemp{j}),bw,'r'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Slot Open [p.u.]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        elseif strcmp(nameTemp{j},'TTD')
            figure;
            set(gca,'xTickLabel',I);subplot()
            bar(eval(nameTemp{j}),bw,'g'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Tooth Tang Depth [mm]');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        else
            figure;
            set(gca,'xTickLabel',I);
            bar(eval(nameTemp{j}),bw,'b'), grid,
            xlim([0 n_mot+1]);set(gca,'xTickLabel',I),set(gca,'xTick',1:1:length(I)),
            legend(nameTemp{j});
            xlabel('Machine Number'); ylabel('Electrical Deg');
            h=gcf();
            if isoctave()
                fig_name=strcat(pathname_res, 'BarChart - ', name_case, '_sort', num2str(pivot_cost), '_', nameTemp{j});
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[pathname_res 'BarChart - ' name_case '_sort' num2str(pivot_cost) '_' nameTemp{j} '.fig']);
            end
        end
    end
    
    % cost functions (already sorted)
    figure()
    for ii=1:OUT.Param.NOBJ
        subplot(OUT.Param.NOBJ,1,ii)
        grid on
        bar(sign(per.objs(ii,1))*COST(I,ii),bw,'b')
        set(gca,'XLim',[0 n_mot+1],'xTick',1:1:length(I),'xTickLabel',I);
        xlabel('Machine Number');
        ylabel(geo0.OBJnames{ii})
        legend(geo0.OBJnames{ii});
    end
end

%% Evaluation of the Progressive distribution of the Pareto Front

% if OUT.eval_type=='MO_GA'
%     Pop=OUT.PSet;
%     Fit=OUT.PFront;
% else
%     Pop=OUT.MatrixPop;
%     Fit=OUT.MatrixFitness(:,:,end);
% end
% legenda={};
% color={'k' 'r' 'g' 'c','m'};
% figure(100);
% ii=1;
% for jk=ceil(linspace(1,size(Fit,3),5))
%     [front,idx] = FastParetoEstimation(Pop(:,:,jk),Fit(:,:,jk));
%     nmot_actual_front=find(idx); %macchine sul fronte vero
%     hold on;
%     plot(Fit(nmot_actual_front,1,jk),Fit(nmot_actual_front,2,jk),'*','Color',color{ii});
%     legenda{ii}=['front ',num2str(jk),' of ',num2str(size(Fit,3))];
%     ii=ii+1;
%
% end
% hold off;grid on; xlabel('Torque [Nm]'); ylabel('ripple [%]'); title('Evolution of the Pareto Front during the optimization process'); legend(legenda);
% saveas(gcf,[pathname '\pareto2x_evolution_during_optimization-' name_case '.fig']);

if OUT.eval_type=='MO_OA'
    if OUT.Param.NOBJ>1
        save([pathname_res '\' filename],'front','idx','-append');
    end
end