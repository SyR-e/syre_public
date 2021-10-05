% Copyright 2021
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

function eval_operatingPointAM(dataIn) %dataIn from GUI

%% Setup
pathname = dataIn.currentpathname;
filename = dataIn.currentfilename;
load([dataIn.currentpathname dataIn.currentfilename]);

RatedCurrent     = dataIn.RatedCurrent;
CurrLoPP         = dataIn.CurrLoPP;
SimulatedCurrent = dataIn.SimulatedCurrent;
GammaPP          = dataIn.GammaPP;
BrPP             = dataIn.BrPP;
TempPP           = dataIn.tempPP;
NumOfRotPosPP    = dataIn.NumOfRotPosPP;
AngularSpanPP    = dataIn.AngularSpanPP;
EvalType         = dataIn.EvalType;

per.EvalSpeed    = dataIn.EvalSpeed;

%Manual input Simulation Data
% eval_type           = 'singt'; %-singt  -coreloss
% per.gamma           = 55; %idq angle  %dataSet.GammaPP;
% per.CurrLoPP        = 1; %Current [pu] %dataSet.CurrLoPP;
% per.EvalSpeed       = 1800; %Evaluation Speed [rpm]
% per.nsim_singt      = 10;%Simulation Points %dataSet.NumOfRotPosPP
% per.delta_sim_singt = 60;%Angular Excursion %dataSet.AngularSpanPP
% per.tempPP          = 80; %permanent magnets temperature
% per.corelossflag    = strcmp(eval_type,'singtIron');
% per.save_fields     = 0; %activate and save fields every # step (0 deactiveted) must be < nsim_singt

eval_type           = EvalType;       %-singt  -coreloss
%per.gamma           = 55;             %idq angle  %dataSet.GammaPP;
%per.CurrLoPP       = CurrLoPP; %Current [pu] %dataSet.CurrLoPP;
per.i0              = RatedCurrent;
per.overload        = CurrLoPP;
%per.EvalSpeed      = 1800; %Evaluation Speed [rpm]
per.nsim_singt      = NumOfRotPosPP;   %Simulation Points %dataSet.NumOfRotPosPP
per.delta_sim_singt = AngularSpanPP;   %Angular Excursion %dataSet.AngularSpanPP
per.tempPP          = TempPP;          %permanent magnets temperature
per.BrPP            = BrPP;
per.corelossflag    = strcmp(eval_type,'singtIron');
per.save_fields     = 0;               %activate and save fields every # step (0 deactiveted) must be < nsim_singt


overload_temp = CurrLoPP;   % current to be simulated
gamma_temp    = GammaPP;    % current phase angle

%Python ansys path
ipypath  = 'C:\Program Files\AnsysEM\AnsysEM20.1\Win64\common\IronPython\';
ipy64exe = 'ipy64.exe';
ipypath  = strcat('"',ipypath,ipy64exe,'" "');

%% Simulation

performance = cell(1,length(overload_temp));
output = cell(1,length(overload_temp));
geometry = cell(1,length(overload_temp));
tempDirName = cell(1,length(overload_temp));
for ii = 1:length(overload_temp)
    performance{ii}          = per;
    performance{ii}.overload = overload_temp(ii);
    performance{ii}.gamma    = gamma_temp(ii);
end

geo.RemoveTMPfile = 'OFF';
for ii = 1:length(overload_temp)
    [geometry{ii},mat,output{ii},tempDirName{ii}] = AMfitness(geo,performance{ii},mat,eval_type,pathname,filename,ipypath);
end

% [~,~,out,~] = AMfitness(geo,per,mat,eval_type,pathname,filename,ipypath);

%% Save output into individual folders
for ii = 1:length(SimulatedCurrent)
    %output data
    %Theta=output{ii}.th;
    T    = output{ii}.T;
    fd   = output{ii}.fd;
    fq   = output{ii}.fq;
    id   = output{ii}.id;
    iq   = output{ii}.iq;
    IPF  = output{ii}.IPF;
    xdeg = output{ii}.xdeg;
    
    if per.corelossflag==1
        CoreLoss          = output{ii}.CoreLoss;
        CoreLossAvg       = output{ii}.CoreLossAvg;
        HysteresisLossAvg = output{ii}.HysteresisLossAvg;
        EddyLossAvg       = output{ii}.EddyLossAvg;
        ExcessLossAvg     = output{ii}.ExcessLossAvg;
        per.savefield     = 0;
    end
    
    %Repeat the 60° results for an entire round of 360°
    if xdeg < 360
        nrep = 360/xdeg;
    else
        nrep = 1;
    end
    
    T     = repmat(T,nrep,1);
    fd    = repmat(fd',nrep,1);
    fq    = repmat(fq',nrep,1);
    id    = repmat(id',nrep,1);
    iq    = repmat(iq',nrep,1);
    IPF   = repmat(IPF',nrep,1);
    Theta = linspace(0,nrep*xdeg,length(T));
    
    
    % plots
    FontSize = 12;
    FontName = 'TimesNewRoman';
    
    %save
    iAmp = round(CurrLoPP*RatedCurrent,2);
    resFolder = [pathname,filename(1:end-4),'_results\FEA results\Ansys_', eval_type , '_', num2str(round(iAmp,2)), 'A\'];
    save([resFolder 'Data.mat'], 'Theta', 'T', 'fd', 'fq', 'IPF');
    
    %Torque
    hT = figure
    figSetting()
    set(hT,'FileName',[resFolder 'Torque_PF.fig'])
    subplot(2,1,1)
    pl = plot(Theta,abs(T)); 
    grid on
    title(['Mean Torque = ' num2str(mean(T))]);
    %set(pl,'LineWidth',[2]);
    xlim([0 Theta(end)]), %ylim([ymin ymax]),
    %set(gca,'FontName',FontName,'FontSize',FontSize);
    ti = 0:60:360; set(gca,'XTick',ti);
    xl = xlabel('$\theta$ - degrees'); %set(xl,'Rotation',[0],'Fontsize',FontSize);
    yl = ylabel('T [Nm]');
    %set(yl,'Rotation',[90],'FontName',FontName,'Fontsize',FontSize);
    
    %IPF
    subplot(2,1,2)
    pl = plot(Theta,IPF);
    title(['Mean IPF = ' num2str(mean(IPF))]);
    grid on
    %set(pl,'LineWidth',2);
    xl = xlabel('$\theta$ - [degrees]'); set(xl,'Rotation',0,'FontName',FontName,'Fontsize',FontSize);
    yl=ylabel('IPF');
    %set(yl,'Rotation',90,'FontName',FontName,'Fontsize',FontSize);
    xlim([0 Theta(end)]);
    %set(gca,'FontName',FontName,'FontSize',FontSize);
    ti = 0:60:xdeg*nrep;  
    set(gca,'XTick',ti);
    
    savePrintFigure(hT)
    
    %Flux
    hdq = figure();
    figSetting()
    set(hdq,'FileName',[resFolder 'Flux_dq.fig'])
    subplot(2,1,1);
    hd1 = plot(Theta,fd);
    title(['Mean $\lambda_d$ = ' num2str(mean(fd))]);
    grid on
    xl_hd = xlabel('$\theta$ [Electrical degrees]');
    %set(xl_hd,'Rotation',0,'FontName',FontName,'Fontsize',FontSize); %,'FontWeight','Bold');
    yl_hd = ylabel('$\lambda_d$ [Wb]');
    %set(yl_hd,'Rotation',90,'FontName',FontName,'Fontsize',FontSize); %,'FontWeight','Bold');
    xlim([0 Theta(end)]);
    %set(gca,'FontSize',FontSize), %,'FontWeight','Bold');
    ti = 0:60:xdeg*nrep;  
    set(gca,'XTick',ti);
    
    hq = subplot(2,1,2);
    hq1 = plot(Theta,fq);
    title(['Mean $\lambda_q$ = ' num2str(mean(fq))]);
    grid on
    %set(hq1,'LineWidth',2);
    xl_hq=xlabel('$\theta$ [Electrical degrees]');
    %set(xl_hq,'Rotation',0,'FontName',FontName,'Fontsize',FontSize); %'FontWeight','Bold');
    yl_hq=ylabel('$\lambda_q$ [Wb]');
    %set(yl_hq,'Rotation',90,'FontName',FontName,'Fontsize',FontSize); %,'FontWeight','Bold');
    xlim([0 Theta(end)]);
    %set(gca,'FontSize',FontSize); %,'FontWeight','Bold');
    ti = 0:60:xdeg*nrep;  
    set(gca,'XTick',ti);
    
    savePrintFigure(hdq)
    
    %CoreLoss
    if per.corelossflag==1
        hl = figure;
        figSetting()
        set(hl,'FileName',[resFolder 'PlossVsTime.fig'])
        pl = plot(Theta,CoreLoss); grid on
        title(['Mean Core Loss = ' num2str(mean(CoreLoss(ceil(6*end/7)+1:end)))]);
        set(pl,'LineWidth',[2]);
        xlim([0 Theta(end)]), %ylim([ymin ymax]),
        set(gca,'FontName',FontName,'FontSize',FontSize);
        ti = 0:60:xdeg*nrep; set(gca,'XTick',ti);
        xl = xlabel('$\theta$ - [degrees]'); set(xl,'Rotation',[0],'Fontsize',FontSize);
        yl = ylabel('Core Loss [W]');
        set(yl,'Rotation',[90],'FontName',FontName,'Fontsize',FontSize);
        savePrintFigure(hl)
        
        Pfe_tot = CoreLossAvg(1); 
        Pfes_tot = CoreLossAvg(2); 
        Pfer_tot = CoreLossAvg(3);
        Pfe_h = HysteresisLossAvg(1); 
        Pfes_h = HysteresisLossAvg(2); 
        Pfer_h = HysteresisLossAvg(3);
        Pfe_c = EddyLossAvg(1);
        Pfes_c = EddyLossAvg(2);
        Pfer_c = EddyLossAvg(3);
        Pfe_ex = ExcessLossAvg(1); 
        Pfes_ex = ExcessLossAvg(2); 
        Pfer_ex = ExcessLossAvg(3);
        
        c = categorical({'Stat Hys','Stat Eddy','Stat Exc','Rot Hys','Rot Eddy','Rot Exc','Total Loss'},'Ordinal',true);
        Pfe = [Pfes_h,Pfes_c,Pfes_ex,Pfer_h,Pfer_c,Pfer_ex,Pfe_tot];
        Pfe(isnan(Pfe)) = 0;
 
        
        % Plot
        br=figure();
        figSetting(10,10,8)
        set(br,'FileName',[resFolder 'Ploss.fig'])
        set(gca,'YLim',[0 Pfe_tot*1.1])
        b = bar(c,Pfe,0.8);
        xtips1 = (b(1).XEndPoints);
        ytips1 = (b(1).YEndPoints);
        labels1 = string(round(b(1).YData));
        text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
        %legend ('Ansys')
        title('Iron Loss')
        ylabel('P [W]')
        save([resFolder 'Data.mat'], 'Pfe');
        savePrintFigure(br)
    end
end