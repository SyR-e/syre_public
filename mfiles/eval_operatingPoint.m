% Copyright 2019
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

function [out,senseOut,newDir] = eval_operatingPoint(dataIn,flagSave)

% simulates single or multiple (id,iq) conditions
% example inputs:
% single condition: CurrLoPP = 1, GammaPP = 45
% multiple points:  CurrLoPP = [1 1.5 2], gamma = [45 45 45]

% Uses matlabpool (parfor)

% Key INPUTs: CurrLoPP: current to be simulated
%             GammaPP: current phase angle
%             BrPP: remanence of all barriers magnets
%             NumOfRotPosPP: # simulated positions
%             AngularSpanPP: angular span of simulation
%=========================================================================

pathname=dataIn.currentpathname;
filemot = strrep(dataIn.currentfilename,'.mat','.fem');
load([dataIn.currentpathname dataIn.currentfilename]);

RatedCurrent      = dataIn.RatedCurrent;
CurrLoPP          = dataIn.CurrLoPP;
GammaPP           = dataIn.GammaPP;
iOffsetPP         = dataIn.SimulatedCurrentOffset;

if(isnan(dataIn.Num3PhaseCircuit))
    if(dataIn.CustomCurrentEnable)
        CurrLoPP = CurrLoPP(1)*ones(1,length(dataIn.CustomCurrentRotorPos));
        GammaPP = GammaPP(1)*ones(1,length(dataIn.CustomCurrentRotorPos));
        iOffsetPP = iOffsetPP(1)*ones(1,length(dataIn.CustomCurrentRotorPos));
    end
end

%SimulatedCurrent = dataIn.SimulatedCurrent;
SimulatedCurrent  = RatedCurrent*CurrLoPP;
BrPP              = dataIn.BrPP;
NumOfRotPosPP     = dataIn.NumOfRotPosPP;
AngularSpanPP     = dataIn.AngularSpanPP;
per.flag3phaseSet = dataIn.Active3PhaseSets;
n3ph              = dataIn.Num3PhaseCircuit;
% iOffsetPP = 0.*ones(1,n3ph);

geo.XFEMMsimulation = dataIn.XFEMMsimulation;

per.EvalSpeed = dataIn.EvalSpeed;
EvalSpeed = per.EvalSpeed*ones(size(CurrLoPP));

clc;

if nargin()==1
    flagSave = 2;
end

% back compatibility EESM
if strcmp(geo.RotType,'EESM')
    geo.axisType = 'SR';
    dataIn.axisType = 'SR';
end

if ~isfield(geo,'axisType')
    if strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype')||strcmp(geo.RotType,'Spoke-type')
        geo.axisType = 'PM';
    else
        geo.axisType = 'SR';
    end
end

if ~strcmp(geo.axisType,dataIn.axisType)
    %geo.axisType = dataIn.axisType;
    if strcmp(dataIn.axisType,'PM')
        geo.th0 = geo.th0 - 90;
    else
        geo.th0 = geo.th0 + 90;
    end
end


eval_type = dataIn.EvalType;

per.overload=CurrLoPP;
per.i0 = RatedCurrent;
per.BrPP=BrPP;
per.iOffset = iOffsetPP;

per.nsim_singt = NumOfRotPosPP;       % # simulated positions
per.delta_sim_singt = AngularSpanPP;  % angular span of simulation
per.if = dataIn.FieldCurrent;         % Field Current EESM
FieldCurrent = dataIn.FieldCurrent;   % For multiple simultions

if ~strcmp(dataSet.TypeOfRotor,'EESM')
    FieldCurrent = zeros(size(CurrLoPP));
end


%custom current
if dataIn.CustomCurrentEnable
    if(isnan(n3ph))
        per.custom_ia         = dataIn.CustomCurrentA;
        per.custom_ib         = dataIn.CustomCurrentB;
        per.custom_ic         = dataIn.CustomCurrentC;
        per.custom_id         = dataIn.CustomCurrentD;
        per.custom_ie         = dataIn.CustomCurrentE;
        per.custom_rotorPos       = dataIn.CustomCurrentRotorPos;
        per.custom_act        = dataIn.CustomCurrentEnable;
    else
        per.custom_ia         = dataIn.CustomCurrentA;
        per.custom_ib         = dataIn.CustomCurrentB;
        per.custom_ic         = dataIn.CustomCurrentC;
        per.custom_time       = dataIn.CustomCurrentTime;
        per.custom_act        = dataIn.CustomCurrentEnable;
    end
else
    per.custom_act = 0;
end

offset = zeros(1,length(CurrLoPP));

% check parallel computing
ppState=parallelComputingCheck();

if ppState==0
    ppState1 = 1;
else
    ppState1 = ppState;
end

flag_singtpp = 0;   %no parallel pool 
% flag_singtpp = 1; % parallel pool 
if ~strcmp(eval_type,'singt')
    flag_singtpp = 0;
end
per.flag_singtpp = flag_singtpp;

% if (strcmp(eval_type,'singtIron') || flag_singtpp == 1)
if (strcmp(eval_type,'singtIron') && flag_singtpp == 1)
    CurrLoPP = CurrLoPP*ones(1,ppState1);
    GammaPP  = GammaPP*ones(1,ppState1);
    iOffsetPP = iOffsetPP*ones(1,ppState1);
    tmp      = linspace(0,AngularSpanPP,ppState1+1);
    offset              = tmp(1,1:end-1);
    if length(offset)>1
        per.delta_sim_singt = offset(2);
    end
    per.nsim_singt      = ceil(per.nsim_singt/ppState1);
    %     offset(1,2:end)     = offset(1,2:end)+offset(1,2)/per.nsim_singt;
end

% single point or array of points simulation
performance = cell(1,length(CurrLoPP));
output = cell(1,length(CurrLoPP));

geometry = cell(1,length(CurrLoPP));
tempDirName = cell(1,length(CurrLoPP));
for ii = 1:length(CurrLoPP)
    performance{ii} = per;
    performance{ii}.overload  = CurrLoPP(ii);
    performance{ii}.gamma     = GammaPP(ii);
    performance{ii}.offset    = offset(ii);
    performance{ii}.if        = FieldCurrent(ii);
    performance{ii}.EvalSpeed = EvalSpeed(ii);
    performance{ii}.iOffset   = iOffsetPP(ii);

    if(isnan(n3ph))
        if(per.custom_act)
            performance{ii}.rotorPos  = ii;
        end
    end
end

% if strcmp(eval_type,'singtIron') && length(CurrLoPP)>1
%     performance{1}.nsim_singt = performance{1}.nsim_singt+1;
% end

nSim = length(CurrLoPP);

if(isnan(n3ph))
    if nSim==1
        performance{1}          = per;
        performance{1}.overload = CurrLoPP;
        performance{1}.gamma    = GammaPP;
        performance{1}.offset   = 0;
        performance{1}.if       = FieldCurrent(1);
        performance{1}.iOffset  = iOffsetPP;
        performance{1}.rotorPos  = 1;
    end
else
    if nSim==n3ph
        nSim = 1;
        performance{1}          = per;
        performance{1}.overload = CurrLoPP;
        performance{1}.gamma    = GammaPP;
        performance{1}.offset   = 0;
        performance{1}.if       = FieldCurrent(1);
        performance{1}.iOffset  = iOffsetPP;
    end
end

geo.RemoveTMPfile = 'OFF';

if (ppState==0 && nSim>4)
    parpool();
    ppState=parallelComputingCheck();
end

fileMotWithPath=[pathname filemot];

geo0 = geo;
mat0 = mat;
% evaluation
if ppState<1
    for ii = 1:nSim
        geoTmp = geo0;
        perTmp = performance{ii};
        matTmp = mat0;
        [~,geometry{ii},~,output{ii},tempDirName{ii}] = FEMMfitness([],geoTmp,perTmp,matTmp,eval_type,fileMotWithPath);
    end
else
    parfor ii = 1:nSim %%%
        geoTmp = geo0;
        perTmp = performance{ii};
        matTmp = mat0;
        [~,geometry{ii},~,output{ii},tempDirName{ii}] = FEMMfitness([],geoTmp,perTmp,matTmp,eval_type,fileMotWithPath);
    end
end


if ((strcmp(eval_type,'singtIron')||strcmp(eval_type,'singmIron'))&&flag_singtpp == 1)

    out.SOL.th    = [];
    out.SOL.id    = [];
    out.SOL.iq    = [];
    out.SOL.fd    = [];
    out.SOL.fq    = [];
    out.SOL.T     = [];
    out.SOL.ia    = [];
    out.SOL.ib    = [];
    out.SOL.ic    = [];
    out.SOL.fa    = [];
    out.SOL.fb    = [];
    out.SOL.fc    = [];
    out.SOL.we    = [];
    out.SOL.wc    = [];
    out.SOL.bs    = [];
    out.SOL.br    = [];
    out.SOL.am    = [];
    out.SOL.pos   = [];
    out.SOL.vol   = [];
    out.SOL.groNo = [];
    per.delta_sim_singt = AngularSpanPP;

    for ii=1:length(output)
        out.SOL.th    = [out.SOL.th output{ii}.SOL.th];
        out.SOL.id    = [out.SOL.id output{ii}.SOL.id];
        out.SOL.iq    = [out.SOL.iq output{ii}.SOL.iq];
        out.SOL.fd    = [out.SOL.fd output{ii}.SOL.fd];
        out.SOL.fq    = [out.SOL.fq output{ii}.SOL.fq];
        out.SOL.T     = [out.SOL.T  output{ii}.SOL.T];
        out.SOL.ia    = [out.SOL.ia output{ii}.SOL.ia];
        out.SOL.ib    = [out.SOL.ib output{ii}.SOL.ib];
        out.SOL.ic    = [out.SOL.ic output{ii}.SOL.ic];
        out.SOL.fa    = [out.SOL.fa output{ii}.SOL.fa];
        out.SOL.fb    = [out.SOL.fb output{ii}.SOL.fb];
        out.SOL.fc    = [out.SOL.fc output{ii}.SOL.fc];
        out.SOL.we    = [out.SOL.we output{ii}.SOL.we];
        out.SOL.wc    = [out.SOL.wc output{ii}.SOL.wc];
        if flag_singtpp == 0
            out.SOL.bs    = [out.SOL.bs; output{ii}.SOL.bs];
            out.SOL.br    = [out.SOL.br; output{ii}.SOL.br];
            out.SOL.am    = [out.SOL.am; output{ii}.SOL.am];
        end
    end
    if flag_singtpp == 0
        out.SOL.pos   = [out.SOL.pos output{1}.SOL.pos];
        out.SOL.vol   = [out.SOL.vol output{1}.SOL.vol];
        out.SOL.groNo = [out.SOL.groNo output{1}.SOL.groNo];
    end
    if flag_singtpp == 0
        [SOL] = evalIronLossFEMM(geo,per,mat,out.SOL,2);
        if isfield(SOL,'psh')
            out.Pfes_h        = sum(sum(SOL.psh))*(2*geo.p/geo.ps);
            out.Pfes_c        = sum(sum(SOL.psc))*(2*geo.p/geo.ps);
            out.Pfer_h        = sum(sum(SOL.prh))*(2*geo.p/geo.ps);
            out.Pfer_c        = sum(sum(SOL.prc))*(2*geo.p/geo.ps);
            out.Ppm           = sum(sum(SOL.ppm))*(2*geo.p/geo.ps);
            out.ppm_no3D      = sum(sum(SOL.ppm_no3D))*(2*geo.p/geo.ps);
            out.ppm_noRFno3D  = sum(sum(SOL.ppm_noRFno3D))*(2*geo.p/geo.ps);
            out.Ppm_breakdown = SOL.ppm_PM*(2*geo.p/geo.ps);
            out.Pfe           = out.Pfes_h + out.Pfes_c + out.Pfer_h + out.Pfer_c;
            out.velDim        = per.EvalSpeed;

            if strcmp(eval_type,'singmIron')
                % remove all the debug data from SOL, to avoid excessive data size
                SOL = rmfield(SOL,'psh');
                SOL = rmfield(SOL,'psc');
                SOL = rmfield(SOL,'prh');
                SOL = rmfield(SOL,'prc');
                SOL = rmfield(SOL,'ppm');
                %SOL = rmfield(SOL,'ppm_RF');
                %SOL = rmfield(SOL,'ppm_noRF');
                SOL = rmfield(SOL,'ppm_PM');
                SOL = rmfield(SOL,'freq');
                SOL = rmfield(SOL,'bs');
                SOL = rmfield(SOL,'br');
                SOL = rmfield(SOL,'am');
                SOL = rmfield(SOL,'Jm');
                SOL = rmfield(SOL,'pos');
                SOL = rmfield(SOL,'vol');
                SOL = rmfield(SOL,'groNo');
                out.SOL = SOL;
            end
        end
    end
    % standard results
    out.id     = mean(out.SOL.id);                                      % [A]
    out.iq     = mean(out.SOL.iq);                                      % [A]
    out.fd     = mean(out.SOL.fd);                                      % [Vs]
    out.fq     = mean(out.SOL.fq);                                      % [Vs]
    out.T      = mean(out.SOL.T);                                       % [Nm]
    out.dT     = std(out.SOL.T);                                        % [Nm]
    out.dTpu   = std(out.SOL.T)/out.T;                                  % [pu]
    out.dTpp   = max(out.SOL.T)-min(out.SOL.T);                             % [Nm]
    out.IPF    = sin(atan2(out.iq,out.id)-atan2(out.fq,out.fd));
    out.We     = mean(out.SOL.we);                                      % [J]
    out.Wc     = mean(out.SOL.wc);                                      % [J]
    %     out.SOL    = SOL;

    % check Torque sign
    if sign(out.T)~=sign(out.fd*out.iq-out.fq*out.id)
        out.T = -out.T;
        out.SOL.T = -out.SOL.T;
    end

    if isfield(out.SOL,'F')
        out.F=mean(out.SOL.F);
    end

    output{1} =  out;
    SimulatedCurrent = SimulatedCurrent(1);
    CurrLoPP = CurrLoPP(1);
    performance{1}.delta_sim_singt = AngularSpanPP;
end

% if flag_singtpp == 1
% %     output{1} =  out;
%     SimulatedCurrent = SimulatedCurrent(1);
%     CurrLoPP = CurrLoPP(1);
%     performance{1}.delta_sim_singt = AngularSpanPP;
% end

% save output into individual folders
for ii = 1:nSim

    geo = geometry{ii};
    out = output{ii};
    per = performance{ii};
    dirName = tempDirName{ii};

    iStr=num2str(SimulatedCurrent(ii),3); iStr = strrep(iStr,'.','A');
    gammaStr=num2str(GammaPP(ii),4); gammaStr = strrep(gammaStr,'.','d');
    if ~contains(gammaStr, 'd')
        gammaStr = [gammaStr 'd'];
    end

    if dataIn.CustomCurrentEnable
        FILENAME = ['T_eval_CustomCurrent_' datestr(now,30)];
    else
        if strcmp(geo.RotType,'EESM')
            FILENAME = ['T_eval_',iStr,'_',gammaStr '_' int2str(floor(dataIn.FieldCurrent)) 'f' int2str(10*(dataIn.FieldCurrent-floor(dataIn.FieldCurrent)))];
        else
            FILENAME = ['T_eval_',iStr,'_',gammaStr '_' int2str(dataIn.tempPP) 'deg'];
        end
        %     FILENAME = [filemot(1:end-4) '_T_eval_',iStr,'_',gammaStr];
    end

    if(isnan(geo.win.n3phase))
            FILENAME = [FILENAME '_n5ph_' char(datetime('now','Format','uuuuMMdd''T''HHmmss'))];
    else
        if sum(per.flag3phaseSet)~=geo.win.n3phase
            FILENAME = [FILENAME '_' mat2str(per.flag3phaseSet)];
        end
    
        if ((length(SimulatedCurrent)==geo.win.n3phase)&&(geo.win.n3phase~=1))
            % FILENAME = [FILENAME '_n3ph_' char(datetime('now','Format','uuuuMMdd''T''HHmmss'))];
            FILENAME = [FILENAME '_n3ph_' mat2str(CurrLoPP*RatedCurrent,3) '_' mat2str(GammaPP,3)];

        end
    end
    
    switch eval_type
        case 'flxdn'
            FILENAME = [FILENAME '_flxdn'];
        case 'izero'
            FILENAME = [FILENAME '_izero'];
        case 'force'
            FILENAME = [FILENAME '_force'];
        case 'singtIron'
            nStr = int2str(per.EvalSpeed);
            nStr = strrep(nStr,'.','rpm');
            if ~strcmpi(nStr,'rpm')
                nStr = [nStr 'rpm'];
            end
            FILENAME = [FILENAME '_' nStr '_ironLoss'];
    end

    resFolder = checkPathSyntax([filemot(1:end-4) '_results\FEA results\']);
    if ~exist([pathname resFolder],'dir')
        mkdir([pathname resFolder]);
    end

    if flagSave~=0
        mkdir([pathname resFolder],FILENAME);
        newDir = checkPathSyntax([pathname resFolder FILENAME '\']);
    end
    
    if flagSave>0
        if isoctave()            %OCT
            file_name1= strcat(newDir,filemot(1:end-4),'_',FILENAME,'.mat');
            save('-mat7-binary', file_name1,'geo','per','mat','out');
            dirIn=checkPathSyntax(strcat(dirName, ['\' filemot]));
            dirDest=strcat(newDir,filemot(1:end-4),'_',FILENAME,'.fem');
            movefile(dirIn, dirDest);
            clear file_name1 dirIn dirDest
        else
            %save([newDir filemot(1:end-4) '_' FILENAME '.mat'],'geo','per','mat','out');
            save([newDir filemot(1:end-4) '_OpPointResults.mat'],'geo','per','mat','out');
            copyfile([dirName filemot],[newDir filemot(1:end-4) '_solved.fem']);
        end
    end

    % plot and save figs
    delta_sim_singt = per.delta_sim_singt;
    
    if flagSave==2
        if(isnan(n3ph))
            if(per.custom_act)            
                if(isscalar(per.custom_rotorPos))
                    plot_singt_5(out,delta_sim_singt,newDir,filemot);
                end
            else
                plot_singt_5(out,delta_sim_singt,newDir,filemot);
            end
        else        
            plot_singt(out,delta_sim_singt,newDir,filemot);
        end

        switch eval_type
            case 'flxdn'
                plot_flxdn_fig(geo,out,newDir,filemot);
                plot_flxdn_gif(geo,out,newDir,filemot);
            case'force'
                plot_force_fig(geo,out,newDir,filemot);
                plot_force_gif(geo,out,newDir,filemot);
                forceOut = elab_singt_toothForce(geo,per,out,newDir);
                plot_forceOut_gif(forceOut,newDir)
                elab_singt_NVH(geo,per,out,forceOut,newDir);
            case 'singtIron'
                plot_singtIron(geo,out,newDir,filemot);
                if exist('plot_ironLoss_geometry.m','file')
                    plot_ironLoss_geometry(geo,per,mat,out,newDir);
                end
        end
    end
    
    if nSim>5
        close all
    end

end

senseOut = [];

% extra figs, if input current is array
if nSim>1

    if(isnan(n3ph))
        id1 = nan(1,length(CurrLoPP));
        iq1 = nan(1,length(CurrLoPP));
        id3 = nan(1,length(CurrLoPP));
        iq3 = nan(1,length(CurrLoPP));
        fd1 = nan(1,length(CurrLoPP));
        fq1 = nan(1,length(CurrLoPP));
        fd3 = nan(1,length(CurrLoPP));
        fq3 = nan(1,length(CurrLoPP));
    else
        id = nan(1,length(CurrLoPP));
        iq = nan(1,length(CurrLoPP));
        fd = nan(1,length(CurrLoPP));
        fq = nan(1,length(CurrLoPP));
    end

    T = nan(1,length(CurrLoPP));
    dTpu = nan(1,length(CurrLoPP));
    dTpp = nan(1,length(CurrLoPP));
    ir = nan(1,length(CurrLoPP));
    fr = nan(1,length(CurrLoPP));

    if(isnan(n3ph))
        if(per.custom_act)
            T = nan(length(CurrLoPP),per.nsim_singt);
        end
    end

    for ii = 1:length(CurrLoPP)
        if(isnan(n3ph))
            id1(ii) = output{ii}.id1;
            iq1(ii) = output{ii}.iq1;
            id3(ii) = output{ii}.id3;
            iq3(ii) = output{ii}.iq3;
            fd1(ii) = output{ii}.fd1;
            fq1(ii) = output{ii}.fq1;
            fd3(ii) = output{ii}.fd3;
            fq3(ii) = output{ii}.fq3;
        else
            id(ii) = output{ii}.id;
            iq(ii) = output{ii}.iq;
            fd(ii) = output{ii}.fd;
            fq(ii) = output{ii}.fq;
        end

        if(isnan(n3ph))
            if(per.custom_act)
                T(ii,:) = output{ii}.T;
                ia(ii,:) = output{ii}.ia;
                ib(ii,:) = output{ii}.ib;
                ic(ii,:) = output{ii}.ic;
                id(ii,:) = output{ii}.id;
                ie(ii,:) = output{ii}.ie;
            else
                T(ii) = output{ii}.T;
            end
        else
            T(ii) = output{ii}.T;
        end

        dTpu(ii) = output{ii}.dTpu;
        dTpp(ii) = output{ii}.dTpp;

        if isfield(output{ii},'if')
            ir(ii) = output{ii}.if;
            fr(ii) = output{ii}.ff;
        end
        OUT(ii) = output{ii};
    end
    %dirPower = [pathname resFolder filemot(1:end-4) '_singT - ' int2str(dataIn.tempPP) 'deg\'];
    if flagSave~=0
        dirPower = checkPathSyntax([pathname resFolder 'senseOut - ' int2str(dataIn.tempPP) 'deg - ' datestr(now,30) '\']);
        mkdir(dirPower);
    end

    if flagSave==2
        x = 1:length(CurrLoPP);

        flagPlot3D = false;
        flagCharger = false;

        if(isnan(n3ph))
            if(per.custom_act)
                flagPlot3D = true;
                flagCharger = true;
            end
        end

        figure();
        if ~isoctave()
            figSetting();
        end
        if(flagCharger)
            if(flagPlot3D)
                x = geo.p*per.custom_rotorPos;
                y = linspace(0,2*pi,per.nsim_singt);
                surf(x,y,T'), grid on, 
                xlabel('$\theta_{rotor}$ [el. deg.]')
                ylabel('${\omega}t$ [rad]')
                zlabel('$T$ [Nm]')
            else
                x = geo.p*per.custom_rotorPos;
                for i=1:length(x)
                    Trms(i) = sqrt(mean(T(i,:,:).^2));
                    Tpk(i) = max(max(abs(T(i,:,:))));
                end
                subplot(2,1,1)
                plot(x,Trms), grid on, 
                xlabel('$\theta_{rotor}$ [el. deg.]')
                ylabel('$T_{RMS}$ [Nm]')
                subplot(2,1,2)
                plot(x,Tpk), grid on, 
                xlabel('$\theta_{rotor}$ [el. deg.]')
                ylabel('$T_{Peak}$ [Nm]')
            end
        else
            subplot(2,1,1)
            plot(x,T,'-x',x,T+0.5*dTpp,'r',x,T-0.5*dTpp,'r'), grid on, ylabel('$T$ [Nm]')
            subplot(2,1,2)
            plot(x,dTpp,'-x'), grid on, ylabel('$\Delta T_{pp}$ [Nm]')
            xlabel('simulation \#')
        end
        h=gcf();
        if isoctave() %OCT
            fig_name=strcat(dirPower, filemot(1:end-4), '_torque_sens');
            hgsave(h,[fig_name]);
        else
            saveas(gcf,[dirPower,filemot(1:end-4),'_torque_sens.fig'])
        end

        if(flagCharger)

            figure()
            if ~isoctave()
                figSetting();
            end

            loss_pu = (ia.^2+ib.^2+ic.^2+id.^2+ie.^2)./(ia.^2+0.5*(ib.^2+id.^2)+0.5*(ic.^2+ie.^2));
            x = geo.p*per.custom_rotorPos;
            y = linspace(0,2*pi,per.nsim_singt);
            surf(x,y,loss_pu'), grid on, 
            xlabel('$\theta_{rotor}$ [el. deg.]')
            ylabel('${\omega}t$ [rad]')
            zlabel('$Joule losses$ [pu]')
        else

            figure()
            if ~isoctave()
                figSetting();
            end

            if(isnan(n3ph))
                subplot(2,1,1)
                plot(x,fd1,'-x',x,fq1,'-x'), grid on, ylabel('[Vs]'), legend('$\lambda_{d1}$','$\lambda_{q1}$'),
                subplot(2,1,2)
                plot(x,abs(sin(atan(iq1./id1)-atan(fq1./fd1))),'-x'), grid on, ylabel('$cos \varphi$')
                xlabel('simulation \#'),
            else
                subplot(2,1,1)
                plot(x,fd,'-x',x,fq,'-x'), grid on, ylabel('[Vs]'), legend('$\lambda_d$','$\lambda_q$'),
                subplot(2,1,2)
                plot(x,abs(sin(atan(iq./id)-atan(fq./fd))),'-x'), grid on, ylabel('$cos \varphi$')
                xlabel('simulation \#'),
            end

            h=gcf();
            if isoctave() %OCT
                fig_name=strcat(dirPower, filemot(1:end-4), '_fdq_IPF_sens');
                hgsave(h,[fig_name]);
            else
                saveas(gcf,[dirPower,filemot(1:end-4),'_fdq_IPF_sens.fig'])
            end
        end

        if(flagCharger)

            figure()
            if ~isoctave()
                figSetting();
            end

            x = geo.p*per.custom_rotorPos;
            y = linspace(0,2*pi,per.nsim_singt);
            surf(x,y,ia'), grid on, 
            xlabel('$\theta_{rotor}$ [el. deg.]')
            ylabel('${\omega}t$ [rad]')
            zlabel('$I_a$ [A]')

            figure()
            if ~isoctave()
                figSetting();
            end

            x = geo.p*per.custom_rotorPos;
            y = linspace(0,2*pi,per.nsim_singt);
            surf(x,y,ib'), grid on, 
            xlabel('$\theta_{rotor}$ [el. deg.]')
            ylabel('${\omega}t$ [rad]')
            zlabel('$I_b$ [A]')

            figure()
            if ~isoctave()
                figSetting();
            end

            x = geo.p*per.custom_rotorPos;
            y = linspace(0,2*pi,per.nsim_singt);
            surf(x,y,ic'), grid on, 
            xlabel('$\theta_{rotor}$ [el. deg.]')
            ylabel('${\omega}t$ [rad]')
            zlabel('$I_c$ [A]')

            figure()
            if ~isoctave()
                figSetting();
            end

            x = geo.p*per.custom_rotorPos;
            y = linspace(0,2*pi,per.nsim_singt);
            surf(x,y,id'), grid on, 
            xlabel('$\theta_{rotor}$ [el. deg.]')
            ylabel('${\omega}t$ [rad]')
            zlabel('$I_d$ [A]')

            figure()
            if ~isoctave()
                figSetting();
            end

            x = geo.p*per.custom_rotorPos;
            y = linspace(0,2*pi,per.nsim_singt);
            surf(x,y,ib'), grid on, 
            xlabel('$\theta_{rotor}$ [el. deg.]')
            ylabel('${\omega}t$ [rad]')
            zlabel('$I_e$ [A]')

        else

            figure()
            if ~isoctave()
                figSetting();
            end

            if(isnan(n3ph))
                subplot(2,1,1)
                plot(x,fd1,'-x','DisplayName','$\lambda_{d1}$');
                plot(x,fq1,'-x','DisplayName','$\lambda_{q1}$');
                ylabel('[Vs]')
                legend('show');
                subplot(2,1,2)
                plot(x,id1,'-x','DisplayName','$i_{d1}$');
                plot(x,iq1,'-x','DisplayName','$i_{q1}$');
                xlabel('simulation \#')
                ylabel('[A]')
            else
                subplot(2,1,1)
                plot(x,fd,'-x','DisplayName','$\lambda_d$');
                plot(x,fq,'-x','DisplayName','$\lambda_q$');
                ylabel('[Vs]')
                legend('show');
                subplot(2,1,2)
                plot(x,id,'-x','DisplayName','$i_d$');
                plot(x,iq,'-x','DisplayName','$i_q$');
                xlabel('simulation \#')
                ylabel('[A]')
            end
        end

        legend('show');
        h=gcf();
        if isoctave() %OCT
            fig_name=strcat(dirPower, filemot(1:end-4), '_fdq_idiq_sens');
            hgsave(h,[fig_name]);
        else
            saveas(gcf,[dirPower,filemot(1:end-4),'_fdq_idiq_sens.fig'])
        end
    end

    if(isnan(n3ph))
        senseOut.id1   = id1;
        senseOut.iq1   = iq1;
        senseOut.id3   = id3;
        senseOut.iq3   = iq3;
        senseOut.fd1   = fd1;
        senseOut.fq1   = fq1;
        senseOut.fd3   = fd3;
        senseOut.fq3   = fq3;
        senseOut.PF   = abs(sin(atan(iq1./id1)-atan(fq1./fd1)));
    else
        senseOut.id   = id;
        senseOut.iq   = iq;
        senseOut.fd   = fd;
        senseOut.fq   = fq;
        senseOut.PF   = abs(sin(atan(iq./id)-atan(fq./fd)));
    end

    senseOut.T    = T;
    senseOut.dTpp = dTpp;
    senseOut.if   = ir;
    senseOut.ff   = fr;
    senseOut.OUT  = OUT;
    if flagSave>0
        save([dirPower,filemot(1:end-4),'_senseResults.mat'],'senseOut');
    end
end

if nargout()==0
    clear out senseOut newDir
end

