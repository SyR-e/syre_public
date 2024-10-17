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

function [out,senseOut,newDir] = eval_operatingPoint(dataIn)

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

RatedCurrent = dataIn.RatedCurrent;
CurrLoPP = dataIn.CurrLoPP;
%SimulatedCurrent = dataIn.SimulatedCurrent;
SimulatedCurrent = RatedCurrent*CurrLoPP;
GammaPP  = dataIn.GammaPP;
BrPP = dataIn.BrPP;
NumOfRotPosPP = dataIn.NumOfRotPosPP;
AngularSpanPP = dataIn.AngularSpanPP;
per.flag3phaseSet = dataIn.Active3PhaseSets;
n3ph = dataIn.Num3PhaseCircuit;


per.EvalSpeed = dataIn.EvalSpeed;

clc;

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

per.nsim_singt = NumOfRotPosPP;       % # simulated positions
per.delta_sim_singt = AngularSpanPP;  % angular span of simulation

%custom current
if dataIn.CustomCurrentEnable
    per.custom_ia         = dataIn.CustomCurrentA;
    per.custom_ib         = dataIn.CustomCurrentB;
    per.custom_ic         = dataIn.CustomCurrentC;
    per.custom_time       = dataIn.CustomCurrentTime;
    per.custom_act        = dataIn.CustomCurrentEnable;
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

if (strcmp(eval_type,'singtIron') || flag_singtpp == 1)
    CurrLoPP = CurrLoPP*ones(1,ppState1);
    GammaPP  = GammaPP*ones(1,ppState1);
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
    performance{ii}.overload = CurrLoPP(ii);
    performance{ii}.gamma    = GammaPP(ii);
    performance{ii}.offset   = offset(ii);
end

% if strcmp(eval_type,'singtIron') && length(CurrLoPP)>1
%     performance{1}.nsim_singt = performance{1}.nsim_singt+1;
% end

nSim = length(CurrLoPP);
if nSim==n3ph
    nSim = 1;
    performance{1}          = per;
    performance{1}.overload = CurrLoPP;
    performance{1}.gamma    = GammaPP;
    performance{1}.offset   = zeros(1,n3ph);
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


if (strcmp(eval_type,'singtIron')||strcmp(eval_type,'singmIron')||flag_singtpp == 1)

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
        FILENAME = ['T_eval_',iStr,'_',gammaStr '_' int2str(dataIn.tempPP) 'deg'];
        %     FILENAME = [filemot(1:end-4) '_T_eval_',iStr,'_',gammaStr];
    end

    if sum(per.flag3phaseSet)~=geo.win.n3phase
        FILENAME = [FILENAME '_' mat2str(per.flag3phaseSet)];
    end

    if length(SimulatedCurrent)==geo.win.n3phase
        FILENAME = [FILENAME '_n3ph_' char(datetime('now','Format','uuuuMMdd''T''HHmmss'))];
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

    resFolder = [filemot(1:end-4) '_results\FEA results\'];
    if ~exist([pathname resFolder],'dir')
        mkdir([pathname resFolder]);
    end

    mkdir([pathname resFolder],FILENAME);
    newDir=[pathname resFolder FILENAME '\'];

    if isoctave()            %OCT
        file_name1= strcat(newDir,filemot(1:end-4),'_',FILENAME,'.mat');
        save('-mat7-binary', file_name1,'geo','per','mat','out');
        dirIn=strcat(dirName, ['\' filemot]);
        dirDest=strcat(newDir,filemot(1:end-4),'_',FILENAME,'.fem');
        movefile(dirIn, dirDest);
        clear file_name1 dirIn dirDest
    else
        %save([newDir filemot(1:end-4) '_' FILENAME '.mat'],'geo','per','mat','out');
        save([newDir filemot(1:end-4) '_OpPointResults.mat'],'geo','per','mat','out');
        copyfile([dirName filemot],[newDir filemot(1:end-4) '_solved.fem']);
    end

    % plot and save figs
    delta_sim_singt = per.delta_sim_singt;

    plot_singt(out,delta_sim_singt,newDir,filemot);

    switch eval_type
        case 'flxdn'
            plot_flxdn_fig(geo,out,newDir,filemot);
            plot_flxdn_gif(geo,out,newDir,filemot);
        case'force'
            plot_force_fig(geo,out,newDir,filemot);
            plot_force_gif(geo,out,newDir,filemot);
            forceOut = elab_singt_toothForce(geo,per,out,newDir);
            plot_forceOut_gif(forceOut,newDir)
        case 'singtIron'
            plot_singtIron(geo,out,newDir,filemot);
    end

end

senseOut = [];

% extra figs, if input current is array
if nSim>1

    id = zeros(1,length(CurrLoPP));
    iq = zeros(1,length(CurrLoPP));
    T = zeros(1,length(CurrLoPP));
    dTpu = zeros(1,length(CurrLoPP));
    dTpp = zeros(1,length(CurrLoPP));
    fd = zeros(1,length(CurrLoPP));
    fq = zeros(1,length(CurrLoPP));

    for ii = 1:length(CurrLoPP)
        id(ii) = output{ii}.id;
        iq(ii) = output{ii}.iq;
        T(ii) = output{ii}.T;
        dTpu(ii) = output{ii}.dTpu;
        dTpp(ii) = output{ii}.dTpp;
        fd(ii) = output{ii}.fd;
        fq(ii) = output{ii}.fq;
    end
    %dirPower = [pathname resFolder filemot(1:end-4) '_singT - ' int2str(dataIn.tempPP) 'deg\'];
    dirPower = [pathname resFolder 'senseOut - ' int2str(dataIn.tempPP) 'deg - ' datestr(now,30) '\'];
    mkdir(dirPower);

    x = 1:length(CurrLoPP);
    figure();
    if ~isoctave()
        figSetting();
    end
    subplot(2,1,1)
    plot(x,T,'-x',x,T+0.5*dTpp,'r',x,T-0.5*dTpp,'r'), grid on, ylabel('$T$ [Nm]')
    subplot(2,1,2)
    plot(x,dTpp,'-x'), grid on, ylabel('$\Delta T_{pp}$ [Nm]')
    xlabel('simulation \#')
    h=gcf();
    if isoctave() %OCT
        fig_name=strcat(dirPower, filemot(1:end-4), '_torque_sens');
        hgsave(h,[fig_name]);
    else
        saveas(gcf,[dirPower,filemot(1:end-4),'_torque_sens.fig'])
    end

    figure()
    if ~isoctave()
        figSetting();
    end
    subplot(2,1,1)
    plot(x,fd,'-x',x,fq,'-x'), grid on, ylabel('[Vs]'), legend('$\lambda_d$','$\lambda_q$'),
    subplot(2,1,2)
    plot(x,abs(sin(atan(iq./id)-atan(fq./fd))),'-x'), grid on, ylabel('$cos \varphi$')
    xlabel('simulation \#'),
    h=gcf();
    if isoctave() %OCT
        fig_name=strcat(dirPower, filemot(1:end-4), '_fdq_IPF_sens');
        hgsave(h,[fig_name]);
    else
        saveas(gcf,[dirPower,filemot(1:end-4),'_fdq_IPF_sens.fig'])
    end

    figure()
    if ~isoctave()
        figSetting();
    end
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
    legend('show');
    h=gcf();
    if isoctave() %OCT
        fig_name=strcat(dirPower, filemot(1:end-4), '_fdq_idiq_sens');
        hgsave(h,[fig_name]);
    else
        saveas(gcf,[dirPower,filemot(1:end-4),'_fdq_idiq_sens.fig'])
    end
    senseOut.id   = id;
    senseOut.iq   = iq;
    senseOut.fd   = fd;
    senseOut.fq   = fq;
    senseOut.T    = T;
    senseOut.dTpp = dTpp;
    senseOut.PF   = abs(sin(atan(iq./id)-atan(fq./fd)));
    save([dirPower,filemot(1:end-4),'_senseResults.mat'],'senseOut');

end

if nargout()==0
    clear out senseOut newDir
end

