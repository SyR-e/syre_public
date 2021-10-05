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

function eval_operatingPoint(dataIn)

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
SimulatedCurrent = dataIn.SimulatedCurrent;
GammaPP  = dataIn.GammaPP;
BrPP = dataIn.BrPP;
NumOfRotPosPP = dataIn.NumOfRotPosPP;
AngularSpanPP = dataIn.AngularSpanPP;
NumGrid = dataIn.NumGrid;

per.EvalSpeed = dataIn.EvalSpeed;

clc;

eval_type = dataIn.EvalType;

per.overload=CurrLoPP;
per.i0 = RatedCurrent;
per.BrPP=BrPP;

per.nsim_singt = NumOfRotPosPP;       % # simulated positions
per.delta_sim_singt = AngularSpanPP;  % angular span of simulation

% single point or array of points simulation
performance = cell(1,length(CurrLoPP));
output = cell(1,length(CurrLoPP));
geometry = cell(1,length(CurrLoPP));
tempDirName = cell(1,length(CurrLoPP));
for ii = 1:length(CurrLoPP)
    performance{ii} = per;
    performance{ii}.overload = CurrLoPP(ii);
    performance{ii}.gamma=GammaPP(ii);
end
geo.RemoveTMPfile = 'OFF';
% check parallel computing
ppState=parallelComputingCheck();
if (ppState==0 && length(CurrLoPP)>4)
    parpool();
    ppState=parallelComputingCheck();
end

fileMotWithPath=[pathname filemot];

geo0=geo;
mat0=mat;
% evaluation
if ppState<1
    for ii = 1:length(CurrLoPP)
        geoTmp = geo0;
        perTmp = performance{ii};
        matTmp = mat0;
        [~,geometry{ii},~,output{ii},tempDirName{ii}] = FEMMfitness([],geoTmp,perTmp,matTmp,eval_type,fileMotWithPath);
    end
else
    parfor ii = 1:length(CurrLoPP) %%%
        geoTmp = geo0;
        perTmp = performance{ii};
        matTmp = mat0;
        [~,geometry{ii},~,output{ii},tempDirName{ii}] = FEMMfitness([],geoTmp,perTmp,matTmp,eval_type,fileMotWithPath);
    end
end

% save output into individual folders
for ii = 1:length(SimulatedCurrent)
    
    geo = geometry{ii};
    out = output{ii};
    per = performance{ii};
    dirName = tempDirName{ii};
    
    iStr=num2str(SimulatedCurrent(ii),3); iStr = strrep(iStr,'.','A');
    gammaStr=num2str(GammaPP(ii),4); gammaStr = strrep(gammaStr,'.','d');
    if ~contains(gammaStr, 'd')
        gammaStr = [gammaStr 'd'];
    end
    
    FILENAME = ['T_eval_',iStr,'_',gammaStr '_' int2str(dataIn.tempPP) 'deg'];
%     FILENAME = [filemot(1:end-4) '_T_eval_',iStr,'_',gammaStr];
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
        save([newDir filemot(1:end-4) '_' FILENAME '.mat'],'geo','per','mat','out');
        copyfile([dirName filemot],[newDir filemot(1:end-4) '_' FILENAME '.fem']);
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
        case 'singtIron'
            plot_singtIron(geo,out,newDir,filemot);
    end
    
end

% extra figs, if input current is array
if length(CurrLoPP)>1
    
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



