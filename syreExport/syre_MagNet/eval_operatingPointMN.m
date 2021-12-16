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

function eval_operatingPointMN(dataIn)

% simulates an existing machine in MN, singt mode

% simulates single or multiple (id,iq) conditions
% example inputs:
% single condition: CurrLoPP = 1, GammaPP = 45
% multiple points:  CurrLoPP = [1 1.5 2], gamma = [45 45 45]

% Open matlabpool manually prior to execution

pathname=dataIn.currentpathname;
filemot= dataIn.currentfilename;
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
% syreRoot = fileparts(which('GUI_Syre.mlapp'));
% current_path = syreRoot;

overload_temp =  CurrLoPP;   % current to be simulated
gamma_temp = GammaPP;        % current phase angle
Br = BrPP;                   % remanence of all barriers magnets

eval_type='singt';

per.overload=CurrLoPP;
per.i0 = RatedCurrent;
per.BrPP=BrPP;

per.nsim_singt = NumOfRotPosPP;       % # simulated positions
per.delta_sim_singt = AngularSpanPP;  % angular span of simulation


% single or multiple simulation
performance = cell(1,length(overload_temp));
output = cell(1,length(overload_temp));
geometry = cell(1,length(overload_temp));
tempDirName = cell(1,length(overload_temp));
for ii = 1:length(overload_temp)
    performance{ii} = per;
    performance{ii}.overload = overload_temp(ii);
    performance{ii}.gamma=gamma_temp(ii);
end

geo.RemoveTMPfile = 'OFF';
for ii = 1:length(overload_temp)
    [geometry{ii},mat,output{ii},tempDirName{ii}] = MNfitness([],geo,performance{ii},mat,eval_type,pathname,filemot);
end

% save output into individual folders
for ii = 1:length(SimulatedCurrent)
    
    geo = geometry{ii};
    out = output{ii};
    per = performance{ii};
    dirName = tempDirName{ii};
    
    iStr=num2str(SimulatedCurrent(ii),3); iStr = strrep(iStr,'.','A');
    gammaStr=num2str(gamma_temp(ii),4); gammaStr = strrep(gammaStr,'.','d');
    if isempty(strfind(gammaStr, 'd'))
        gammaStr = [gammaStr 'd'];
    end
    nStr = int2str(per.EvalSpeed);
    nStr = strrep(nStr,'.','rpm');
    if ~strcmpi(nStr,'rpm')
        nStr = [nStr 'rpm'];
    end
    
    resFolder = [filemot(1:end-4) '_results\FEA results\'];
    if ~exist([pathname resFolder],'dir')
        mkdir([pathname resFolder]);
    end
    
    FILENAME = ['T_eval_',iStr,'_',gammaStr '_' int2str(dataIn.tempPP) 'deg' '_' nStr '_MN'];
    
    mkdir([pathname resFolder],FILENAME);
    newDir=[pathname resFolder FILENAME '\'];
    
    %     copyfile([dirName filemot],[newDir filemot]);
    %     save([newDir filemot],'geo','per','mat','out','-append');
    save([newDir filemot],'geo','per','mat','out');
    copyfile([dirName strrep(filemot,'.mat','.mn')],[newDir strrep(filemot,'.mat','.mn')]);
    
    % plot and save figs
    delta_sim_singt = per.delta_sim_singt;
    plot_singt(out,delta_sim_singt,newDir,filemot);
    if delta_sim_singt==360
        plot_singtIron(geo,out,newDir,filemot);
    end
    
end

% extra figs, if input current is array
if length(overload_temp)>1
    
    id = zeros(1,length(overload_temp));
    iq = zeros(1,length(overload_temp));
    T = zeros(1,length(overload_temp));
    dTpu = zeros(1,length(overload_temp));
    dTpp = zeros(1,length(overload_temp));
    fd = zeros(1,length(overload_temp));
    fq = zeros(1,length(overload_temp));
    
    for ii = 1:length(overload_temp)
        id(ii) = output{ii}.id;
        iq(ii) = output{ii}.iq;
        T(ii) = output{ii}.T;
        dTpu(ii) = output{ii}.dTpu;
        dTpp(ii) = output{ii}.dTpp;
        fd(ii) = output{ii}.fd;
        fq(ii) = output{ii}.fq;
    end
    dirPower=[pathname resFolder filemot(1:end-4) 'singT_MN\'];
    mkdir(dirPower);
    
    x = 1:length(overload_temp);
    figure(10), subplot(2,1,1)
    plot(x,T,'-x',x,T+0.5*dTpp,'r',x,T-0.5*dTpp,'r'), grid on, ylabel('torque [Nm]')
    subplot(2,1,2)
    plot(x,dTpp,'-x'), grid on, ylabel('torque ripple pk-pk [Nm]')
    xlabel('simulation #'),
    h=gcf();
    if isoctave() %OCT
        fig_name=strcat(dirPower, filemot(1:end-4), '_torque_sens');
        hgsave(h,[fig_name]);
    else
        saveas(gcf,[dirPower,filemot(1:end-4),'_torque_sens.fig'])
    end
    
    figure(11), subplot(2,1,1)
    plot(x,fd,'-x',x,fq,'-x'), grid on, ylabel('[Vs]'), legend('\lambda_d','\lambda_q'),
    subplot(2,1,2)
    plot(x,sin(atan(iq./id)-atan(fq./fd)),'-x'), grid on, ylabel('IPF')
    xlabel('simulation #'),
    h=gcf();
    if isoctave() %OCT
        fig_name=strcat(dirPower, filemot(1:end-4), '_fdq_IPF_sens');
        hgsave(h,[fig_name]);
    else
        saveas(gcf,[dirPower,filemot(1:end-4),'_fdq_IPF_sens.fig'])
    end
end




