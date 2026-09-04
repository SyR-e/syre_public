%  Copyright 2024
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

function eval_operatingPointJMAG(dataIn)

% Open matlabpool manually prior to execution
% dataIn.currentpathname='D:\syre_maedeh\syreExport\syre_JMAG\';
% dataIn.currentfilename='syreDefaultMotor.mat';%syreDefaultMotor,THOR,SPM1

pathname=dataIn.currentpathname;
filemot= dataIn.currentfilename;
load([pathname filemot]);

% geo.RotType = 'Hybrid';

n3ph = dataIn.Num3PhaseCircuit;
NumOfRotPosPP_     = dataIn.NumOfRotPosPP;

flagCharger = false;
if(isinf(NumOfRotPosPP_))
    flagCharger = true;
    NumOfRotPosPP_ = length(dataIn.CustomCurrentRotorPos);
end

if(flagCharger)
    nSim = length(dataIn.CustomCurrentRotorPos);
else
    nSim = 1;
end

if ~isfield(geo,'axisType')
    if strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype')
        geo.axisType = 'PM';
    else
        geo.axisType = 'SR';
    end
end

if ~strcmp(geo.axisType,dataIn.axisType)
    %geo.axisType = dataIn.axisType;
    if strcmp(dataIn.axisType,'PM')
        geo.th0 = geo.th0-90;
    else
        geo.th0 = geo.th0+90;
    end
end

%custom current
if dataIn.CustomCurrentEnable
    if(flagCharger)
        if(isnan(n3ph))
            per.custom_ia         = dataIn.CustomCurrentA;
            per.custom_ib         = dataIn.CustomCurrentB;
            per.custom_ic         = dataIn.CustomCurrentC;
            per.custom_id         = dataIn.CustomCurrentD;
            per.custom_ie         = dataIn.CustomCurrentE;
        else
            if(n3ph==2)
                per.custom_ia         = dataIn.CustomCurrentA;
                per.custom_ib         = dataIn.CustomCurrentB;
                per.custom_ic         = dataIn.CustomCurrentC;
                per.custom_id         = dataIn.CustomCurrentD;
                per.custom_ie         = dataIn.CustomCurrentE;
                per.custom_if         = dataIn.CustomCurrentF;
            end
            per.custom_win_delta      = dataIn.CustomWinDelta;
        end
        per.custom_Amp         = dataIn.CustomCurrentAmp;
        per.custom_Ph         = dataIn.CustomCurrentPh;
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

% '------------------------------------------------------------------------
%% 'Motor Drive definitions:
% Drive_Type = 'Current';
% '------------------------------------------------------------------------
per.overload = dataIn.CurrLoPP;
per.i0 = dataIn.RatedCurrent;
per.BrPP = dataIn.BrPP;
per.EvalSpeed = dataIn.EvalSpeed;
gamma_temp = dataIn.GammaPP;%Phase Agle of Drive
% '------------------------------------------------------------------------

%% 'Study property settings: Step Control and Resolution
% 'Resolution or number of Step division per cycle of calculations
NumOfRotPosPP = NumOfRotPosPP_; 
%'#Number of electric periods to be calculed
AngularSpanPP = dataIn.AngularSpanPP;
% '------------------------------------------------------------------------
per.nsim_singt = NumOfRotPosPP;       % # simulated positions
per.delta_sim_singt = AngularSpanPP;  % angular span of simulation
% '------------------------------------------------------------------------
% ' Parallel Computing: (--> SMP is not supported in hysteresis analysis or when using hysteresis loop data)
% 0 or 'NonParallel' or Off   % 1 or 'SMP' or On    % 2 or 'MPP' or DMP
% useCPU=1;
% '#Number of CPU for multiprocessing calculation (even number > 0)% For SMP: nCPU <= 36% For MPP: nCPU <= 512
% nCPU=8; 
%% '-----------------------------------------------------------------------
performance = per;
performance.gamma=gamma_temp;

for ii = 1:nSim

    % disp(strcat(num2str(ii),'/',num2str(nSim)))
    geoTmp = geo;
    perTmp = performance;
    matTmp = mat;

    if(flagCharger)
        perTmp.rotorPos = ii;
    else
        perTmp.rotorPos = NaN;
    end

    % [geometry,mat,output,tempDirName] = JMAGfitness([],geo,performance,mat,pathname,filemot);
    [geometry{ii},~,output{ii},tempDirName] = JMAGfitness([],geoTmp,perTmp,matTmp,pathname,filemot);
end

% save output into individual folders
geo = geometry{1};
out = output{1};
per = performance;
dirName = tempDirName;

iStr=num2str(dataIn.SimulatedCurrent,3); iStr = strrep(iStr,'.','A');
gammaStr=num2str(gamma_temp,4); gammaStr = strrep(gammaStr,'.','d');
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

FILENAME = ['T_eval_',iStr,'_',gammaStr '_' int2str(dataIn.tempPP) 'deg' '_' nStr '_JMAG'];

mkdir([pathname resFolder],FILENAME);
newDir=[pathname resFolder FILENAME '\'];

save([newDir filemot],'geo','per','mat','out');
copyfile(fullfile(dirName, strcat(strrep(filemot,'.mat','.jmag'),'.jproj')),fullfile (newDir, strcat(strrep(filemot,'.mat','.jmag'),'.jproj'))); % copy .jproj in the temporary folder

% plot and save figs
delta_sim_singt = 360;

if nSim>1
    if(flagCharger)
       plot_singt_chargerJMAG(geo.p, nSim, output, n3ph, per.custom_rotorPos);
    end
else
    if(flagCharger)
        if(isnan(n3ph))
            plot_singt_5(out,delta_sim_singt,newDir,filemot);
        else
            plot_singt_3n(n3ph, out,delta_sim_singt,newDir,filemot);
        end
    else
        plot_singt(out,delta_sim_singt,newDir,filemot);
    end
    
    if delta_sim_singt==360
        plot_singtIron(geo,out,newDir,filemot);
    end
end

