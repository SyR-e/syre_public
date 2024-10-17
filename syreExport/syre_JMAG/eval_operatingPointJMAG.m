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
NumOfRotPosPP = dataIn.NumOfRotPosPP; 
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

[geometry,mat,output,tempDirName] = JMAGfitness([],geo,performance,mat,pathname,filemot);

% save output into individual folders
geo = geometry;
out = output;
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
plot_singt(out,delta_sim_singt,newDir,filemot);
if delta_sim_singt==360
    plot_singtIron(geo,out,newDir,filemot);
end
end

