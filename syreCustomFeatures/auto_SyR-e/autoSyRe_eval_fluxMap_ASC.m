% Copyright 2023
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

function autoSyRe_eval_fluxMap_ASC(dataSet,mat)

%% Input
prompt = {'Rotor angular excursion [elt °]','Number of rotor positions:','Phase current [A]:','PM temperature [° C]:','Number of points in q-axis:'};
dlgtitle = 'Input';
dims = [1 40];
definput = {num2str(60),'12',num2str(round(dataSet.SimulatedCurrent)),num2str(dataSet.tempPP),num2str(dataSet.NumGrid)};
answer = inputdlg(prompt,dlgtitle,dims,definput);

dataSet.AngularSpanPP    = str2double(answer{1});
dataSet.NumOfRotPosPP    = str2double(answer{2});
dataSet.SimulatedCurrent = str2double(answer{3});
dataSet.tempPP           = str2double(answer{4});
dataSet.NumGrid          = str2double(answer{5});

dataSet.MapQuadrants = 2;
dataSet.CurrLoPP = dataSet.SimulatedCurrent/dataSet.RatedCurrent;

dataSet.BrPP = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,dataSet.tempPP);


if strcmp(dataSet.axisType,'SR')
    dataSet.GammaPP = 0;
else
    dataSet.GammaPP = 90;
end


%% 
[pkSCout] = eval_peakShortCircuitCurrent(dataSet);
dataSet.EvalType = 'singm';
dataSet.SimulatedCurrent_HWC = abs(pkSCout.idq);

dataSet.flagHWCMap = 1;


eval_fluxMap(dataSet)

