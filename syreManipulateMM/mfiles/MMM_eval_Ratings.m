% Copyright 2020
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [motorModel] = MMM_eval_Ratings(motorModel,flagPlot)

if nargin==1
    flagPlot=1;
end

i0    = motorModel.data.i0;
Vdc   = motorModel.data.Vdc;

Plim{1} = OpLimEval(motorModel,i0,Vdc);

if flagPlot
    OpLimPlot(Plim,i0,motorModel);
end

motorModel.data.n0 = round(Plim{1}.n_A,0);
motorModel.data.T0 = round(Plim{1}.T_A,2);
motorModel.data.P0 = round(Plim{1}.T_A*Plim{1}.n_A*pi/30,2);



