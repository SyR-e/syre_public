% Copyright 2022
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

function [Inductance] = MMM_eval_appInductanceMap(motorModel)

% Computation of the apparent inductances

Id = motorModel.FluxMap_dq.Id;
Iq = motorModel.FluxMap_dq.Iq;
Fd = motorModel.FluxMap_dq.Fd;
Fq = motorModel.FluxMap_dq.Fq;

axisType = motorModel.data.axisType;

switch axisType
    case 'PM'
        Fm = interp2(Id,Iq,Fd,zeros(size(Id)),Iq);
        Ld = (Fd-Fm)./Id;
        Lq = Fq./Iq;
    case 'SR'
        Fm = interp2(Id,Iq,Fq,Id,zeros(size(Iq)));
        Ld = Fd./Id;
        Lq = (Fq-Fm)./Iq;
end

% output data
Inductance.Id = Id;
Inductance.Iq = Iq;
Inductance.Ld = Ld;
Inductance.Lq = Lq;
Inductance.Fm = Fm;
end