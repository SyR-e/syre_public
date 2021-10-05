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

function [Inductance] = MMM_eval_inductanceMap(motorModel)

% Computation of the incremental inductances

Id = motorModel.FluxMap_dq.Id;
Iq = motorModel.FluxMap_dq.Iq;
Fd = motorModel.FluxMap_dq.Fd;
Fq = motorModel.FluxMap_dq.Fq;

[dFdd,dFdq] = gradient(Fd);
[dFqd,dFqq] = gradient(Fq);
[dIdd,~]    = gradient(Id);
[~,dIqq]    = gradient(Iq);

Ldd = dFdd./dIdd;
Ldq = dFdq./dIqq;
Lqd = dFqd./dIdd;
Lqq = dFqq./dIqq;

% output data
Inductance.Id = Id;
Inductance.Iq = Iq;
Inductance.Ldd = Ldd;
Inductance.Ldq = Ldq;
Inductance.Lqd = Lqd;
Inductance.Lqq = Lqq;
end