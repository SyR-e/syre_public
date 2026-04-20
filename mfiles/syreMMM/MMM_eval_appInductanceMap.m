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

axisType  = motorModel.data.axisType;
motorType = motorModel.data.motorType;
p         = motorModel.data.p;

switch axisType
    case 'PM'
        if strcmp(motorType, 'VFM')
            MS = motorModel.FluxMap_dq.MS;            
            Fdm = zeros(1,size(MS,3));
            Fqm = zeros(1,size(MS,3));
            FLd = zeros(size(Fd));
            FLq = zeros(size(Fd));

            Fdm = interp3(Id, Iq, MS, Fd, zeros(size(Id)), Iq, MS);
            Fqm = interp3(Id, Iq, MS, Fq, Id, zeros(size(Iq)), MS);
            FLd = Fd - Fdm;
            FLq = Fq - Fqm;

            Ld = FLd./Id;
            Lq = FLq./Iq;
        else
            Fm = interp2(Id,Iq,Fd,zeros(size(Id)),Iq);
            Ld = (Fd-Fm)./Id;
            Lq = Fq./Iq;
        end
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
if strcmp(motorType, 'VFM')
    Inductance.MS = MS; 
    Inductance.Fdm = Fdm;
    Inductance.Fqm = Fqm;
    Inductance.Tm = 3/2*p*(Fdm.*Iq - Fqm.*Id);
    Inductance.Tr = 3/2*p*(Ld-Lq).*Id.*Iq;

else
    Inductance.Fm = Fm;
    Inductance.Tm = 3/2*p*Fm.*Iq;
    Inductance.Tr = 3/2*p*(Ld-Lq).*Id.*Iq;
end



end