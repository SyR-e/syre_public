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

function [F_map,F_IM] = elab_Fmap_IM(F_map,geo,per)

[nR,nC] = size(F_map.Id);
F_IM = F_map.IM;
F_bar = F_map.bar;

% aggiunta induttanza di testata statore
F_map.Fd = F_map.Fd+per.Lend*F_map.Id;
F_map.Fq = F_map.Fq+per.Lend*F_map.Iq;
F_map.T  = 3/2*geo.p*(F_map.Fd.*F_map.Iq-F_map.Fq.*F_map.Id);

% filtro kr=0
F_IM.kr(F_IM.kr==0) = NaN;

% calcolo della resistenza di rotore tramite potenze
F_IM.Pbar = zeros([nR,nC]);
F_IM.Pring = zeros([nR,nC]);


for rr=1:nR
    for cc=1:nC
        F_IM.Pbar(rr,cc)  = sum(F_bar.V{rr,cc}.*F_bar.I{rr,cc})*2*geo.p/geo.ps;
        F_IM.Pring(rr,cc) = sum((2*per.IM.Rring)*(F_bar.I{rr,cc}*geo.IM.k).^2)*2*geo.p/geo.ps;
    end
end

F_IM.Pr = F_IM.Pbar+F_IM.Pring;
F_IM.Rr = F_IM.Pr./(3/2*F_IM.Ir.^2);
F_IM.Rs = per.Rs*ones([nR,nC]);

% Inverse-Gamma equivalent circuit parameters
F_IM.Ls    = F_map.Fd./F_map.Id;            % stator inductance (Ld)
F_IM.Lt    = F_map.Fq./F_map.Iq;            % transient inductance, sigma*Ls (Lq)
F_IM.sigma = F_IM.Lt./F_IM.Ls;              % leakage factor (Lq/Ld, inverse anisotropy factor)
F_IM.Lm    = (F_IM.Ls-F_IM.Lt)./F_IM.kr;    % magnetization inductance
F_IM.Lr    = F_IM.Lm./F_IM.kr;              % rotor inductance
F_IM.Lls   = F_IM.Ls-F_IM.Lm;               % stator leakage inductance
F_IM.Llr   = F_IM.Lr-F_IM.Lm;               % rotor leakage inductance
F_IM.ks    = F_IM.Lm./F_IM.Ls;              % stator coupling factor

F_IM.tr    = F_IM.Lr./F_IM.Rr;              % rotor time constant [s]
F_IM.wslip = F_map.Iq./F_map.Id./F_IM.tr;   % slip pulsation, rotor pulsation [rad/s]

F_map.IM = F_IM;



