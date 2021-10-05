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

function [Plim] = OpLimEval(motorModel,Imax,Vdc)

% load data
Id  = motorModel.FluxMap_dq.Id;
Iq  = motorModel.FluxMap_dq.Iq;
Fd  = motorModel.FluxMap_dq.Fd;
Fq  = motorModel.FluxMap_dq.Fq;
T   = motorModel.FluxMap_dq.T;
IPF = sin(atan2(Iq,Id)-atan2(Fq,Fd));

axisType  = motorModel.data.axisType;
motorType = motorModel.data.motorType;
nmax      = motorModel.data.nmax;
p         = motorModel.data.p;
Rs        = motorModel.data.Rs;

Vmax = Vdc/sqrt(3);

id_MTPA = motorModel.controlTrajectories.MTPA.id;
iq_MTPA = motorModel.controlTrajectories.MTPA.iq;
T_MTPA  = motorModel.controlTrajectories.MTPA.T;

I_MTPA = abs(id_MTPA+1i*iq_MTPA);

id_MTPV = motorModel.controlTrajectories.MTPV.id;
iq_MTPV = motorModel.controlTrajectories.MTPV.iq;
T_MTPV  = motorModel.controlTrajectories.MTPV.T;

I_MTPV = abs(id_MTPV+1i*iq_MTPV);

I = abs(Id+1i*Iq);
F = abs(Fd+1i*Fq);

if strcmp(motorType,'PM')
    ich = min(I_MTPV);
    if isempty(ich)
        ich = inf;
    end
else
    ich=0;
end

% Important points:
% A --> MTPA @ I=imax
% B --> start of MTPV @ I=Imax
% C --> characteristic current

id_A = interp1(I_MTPA,id_MTPA,Imax);
iq_A = interp1(I_MTPA,iq_MTPA,Imax);

if ich>Imax
    if strcmp(axisType,'SR')
        % no MTPV, or no MTPV and Imax crossing
        id_B = 0;
        iq_B = Imax;
    else
        id_B = -Imax;
        iq_B = 0;
    end
    id_C = id_B;
    iq_C = iq_B;
else
    id_B = interp1(I_MTPV,id_MTPV,Imax);
    iq_B = interp1(I_MTPV,iq_MTPV,Imax);
    if strcmp(axisType,'SR')
        id_C = 0;
        iq_C = Imax;
    else
        id_C = -Imax;
        iq_C = 0;
    end
end

% A-B segment
c = contourc(unique(Id),unique(Iq),I,[Imax Imax]);
id_Imax = (c(1,2:end));
iq_Imax = (c(2,2:end));

if strcmp(axisType,'SR')
    % A --> B
    id_AB = id_Imax((iq_Imax >= iq_A) & (iq_Imax <= iq_B));
    iq_AB = iq_Imax((iq_Imax >= iq_A) & (iq_Imax <= iq_B)); %avoid interp errors
    % B --> C
    if (ich<=Imax)
        id_BC = linspace(id_B,id_C,50);
        id_BC = id_BC(1:end-1);
        iq_BC = interp1(id_MTPV,iq_MTPV,id_BC,'linear','extrap');
    else
        id_BC = [id_B id_C];
        iq_BC = [iq_B iq_C];
    end
else
    % SPM style axes
    id_AB = id_Imax((iq_Imax >= iq_B) & (iq_Imax <= iq_A));
    iq_AB = iq_Imax((iq_Imax >= iq_B) & (iq_Imax <= iq_A)); %avoid interp errors
    % tratto B --> C
    if (ich<=Imax)
        iq_BC = linspace(iq_B,iq_C,50);
        iq_BC = iq_BC(1:end-1);
        id_BC = interp1(iq_MTPV,id_MTPV,iq_BC,'linear','extrap');
    else
        id_BC = [id_B id_C];
        iq_BC = [iq_B iq_C];
    end
end

% PROFILO COPPIA - POT MAX
fd_AB = interp2(Id,Iq,Fd,id_AB,iq_AB,'cubic');
fq_AB = interp2(Id,Iq,Fq,id_AB,iq_AB,'cubic');
F_AB = fd_AB + 1i*fq_AB;
if abs(F_AB(end)) > abs(F_AB(1))
    F_AB  = fliplr(F_AB);
    id_AB = fliplr(id_AB);
    iq_AB = fliplr(iq_AB);
    fd_AB = fliplr(fd_AB);
    fq_AB = fliplr(fq_AB);
end
fd_BC = interp2(Id,Iq,Fd,id_BC,iq_BC,'cubic');
fq_BC = interp2(Id,Iq,Fq,id_BC,iq_BC,'cubic');
F_BC = fd_BC + 1i*fq_BC;

T_AB = interp2(Id,Iq,T,id_AB,iq_AB,'cubic');
I_AB = (id_AB+1i*iq_AB);
% voltage equation (Vdq)max = R*idq + j*w*fdq --> w
b = Rs*(fd_AB.*iq_AB-fq_AB.*id_AB)./abs(F_AB).^2;
c = (Rs^2*abs(I_AB).^2-Vmax^2)./abs(F_AB).^2;
w_AB = -b + sqrt(b.^2-c);
V_AB = Rs*I_AB + 1i.*w_AB.*F_AB;

T_BC = interp2(Id,Iq,T,id_BC,iq_BC,'cubic');
I_BC = (id_BC+1i*iq_BC);
% voltage equation (Vdq)max = R*idq + j*w*fdq --> w
b = Rs*(fd_BC.*iq_BC-fq_BC.*id_BC)./(fd_BC.^2+fq_BC.^2);
c = (Rs^2*(id_BC.^2+iq_BC.^2)-Vmax^2)./(fd_BC.^2+fq_BC.^2);
w_BC = -b + sqrt(b.^2-c);
V_BC = Rs*I_BC + 1i.*w_BC.*F_BC;

w_A = w_AB(1);     % rated pulsation [rad/s elt]
w_B = w_AB(end);     
w_C = w_BC(end);     

n_A = w_A*30/pi/p;  % rated speed [rpm]
n_B = w_B*30/pi/p;  
n_C = w_C*30/pi/p;  

% low speed values (w < w1, MTPA)
w_0A = linspace(0,w_A,20);
w_0A = w_0A(1:end-1);

T_A = interp2(Id,Iq,T,id_A,iq_A,'cubic');   % rated torque
T_B = interp2(Id,Iq,T,id_B,iq_B,'cubic');
T_C = interp2(Id,Iq,T,id_C,iq_C,'cubic');

T_0A = ones(size(w_0A)) * T_A;
V_0A = (Rs*(id_A + 1i*iq_A) + 1i.*w_0A.*(fd_AB(1) + 1i*fq_AB(1)));

% merge 0A - AB - BC
id_max = [id_A*ones(size(w_0A)) id_AB id_BC];
iq_max = [iq_A*ones(size(w_0A)) iq_AB iq_BC];
I = id_max + 1i*iq_max;
fd_max = [fd_AB(1)*ones(size(w_0A)) fd_AB fd_BC];
fq_max = [fq_AB(1)*ones(size(w_0A)) fq_AB fq_BC];
F = fd_max + 1i*fq_max;

% mechanical speed evaluation (IM only)
% synchronous speed
if exist('Wslip','var')
    wslip = interp2(Id,Iq,Wslip,id_max,iq_max,'cubic');
    [~,last_number] = find(not(isnan(wslip)),1,'last');
    wslip(isnan(wslip)) = wslip(last_number);
    temp_coeff=1;
    wslip = wslip * temp_coeff;
else
    wslip = id_max * 0;
end

w = [w_0A w_AB w_BC];
wr = (w - wslip)/p;

T = [T_0A T_AB T_BC];
P = T.*wr;
V = [V_0A V_AB V_BC];
PF = cos(angle(V) - angle(I));

% save Plim points
Plim.wr     = wr;
Plim.n      = wr * 30/pi;
Plim.w      = w;
Plim.P      = P;
Plim.V      = abs(V);
Plim.I      = abs(I);
Plim.T      = T;
Plim.F      = abs(F);
Plim.fd_max = fd_max;
Plim.fq_max = fq_max;
Plim.id_max = id_max;
Plim.iq_max = iq_max;
Plim.PF     = PF;

Plim.id_A   = id_A;
Plim.iq_A   = iq_A;
Plim.F_A    = abs(F_AB(1));
Plim.T_A    = T_A;
Plim.n_A    = n_A;
Plim.id_B   = id_B;
Plim.iq_B   = iq_B;
Plim.F_B    = abs(F_AB(end));
Plim.T_B    = T_B;
Plim.n_B    = n_B;
Plim.id_C   = id_C;
Plim.iq_C   = iq_C;
Plim.F_C    = abs(F_BC(end));
Plim.T_C    = T_C;
Plim.n_C    = n_C;

% Plim.id_AB  = id_AB;
% Plim.iq_AB  = iq_AB;
% Plim.id_BC  = id_BC;
% Plim.iq_BC  = iq_BC;

% Plim.w_BC = w_BC;
% Plim.I_BC = I_BC;
% Plim.V_BC = V_BC;
% Plim.P_BC = P_BC;
% Plim.T_BC = T_BC;








