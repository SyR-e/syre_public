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
Id  = motorModel.fdfq.Id;
Iq  = motorModel.fdfq.Iq;
Fd  = motorModel.fdfq.Fd;
Fq  = motorModel.fdfq.Fq;
T   = motorModel.fdfq.T;
IPF = sin(atan2(Iq,Id)-atan2(Fq,Fd));

axisType = motorModel.data.axisType;
nmax     = motorModel.data.nmax;
p        = motorModel.data.p;

Vmax = Vdc/sqrt(3);

id_MTPA = motorModel.AOA.MTPA.id;
iq_MTPA = motorModel.AOA.MTPA.iq;
T_MTPA  = motorModel.AOA.MTPA.T;

I_MTPA = abs(id_MTPA+j*iq_MTPA);

id_MTPV = motorModel.AOA.MTPV.id;
iq_MTPV = motorModel.AOA.MTPV.iq;
T_MTPV  = motorModel.AOA.MTPV.T;

I_MTPV = abs(id_MTPV+j*iq_MTPV);

I = abs(Id+j*Iq);
F = abs(Fd+j*Fq);

ich = min(I_MTPV);
if isempty(ich)
    ich = inf;
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

F_A = interp2(Id,Iq,F,id_A,iq_A,'cubic');   % rated flux
T_A = interp2(Id,Iq,T,id_A,iq_A,'cubic');   % rated torque

w_A = Vmax/F_A;     % rated pulsation [rad/s elt]
n_A = w_A*30/pi/p;  % rated speed [rpm]

F_B = interp2(Id,Iq,F,id_B,iq_B,'cubic');
T_B = interp2(Id,Iq,T,id_B,iq_B,'cubic');

w_B = Vmax/F_B;
n_B = w_B*30/pi/p;

F_C = interp2(Id,Iq,F,id_C,iq_C,'cubic');
T_C = interp2(Id,Iq,T,id_C,iq_C,'cubic');

w_C = Vmax/F_C;
n_C = w_C*30/pi/p;

% tratto A-B
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
        %iq_BC = polyval(p_KvMax_i,id_BC);
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
        %id_BC = polyval(p_KvMax_i,iq_BC);
     %Definisco id_BC, prima non lo era - rev.Gallo
        id_BC = interp1(iq_MTPV,id_MTPV,iq_BC,'linear','extrap');
    else
        id_BC = [id_B id_C];
        iq_BC = [iq_B iq_C];
    end
end

% PROFILO COPPIA - POT MAX
F_AB = interp2(Id,Iq,F,id_AB,iq_AB,'cubic');

if F_AB(end) > F_AB(1)
    F_AB  = fliplr(F_AB);
    id_AB = fliplr(id_AB);
    iq_AB = fliplr(iq_AB);
end
w_AB = Vmax ./ F_AB;
T_AB = interp2(Id,Iq,T,id_AB,iq_AB,'cubic');
V_AB = w_AB .* F_AB;
I_AB = abs(id_AB+j*iq_AB);

F_BC = interp2(Id,Iq,F,id_BC,iq_BC,'cubic');
w_BC = Vmax ./ F_BC;
T_BC = interp2(Id,Iq,T,id_BC,iq_BC,'cubic');
V_BC = w_BC .* F_BC;
I_BC = abs(id_BC+j*iq_BC);

% low speed values (w < w1)
w_0A = linspace(0,w_A,20);
T_0A = ones(size(w_0A)) * T_A;
V_0A = F_A * w_0A;
I_0A = ones(size(w_0A)) * Imax;

% limiti iq
id_max = [id_A*ones(size(w_0A)) id_AB id_BC];
iq_max = [iq_A*ones(size(w_0A)) iq_AB iq_BC];

F_max = [F_A*ones(size(w_0A)) F_AB F_BC];
iq_min = zeros(size(iq_max));

% a pieno carico
fd_AB = interp2(Id,Iq,Fd,id_AB,iq_AB,'cubic');
fd_BC = interp2(Id,Iq,Fd,id_BC,iq_BC,'cubic');
fq_AB = interp2(Id,Iq,Fq,id_AB,iq_AB,'cubic');
fq_BC = interp2(Id,Iq,Fq,id_BC,iq_BC,'cubic');
% %  a vuoto
% fq_AB_0 = interp2(id,iq,Fq,id_AB,0,'cubic');
% fq_BC_0 = interp2(id,iq,Fq,id_BC,0,'cubic');
% % fq_AB_0 = interp2(id,iq,Fd,id_AB,0);
% % fq_BC_0 = interp2(id,iq,Fd,id_BC,0);
% a vuoto

fd_max = [fd_AB(1)*ones(size(w_0A)) fd_AB fd_BC];
fq_max = [fq_AB(1)*ones(size(w_0A)) fq_AB fq_BC];

% fq_0 = [fq_AB_0(1)*ones(size(w_0A)) fq_AB_0 fq_BC_0];

% F0 = abs(fd_max + j*fq_0);
% temp = abs(fd_max + j*fq_max);
% iq_min(lm > F_max) = iq_max(lm > F_max);

% mechanical speed evaluation (IM only)
% synchronous speed
if exist('Wslip','var')
    wslip = interp2(Id,Iq,Wslip,id_max,iq_max,'cubic');
    [a,last_number] = find(not(isnan(wslip)),1,'last');
    wslip(isnan(wslip)) = wslip(last_number);
    %     rot_temperature = 20;
    %     ref_temperature = 100;
    %temp_coeff = (234.5 + rot_temperature)/(234.5 + Rr_temp);
    temp_coeff=1;
    wslip = wslip * temp_coeff;
else
    wslip = id_max * 0;
end

w = [w_0A w_AB w_BC];
wr = w - wslip;

% potenza meccanica corretta (a parte Pfe)
P = [T_0A T_AB T_BC] .* wr/p;
% calcolo sbagliato, manca wslip .. lascio per compatibilità con versioni prec
if ~exist('Wslip','var')
    P_AB = T_AB .* w_AB/p;
    P_BC = T_BC .* w_BC/p;
    P_0A = T_0A .* w_0A/p;
    % Tmax = T_A;
    % Pmax = max([P_AB P_BC])
    
%     PFlim = P ./ (3/2 * n3phase * [V_0A V_AB V_BC] .* [I_0A I_AB I_BC]); %AS
    PF = interp2(Id,Iq,IPF,id_max,iq_max);
else
    w_0A=w_0A-wslip(1:length(w_0A));
    w_AB=w_AB-wslip(length(w_0A)+1:length([w_0A w_AB]));
    w_BC=w_BC-wslip(length([w_0A w_AB])+1:end);
    
    P_AB = T_AB .* w_AB/p;
    P_BC = T_BC .* w_BC/p;
    P_0A = T_0A .* w_0A/p;
    
    PF = interp2(Id,Iq,IPF,id_max,iq_max);
end

% save Plim points
Plim.wr     = wr;
Plim.n      = wr * 30/pi/p;
Plim.w      = w;
Plim.P      = [P_0A P_AB P_BC];
Plim.V      = [V_0A V_AB V_BC];
Plim.I      = [I_0A I_AB I_BC];
Plim.T      = [T_0A T_AB T_BC];
Plim.F      = F_max;
Plim.fd_max = fd_max;
Plim.fq_max = fq_max;
Plim.id_max = id_max;
Plim.iq_max = iq_max;
Plim.PF     = PF;

Plim.id_A   = id_A;
Plim.iq_A   = iq_A;
Plim.F_A    = F_A;
Plim.T_A    = T_A;
Plim.n_A    = n_A;
Plim.id_B   = id_B;
Plim.iq_B   = iq_B;
Plim.F_B    = F_B;
Plim.T_B    = T_B;
Plim.n_B    = n_B;
Plim.id_C   = id_C;
Plim.iq_C   = iq_C;
Plim.F_C    = F_C;
Plim.T_C    = T_C;
Plim.n_C    = n_C;

Plim.id_AB  = id_AB;
Plim.iq_AB  = iq_AB;
Plim.id_BC  = id_BC;
Plim.iq_BC  = iq_BC;

Plim.w_BC = w_BC;
Plim.I_BC = I_BC;
Plim.V_BC = V_BC;
Plim.P_BC = P_BC;
Plim.T_BC = T_BC;








