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

axisType  = motorModel.data.axisType;
motorType = motorModel.data.motorType;
nmax      = motorModel.data.nmax;
p         = motorModel.data.p;
Rs        = motorModel.data.Rs;

Vmax = Vdc/sqrt(3);

if strcmp(motorType,'IM')
    IM = motorModel.FluxMap_dq.IM;
end

% id_MTPA = motorModel.controlTrajectories.MTPA.id;
% iq_MTPA = motorModel.controlTrajectories.MTPA.iq;
% T_MTPA  = motorModel.controlTrajectories.MTPA.T;

MTPA = motorModel.controlTrajectories.MTPA;
MTPV = motorModel.controlTrajectories.MTPV;

% I_MTPA = abs(id_MTPA+1i*iq_MTPA);

% id_MTPV = motorModel.controlTrajectories.MTPV.id;
% iq_MTPV = motorModel.controlTrajectories.MTPV.iq;
% T_MTPV  = motorModel.controlTrajectories.MTPV.T;

% I_MTPV = abs(id_MTPV+1i*iq_MTPV);

% I = abs(Id+1i*Iq);
% F = abs(Fd+1i*Fq);

if strcmp(motorType,'PM')
    ich = min(abs(MTPV.id+j*MTPV.iq));
    if isempty(ich)
        ich = inf;
    end
else
    ich=0;
end

% a = (Fd.^2+Fq.^2);
% b = 2*Rs*(Fd.*Iq-Fq.*Id);
% c = Rs^2*(Id.^2+Iq.^2)-Vmax^2;
% wLim = (-b+(b.^2-4.*a.*c).^0.5)./(2*a); % electric pulsation @ Vmax [rad/s]

wLim = calcLimitPulsation(Id,Iq,Fd,Fq,Rs,Vmax); % electric pulsation @ Vmax [rad/s]
wLim(isnan(wLim))=max(wLim,[],'all')*1.1;



% Important points:
% A --> MTPA @ I=imax
% B --> start of MTPV @ I=Imax
% C --> characteristic current
% M --> max speed point

id_A = interp1(abs(MTPA.id+j*MTPA.iq),MTPA.id,Imax);
iq_A = interp1(abs(MTPA.id+j*MTPA.iq),MTPA.iq,Imax);

if ich>Imax
    if strcmp(axisType,'SR')
        % no MTPV, or no MTPV and Imax crossing
        id_B = 0;
        iq_B = Imax;
    else
        id_B = -Imax;
        iq_B = 0;
    end
    id_C = [];
    iq_C = [];
else
    id_B = interp1(abs(MTPV.id+j*MTPV.iq),MTPV.id,Imax);
    iq_B = interp1(abs(MTPV.id+j*MTPV.iq),MTPV.iq,Imax);
    if strcmp(axisType,'SR')
        id_C = 0;
        iq_C = ich;
    else
        id_C = -ich;
        iq_C = 0;
    end
end

% other quantities
T_A  = interp2(Id,Iq,T,id_A,iq_A);
fd_A = interp2(Id,Iq,Fd,id_A,iq_A);
fq_A = interp2(Id,Iq,Fq,id_A,iq_A);
F_A  = interp2(Id,Iq,abs(Fd+j*Fq),id_A,iq_A);
%w_A  = interp2(Id,Iq,wLim,id_A,iq_A);
w_A = calcLimitPulsation(id_A,iq_A,fd_A,fq_A,Rs,Vmax);

T_B  = interp2(Id,Iq,T,id_B,iq_B);
fd_B = interp2(Id,Iq,Fd,id_B,iq_B);
fq_B = interp2(Id,Iq,Fq,id_B,iq_B);
F_B  = interp2(Id,Iq,abs(Fd+j*Fq),id_B,iq_B);
%w_B  = interp2(Id,Iq,wLim,id_B,iq_B);
w_B = calcLimitPulsation(id_B,iq_B,fd_B,fq_B,Rs,Vmax);

T_C  = interp2(Id,Iq,T,id_C,iq_C);
fd_C = interp2(Id,Iq,Fd,id_C,iq_C);
fq_C = interp2(Id,Iq,Fq,id_C,iq_C);
F_C  = interp2(Id,Iq,abs(Fd+j*Fq),id_C,iq_C);
%w_C  = interp2(Id,Iq,wLim,id_C,iq_C);
w_C = calcLimitPulsation(id_C,iq_C,fd_C,fq_C,Rs,Vmax);



if strcmp(motorType,'IM')
    wslip_A = interp2(Id,Iq,IM.wslip,id_A,iq_A,'cubic');
    wslip_B = interp2(Id,Iq,IM.wslip,id_B,iq_B,'cubic');
    wslip_C = interp2(Id,Iq,IM.wslip,id_C,iq_C,'cubic');
    wrLim = (wLim-IM.wslip)/p;
else
    wslip_A = 0;
    wslip_B = 0;
    wslip_C = 0;
    wrLim = wLim/p;
end

wr_A = (w_A-wslip_A)/p;
wr_B = (w_B-wslip_B)/p;
wr_C = (w_C-wslip_C)/p;


nLim  = wrLim*30/pi;

n_A = wr_A*30/pi;
n_B = wr_B*30/pi;
n_C = wr_C*30/pi;

%% Profiles

% A-B segment
gamma_A = atan2(iq_A,id_A);
gamma_B = atan2(iq_B,id_B);
fi = linspace(gamma_A,gamma_B,201);
id_AB = Imax*cos(fi);
iq_AB = Imax*sin(fi);

if strcmp(axisType,'SR')
    % B --> C
    if (ich<=Imax)
        id_BC = linspace(id_B,id_C,101);
        iq_BC = interp1(MTPV.id,MTPV.iq,id_BC,'linear','extrap');
    else
        id_BC = [id_B id_B];
        iq_BC = [iq_B iq_B];
    end
else
    % SPM style axes
    % B --> C
    if (ich<=Imax)
        iq_BC = linspace(iq_B,iq_C,101);
        id_BC = interp1(MTPV.iq,MTPV.id,iq_BC,'linear','extrap');
    else
        id_BC = [id_B id_B];
        iq_BC = [iq_B iq_B];
    end
end

% Profiles
fd_AB = interp2(Id,Iq,Fd,id_AB,iq_AB);
fq_AB = interp2(Id,Iq,Fq,id_AB,iq_AB);
w_AB  = calcLimitPulsation(id_AB,iq_AB,fd_AB,fq_AB,Rs,Vmax);

fd_BC = interp2(Id,Iq,Fd,id_BC,iq_BC);
fq_BC = interp2(Id,Iq,Fq,id_BC,iq_BC);
w_BC  = calcLimitPulsation(id_BC,iq_BC,fd_BC,fq_BC,Rs,Vmax);

% n_AB = interp2(Id,Iq,nLim,id_AB,iq_AB,'linear');
% n_BC = interp2(Id,Iq,nLim,id_BC,iq_BC,'linear');

id = [id_A id_A id_AB(2:end-1) id_B id_BC(2:end-1) id_C];
iq = [iq_A iq_A iq_AB(2:end-1) iq_B iq_BC(2:end-1) iq_C];
w  = [0    w_A  w_AB(2:end-1)  w_B  w_BC(2:end-1)  w_C];
if isnan(w(end))
    w(end) = w(end-1)*1.1;
end


if strcmp(motorType,'IM')
    wSlip = interp2(Id,Iq,IM.wslip,id,iq,'cubic');
    wr = (w-wSlip)/p;
else
    wr = w/p;
end

n = wr*30/pi;

% filt NaN in speed (problem with IM, if the axis are not included in the map)

indexFilt = ~isnan(n);
id = id(indexFilt);
iq = iq(indexFilt);
w  = w(indexFilt);
n  = n(indexFilt);

if nmax<n(end)
    id_M = interp1(n,id,nmax);
    iq_M = interp1(n,iq,nmax);
    n_M  = nmax;
else
    id_M = id(end);
    iq_M = iq(end);
    n_M  = n(end);
end

% add the M point, if needed
flagAdd=1;
if n_A==n_M
    flagAdd=0;
end
if n_B==n_M
    flagAdd=0;
end
if ~isempty(n_C)
    if ~isnan(n_C)
        if n_C==n_M
            flagAdd=0;
        end
    end
end

% if((n_M~=n_A)&&(n_M~=n_B)&&(n_M~=n_C))
%     index = length(n(n<n_M));
%     n  = [n(1:index) n_M n(index+1:end)];
%     id = [id(1:index) id_M id(index+1:end)];
%     iq = [iq(1:index) iq_M iq(index+1:end)];
% end

if flagAdd
    index = length(n(n<n_M));
    n  = [n(1:index) n_M n(index+1:end)];
    id = [id(1:index) id_M id(index+1:end)];
    iq = [iq(1:index) iq_M iq(index+1:end)];
end


T_M  = interp2(Id,Iq,T,id_M,iq_M);
fd_M = interp2(Id,Iq,Fd,id_M,iq_M);
fq_M = interp2(Id,Iq,Fq,id_M,iq_M);
F_M  = interp2(Id,Iq,abs(Fd+j*Fq),id_M,iq_M);
w_M  = interp2(Id,Iq,wLim,id_M,iq_M);

index = n<=n_M;
n = n(index);
id = id(index);
iq = iq(index);


% other quantities profiles
fd = interp2(Id,Iq,Fd,id,iq);
fq = interp2(Id,Iq,Fq,id,iq);
F  = abs(fd+j*fq);
I  = abs(id+j*iq);
T  = interp2(Id,Iq,T,id,iq);
wr = n*pi/30;

if strcmp(motorType,'IM')
    wSlip = interp2(Id,Iq,IM.wslip,id,iq);
    w = p*wr+wSlip;
else
    w = p*wr;
end

ed  = -w.*fq;
eq  = +w.*fd;
E   = abs(ed+j*eq);  % back-emf [Vpk]
vd  = Rs*id+ed;
vq  = Rs*iq+eq;
V   = abs(vd+j*vq);  % phase voltage [Vpk]
IPF = sin(atan2(iq,id)-atan2(fq,fd));
PF  = cos(atan2(vq,vd)-atan2(iq,id));

P = T.*wr;

% save Plim points
Plim.wr  = wr;
Plim.n   = n;
Plim.w   = w;
Plim.P   = P;
Plim.V   = V;
Plim.E   = E;
Plim.I   = I;
Plim.T   = T;
Plim.F   = F;
Plim.fd  = fd;
Plim.fq  = fq;
Plim.id  = id;
Plim.iq  = iq;
Plim.PF  = PF;
Plim.IPF = IPF;
Plim.ed  = ed;
Plim.eq  = eq;
Plim.vd  = vd;
Plim.vq  = vq;


Plim.id_A = id_A;
Plim.iq_A = iq_A;
Plim.T_A  = T_A;
Plim.n_A  = n_A;

Plim.id_B = id_B;
Plim.iq_B = iq_B;
Plim.T_B  = T_B;
Plim.n_B  = n_B;

Plim.id_C = id_C;
Plim.iq_C = iq_C;
Plim.T_C  = T_C;
Plim.n_C  = n_C;

Plim.id_M = id_M;
Plim.iq_M = iq_M;
Plim.T_M  = T_M;
Plim.n_M  = n_M;


