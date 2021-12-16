% Copyright 2021
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

function [out] = calcTnPoint(motorModel,Tref,nref)

% fdfq   = motorModel.FluxMap_dq;
% TwData = motorModel.TnSetup;

Id  = motorModel.FluxMap_dq.Id;
Iq  = motorModel.FluxMap_dq.Iq;
Fd  = motorModel.FluxMap_dq.Fd;
Fq  = motorModel.FluxMap_dq.Fq;
Tem = motorModel.FluxMap_dq.T;

p         = motorModel.data.p;
axisType  = motorModel.data.axisType;
motorType = motorModel.data.motorType;
temp0     = motorModel.data.tempCu;
Rs0       = motorModel.data.Rs;
n3phase   = motorModel.data.n3phase;
Imax      = motorModel.data.Imax;
Vdc       = motorModel.data.Vdc;
l         = motorModel.data.l;
lend      = motorModel.data.lend;

temp             = motorModel.TnSetup.temperature;
IronLossFlag     = motorModel.TnSetup.IronLossFlag;
PMLossFlag       = motorModel.TnSetup.PMLossFlag;
PMLossFactor     = motorModel.TnSetup.PMLossFactor;
IronLossFactor   = motorModel.TnSetup.IronLossFactor;
SkinEffectFlag   = motorModel.TnSetup.SkinEffectFlag;
SkinEffectMethod = motorModel.TnSetup.SkinEffectMethod;
MechLoss         = motorModel.TnSetup.MechLoss;
Control          = motorModel.TnSetup.Control;


I = Id+j*Iq;
F = Fd+j*Fq;


% 1) stator frequency evaluation
if strcmp(motorType,'IM')
    Wslip=IM.Wslip*(1+0.004*(temp-temp0));
    FreqElet=(sign(Tref)*Wslip+nref*p*pi/30)/(2*pi);
else
    FreqElet=(nref*p*pi/30)/(2*pi);
end

% 2) Iron loss evaluation
if strcmp(IronLossFlag,'Yes')
    [~,Pfes_h,Pfes_c,Pfer_h,Pfer_c,Ppm] = calcIronLoss(motorModel.IronPMLossMap_dq,motorModel.FluxMap_dq,FreqElet);
    if strcmp(PMLossFlag,'Yes')
        Ppm = Ppm*PMLossFactor;
    else
        Ppm = Ppm*0;
    end
else
    Pfes_h = zeros(size(Id));
    Pfes_c = zeros(size(Id));
    Pfer_h = zeros(size(Id));
    Pfer_c = zeros(size(Id));
    Ppm    = zeros(size(Id));
end
Pfes = IronLossFactor*(Pfes_h+Pfes_c);
Pfer = IronLossFactor*(Pfer_h+Pfer_c);
Pfe  = (Pfes+Pfer+Ppm);

% 3) Joule rotor loss (IM only)
if strcmp(motorType,'IM')
    Prot=3/2*IM.RR.*IM.Ir.^2;
    Prot=Prot*(1+0.004*(temp-temp0));
else
    Prot=zeros(size(Id));
end

% 4) Back-emf computation
if Tref>=0
    Vind=1j*2*pi*FreqElet.*F;
else
    if strcmp(motorType,'PM')
        if strcmp(axisType,'SR')
            Vind=1j*2*pi*FreqElet.*(-j.*conj(j*F));
        else
            Vind=1j*2*pi*FreqElet.*conj(F);
        end
    else
        Vind=1j*2*pi*FreqElet.*conj(F);
    end
end

% 5) Current component representing Fe and PM loss and total current
Ife=2/3*Pfe./conj(Vind);
Ife(Pfe==0)=0;

if Tref>=0
    Io=I+Ife;
else
    if strcmp(motorType,'PM')
        if strcmp(axisType,'SR') % invert the d-axis
            Io=-j*conj(j*I)+Ife;
        else    %invert q-axis
            Io=conj(I)+Ife;
        end
    else    % invert q-axis
        Io=conj(I)+Ife;
    end
end

% 6) Phase resistance computation (with temperature and skin effect)
if strcmp(SkinEffectFlag,'Yes')
    kAC = calcSkinEffect(motorModel.acLossFactor,FreqElet,temp,SkinEffectMethod);
else
    kAC = 1;
end
Rs  = Rs0.*(kAC*l/(lend+l)+lend/(lend+l)).*(1+0.004*(temp-temp0)).*ones(size(Id));

% 7) Voltage computation
Vof = Vind+Rs.*Io;  % phase voltage
Voc = Vof*sqrt(3);  % line voltage

cosfi = cos(angle(Io)-angle(Vof));

% 8) Mechanical loss computation
Pmech = polyval(MechLoss,abs(nref));

% 9) Total loss map computation
Ploss = Pfe+Prot+3/2*Rs*n3phase.*abs(Io).^2+Pmech;

% 10) Voltage and current limits
Io_m = abs(Io);
Voc_m = abs(Voc);

lim = ones(size(Id));
lim(Voc_m>Vdc) = NaN;
lim(Io_m>Imax) = NaN;

Ploss = Ploss.*lim;

% 11) Mechanical torque computation
if Tref>=0
    T = Tem-(Pmech)./(nref*p*pi/30);
else
    T = Tem+(Pmech)./(nref*p*pi/30);
end

T(isnan(T)) = Tem(isnan(T)); % avoid error for zero speed

if abs(Tref)<=max(max(T))
    % Extract (id,iq) curve @ T=cost
    c = contourc(unique(Id),unique(Iq),T,abs(Tref*[1 1]));
    idIso = c(1,2:end);
    iqIso = c(2,2:end);
    if strcmp(Control,'Maximum efficiency')
        PlossIso = interp2(Id,Iq,Ploss,idIso,iqIso);
    elseif strcmp(Control,'MTPA')
        PlossIso = interp2(Id,Iq,3/2*Rs.*Io_m.^2.*lim,idIso,iqIso);
    end
    idIso    = idIso(~isnan(PlossIso));
    iqIso    = iqIso(~isnan(PlossIso));
    PlossIso = PlossIso(~isnan(PlossIso));
    [~,index] = min(PlossIso);
    if isempty(index)
        id = NaN;
        iq = NaN;
        limIso = NaN;
    else
        id = idIso(index);
        iq = iqIso(index);
        limIso = 1;
    end
else
    limIso = NaN;
end

% update figure and matrices
if ((Tref==0)&&(nref==0))
    IdMin = 0;
    IqMin = 0;
    FdMin = +interp2(Id,Iq,Fd,IdMin,IqMin);
    FqMin = +interp2(Id,Iq,Fq,IdMin,IqMin);

    out.T     = Tref;
    out.Id    = IdMin;
    out.Iq    = IqMin;
    out.Fd    = FdMin;
    out.Fq    = FqMin;
    out.Tem   = 0;
    out.Vo    = 0;
    out.Io    = 0;
    out.Im    = IdMin+j*IqMin;
    out.Iph   = 0;
    out.Vph   = 0;
    out.PF    = NaN;
    out.P     = Tref*nref*pi/30;
    out.Ploss = 0;
    out.Pjs   = 0;
    out.PjDC  = 0;
    out.PjAC  = out.Pjs-out.PjDC;
    out.Pfe   = 0;
    out.Pfes  = 0;
    out.Pfer  = 0;
    out.Ppm   = 0;
    out.Pjr   = 0;
    out.Pmech = 0;
    out.Eo    = 0;
    out.Ife   = 0;
    if strcmp(motorType,'IM')
        out.slip  = 0;
        out.Ir    = 0;
    else
        out.slip  = NaN;
        out.Ir    = NaN;
    end
    out.Rs    = interp2(Id,Iq,Rs,id,iq);
    out.eff   = 0;
elseif ~isnan(limIso)
    if Tref<0
        if strcmp(motorType,'PM')
            if strcmp(axisType,'SR') % invert d current for generator (PM on -q)
                IdMin = -id;
                FdMin = -interp2(Id,Iq,Fd,id,iq);
                IqMin = +iq;
                FqMin = +interp2(Id,Iq,Fq,id,iq);
            else % invert q current for generator (PM on d axis)
                IdMin = +id;
                FdMin = +interp2(Id,Iq,Fd,id,iq);
                IqMin = -iq;
                FqMin = -interp2(Id,Iq,Fq,id,iq);
            end
            IqMin = iq;
            FqMin = interp2(Id,Iq,Fq,id,iq);
        else    % SyR motor: invert q current for generator
            IdMin = +id;
            FdMin = +interp2(Id,Iq,Fd,id,iq);
            IqMin = -iq;
            FqMin = -interp2(Id,Iq,Fq,id,iq);
        end
    else
        IdMin = +id;
        FdMin = +interp2(Id,Iq,Fd,id,iq);
        IqMin = +iq;
        FqMin = +interp2(Id,Iq,Fq,id,iq);
    end
    out.T     = Tref;
    out.Id    = IdMin;
    out.Iq    = IqMin;
    out.Fd    = FdMin;
    out.Fq    = FqMin;
    out.Tem   = interp2(Id,Iq,Tem,id,iq);
    out.Vo    = interp2(Id,Iq,Voc_m,id,iq);
    out.Io    = interp2(Id,Iq,Io_m,id,iq);
    out.Im    = IdMin+j*IqMin;
    out.Iph   = interp2(Id,Iq,Io,id,iq);
    out.Vph   = interp2(Id,Iq,Vof,id,iq);
    out.PF    = interp2(Id,Iq,cosfi,id,iq);
    out.P     = Tref*nref*pi/30;
    out.Ploss = interp2(Id,Iq,Ploss,id,iq);
    out.Pjs   = interp2(Id,Iq,3/2*Rs*n3phase.*Io_m.^2,id,iq);
    out.PjDC  = interp2(Id,Iq,3/2*Rs*n3phase.*Io_m.^2./(kAC*l/(lend+l)+lend/(lend+l)),id,iq);
    out.PjAC  = out.Pjs-out.PjDC;
    out.Pfe   = interp2(Id,Iq,Pfe,id,iq);
    out.Pfes  = interp2(Id,Iq,Pfes,id,iq);
    out.Pfer  = interp2(Id,Iq,Pfer,id,iq);
    out.Ppm   = interp2(Id,Iq,Ppm,id,iq);
    out.Pjr   = interp2(Id,Iq,Prot,id,iq);
    out.Pmech = Pmech;
    out.Eo    = interp2(Id,Iq,Vind,id,iq);
    out.Ife   = interp2(Id,Iq,Ife,id,iq);
    if strcmp(motorType,'IM')
        out.slip  = interp2(Id,Iq,Wslip./(2*pi*FreqElet),id,iq);
        out.Ir    = interp2(Id,Iq,IM.Ir,id,iq);
    else
        out.slip  = NaN;
        out.Ir    = NaN;
    end
    out.Rs    = interp2(Id,Iq,Rs,id,iq);
    if Tref>0
        out.eff = out.P/(out.P+out.Ploss);
    else
        out.eff = -out.P/(-out.P+out.Ploss);
    end
else
    out.T     = NaN;
    out.Id    = NaN;
    out.Iq    = NaN;
    out.Fd    = NaN;
    out.Fq    = NaN;
    out.Tem   = NaN;
    out.Vo    = NaN;
    out.Io    = NaN;
    out.Im    = NaN;
    out.Iph   = NaN;
    out.Vph   = NaN;
    out.PF    = NaN;
    out.P     = NaN;
    out.Ploss = NaN;
    out.Pjs   = NaN;
    out.PjDC  = NaN;
    out.PjAC  = NaN;
    out.Pfe   = NaN;
    out.Pfes  = NaN;
    out.Pfer  = NaN;
    out.Ppm   = NaN;
    out.Pjr   = NaN;
    out.Pmech = NaN;
    out.Eo    = NaN;
    out.Ife   = NaN;
    out.slip  = NaN;
    out.Ir    = NaN;
    out.Rs    = NaN;
    out.eff   = NaN;
end