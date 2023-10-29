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

function [SOL] = simulate_FOC_IM(geo,per,mat,eval_type,pathname,filename)


p       = geo.p;
Nbars   = geo.IM.Nbars;
Nbob    = geo.win.Nbob;
offset  = geo.IM.offset;
th      = geo.IM.th;
thR     = geo.IM.thR;
kturns  = geo.IM.kturns;
gamma   = per.gamma;
n3phase = geo.win.n3phase;
ps      = geo.ps;
NbarSim = Nbars/(2*p)*ps;

if isfield(per,'flag3phaseSet')
    flag3phSet = per.flag3phaseSet;
else
    flag3phSet = ones(1,n3phase);
end

% evaluation of the stator current (dq) and stator current FEMM update
iAmp     = per.overload*per.i0;
iAmpCoil = iAmp*Nbob;
id       = iAmpCoil*cos(gamma*pi/180);
iq       = iAmpCoil*sin(gamma*pi/180);



openfemm(1)
opendocument([pathname filename])

for ik=0:(n3phase-1)
    i123 = dq2abc(id,iq,th(ik+1)*pi/180);      % each 3phase set has its own offset angle
    i_tmp((3*ik)+1,1) = (i123(1))*flag3phSet(ik+1);
    i_tmp((3*ik)+2,1) = (i123(2))*flag3phSet(ik+1);
    i_tmp((3*ik)+3,1) = (i123(3))*flag3phSet(ik+1);

    phase_name{3*ik+1} = strcat('fase',num2str(3*ik+1));
    phase_name{3*ik+2} = strcat('fase',num2str(3*ik+2));
    phase_name{3*ik+3} = strcat('fase',num2str(3*ik+3));
    mi_modifycircprop(phase_name{3*ik+1}, 1,i_tmp((3*ik)+1,1));
    mi_modifycircprop(phase_name{3*ik+2}, 1,i_tmp((3*ik)+2,1));
    mi_modifycircprop(phase_name{3*ik+3}, 1,i_tmp((3*ik)+3,1));
end


% i123 = dq2abc(id,iq,th*pi/180);
% Ni1  = i123(1);
% Ni2  = i123(2);
% Ni3  = i123(3);

if iq~= 0
    maxIter = 10; %% 10 is nosense .. maximum number of attempts is 3
else
    maxIter = 1;
end

% % define phase names
% phase_name = cell(n3phase*3,1);
% for jj=0:(n3phase-1)
%     phase_name{3*jj+1}=strcat('fase',num2str(3*jj+1));
%     phase_name{3*jj+2}=strcat('fase',num2str(3*jj+2));
%     phase_name{3*jj+3}=strcat('fase',num2str(3*jj+3));
% end

% modify rotor bar temperature
sigmaBar = mat.BarCond.sigma*1/(1+mat.BarCond.alpha*(per.tempPP-20));
mi_modifymaterial(mat.BarCond.MatName,5,sigmaBar);

done=0;
ii=1;
while ~done
    switch ii
        case 1
            iRA = -0.95 * iq/Nbob;
            iR = iRA;
        case 2
            iRB = -1.05 * iq/Nbob;
            iR = iRB;
        otherwise
            fqB;
            iR = iRA - fqA * (iRB - iRA)/(fqB - fqA);
    end
    
    % rotor bar currents
    ibar = kturns(1) * dq2bar(0,iR,(th(1)-thR-offset(1))*pi/180,Nbars/p);

    % rotor current in FEMM
    
    for jj  = 1:NbarSim
        circ_name = ['bar' num2str(jj)];
        mi_modifycircprop(circ_name,1,ibar(jj));
    end

    mi_analyze(1);
    mi_loadsolution;
    %     out = post_proc(geo,0);

    % evaluate stator flux
    temp_out = mo_getcircuitproperties('fase1');
    f1 = temp_out(3)*2*p/ps;
    temp_out = mo_getcircuitproperties('fase2');
    f2 = temp_out(3)*2*p/ps;
    temp_out = mo_getcircuitproperties('fase3');
    f3 = temp_out(3)*2*p/ps;
    fdq = abc2dq(f1,f2,f3,th*pi/180);

    % evaluate rotor flux
    Fbar = zeros(1,NbarSim);
    Ibar = zeros(1,NbarSim);
    Vbar = zeros(1,NbarSim);
    for jj  = 1:NbarSim
        circ_name = ['bar' num2str(jj)];
        temp_out = mo_getcircuitproperties(circ_name); % [current,voltage,flux]
        Fbar(jj) = temp_out(3);
        Ibar(jj) = temp_out(1);
        Vbar(jj) = temp_out(2);
    end
    % aggiunta flussi di anello di corto circuito (end ring)
    FbarTot=Fbar+geo.IM.k^2*(2*per.IM.Lring)*Ibar;
    if rem(ps,2)==1
        FbarTot = [FbarTot -FbarTot];
    end
    % flussi dq rotore (calcolo diretto da barre)
    fdqR = kturns(1)*Nbars/3*bar2dq(FbarTot',thR*pi/180,Nbars/p);

    if iq==0
        kr=0;
    else
        kr = -iR/(iq/Nbob);
    end


    switch ii %find correct iR by trial and error
        case 1
            fqA = fdqR(2);
        case 2
            fqB = fdqR(2);
        otherwise
            if abs(fqB) < abs(fqA)
                iRA = iRB;
                fqA = fqB;
            end
            iRB = iR;
            fqB = fdqR(2);
    end

    if ii==maxIter
        done=1;
    end

    % stop if fqR = 0
    if abs(fdqR(2))<(0.01*abs(fdq(2)))
        done=1;
    end

    ii = ii+1;
end


T = -mo_gapintegral('AGap',0);

mo_groupselectblock();
we = mo_blockintegral(2)*2*p/ps;
wc = mo_blockintegral(17)*2*p/ps;
mo_clearblock();

mo_close, mi_close
closefemm



%sol = [Id,Iq,fdq(1),fdq(2),iR,fdqR(1),fdqR(2),T,kr,ii];
SOL.th = 0;
SOL.id = id/Nbob;
SOL.iq = iq/Nbob;
SOL.fd = fdq(1)*Nbob;
SOL.fq = fdq(2)*Nbob;
SOL.T  = T;
SOL.we = we;
SOL.wc = wc;

SOL.IM.ir      = iR;
SOL.IM.fdr     = fdqR(1);
SOL.IM.fqr     = fdqR(2);
SOL.IM.kr      = kr;
SOL.IM.Ibar    = Ibar;
SOL.IM.Vbar    = Vbar;
SOL.IM.Fbar    = Fbar;
SOL.IM.FbarTot = FbarTot;



% cd(currentDir);