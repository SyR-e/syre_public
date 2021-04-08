% Copyright 2018
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [SOL] = simulate_xdeg_MN(geo,per,eval_type,pathname,filename)

h = OpenMagnet(1);

switch eval_type
    case 'singt'
        nsim = per.nsim_singt;
        xdeg = per.delta_sim_singt;
        gamma = per.gamma;
    case 'singm'
        nsim = per.nsim_singt;
        xdeg = per.delta_sim_singt;
        gamma = per.gamma;
    otherwise
        error('MagNet Simulations not available for the selected evaluation type!');
end

th0 = geo.th0;
p   = geo.p;
r  = geo.r;
gap = geo.g;
pc  = 360/(6*geo.q*p)/2;
ps  = geo.ps;
n3ph = geo.win.n3phase; %AS number of 3-phase circuits
Q = 6*geo.q*geo.p*geo.win.n3phase;
skew_angle = 0;
N_parallel = 1;
N_cond = geo.win.Ns/geo.p/(geo.q)/size(geo.win.avv,1);
Br = per.BrPP;
tempPP = per.tempPP;
q = geo.q;
gamma_ = gamma*pi/180;

% evaluation of the phase current values for all positions to be simulated
iAmp = per.overload*per.i0;
% iAmp = per.overload*calc_io(geo,per);

id = iAmp * cos(gamma * pi/180);
iq = iAmp * sin(gamma * pi/180);
Imod = abs(id + 1i* iq);       % Modulo corrente (A)

DocSV = invoke(h.magnetHandler, 'openDocument',[pathname,filename(1:end-4) '.mn']);
Doc = invoke(h.magnetHandler, 'getDocument');
View = invoke(Doc, 'getCurrentView');
Solution = invoke(Doc, 'getSolution');
calculator = invoke (h.magnetHandler, 'getStackCalculator');

% set the operating parameters
if (Br==0)
    n_mag_simulati=0;
else
    %n_mag_simulati=size(geo.BLKLABELS.rotore.BarName,1);
    tmp=geo.BLKLABELS.rotore.xy(:,3);
    n_mag_simulati=length(tmp(tmp==6));
end

% type of solution (order)
% Polynomial Order = 1 (3 field points per triangle)
% Command = ['Call getDocument().setParameter("", "PolynomialOrder", "1", infoNumberParameter)'];
% invoke(MN6, 'processCommand', Command);

% SourcesOnAtTransientStart: smooth field at t = 0
Command = ['Call getDocument().setParameter("", "SourcesOnAtTransientStart", "Yes", infoStringParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

% 09 03 2010 - MeshAllowUnconstrainedHoles
Command = ['Call getDocument().setParameter("", "MeshAllowUnconstrainedHoles", "Yes", infoStringParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

% Model temperature
Command = ['Call getDocument().setParameter("", "Temperature", "' num2str(tempPP) '", infoNumberParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

% set transient options (time0 etc..)
time_end = 1000 * xdeg/(360*p)*60/per.EvalSpeed;   % ms
invoke(Doc,'setFixedIntervalTimeSteps',0, time_end, time_end/nsim);

% Iron loss time interval - global
Command = ['Call getDocument().setParameter("", "TransientAveragePowerLossStartTime", "0", infoNumberParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

Command = ['Call getDocument().setParameter("", "TransientAveragePowerLossEndTime", "' num2str(time_end) ' %ms", infoNumberParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

% Iron loss time interval - rotor
% time_start_pferot = time_end - time_end/Q;
%% temp -- macchina vagati
% time_start_pferot = time_end - q/Q*time_end;
time_start_pferot = 0;

Command = ['Call getDocument().setParameter("rotor", "TransientAveragePowerLossStartTime", "' num2str(time_start_pferot) ' %ms", infoNumberParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

Command = ['Call getDocument().setParameter("rotor", "TransientAveragePowerLossEndTime", "' num2str(time_end) ' %ms", infoNumberParameter)'];
invoke(h.magnetHandler, 'processCommand', Command);

Command = ['Call getDocument().setMotionSourceType("Motion#1", infoVelocityDriven)'];
invoke(h.magnetHandler, 'processCommand', Command);
%%

% Motion#1: transient angle Vs time
% start time = 0
Command = 'REDIM ArrayOfValues1(1)';
invoke(h.magnetHandler, 'processCommand', Command);
Command = 'ArrayOfValues1(0)= 0';
invoke(h.magnetHandler, 'processCommand', Command);
% end time == sim_angle (elt deg)
Command = ['ArrayOfValues1(1)= ' num2str(time_end)];
invoke(h.magnetHandler, 'processCommand', Command);

% start angle = Cas.skew_angle
Command = 'REDIM ArrayOfValues2(1)';
invoke(h.magnetHandler, 'processCommand', Command);
Command = ['ArrayOfValues2(0)= ' num2str(skew_angle)]; % Skew_New
invoke(h.magnetHandler, 'processCommand', Command);
% end angle = start + sim_angle
Command = ['ArrayOfValues2(1)= ' num2str(xdeg/p + skew_angle)];   %sim_angle elettrici % Skew_New
invoke(h.magnetHandler, 'processCommand', Command);
Command = 'Call getDocument().setMotionPositionVsTime("Motion#1", ArrayOfValues1, ArrayOfValues2)';
invoke(h.magnetHandler, 'processCommand', Command);

%%
phase_index = {'U','V','W'};
if exist('th0')
    % casi successivi feb 2010 (+90 perch? Waveform in Magnet\coil ? un sin)
    phase_angle = [0 -120 120] + (gamma_*180/pi + th0(1) + 90);
else
    % casi precedenti feb 2010
    phase_angle = [0 -120 120] - 180/pi*(pi-gamma_);
end

for nn=1:n3ph
    for ii = 1:3
        Command=['Call getDocument().setCoilSourceType("',phase_index{ii},int2str(nn),'", infoCurrentDriven)'];
        invoke(h.magnetHandler, 'processCommand', Command);
        
        Command = 'REDIM ArrayOfValues(5)';
        invoke(h.magnetHandler, 'processCommand', Command);
        Command = 'ArrayOfValues(0)= 0';
        invoke(h.magnetHandler, 'processCommand', Command);
        Command = ['ArrayOfValues(1)= ' num2str(Imod/N_parallel)];
        invoke(h.magnetHandler, 'processCommand', Command);
        Command = ['ArrayOfValues(2)= ' num2str(per.EvalSpeed/60*p)];
        invoke(h.magnetHandler, 'processCommand', Command);
        Command = 'ArrayOfValues(3)= 0';
        invoke(h.magnetHandler, 'processCommand', Command);
        Command = 'ArrayOfValues(4)= 0';
        invoke(h.magnetHandler, 'processCommand', Command);
        Command = ['ArrayOfValues(5)= ' num2str(phase_angle(ii)+th0(nn)-th0(1))];
        invoke(h.magnetHandler, 'processCommand', Command);
        Command = ['Call getDocument.setSourceWaveform("' phase_index{ii} int2str(nn) '","SIN", ArrayOfValues)'];
        invoke(h.magnetHandler, 'processCommand', Command);
        % Command = ['getDocument.setParameter("' index{i} '", "WaveFormValues", "[0, ' num2str(Imod) ', ' num2str(n/60*p) ', 0, 0, 90]", infoArrayParameter)'];
        % invoke(h.magnetHandler, 'processCommand', Command);
        
        Command = ['Call getDocument().setCoilNumberOfTurns("U' int2str(nn) '", ' num2str(N_cond) ')'];
        invoke(h.magnetHandler, 'processCommand', Command);
        
        Command = ['Call getDocument().setCoilNumberOfTurns("V' int2str(nn) '", ' num2str(N_cond) ')'];
        invoke(h.magnetHandler, 'processCommand', Command);
        
        Command = ['Call getDocument().setCoilNumberOfTurns("W' int2str(nn) '", ' num2str(N_cond) ')'];
        invoke(h.magnetHandler, 'processCommand', Command);
        
    end
end

% solve
invoke(Doc, 'solveTransient2dWithMotion');

%% Post processing
NTime       = invoke(Solution,'getFieldSolutionTimeInstants',1);  % read simulation time values
Command     = 'TimeInstants = getDocument().getSolution().getFieldSolutionTimeInstants(1,time_instants)';
invoke(h.magnetHandler, 'processCommand', Command);

time = [];
for ij = 1:NTime
    invoke(h.magnetHandler, 'processCommand', ['Call setVariant(0, time_instants(' num2str(ij-1) '))']);
    time(ij) = invoke(h.magnetHandler, 'getVariant', 0);
end

Nbodies = invoke(Doc,'getNumberOfBodies',1);
invoke(h.magnetHandler, 'processCommand', 'REDIM ArrayOfValues(0)');

Flux_1sim = [];
Torque_1sim = [];
pm_loss = zeros(n_mag_simulati,NTime);

Ncoils = invoke(Doc,'getNumberOfCoils');

Nbarrier=0;
rotor_barrier_loss=zeros(NTime,Nbarrier);

for ii=1:NTime
    % flux, currents, resitance
    Flux_1sim_1ph = [];
    Torque_r = [];
    
    % problem ID
    Command = 'ReDim ProblemID(1)';
    invoke(h.magnetHandler, 'processCommand', Command);
    Command = 'ProblemID(0) = 1';
    invoke(h.magnetHandler, 'processCommand', Command);
    Command = ['ProblemID(1) = ' num2str(time(ii) * 1.001)];
    invoke(h.magnetHandler, 'processCommand', Command);
    
    for jj=1:Ncoils
        Command = ['Call getDocument.getSolution.getFluxLinkageThroughCoil(ProblemID, ' num2str(jj) ', magnitude, phase)'];
        invoke(h.magnetHandler, 'processCommand', Command);
        invoke(h.magnetHandler, 'processCommand', 'Call setVariant(0, magnitude)');
        x = invoke(h.magnetHandler, 'getVariant', 0);
        Flux_1sim_1ph = [Flux_1sim_1ph x];
    end
    
    % torque
    for jjj=1:Nbodies
        Command = ['Call getDocument.getSolution.getTorqueAboutOriginOnBody(ProblemID, ' num2str(jjj) ', torque_x, torque_y, torque_z)'];
        invoke(h.magnetHandler, 'processCommand', Command);
        invoke(h.magnetHandler, 'processCommand', 'Call setVariant(0, torque_z)');
        torque = invoke(h.magnetHandler, 'getVariant', 0);
        Torque_r = [Torque_r torque];
    end
    
    % load iron and pm losses (only for 360 deg simulation)
    if (ii == 1) && (xdeg == 360)
        
        % stator: hysteresis and classical
        %         h_loss_s = zeros(1,Mac.Qs);
        %         c_loss_s = zeros(1,Mac.Qs);
        
        for jjj = 1:geo.Qs
            Command = ['IronLoss = getDocument.getSolution.getIronLossInComponent (ProblemID,"statore_' num2str(jjj) '" , Losses)'];
            invoke(h.magnetHandler, 'processCommand', Command);
            invoke(h.magnetHandler, 'processCommand', 'Call setVariant(0, Losses)');
            test0 = invoke(h.magnetHandler, 'getVariant', 0);
            test0 = cell2mat(test0);
            h_loss_s(jjj) = test0(1);
            c_loss_s(jjj) = test0(2);
        end
        % rotor
        Command = ['IronLoss = getDocument.getSolution.getIronLossInComponent (ProblemID,"rotor" , Losses)'];
        invoke(h.magnetHandler, 'processCommand', Command);
        invoke(h.magnetHandler, 'processCommand', 'Call setVariant(0, Losses)');
        test0 = invoke(h.magnetHandler, 'getVariant', 0);
        test0 = cell2mat(test0);
        h_loss_r = test0(1);
        c_loss_r = test0(2);
    end
    % PM loss
    for jjj = 1:n_mag_simulati
        Command = ['PowerLoss=getDocument().getSolution().getOhmicLossInConductor(ProblemID,"Magnet_Bar_' num2str(jjj) '")'];
        invoke(h.magnetHandler, 'processCommand', Command);
        invoke(h.magnetHandler, 'processCommand', 'Call setVariant(0, PowerLoss)');
        pm_loss(jjj,ii) = invoke(h.magnetHandler, 'getVariant', 0);
    end
    Flux_1sim = [Flux_1sim;Flux_1sim_1ph];
    Torque_1sim = [Torque_1sim;Torque_r];
    ii;
    % Rotor conductor Losses:
    if (Nbarrier==0)
        jjj=1;
        rotor_barrier_loss(ii,jjj) = 0;
    else
        jk=1;
        for jjj=1:Nbarrier
            for jj=1:geo.ps
                Command = ['PowerLoss=getDocument().getSolution().getOhmicLossInConductor(ProblemID,"barrier#' num2str(jjj),'_',num2str(jj) '")'];
                invoke(h.magnetHandler, 'processCommand', Command);
                invoke(h.magnetHandler, 'processCommand', 'Call setVariant(0, PowerLoss)');
                rotor_barrier_loss(ii,jk) = invoke(h.magnetHandler, 'getVariant', 0);
                jk=jk+1;
            end
        end
    end
end
Torque_1sim = Torque_1sim * Q /geo.Qs;
Flux_1sim = Flux_1sim * Q /geo.Qs / N_parallel;
pm_loss = pm_loss * Q /geo.Qs;
Ppm = mean(sum(pm_loss));
if (xdeg == 360)
    Period_mot=gcd(geo.Qs,geo.p);
    if (geo.Qs==1)
        Pfes_h = mean(h_loss_s) * Q /geo.Qs;
        Pfes_c = mean(c_loss_s) * Q /geo.Qs;
    else
        Pfes_h = mean(h_loss_s) * Q;
        Pfes_c = mean(c_loss_s) * Q;
    end
    Pfer_h = h_loss_r / geo.ps * 2 * geo.p;
    Pfer_c = c_loss_r / geo.ps * 2 * geo.p;
    PjrBar=sum(mean(rotor_barrier_loss))*2 * geo.p/geo.ps;
end


%% FluxElab (RunS_FluxElab_MN_2)
% PostProcess Flux_UVW --> Flux_dq

% phase flux
% time sequence is POSITIVE
% physical sequence is POSITIVE

% coil positive is towards airgap
%
% % eliminates time = 0
% Fluxa = Flux_1sim(:,1);
% Fluxb = Flux_1sim(:,2);
% Fluxc = Flux_1sim(:,3);

% rotor mechanical position (RunS_SetParCase_new)
% t = 0, fase corrente in base a gamma, posizione rotore 0
pos_0 = 0;
pos_Mov = time/1000 * per.EvalSpeed * pi/30 + pos_0;  % rotor position (rad mec)
pos_Mov = pos_Mov(1:end)';

%% posizione del rotore (rad elt)
if exist('th0')
    % casi successivi feb 2010 (+90 perch� Waveform in Magnet\coil � un sin)
    theta = th0(1) * pi/180 + pos_Mov * geo.p;
else
    % casi precedenti feb 2010
    % d_axis Vs alpha_axis mechanical position
    theta = pi/2 + pos_Mov * geo.p;  % mec rad
end

Fluxa = zeros(size(theta,1),n3ph);
Fluxb = zeros(size(theta,1),n3ph);
Fluxc = zeros(size(theta,1),n3ph);
Fluxd = zeros(size(theta,1),n3ph);
Fluxq = zeros(size(theta,1),n3ph);

for ii=1:n3ph
    theta = theta+(th0(ii)-th0(1))*pi/180;
    
    Fluxa(:,ii) = Flux_1sim(:,1+3*(ii-1));
    Fluxb(:,ii) = Flux_1sim(:,2+3*(ii-1));
    Fluxc(:,ii) = Flux_1sim(:,3+3*(ii-1));
    
    Fluxd(:,ii) = 2/3*(+(Fluxa(:,ii)-0.5*Fluxb(:,ii)-0.5*Fluxc(:,ii)).*cos(theta)+sqrt(3)/2*(Fluxb(:,ii)-Fluxc(:,ii)).*sin(theta));
    Fluxq(:,ii) = 2/3*(-(Fluxa(:,ii)-0.5*Fluxb(:,ii)-0.5*Fluxc(:,ii)).*sin(theta)+sqrt(3)/2*(Fluxb(:,ii)-Fluxc(:,ii)).*cos(theta));
end

Fluxd = Fluxd(2:end,:);
Fluxd = [Fluxd(NTime-1,:);Fluxd];
Fluxq = Fluxq(2:end,:);
Fluxq = [Fluxq(NTime-1,:);Fluxq];

Fluxd = Fluxd(1:(end-1),:);
Fluxq = Fluxq(1:(end-1),:);

Fluxd = mean(Fluxd,2);
Fluxq = mean(Fluxq,2);

%% Voltage Elab (RunS_VoltageElab_MN)
% dq e.m.f.
Ed = per.EvalSpeed * geo.p * pi/30 * (-Fluxq);
Eq = per.EvalSpeed * geo.p * pi/30 * ( Fluxd);

% Ed_g = per.EvalSpeed * geo.p * pi/30 * ( Fluxq);
% Eq_g = per.EvalSpeed * geo.p * pi/30 * ( Fluxd);

pos_EMF = (pos_Mov(2:end) + pos_Mov(1:end-1))/2;
dT = (time(2) - time(1))/1000;
% Ea_ =  (Fluxa(2:end)-Fluxa(1:end-1))/dT;
% Eb_ =  (Fluxb(2:end)-Fluxb(1:end-1))/dT;
% Ec_ =  (Fluxc(2:end)-Fluxc(1:end-1))/dT;

% Ea = interp1(pos_EMF,Ea_,pos_Mov,'pchip');
% Eb = interp1(pos_EMF,Eb_,pos_Mov,'pchip');
% Ec = interp1(pos_EMF,Ec_,pos_Mov,'pchip');

% Ea = Ea(1:(end-1),:);
% Eb = Eb(1:(end-1),:);
% Ec = Ec(1:(end-1),:);


Ea_ =  (Fluxa(2:end,:)-Fluxa(1:end-1,:))/dT;
Eb_ =  (Fluxb(2:end,:)-Fluxb(1:end-1,:))/dT;
Ec_ =  (Fluxc(2:end,:)-Fluxc(1:end-1,:))/dT;

Ea = zeros(length(time),n3ph);
Eb = zeros(length(time),n3ph);
Ec = zeros(length(time),n3ph);

for ii=1:n3ph
    Ea(:,ii) = interp1(pos_EMF,Ea_(:,1),pos_Mov,'pchip');
    Eb(:,ii) = interp1(pos_EMF,Eb_(:,1),pos_Mov,'pchip');
    Ec(:,ii) = interp1(pos_EMF,Ec_(:,1),pos_Mov,'pchip');
end

Ea = Ea(1:(end-1),:);
Eb = Eb(1:(end-1),:);
Ec = Ec(1:(end-1),:);

Torque_0 = -mean(Torque_1sim(:,1));
tmp1 = id*ones(size(time,2),1);
tmp2 = iq*ones(size(time,2),1);

%per avere il grafico Torque coincidente come in FEMM devo:
%parto dalla seconda riga, metto la sesta riga nella prima riga,elimino
%l'ultima riga
Torque_1sim = Torque_1sim((2:(end)),:);
Torque_1sim = [Torque_1sim(NTime-1,:);Torque_1sim];
Torque_1sim = Torque_1sim(1:(end-1),:);

%save the model (.mn)
Command=['Call getDocument().save("',[pathname,filename(1:end-4) '.mn'],'")'];
% Command=['Call getDocument().save("',[CaseDir,filename(1:end-4) '.mn'],'")'];
invoke(h.magnetHandler, 'processCommand', Command);

% output
tempo   = time(2:end);         %ms
th = (geo.th0(1) * pi/180 + tempo/1000 * per.EvalSpeed * pi/30 * geo.p);

if (xdeg == 360)
    SOL.Pfes_h = Pfes_h;
    SOL.Pfes_c = Pfes_c;
    SOL.Pfer_h = Pfer_h;
    SOL.Pfer_c = Pfer_c;
    SOL.PjrBar = PjrBar;
    SOL.Ppm    = Ppm;
end
SOL.th = (th)';
SOL.id = tmp1(2:end)';
SOL.iq = tmp2(2:end)';
SOL.fd = Fluxd';
SOL.fq = Fluxq';
SOL.T = (mean([-Torque_1sim(:,1) Torque_1sim(:,2)],2))';

CloseMagnet(h);






