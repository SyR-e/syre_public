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

function [F_map,F_IM,F_bar]=elab_IM_Map(pathname,filename)

%% input
if nargin()<2
    [filename, pathname,~] = uigetfile('*.mat', 'Pick an output-file');
end

load([pathname filename])

%Rend=1.48241663779955e-06;


%% initialization
[nR,nC]=size(OUT);

F_map.Id  = zeros(nR,nC);
F_map.Iq  = zeros(nR,nC);
F_map.Fd  = zeros(nR,nC);
F_map.Fq  = zeros(nR,nC);
F_map.T   = zeros(nR,nC);
F_map.FdR = zeros(nR,nC);
F_map.FqR = zeros(nR,nC);
F_map.Ir  = zeros(nR,nC);

F_bar.I = {};
F_bar.V = {};
F_bar.F = {};
F_bar.Ftot = {};

F_IM.KR = zeros(nR,nC);
F_IM.tr = zeros(nR,nC);

%% reading results
for rr=1:nR
    for cc=1:nC
      F_map.Id(rr,cc)   = OUT{rr,cc}.Id;
      F_map.Iq(rr,cc)   = OUT{rr,cc}.Iq;
      F_map.Fd(rr,cc)   = OUT{rr,cc}.Fd;
      F_map.Fq(rr,cc)   = OUT{rr,cc}.Fq;
      F_map.Tist(rr,cc) = OUT{rr,cc}.T;
      F_map.Ir(rr,cc)   = OUT{rr,cc}.Ir;
      F_map.FdR(rr,cc)  = OUT{rr,cc}.Fdr;
      F_map.FqR(rr,cc)  = OUT{rr,cc}.Fqr;
            
      F_bar.I{rr,cc}    = OUT{rr,cc}.Ibar;
      F_bar.V{rr,cc}    = OUT{rr,cc}.Vbar;
      F_bar.F{rr,cc}    = OUT{rr,cc}.Fbar;
      F_bar.Ftot{rr,cc} = OUT{rr,cc}.FbarTot;
      
      F_IM.KR(rr,cc)    = OUT{rr,cc}.kr;
      
      
      % 2017 08 21 - calcolo trasformazioni bar-dq per calcolo costante di tempo
      Vr = geo.IM.kturns * geo.IM.Nbars/3 * bar2dq([F_bar.V{rr,cc} -F_bar.V{rr,cc}]',geo.IM.thR*pi/180,geo.IM.Nbars/geo.p);
      F_map.Vr(rr,cc)=abs(Vr(1)+j*Vr(2));
      
    end
end

% aggiunta dell'induttanza di testata
F_map.Fd = F_map.Fd+per.Lend*F_map.Id;
F_map.Fq = F_map.Fq+per.Lend*F_map.Iq;
F_map.T  = 3/2*geo.p*(F_map.Fd.*F_map.Iq-F_map.Fq.*F_map.Id);
% F_map.T=3/2*geoIM.p*(F_map.FdR.*F_map.Ir);

F_IM.KR(F_IM.KR==0) = NaN;

% 2017 09 21 - calcolo RR tramite potenze
F_map.Pbar  = zeros(nR,nC);
F_map.Pring = zeros(nR,nC);
for rr=1:nR
    for cc=1:nC
        F_map.Pbar(rr,cc)  = sum(F_bar.V{rr,cc}.*F_bar.I{rr,cc})*2*geo.p/geo.ps;
        F_map.Pring(rr,cc) = sum((2*per.IM.Rring)*(F_bar.I{rr,cc}*geo.IM.k).^2)*2*geo.p/geo.ps;
    end
end

F_map.Prot = F_map.Pbar+F_map.Pring;
F_IM.RRp   = F_map.Prot./(3/2*F_map.Ir.^2);

F_IM.RRv   = F_map.Vr./abs(F_map.Ir);

F_IM.RR     = F_IM.RRp;
F_IM.RS     = per.Rs*ones(size(F_map.Id));
F_IM.RRbar  = F_map.Pbar./(3/2*F_map.Ir.^2);
F_IM.RRring = F_IM.RR-F_IM.RRbar;

% Inverse-Gamma equivalent circuit parameters
F_IM.LS    = F_map.Fd./F_map.Id;            % stator inductance (Ld)
F_IM.Lt    = F_map.Fq./F_map.Iq;            % transient inductance, sigma*Ls (Lq)
F_IM.sigma = F_IM.Lt./F_IM.LS;              % leakage factor (Lq/Ld, inverse anisotropy factor)
F_IM.LM    = (F_IM.LS-F_IM.Lt)./F_IM.KR;    % magnetization inductance
F_IM.LR    = F_IM.LM./F_IM.KR;              % rotor inductance
F_IM.Lls   = F_IM.LS-F_IM.LM;               % stator leakage inductance
F_IM.Llr   = F_IM.LR-F_IM.LM;               % rotor leakage inductance
F_IM.KS    = F_IM.LM./F_IM.LS;              % stator coupling factor

F_IM.tr    = F_IM.LR./F_IM.RR;              % rotor time constant [s]
F_IM.wslip = F_map.Iq./F_map.Id./F_IM.tr;   % slip pulsation, rotor pulsation [rad/s]

save([pathname 'F_map.mat'],'F_map','F_IM','F_bar')

%% plot

if rr>1 && cc>1
    figure()
    figSetting()
    plot(F_map.Id',F_map.Fd','DisplayName','$\lambda_d$')
    plot(F_map.Iq,F_map.Fq,'DisplayName','$\lambda_q$')
    xlabel('i - A'), ylabel('[Vs]')
    legend('show')
    saveas(gcf,[pathname 'fluxVScurrent.fig'])
    
    figure()
    figSetting()
    surf(F_map.Id,F_map.Iq,F_map.Fd)
    view(3)
    xlabel('$i_{s,d}$ [A]')
    ylabel('$i_{s,q}$ [A]')
    zlabel('$\lambda_{s,d}$ [Vs]');
    saveas(gcf,[pathname 'FluxSD.fig'])
    
    figure()
    figSetting()
    surf(F_map.Id,F_map.Iq,F_map.Fq)
    view(3)
    xlabel('$i_{s,d}$ [A]')
    ylabel('$i_{s,q}$ [A]')
    zlabel('$\lambda_{s,q}$ [Vs]');
    saveas(gcf,[pathname 'FluxSQ.fig'])
    
    figure()
    figSetting()
    contour(F_map.Id,F_map.Iq,F_map.Ir,'ShowText','on');
    xlabel('$i_{s,d}$ [A]')
    ylabel('$i_{s,q}$ [A]')
    title('$i_r$ [A]');
    saveas(gcf,[pathname 'CurrentR.fig'])
    
    figure()
    figSetting()
    surf(F_map.Id,F_map.Iq,F_map.FdR)
    view(3)
    xlabel('$i_{s,d}$ [A]')
    ylabel('$i_{s,q}$ [A]')
    zlabel('$\lambda_r$ [Vs]')
    saveas(gcf,[pathname 'FluxR.fig'])

    figure()
    figSetting()
    surf(F_map.Id,F_map.Iq,F_map.FqR)
    view(3)
    xlabel('$i_{s,d}$ [A]')
    ylabel('$i_{s,q}$ [A]')
    zlabel('$\lambda_{r,q}$ [Vs]')
    saveas(gcf,[pathname 'FluxRQ.fig'])
    
    figure()
    figSetting()
    view(3)
    surf(F_map.Id,F_map.Iq,abs(F_map.T))
    xlabel('$i_{s,d}$ [A]')
    ylabel('$i_{s,q}$ [A]')
    zlabel('$T$ [Nm]')
    saveas(gcf,[pathname 'Torque.fig'])
    
    figure()
    figSetting()
    view(3)
    surf(F_map.Id,F_map.Iq,F_IM.KR)
    xlabel('$i_{s,d}$ [A]')
    ylabel('$i_{s,q}$ [A]')
    zlabel('$k_r$')
    saveas(gcf,[pathname 'KR.fig'])
    
    figure()
    figSetting()
    view(3)
    surf(F_map.Id,F_map.Iq,F_IM.Lls)
    xlabel('$i_{s,d}$ [A]')
    ylabel('$i_{s,q}$ [A]')
    zlabel('$L_{\sigma,s}$ [H]')
    saveas(gcf,[pathname 'InductanceLeakageS.fig'])
    
    figure()
    figSetting()
    view(3)
    surf(F_map.Id,F_map.Iq,F_IM.Llr)
    xlabel('$i_{s,d}$ [A]')
    ylabel('$i_{s,q}$ [A]')
    zlabel('$L_{\sigma,r}$ [H]')
    saveas(gcf,[pathname 'InductanceLeakageR.fig'])
    
    figure()
    figSetting()
    view(3)
    surf(F_map.Id,F_map.Iq,F_IM.tr)
    xlabel('$i_{s,d}$ [A]')
    ylabel('$i_{s,q}$ [A]')
    zlabel('$\tau_r$ [s]')
    saveas(gcf,[pathname 'TauR.fig'])
    
    figure()
    figSetting()
    view(3)
    surf(F_map.Id,F_map.Iq,F_IM.wslip)
    xlabel('$i_{s,d}$ [A]')
    ylabel('$i_{s,q}$ [A]')
    zlabel('$\omega_{slip}$ [rad/s]')
    saveas(gcf,[pathname 'OmegaSlip.fig'])
    
    figure()
    figSetting()
    view(3)
    mesh(F_map.Id,F_map.Iq,F_map.T,'EdgeColor','b','DisplayName','mean')
    mesh(F_map.Id,F_map.Iq,-F_map.Tist,'EdgeColor','r','DisplayName','FEMM')
    xlabel('$i_{s,d}$ [A]')
    ylabel('$i_{s,q}$ [A]')
    zlabel('$T$ [Nm]')
    legend('show','Location','northeastoutside')
    saveas(gcf,[pathname 'TorqueComparison.fig'])
    
    figure()
    figSetting()
    surf(F_map.Id,F_map.Iq,atan2(F_map.FqR,F_map.FdR)*180/pi)
    view(3)
    xlabel('$i_{s,d}$ [A]')
    ylabel('$i_{s,q}$ [A]')
    zlabel('$\delta$ [$^\circ$]')
    saveas(gcf,[pathname 'Delta.fig'])
    
    figure()
    figSetting()
    surf(F_map.Id,F_map.Iq,F_IM.Lt)
    view(3)
    xlabel('$i_{s,d}$ [A]')
    ylabel('$i_{s,q}$ [A]')
    zlabel('$\sigma\cdot L_s$')
    saveas(gcf,[pathname 'SigmaLs.fig'])
    
end

% Interpolate over finer grid
n = 256;

Id = linspace(round(min(min(F_map.Id))),round(max(max(F_map.Id))),n);
Iq = linspace(round(min(min(F_map.Iq))),round(max(max(F_map.Iq))),n);
[Id,Iq] = meshgrid(Id,Iq);

Fd = interp2(F_map.Id,F_map.Iq,F_map.Fd,Id,Iq,'spline');
Fq = interp2(F_map.Id,F_map.Iq,F_map.Fq,Id,Iq,'spline');
T  = interp2(F_map.Id,F_map.Iq,F_map.T,Id,Iq,'spline');

IM.FdR   = interp2(F_map.Id,F_map.Iq,F_map.FdR,Id,Iq,'spline');
IM.FqR   = interp2(F_map.Id,F_map.Iq,F_map.FqR,Id,Iq,'spline');
IM.Ir    = interp2(F_map.Id,F_map.Iq,F_map.Ir,Id,Iq,'spline');
IM.tr    = interp2(F_map.Id,F_map.Iq,F_IM.tr,Id,Iq,'spline');
IM.KR    = interp2(F_map.Id,F_map.Iq,F_IM.KR,Id,Iq,'spline');
IM.KS    = interp2(F_map.Id,F_map.Iq,F_IM.KS,Id,Iq,'spline');
IM.LM    = interp2(F_map.Id,F_map.Iq,F_IM.LM,Id,Iq,'spline');
IM.LS    = interp2(F_map.Id,F_map.Iq,F_IM.LS,Id,Iq,'spline');
IM.LR    = interp2(F_map.Id,F_map.Iq,F_IM.LR,Id,Iq,'spline');
IM.Lls   = interp2(F_map.Id,F_map.Iq,F_IM.Lls,Id,Iq,'spline');
IM.Llr   = interp2(F_map.Id,F_map.Iq,F_IM.Llr,Id,Iq,'spline');
IM.Lt    = interp2(F_map.Id,F_map.Iq,F_IM.Lt,Id,Iq,'spline');
IM.RR    = interp2(F_map.Id,F_map.Iq,F_IM.RR,Id,Iq,'spline');
IM.RS    = interp2(F_map.Id,F_map.Iq,F_IM.RS,Id,Iq,'spline');
IM.RRbar = interp2(F_map.Id,F_map.Iq,F_IM.RRbar,Id,Iq,'spline');
IM.Wslip = Iq./Id./IM.tr;

save ([pathname 'fdfq_idiq_n' num2str(n) '.mat'],'Fd','Fq','Id','Iq','T','IM');


% build_fdfq_idiq_IM(pathname,'F_map.mat',256);

