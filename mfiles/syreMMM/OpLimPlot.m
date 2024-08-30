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

function [hfig,pathnameOut] = OpLimPlot(Plim,Ivect,motorModel)

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = [motName '_results\MMM results\' 'OpLim_' datestr(now,30) '\'];
pathnameOut = [pathname resFolder];

nmax = motorModel.data.nmax;
% p    = motorModel.data.p;

Id = motorModel.FluxMap_dq.Id;
Iq = motorModel.FluxMap_dq.Iq;
% Fd = motorModel.FluxMap_dq.Fd;
% Fq = motorModel.FluxMap_dq.Fq;
T  = motorModel.FluxMap_dq.T;
% PF = sin(atan2(Iq,Id)-atan2(Fq,Fd));

MTPA  = motorModel.controlTrajectories.MTPA;
MTPV  = motorModel.controlTrajectories.MTPV;
MPFPA = motorModel.controlTrajectories.MPFPA;


% create figures
for ii=1:6
    hfig(ii) = figure();
    figSetting()
    hax(ii) = axes('OuterPosition',[0 0 1 1],...
        'XLim',[0 motorModel.data.nmax]);
    xlabel('$n$ (rpm)')
    switch ii
        case 1
            ylabel('$T$ (Nm)')
            title('Torque')
            set(hfig(1),'FileName',[pathname resFolder 'Torque.fig']);
        case 2
            ylabel('$P$ (W)')
            title('Power')
            set(hfig(2),'FileName',[pathname resFolder 'Power.fig']);
        case 3
            ylabel('$V_{l}$ (Vpk)')
            title('Peak line voltage')
            set(hfig(3),'FileName',[pathname resFolder ,'Voltage.fig'])
        case 4
            ylabel('$I_{ph}$ (Apk)')
            title('Peak phase current')
            set(hfig(4),'FileName',[pathname resFolder 'Current.fig'])
        case 5
            ylabel('$cos \varphi$')
            title('Power factor')
            set(hfig(5),'FileName',[pathname resFolder 'PowerFactor.fig'])
        case 6
            ylabel('$\lambda$ (Vs)')
            title('Flux linkage')
            set(hfig(6),'FileName',[pathname resFolder 'FluxLinkage.fig'])
        case 7
    end
end

hfig(7) = figure();
figSetting()
hax(7) = axes('OuterPosition',[0 0 1 1],...
    'DataAspectRatio',[1 1 1],...
    'XLim',[min(Id,[],'all') max(Id,[],'all')],...
    'YLim',[min(Iq,[],'all') max(Iq,[],'all')]);
xlabel('$i_d$ (A)')
ylabel('$i_q$ (A)')
set(hfig(7),'FileName',[pathname resFolder 'DQplane.fig'])
[c,h] = contour(Id,Iq,abs(Id+j*Iq),'-k','DisplayName','$I$');
clabel(c,h)
[c,h] = contour(Id,Iq,T,'-b','DisplayName','$T$');
clabel(c,h)

plot(MTPA.id,MTPA.iq,'-k','DisplayName','MTPA')
plot(MTPV.id,MTPV.iq,':k','DisplayName','MTPV')
% plot(MPFPA.id,MPFPA.iq,'--k','DisplayName','MPFPA')

for ii=1:length(Plim)
    % curves
    for jj=1:length(hax)
        set(hax(jj),'ColorOrderIndex',ii);
    end
    iName = ['$I=' num2str(round(Ivect(ii),2)) '$ A'];
    index = Plim{ii}.n<=nmax;
    plot(hax(1),Plim{ii}.n(index),Plim{ii}.T(index),'DisplayName',iName);
    plot(hax(2),Plim{ii}.n(index),Plim{ii}.P(index),'DisplayName',iName);
    plot(hax(3),Plim{ii}.n(index),Plim{ii}.V(index)*sqrt(3),'DisplayName',iName);
    set(hax(3),'ColorOrderIndex',ii);
    plot(hax(3),Plim{ii}.n(index),Plim{ii}.E(index)*sqrt(3),'--','DisplayName',iName);
    plot(hax(4),Plim{ii}.n(index),Plim{ii}.I(index),'DisplayName',iName);
    plot(hax(5),Plim{ii}.n(index),Plim{ii}.PF(index),'DisplayName',iName);
    set(hax(5),'ColorOrderIndex',ii);
    plot(hax(5),Plim{ii}.n(index),Plim{ii}.IPF(index),'--','DisplayName',iName);
    plot(hax(6),Plim{ii}.n(index),Plim{ii}.F(index),'DisplayName',iName);
    plot(hax(7),Plim{ii}.id(index),Plim{ii}.iq(index),'DisplayName',iName);
    
    % Point A
    for jj=1:length(hax)
        set(hax(jj),'ColorOrderIndex',ii);
    end
    plot(hax(1),Plim{ii}.n_A,Plim{ii}.T_A,'o','DisplayName','A');
    plot(hax(2),Plim{ii}.n_A,Plim{ii}.T_A*Plim{ii}.n_A*pi/30,'o','DisplayName','A')
    plot(hax(3),Plim{ii}.n_A,interp1(Plim{ii}.n,Plim{ii}.V,Plim{ii}.n_A)*sqrt(3),'o','DisplayName','A')
    set(hax(3),'ColorOrderIndex',ii);
    plot(hax(3),Plim{ii}.n_A,interp1(Plim{ii}.n,Plim{ii}.E,Plim{ii}.n_A)*sqrt(3),'o','DisplayName','A')
    plot(hax(4),Plim{ii}.n_A,abs(Plim{ii}.id_A+j*Plim{ii}.iq_A),'o','DisplayName','A')
    plot(hax(5),Plim{ii}.n_A,interp1(Plim{ii}.n,Plim{ii}.PF,Plim{ii}.n_A),'o','DisplayName','A')
    set(hax(5),'ColorOrderIndex',ii);
    plot(hax(5),Plim{ii}.n_A,interp1(Plim{ii}.n,Plim{ii}.IPF,Plim{ii}.n_A),'o','DisplayName','A')
    plot(hax(6),Plim{ii}.n_A,interp1(Plim{ii}.n,Plim{ii}.F,Plim{ii}.n_A),'o','DisplayName','A')
    plot(hax(7),Plim{ii}.id_A,Plim{ii}.iq_A,'o','DisplayName','A')

    % Point B
    if Plim{ii}.n_B<Plim{ii}.n_M
        for jj=1:length(hax)
            set(hax(jj),'ColorOrderIndex',ii);
        end
        plot(hax(1),Plim{ii}.n_B,Plim{ii}.T_B,'d','DisplayName','B');
        plot(hax(2),Plim{ii}.n_B,Plim{ii}.T_B*Plim{ii}.n_B*pi/30,'d','DisplayName','B')
        plot(hax(3),Plim{ii}.n_B,interp1(Plim{ii}.n,Plim{ii}.V,Plim{ii}.n_B)*sqrt(3),'d','DisplayName','B')
        set(hax(3),'ColorOrderIndex',ii);
        plot(hax(3),Plim{ii}.n_B,interp1(Plim{ii}.n,Plim{ii}.E,Plim{ii}.n_B)*sqrt(3),'d','DisplayName','B')
        plot(hax(4),Plim{ii}.n_B,abs(Plim{ii}.id_B+j*Plim{ii}.iq_B),'d','DisplayName','B')
        plot(hax(5),Plim{ii}.n_B,interp1(Plim{ii}.n,Plim{ii}.PF,Plim{ii}.n_B),'d','DisplayName','B')
        set(hax(5),'ColorOrderIndex',ii);
        plot(hax(5),Plim{ii}.n_B,interp1(Plim{ii}.n,Plim{ii}.IPF,Plim{ii}.n_B),'d','DisplayName','B')
        plot(hax(6),Plim{ii}.n_B,interp1(Plim{ii}.n,Plim{ii}.F,Plim{ii}.n_B),'d','DisplayName','B')
        plot(hax(7),Plim{ii}.id_B,Plim{ii}.iq_B,'d','DisplayName','B')
    end

    % Max speed
    for jj=1:length(hax)
        set(hax(jj),'ColorOrderIndex',ii);
    end
    plot(hax(1),Plim{ii}.n_M,Plim{ii}.T_M,'p','DisplayName','M');
    plot(hax(2),Plim{ii}.n_M,Plim{ii}.T_M*Plim{ii}.n_M*pi/30,'p','DisplayName','M')
    plot(hax(3),Plim{ii}.n_M,interp1(Plim{ii}.n,Plim{ii}.V,Plim{ii}.n_M)*sqrt(3),'p','DisplayName','M')
    set(hax(3),'ColorOrderIndex',ii);
    plot(hax(3),Plim{ii}.n_M,interp1(Plim{ii}.n,Plim{ii}.E,Plim{ii}.n_M)*sqrt(3),'p','DisplayName','M')
    plot(hax(4),Plim{ii}.n_M,abs(Plim{ii}.id_M+j*Plim{ii}.iq_M),'p','DisplayName','M')
    plot(hax(5),Plim{ii}.n_M,interp1(Plim{ii}.n,Plim{ii}.PF,Plim{ii}.n_M),'p','DisplayName','M')
    set(hax(5),'ColorOrderIndex',ii);
    plot(hax(5),Plim{ii}.n_M,interp1(Plim{ii}.n,Plim{ii}.IPF,Plim{ii}.n_M),'p','DisplayName','M')
    plot(hax(6),Plim{ii}.n_M,interp1(Plim{ii}.n,Plim{ii}.F,Plim{ii}.n_M),'p','DisplayName','M')
    plot(hax(7),Plim{ii}.id_M,Plim{ii}.iq_M,'p','DisplayName','M')
end

for ii=1:length(hfig)
    tmp = get(hfig(ii),'FileName');
    [~,name,~] = fileparts(tmp);
    set(hfig(ii),'Name',name);
end
