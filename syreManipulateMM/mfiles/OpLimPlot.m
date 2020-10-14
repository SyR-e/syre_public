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
p    = motorModel.data.p;

Id = motorModel.fdfq.Id;
Iq = motorModel.fdfq.Iq;
Fd = motorModel.fdfq.Fd;
Fq = motorModel.fdfq.Fq;
T  = motorModel.fdfq.T;
PF = sin(atan2(Iq,Id)-atan2(Fq,Fd));

MTPA = motorModel.AOA.MTPA;
MTPV = motorModel.AOA.MTPV;


% create figures
for ii=1:5
    hfig(ii) = figure();
    figSetting()
    hax(ii) = axes('OuterPosition',[0 0 1 1],...
        'XLim',[0 motorModel.data.nmax]);
    xlabel('$n$ [$rpm$]')
    switch ii
        case 1
            ylabel('$T$ [$Nm$]')
            title('Torque')
            set(hfig(1),'FileName',[pathname resFolder 'Torque.fig']);
        case 2
            ylabel('$P$ [$W$]')
            title('Power')
            set(hfig(2),'FileName',[pathname resFolder 'Power.fig']);
        case 3
            ylabel('$V_{l}$ [$V$]')
            title('Peak line voltage')
            set(hfig(3),'FileName',[pathname resFolder ,'Voltage.fig'])
        case 4
            ylabel('$I_{ph}$ [$A$]')
            title('Peak phase current')
            set(hfig(4),'FileName',[pathname resFolder 'Current.fig'])
        case 5
            ylabel('$cos \varphi$')
            title('Power factor')
            set(hfig(5),'FileName',[pathname resFolder 'PowerFactor.fig'])
    end
end

hfig(6) = figure();
figSetting()
hax(6) = axes('OuterPosition',[0 0 1 1],...
    'DataAspectRatio',[1 1 1],...
    'XLim',[min(Id,[],'all') max(Id,[],'all')],...
    'YLim',[min(Iq,[],'all') max(Iq,[],'all')]);
xlabel('$i_d$ [$A$]')
ylabel('$i_q$ [$A$]')
set(hfig(6),'FileName',[pathname resFolder 'DQplane.fig'])
[c,h] = contour(Id,Iq,abs(Id+j*Iq),'-k','DisplayName','$I$');
clabel(c,h)
[c,h] = contour(Id,Iq,T,'-b','DisplayName','$T$');
clabel(c,h)
% [c,h] = contour(Id,Iq,abs(Fd+j*Fq),'-r','DisplayName','$\lambda$');
% clabel(c,h)
plot(MTPA.id,MTPA.iq,'--k','DisplayName','MTPA')
plot(MTPV.id,MTPV.iq,':k','DisplayName','MTPV')

for ii=1:length(Plim)
    % curves
    for jj=1:length(hax)
        set(hax(jj),'ColorOrderIndex',ii);
    end
    iName = ['$I=' num2str(round(Ivect(ii),2)) '\,A$'];
    plot(hax(1),Plim{ii}.n,Plim{ii}.T,'DisplayName',iName);
    plot(hax(2),Plim{ii}.n,Plim{ii}.P,'DisplayName',iName);
    plot(hax(3),Plim{ii}.n,Plim{ii}.V*sqrt(3),'DisplayName',iName);
    plot(hax(4),Plim{ii}.n,Plim{ii}.I,'DisplayName',iName);
    plot(hax(5),Plim{ii}.n,Plim{ii}.PF,'DisplayName',iName);
    plot(hax(6),Plim{ii}.id_max,Plim{ii}.iq_max,'DisplayName',iName);
    
    % Point A
    for jj=1:length(hax)
        set(hax(jj),'ColorOrderIndex',ii);
    end
    plot(hax(1),Plim{ii}.n_A,Plim{ii}.T_A,'o','DisplayName','A');
    plot(hax(2),Plim{ii}.n_A,Plim{ii}.T_A*Plim{ii}.n_A*pi/30,'o','DisplayName','A')
    plot(hax(3),Plim{ii}.n_A,Plim{ii}.F_A*Plim{ii}.n_A*p*pi/30*sqrt(3),'o','DisplayName','A')
    plot(hax(4),Plim{ii}.n_A,abs(Plim{ii}.id_A+j*Plim{ii}.iq_A),'o','DisplayName','A')
    plot(hax(5),Plim{ii}.n_A,interp2(Id,Iq,PF,Plim{ii}.id_A,Plim{ii}.iq_A),'o','DisplayName','A')
    plot(hax(6),Plim{ii}.id_A,Plim{ii}.iq_A,'o','DisplayName','A')
    
    % Point B
    for jj=1:length(hax)
        set(hax(jj),'ColorOrderIndex',ii);
    end
    plot(hax(1),Plim{ii}.n_B,Plim{ii}.T_B,'d','DisplayName','B');
    plot(hax(2),Plim{ii}.n_B,Plim{ii}.T_B*Plim{ii}.n_B*pi/30,'d','DisplayName','B')
    plot(hax(3),Plim{ii}.n_B,Plim{ii}.F_B*Plim{ii}.n_B*p*pi/30*sqrt(3),'d','DisplayName','B')
    plot(hax(4),Plim{ii}.n_B,abs(Plim{ii}.id_B+j*Plim{ii}.iq_B),'d','DisplayName','B')
    plot(hax(5),Plim{ii}.n_B,interp2(Id,Iq,PF,Plim{ii}.id_B,Plim{ii}.iq_B),'d','DisplayName','B')
    plot(hax(6),Plim{ii}.id_B,Plim{ii}.iq_B,'d','DisplayName','B')
    
%     % Point C
%     for jj=1:length(hax)
%         set(hax(jj),'ColorOrderIndex',ii);
%     end
%     plot(hax(1),Plim{ii}.n_C,Plim{ii}.T_C,'^','DisplayName','C');
%     plot(hax(2),Plim{ii}.n_C,Plim{ii}.T_C*Plim{ii}.n_C*pi/30,'^','DisplayName','C')
%     plot(hax(3),Plim{ii}.n_C,Plim{ii}.F_C*Plim{ii}.n_C*p*pi/30*sqrt(3),'^','DisplayName','C')
%     plot(hax(4),Plim{ii}.n_C,abs(Plim{ii}.id_C+j*Plim{ii}.iq_C),'^','DisplayName','C')
%     plot(hax(5),Plim{ii}.n_C,interp2(Id,Iq,PF,Plim{ii}.id_C,Plim{ii}.iq_C),'^','DisplayName','C')
%     plot(hax(6),Plim{ii}.id_C,Plim{ii}.iq_C,'^','DisplayName','C')
    
    % Max speed
    if max(Plim{ii}.n)>=nmax
        for jj=1:length(hax)
            set(hax(jj),'ColorOrderIndex',ii);
        end
        nVect  = Plim{ii}.n(~isnan(Plim{ii}.n));
        idVect = Plim{ii}.id_max(~isnan(Plim{ii}.n));
        iqVect = Plim{ii}.iq_max(~isnan(Plim{ii}.n));
        [~,index,~] = unique(nVect);
        nVect  = nVect(index);
        idVect = idVect(index);
        iqVect = iqVect(index);
        id_M = interp1(nVect,idVect,nmax);
        iq_M = interp1(nVect,iqVect,nmax);
        
        % id_M = interp1(Plim{ii}.n(~isnan(Plim{ii}.n)),Plim{ii}.id_max(~isnan(Plim{ii}.n)),nmax);
        % iq_M = interp1(Plim{ii}.n(~isnan(Plim{ii}.n)),Plim{ii}.iq_max(~isnan(Plim{ii}.n)),nmax);
        
        plot(hax(1),nmax,interp2(Id,Iq,T,id_M,iq_M),'p','DisplayName','M');
        plot(hax(2),nmax,interp2(Id,Iq,T,id_M,iq_M)*nmax*pi/30,'p','DisplayName','M');
        plot(hax(3),nmax,interp2(Id,Iq,abs(Fd+j*Fq),id_M,iq_M)*nmax*p*pi/30*sqrt(3),'p','DisplayName','M')
        plot(hax(4),nmax,abs(id_M+j*iq_M),'p','DisplayName','M')
        plot(hax(5),nmax,interp2(Id,Iq,PF,id_M,iq_M),'p','DisplayName','M')
        plot(hax(6),id_M,iq_M,'p','DisplayName','M')
    end
end

% legend
% for ii=1:length(hax)
%     hleg(ii) = legend(hax(ii),'show','Location','best');
% end

