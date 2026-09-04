% Copyright 2014
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

function plot_singt_5(out,delta_sim_singt,newDir,filemot)
% single working point has been simulated

if delta_sim_singt<=360
    nRep = 360/delta_sim_singt; % number of repetition needed
else
    nRep = 1;
end
T  = repmat(out.SOL.T,1,nRep);    % last point added for plot
fd1 = repmat(out.SOL.fd1,1,nRep) ;  % last point added for plot
fq1 = repmat(out.SOL.fq1,1,nRep) ;  % last point added for plot
fd3 = repmat(out.SOL.fd3,1,nRep) ;  % last point added for plot
fq3 = repmat(out.SOL.fq3,1,nRep) ;  % last point added for plot

if(isfield(out.SOL,'parms'))
    for i=1:length(out.SOL.parms)
        phi_pm(i,:) = out.SOL.parms{i}.Phi_pm;
    end
end

for ii = 1:height(out.SOL.fa)

    id1 = repmat(out.SOL.id1,1,nRep) ;  % last point added for plot
    iq1 = repmat(out.SOL.iq1,1,nRep) ;  % last point added for plot
    id3 = repmat(out.SOL.id3,1,nRep) ;  % last point added for plot
    iq3 = repmat(out.SOL.iq3,1,nRep) ;  % last point added for plot
    i0 = repmat(out.SOL.i0,1,nRep) ;  % last point added for plot

    if isfield(out.SOL,'ia')
        if nRep==1
            ia(ii,:) = out.SOL.ia(ii,:);
            ib(ii,:) = out.SOL.ib(ii,:) ;
            ic(ii,:) = out.SOL.ic(ii,:);
            id(ii,:) = out.SOL.id(ii,:);
            ie(ii,:) = out.SOL.ie(ii,:);
    
            fa(ii,:) = out.SOL.fa(ii,:);
            fb(ii,:) = out.SOL.fb(ii,:);
            fc(ii,:) = out.SOL.fc(ii,:);
            fd(ii,:) = out.SOL.fd(ii,:);
            fe(ii,:) = out.SOL.fe(ii,:);
        else
            iph = phaseQuantityDecoding_5(out.SOL.ia(ii,:),out.SOL.ib(ii,:),out.SOL.ic(ii,:),out.SOL.id(ii,:),out.SOL.ie(ii,:),delta_sim_singt);
            ia(ii,:) = iph.a;
            ib(ii,:) = iph.b ;
            ic(ii,:) = iph.c ;
            id(ii,:) = iph.d ;
            ie(ii,:) = iph.e ;
    
            fph = phaseQuantityDecoding_5(out.SOL.fa(ii,:),out.SOL.fb(ii,:),out.SOL.fc(ii,:),out.SOL.fd(ii,:),out.SOL.fe(ii,:),delta_sim_singt);
            fa(ii,:) = fph.a ;
            fb(ii,:) = fph.b ;
            fc(ii,:) = fph.c ;
            fd(ii,:) = fph.d ;
            fe(ii,:) = fph.e ;
        end
    else
        ia = NaN;
        ib = NaN;
        ic = NaN;
        id = NaN;
        ie = NaN;
        fa = NaN;
        fb = NaN;
        fc = NaN;
        fd = NaN;
        fe = NaN;        
    end
end

% T  = [repmat(out.SOL.T,1,nRep) out.SOL.T(1)];    % last point added for plot
% fd = [repmat(out.SOL.fd,1,nRep) out.SOL.fd(1)];  % last point added for plot
% fq = [repmat(out.SOL.fq,1,nRep) out.SOL.fq(1)];  % last point added for plot
% id = [repmat(out.SOL.id,1,nRep) out.SOL.id(1)];  % last point added for plot
% iq = [repmat(out.SOL.iq,1,nRep) out.SOL.iq(1)];  % last point added for plot
% 
% if isfield(out.SOL,'ia')
%     if nRep==1
%         ia = [out.SOL.ia out.SOL.ia(:,1)];
%         ib = [out.SOL.ib out.SOL.ib(:,1)];
%         ic = [out.SOL.ic out.SOL.ic(:,1)];
% 
%         fa = [out.SOL.fa out.SOL.fa(:,1)];
%         fb = [out.SOL.fb out.SOL.fb(:,1)];       rimosso da simodelsa
%         fc = [out.SOL.fc out.SOL.fc(:,1)];
%     else
%         iph = phaseQuantityDecoding(out.SOL.ia,out.SOL.ib,out.SOL.ic,delta_sim_singt);
%         ia = [iph.a iph.a(:,1)];
%         ib = [iph.b iph.b(:,1)];
%         ic = [iph.c iph.c(:,1)];
% 
%         fph = phaseQuantityDecoding(out.SOL.fa,out.SOL.fb,out.SOL.fc,delta_sim_singt);
%         fa = [fph.a fph.a(:,1)];
%         fb = [fph.b fph.b(:,1)];
%         fc = [fph.c fph.c(:,1)];
%     end
% else
%     ia = NaN;
%     ib = NaN;
%     ic = NaN;
% 
%     fa = NaN;
%     fb = NaN;
%     fc = NaN;
% end


if(isfield(out.SOL,'tpwm'))
    th = (2*pi*50)*out.SOL.tpwm*(180/pi);
else
    % th = linspace(0,360,length(T));
    if length(out.SOL.th)>1
        dth = out.SOL.th(2)-out.SOL.th(1);
        th = 0:dth:dth*(length(T)-1);
    else
        th = 0:60:360;
    end
end

th = th*pi/180;

hfig(1) = figure();
figSetting()
phi_pm_mean = mean(phi_pm - phi_pm(:,5));
G = [0 1*0 0 0 0; -1 0 0 0 0; 0 0 0 3*0 0; 0 0 -3 0 0;  0 0 0 0 0];
idq = [id1; iq1; id3; iq3; i0];
tdq = (G*max(phi_pm-phi_pm(:,5))').*idq;
plot(th, tdq)
hold on
plot(th, sum(tdq))
legend('d1','q1','d3','q3','0','sum')

hfig(1) = figure();
figSetting()
set(gca,'XLim',[0 2*pi],'XTick',0:pi/3:2*pi);
xlabel('$\omega t$ (rad)')
ylabel('Current (A)')
xtickformat('%.0f')
plot(th,iq1,th,iq3);

hfig(2) = figure();
figSetting()
set(gca,'XLim',[0 2*pi],'XTick',0:pi/3:2*pi,'YTick',-5:2:5);
xlabel('$\omega t$ (rad)')
ylabel('Torque (Nm)')
xtickformat('%.0f')
plot(th,T);

gamma = atan2(iq1,id1);
delta = atan2(fq1(1,:),fd1(1,:));
IPF = sin(gamma-delta);

hfig(1) = figure();
figSetting()
subplot(2,1,1)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ (elt deg)')
ylabel('(Nm)')
title(['Mean Torque = ' num2str(mean(out.T)) ' Nm'])
plot(th,T);
subplot(2,1,2)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ (elt deg)')
ylabel('IPF')
% title(['Mean IPF = ' num2str(out.IPF)])
title(['Mean IPF = ' num2str(mean(IPF))])
plot(th,IPF);
if isoctave()
    fig_name=strcat(newDir, filemot(1:end-4), '_T_gamma');
    hgsave(hfig(1),[fig_name]);
else
    % saveas(hfig(1),[newDir filemot(1:end-4) '_T_gamma']);
end

hfig(2) = figure();
figSetting()
subplot(2,1,1)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ (elt deg)')
ylabel('$\lambda_d$ (Vs)')
title(['Mean $\lambda_d$ = ' num2str(mean(fd1(1,:))) ' Vs'])
plot(th,fd1(1,:));
subplot(2,1,2)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ (elt deg)')
ylabel('$\lambda_q$ (Vs)')
title(['Mean $\lambda_q$ = ' num2str(mean(fq1(1,:))) ' Vs'])
plot(th,fq1(1,:));
if isoctave()
    fig_name=strcat(newDir, filemot(1:end-4), '_plot_flux');
    hgsave(hfig(2),[fig_name]);
else
    % saveas(hfig(2),[newDir filemot(1:end-4) '_plot_flux']);
end

hfig(3) = figure();
figSetting()
subplot(2,1,1)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ (elt deg)')
ylabel('$\lambda_{d3}$ (Vs)')
title(['Mean $\lambda_{d3}$ = ' num2str(mean(fd3(1,:))) ' Vs'])
plot(th,fd3(1,:));
subplot(2,1,2)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ (elt deg)')
ylabel('$\lambda_{q3}$ (Vs)')
title(['Mean $\lambda_{q3}$ = ' num2str(mean(fq3(1,:))) ' Vs'])
plot(th,fq3(1,:));
% if isoctave()
%     fig_name=strcat(newDir, filemot(1:end-4), '_plot_flux');
%     hgsave(hfig(2),[fig_name]);
% else
%     saveas(hfig(2),[newDir filemot(1:end-4) '_plot_flux']);
% end

hfig(4) = figure();
figSetting()
subplot(2,1,1)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ (elt deg)')
ylabel('$i_d$ (A)')
title(['Mean $i_d$ = ' num2str(mean(id1(1,:))) ' A'])
plot(th,id1(1,:));
subplot(2,1,2)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ (elt deg)')
ylabel('$i_q$ (A)')
title(['Mean $i_q$ = ' num2str(mean(iq1(1,:))) ' A'])
plot(th,iq1(1,:));
% if isoctave()
%     fig_name=strcat(newDir, filemot(1:end-4), '_plot_flux');
%     hgsave(hfig(2),[fig_name]);
% else
%     saveas(hfig(2),[newDir filemot(1:end-4) '_plot_flux']);
% end

hfig(5) = figure();
figSetting()
subplot(2,1,1)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ (elt deg)')
ylabel('$i_{d3}$ (Vs)')
title(['Mean $i_{d3}$ = ' num2str(mean(id3(1,:))) ' A'])
plot(th,id3(1,:));
subplot(2,1,2)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ (elt deg)')
ylabel('$i_{q3}$ (Vs)')
title(['Mean $i_{q3}$ = ' num2str(mean(iq3(1,:))) ' A'])
plot(th,iq3(1,:));
% if isoctave()
%     fig_name=strcat(newDir, filemot(1:end-4), '_plot_flux');
%     hgsave(hfig(2),[fig_name]);
% else
%     saveas(hfig(2),[newDir filemot(1:end-4) '_plot_flux']);
% end

if ~sum(isnan(fa))
    hfig(6) = figure();
    figSetting()
    subplot(2,1,1)
    set(gca,'XLim',[0 360],'XTick',0:60:360);
    xlabel('$\theta$ (elt deg)')
    ylabel('$\lambda_{abc}$ (Vs)')
    title(['Phase flux linkages'])
    plot(th,fa);
    plot(th,fb);
    plot(th,fc);
    plot(th,fd);
    plot(th,fe);
    subplot(2,1,2)
    set(gca,'XLim',[0 360],'XTick',0:60:360);
    xlabel('$\theta$ (elt deg)')
    ylabel('$i_{abc}$ (A)')
    title(['Phase currents'])
    plot(th,ia,'DisplayName','$i_a$');
    plot(th,ib,'DisplayName','$i_b$');
    plot(th,ic,'DisplayName','$i_c$');
    plot(th,id,'DisplayName','$i_d$');
    plot(th,ie,'DisplayName','$i_e$');
    if length(ia(:,1))<2
        legend show
    end
    if isoctave()
        fig_name=strcat(newDir, filemot(1:end-4), '_plot_phase');
        hgsave(hfig(3),[fig_name]);
    else
        % saveas(hfig(3),[newDir filemot(1:end-4) '_plot_phase']);
    end
end

% if(isfield(out.SOL,'parms'))
%     hfig(7) = figure();
%     figSetting()
%     subplot(2,1,1)
%     set(gca,'XLim',[0 360],'XTick',0:60:360);
%     xlabel('$\theta$ (elt deg)')
%     ylabel('d-axis PM Flux (Vs)')
%     % title(['Mean PM Flux = ' num2str(mean(out.T)) ' Nm'])
%     plot(th,phi_pm(:,[1 3]));
%     subplot(2,1,2)
%     set(gca,'XLim',[0 360],'XTick',0:60:360);
%     xlabel('$\theta$ (elt deg)')
%     ylabel('d3-axis PM Flux (Vs)')
%     % title(['Mean IPF = ' num2str(mean(IPF))])
%     plot(th,phi_pm(:,[2 4]));
% end
% if isoctave()
%     fig_name=strcat(newDir, filemot(1:end-4), '_T_gamma');
%     hgsave(hfig(1),[fig_name]);
% else
%     % saveas(hfig(1),[newDir filemot(1:end-4) '_T_gamma']);
% end