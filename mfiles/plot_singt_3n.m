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

function plot_singt_3n(n3ph, out,delta_sim_singt,newDir,filemot)
% single working point has been simulated

if delta_sim_singt<=360
    nRep = 360/delta_sim_singt; % number of repetition needed
else
    nRep = 1;
end
T  = repmat(out.SOL.T,1,nRep);    % last point added for plot

for i=1:n3ph
    fd(:,i) = repmat(out.SOL.fd(:,i),1,nRep) ;  % last point added for plot
    fq(:,i) = repmat(out.SOL.fq(:,i),1,nRep) ;  % last point added for plot
end

for ii = 1:height(out.SOL.fd)
    for i=1:n3ph
        id(:,i) = repmat(out.SOL.id(:,i),1,nRep) ;  % last point added for plot
        iq(:,i) = repmat(out.SOL.iq(:,i),1,nRep) ;  % last point added for plot
    end

    if isfield(out.SOL,'ia')
        if nRep==1
            ia(ii,:) = out.SOL.ia(ii,:);
            ib(ii,:) = out.SOL.ib(ii,:);
            ic(ii,:) = out.SOL.ic(ii,:);
    
            fa(ii,:) = out.SOL.fa(ii,:);
            fb(ii,:) = out.SOL.fb(ii,:);
            fc(ii,:) = out.SOL.fc(ii,:);
        else
            iph = phaseQuantityDecoding_6(out.SOL.ia(ii,:),out.SOL.ib(ii,:),out.SOL.ic(ii,:),out.SOL.id(ii,:),out.SOL.ie(ii,:),out.SOL.if(ii,:),delta_sim_singt);
            ia(ii,:) = iph.a;
            ib(ii,:) = iph.b ;
            ic(ii,:) = iph.c ;
            id(ii,:) = iph.d ;
            ie(ii,:) = iph.e ;
            lf(ii,:) = iph.f ;
    
            fph = phaseQuantityDecoding_6(out.SOL.fa(ii,:),out.SOL.fb(ii,:),out.SOL.fc(ii,:),out.SOL.fd(ii,:),out.SOL.fe(ii,:),out.SOL.ff(ii,:),delta_sim_singt);
            fa(ii,:) = fph.a ;
            fb(ii,:) = fph.b ;
            fc(ii,:) = fph.c ;
            fd(ii,:) = fph.d ;
            fe(ii,:) = fph.e ;
            ff(ii,:) = fph.f ;
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

gamma = atan2(iq(:,1),id(:,1));
delta = atan2(fq(:,1),fd(:,1));
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
ylabel('$\lambda_da$ (Vs)')
title(['Mean $\lambda_d$ = ' num2str(mean(fd(:,1))) ' Vs'])
plot(th,fd(:,1));
subplot(2,1,2)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ (elt deg)')
ylabel('$\lambda_qa$ (Vs)')
title(['Mean $\lambda_q$ = ' num2str(mean(fq(:,1))) ' Vs'])
plot(th,fq(:,1));
if isoctave()
    fig_name=strcat(newDir, filemot(1:end-4), '_plot_flux');
    hgsave(hfig(2),[fig_name]);
else
    % saveas(hfig(2),[newDir filemot(1:end-4) '_plot_flux']);
end

% hfig(3) = figure();
% figSetting()
% subplot(2,1,1)
% set(gca,'XLim',[0 360],'XTick',0:60:360);
% xlabel('$\theta$ (elt deg)')
% ylabel('$\lambda_{db}$ (Vs)')
% title(['Mean $\lambda_{d3}$ = ' num2str(mean(fdb(1,:))) ' Vs'])
% plot(th,fdb(1,:));
% subplot(2,1,2)
% set(gca,'XLim',[0 360],'XTick',0:60:360);
% xlabel('$\theta$ (elt deg)')
% ylabel('$\lambda_{qb}$ (Vs)')
% title(['Mean $\lambda_{q3}$ = ' num2str(mean(fqb(1,:))) ' Vs'])
% plot(th,fqb(1,:));
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
ylabel('$i_da$ (A)')
title(['Mean $i_d$ = ' num2str(mean(id(:,1))) ' A'])
plot(th,id(:,1));
subplot(2,1,2)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ (elt deg)')
ylabel('$i_qa$ (A)')
title(['Mean $i_q$ = ' num2str(mean(iq(:,1))) ' A'])
plot(th,iq(:,1));
% if isoctave()
%     fig_name=strcat(newDir, filemot(1:end-4), '_plot_flux');
%     hgsave(hfig(2),[fig_name]);
% else
%     saveas(hfig(2),[newDir filemot(1:end-4) '_plot_flux']);
% end

% hfig(5) = figure();
% figSetting()
% subplot(2,1,1)
% set(gca,'XLim',[0 360],'XTick',0:60:360);
% xlabel('$\theta$ (elt deg)')
% ylabel('$i_{db}$ (Vs)')
% title(['Mean $i_{d3}$ = ' num2str(mean(idb(1,:))) ' A'])
% plot(th,idb(1,:));
% subplot(2,1,2)
% set(gca,'XLim',[0 360],'XTick',0:60:360);
% xlabel('$\theta$ (elt deg)')
% ylabel('$i_{qb}$ (Vs)')
% title(['Mean $i_{q3}$ = ' num2str(mean(iqb(1,:))) ' A'])
% plot(th,iqb(1,:));
% % if isoctave()
% %     fig_name=strcat(newDir, filemot(1:end-4), '_plot_flux');
% %     hgsave(hfig(2),[fig_name]);
% % else
% %     saveas(hfig(2),[newDir filemot(1:end-4) '_plot_flux']);
% % end

% if ~sum(isnan(fd))
    hfig(6) = figure();
    figSetting()
    subplot(2,1,1)
    set(gca,'XLim',[0 360],'XTick',0:60:360);
    xlabel('$\theta$ (elt deg)')
    ylabel('$\lambda_{abc}$ (Vs)')
    title(['Phase flux linkages'])
    for i=1:n3ph
        plot(th,fa(:,i));
        plot(th,fb(:,i));
        plot(th,fc(:,i));
    end
    subplot(2,1,2)
    set(gca,'XLim',[0 360],'XTick',0:60:360);
    xlabel('$\theta$ (elt deg)')
    ylabel('$i_{abc}$ (A)')
    title(['Phase currents'])
    for i=1:n3ph
        plot(th,ia(:,i),'DisplayName','$i_a$');
        plot(th,ib(:,i),'DisplayName','$i_b$');
        plot(th,ic(:,i),'DisplayName','$i_c$');
    end

    if length(ia(:,1))<2
        legend show
    end
    if isoctave()
        fig_name=strcat(newDir, filemot(1:end-4), '_plot_phase');
        hgsave(hfig(3),[fig_name]);
    else
        % saveas(hfig(3),[newDir filemot(1:end-4) '_plot_phase']);
    end
% end

