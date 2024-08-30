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

function plot_singt(out,delta_sim_singt,newDir,filemot)
% single working point has been simulated

if delta_sim_singt<=360
    nRep = 360/delta_sim_singt; % number of repetition needed
else
    nRep = 1;
end

T  = [repmat(out.SOL.T,1,nRep) out.SOL.T(1)];    % last point added for plot
fd = [repmat(out.SOL.fd,1,nRep) out.SOL.fd(1)];  % last point added for plot
fq = [repmat(out.SOL.fq,1,nRep) out.SOL.fq(1)];  % last point added for plot
id = [repmat(out.SOL.id,1,nRep) out.SOL.id(1)];  % last point added for plot
iq = [repmat(out.SOL.iq,1,nRep) out.SOL.iq(1)];  % last point added for plot

if isfield(out.SOL,'ia')
    if nRep==1
        ia = [out.SOL.ia out.SOL.ia(:,1)];
        ib = [out.SOL.ib out.SOL.ib(:,1)];
        ic = [out.SOL.ic out.SOL.ic(:,1)];

        fa = [out.SOL.fa out.SOL.fa(:,1)];
        fb = [out.SOL.fb out.SOL.fb(:,1)];
        fc = [out.SOL.fc out.SOL.fc(:,1)];
    else
        iph = phaseQuantityDecoding(out.SOL.ia,out.SOL.ib,out.SOL.ic,delta_sim_singt);
        ia = [iph.a iph.a(:,1)];
        ib = [iph.b iph.b(:,1)];
        ic = [iph.c iph.c(:,1)];

        fph = phaseQuantityDecoding(out.SOL.fa,out.SOL.fb,out.SOL.fc,delta_sim_singt);
        fa = [fph.a fph.a(:,1)];
        fb = [fph.b fph.b(:,1)];
        fc = [fph.c fph.c(:,1)];
    end
else
    ia = NaN;
    ib = NaN;
    ic = NaN;

    fa = NaN;
    fb = NaN;
    fc = NaN;
end

% th = linspace(0,360,length(T));
if length(out.SOL.th)>1
    dth = out.SOL.th(2)-out.SOL.th(1);
    th = 0:dth:dth*(length(T)-1);
else
    th = 0:60:360;
end

gamma = atan2(iq,id);
delta = atan2(fq,fd);
IPF = sin(gamma-delta);

hfig(1) = figure();
figSetting()
subplot(2,1,1)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ (elt deg)')
ylabel('(Nm)')
title(['Mean Torque = ' num2str(out.T) ' Nm'])
plot(th,T);
subplot(2,1,2)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ (elt deg)')
ylabel('IPF')
title(['Mean IPF = ' num2str(out.IPF)])
plot(th,IPF);
if isoctave()
    fig_name=strcat(newDir, filemot(1:end-4), '_T_gamma');
    hgsave(hfig(1),[fig_name]);
else
    saveas(hfig(1),[newDir filemot(1:end-4) '_T_gamma']);
end

hfig(2) = figure();
figSetting()
subplot(2,1,1)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ (elt deg)')
ylabel('$\lambda_d$ (Vs)')
title(['Mean $\lambda_d$ = ' num2str(out.fd) ' Vs'])
plot(th,fd);
subplot(2,1,2)
set(gca,'XLim',[0 360],'XTick',0:60:360);
xlabel('$\theta$ (elt deg)')
ylabel('$\lambda_q$ (Vs)')
title(['Mean $\lambda_q$ = ' num2str(out.fq) ' Vs'])
plot(th,fq);
if isoctave()
    fig_name=strcat(newDir, filemot(1:end-4), '_plot_flux');
    hgsave(hfig(2),[fig_name]);
else
    saveas(hfig(2),[newDir filemot(1:end-4) '_plot_flux']);
end

if ~sum(isnan(fa))
    hfig(3) = figure();
    figSetting()
    subplot(2,1,1)
    set(gca,'XLim',[0 360],'XTick',0:60:360);
    xlabel('$\theta$ (elt deg)')
    ylabel('$\lambda_{abc}$ (Vs)')
    title(['Phase flux linkages'])
    plot(th,fa);
    plot(th,fb);
    plot(th,fc);
    subplot(2,1,2)
    set(gca,'XLim',[0 360],'XTick',0:60:360);
    xlabel('$\theta$ (elt deg)')
    ylabel('$i_{abc}$ (A)')
    title(['Phase currents'])
    plot(th,ia,'DisplayName','$i_a$');
    plot(th,ib,'DisplayName','$i_b$');
    plot(th,ic,'DisplayName','$i_c$');
    if length(ia(:,1))<2
        legend show
    end
    if isoctave()
        fig_name=strcat(newDir, filemot(1:end-4), '_plot_phase');
        hgsave(hfig(3),[fig_name]);
    else
        saveas(hfig(3),[newDir filemot(1:end-4) '_plot_phase']);
    end
end
