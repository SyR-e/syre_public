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

function plot_singt(out,klength,kturns,delta_sim_singt,newDir,filemot)
% single working point has been simulated

t60 = out.SOL.T'* klength;
t60 = [t60;t60(1)];
t = repeat_n(t60',360/delta_sim_singt);

f_d = out.SOL.fd'*klength*kturns;
fd60=[f_d;f_d(1)];
fd=repeat_n(fd60',360/delta_sim_singt);
f_q = out.SOL.fq'*klength*kturns;
fq60=[f_q;f_q(1)];
fq=repeat_n(fq60',360/delta_sim_singt);
gamma = mean(atan2(-out.SOL.id',out.SOL.iq')) * 180/pi;

delta = atan2(fq,fd) * 180/pi;

IPF = cosd(delta-gamma);
th = linspace(0,360,length(t));

% plots
FontSize = 12;
FontName = 'TimesNewRoman';

figure()
subplot(2,1,1)
pl = plot(th,abs(t)); grid on
title(['Mean Torque = ' num2str(mean(t))]);
set(pl,'LineWidth',[2]);
xlim([0 360]), %ylim([ymin ymax]),
set(gca,'FontName',FontName,'FontSize',FontSize);
ti = 0:60:360; set(gca,'XTick',ti);
xl = xlabel('\theta - degrees'); set(xl,'Rotation',[0],'Fontsize',FontSize);
yl = ylabel('Nm');
set(yl,'Rotation',[90],'FontName',FontName,'Fontsize',FontSize);

%figure();
subplot(2,1,2)
pl = plot(th,IPF);
title(['Mean IPF = ' num2str(mean(IPF))]);
grid on
set(pl,'LineWidth',2);
xl = xlabel('\theta - degrees'); set(xl,'Rotation',0,'FontName',FontName,'Fontsize',FontSize);
yl=ylabel('IPF');
set(yl,'Rotation',90,'FontName',FontName,'Fontsize',FontSize);
xlim([0 360]);
set(gca,'FontName',FontName,'FontSize',FontSize);
ti = 0:60:360;
set(gca,'XTick',ti);

h=gcf();       %OCT
if isoctave()
    fig_name=strcat(newDir, filemot(1:end-4), '_T_gamma');
    hgsave(h,[fig_name]);
else
    saveas(h,[newDir filemot(1:end-4) '_T_gamma']);
end

ymin = round(( min(min(f_d),min(f_q)))*1000)/1000;
ymax = round(( max(max(f_d),max(f_q)))*1000)/1000;

if ymax>ymin
    ymax=1.2*ymax;
    ymin=0.8*ymin;
else
    ymax=1.2*ymin;
    ymin=0.8*ymax;
end

hdq=figure();
subplot(2,1,1);
hd1=plot(th,fd);
title(['Mean \lambda_d = ' num2str(mean(fd))]);
grid on
set(hd1,'LineWidth',2);
xl_hd=xlabel('\theta [Electrical degrees]');
set(xl_hd,'Rotation',0,'FontName',FontName,'Fontsize',FontSize); %,'FontWeight','Bold');
yl_hd=ylabel('\lambda_d [Wb]');
set(yl_hd,'Rotation',90,'FontName',FontName,'Fontsize',FontSize); %,'FontWeight','Bold');
xlim([0 360]);
set(gca,'FontSize',FontSize), %,'FontWeight','Bold');
ti = 0:60:360;
set(gca,'XTick',ti);

hq = subplot(2,1,2);
hq1 = plot(th,fq);
title(['Mean \lambda_q = ' num2str(mean(fq))]);
grid on
set(hq1,'LineWidth',2);
xl_hq=xlabel('\theta [Electrical degrees]');
set(xl_hq,'Rotation',0,'FontName',FontName,'Fontsize',FontSize); %'FontWeight','Bold');
yl_hq=ylabel('\lambda_q [Wb]');
set(yl_hq,'Rotation',90,'FontName',FontName,'Fontsize',FontSize); %,'FontWeight','Bold');
xlim([0 360]);
set(gca,'FontSize',FontSize); %,'FontWeight','Bold');
ti = 0:60:360;
set(gca,'XTick',ti);

if isoctave() %OCT
    fig_name=strcat(newDir, filemot(1:end-4), '_plot_Flux');
    hgsave(hdq,[fig_name]);
else
    saveas(hdq,[newDir filemot(1:end-4) '_plot_Flux']);
end
