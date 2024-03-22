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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_singt_vsTime(filename,pathname)

%close all
if nargin()<2
    [filename, pathname, ~] = uigetfile([cd '\.mat'], 'LOAD DATA');
end

load([pathname filename]);

if exist('per','var')
    n               = per.EvalSpeed;
    p               = geo.p;
    Rs              = per.Rs;
    delta_sim_singt = per.delta_sim_singt;
else
    n               = motorModel.WaveformSetup.EvalSpeed;
    p               = motorModel.data.p;
    Rs              = motorModel.data.Rs;
    delta_sim_singt = 360;
end

w   = n*pi/30*p;
if w==0
    n=3000;
    w   = n*pi/30*p;
    warning(['speed equal to zero, set to ' int2str(w) ' rad/s'])
end

nRep = 360/delta_sim_singt; % number of repetition needed

T  = [repmat(out.SOL.T,1,nRep)];   % last point added for plot
fd = [repmat(out.SOL.fd,1,nRep)];  % last point added for plot
fq = [repmat(out.SOL.fq,1,nRep)];  % last point added for plot
id = [repmat(out.SOL.id,1,nRep)];  % last point added for plot
iq = [repmat(out.SOL.iq,1,nRep)];  % last point added for plot

iph = phaseQuantityDecoding(out.SOL.ia,out.SOL.ib,out.SOL.ic,delta_sim_singt);
ia = [iph.a];
ib = [iph.b];
ic = [iph.c];

fph = phaseQuantityDecoding(out.SOL.fa,out.SOL.fb,out.SOL.fc,delta_sim_singt);
fa = [fph.a];
fb = [fph.b];
fc = [fph.c];

th = linspace(0,360,length(T)+1)+out.SOL.th(1);
th = th(1:end-1);
dt = 2*pi/length(th)/w;    %[s]
time = [0:1:length(th)-1]*dt; %[s]

tmpTime = [-dt time time(end)+dt];

va = diff([fa(:,end) fa fa(:,1)]')'/dt;
vb = diff([fb(:,end) fb fb(:,1)]')'/dt;
vc = diff([fc(:,end) fc fc(:,1)]')'/dt;

vd = diff([fd(:,end) fd fd(:,1)]')'/dt-w*[fq fq(:,1)];
vq = diff([fq(:,end) fq fq(:,1)]')'/dt+w*[fd fd(:,1)];

tmpTime = tmpTime(1:end-1)+dt/2;

va = interp1(tmpTime,va,time)+Rs*ia;
vb = interp1(tmpTime,vb,time)+Rs*ib;
vc = interp1(tmpTime,vc,time)+Rs*ic;

vd = interp1(tmpTime,vd,time)+Rs*id;
vq = interp1(tmpTime,vq,time)+Rs*iq;

wf.th_e = th;
wf.th_m = (th-th(1))/p;
wf.t    = time;
wf.id   = id;
wf.iq   = iq;
wf.fd   = fd;
wf.fq   = fq;
wf.vd   = vd;
wf.vq   = vq;
wf.T    = T;
wf.n    = n;
wf.ia   = ia;
wf.ib   = ib;
wf.ic   = ic;
wf.fa   = fa;
wf.fb   = fb;
wf.fc   = fc;
wf.va   = va;
wf.vb   = vb;
wf.vc   = vc;

hfig(1) = figure();
figSetting(12,18)
set(hfig(1),'FileName',[pathname 'timeWaveform_phase.fig']);
hax(1) = axes('OuterPosition',[0,2/3,1,1/3]);
hax(2) = axes('OuterPosition',[0 1/3 1 1/3]);
hax(3) = axes('OuterPosition',[0 0/3 1 1/3]);
for ii=1:3
    xlabel(hax(ii),'$t$ [s]');
    set(hax(ii),'XLim',[0 time(end)+dt]);
    plot(hax(ii),[0 time(end)+dt],[0 0],'-k','lineWidth',1,'HandleVisibility','off');
    set(hax(ii),'ColorOrderIndex',1)
end
ylabel(hax(1),'$i_{abc}$ [A]')
ylabel(hax(2),'$\lambda_{abc}$ [Vs]')
ylabel(hax(3),'$v_{abc}$ [V]')

plot(hax(1),time,ia);
plot(hax(1),time,ib);
plot(hax(1),time,ic);

plot(hax(2),time,fa);
plot(hax(2),time,fb);
plot(hax(2),time,fc);

plot(hax(3),time,va);
plot(hax(3),time,vb);
plot(hax(3),time,vc);

clear hax

hfig(2) = figure();
figSetting(12,18);
set(hfig(2),'FileName',[pathname 'timeWaveform_dq.fig']);
hax(1) = axes('OuterPosition',[0,2/3,1,1/3]);
hax(2) = axes('OuterPosition',[0 1/3 1 1/3]);
hax(3) = axes('OuterPosition',[0 0/3 1 1/3]);
for ii=1:3
    xlabel(hax(ii),'$t$ [s]');
    set(hax(ii),'XLim',[0 time(end)+dt]);
    plot(hax(ii),[0 time(end)+dt],[0 0],'-k','lineWidth',1,'HandleVisibility','off');
    legend(hax(ii),'show','Location','best');
    set(hax(ii),'ColorOrderIndex',1)
end
ylabel(hax(1),'$i_{dq}$ [A]')
ylabel(hax(2),'$\lambda_{dq}$ [Vs]')
ylabel(hax(3),'$v_{dq}$ [V]')

plot(hax(1),time,id,'DisplayName','$i_d$');
plot(hax(1),time,iq,'DisplayName','$i_q$');

plot(hax(2),time,fd,'DisplayName','$\lambda_d$');
plot(hax(2),time,fq,'DisplayName','$\lambda_q$');

plot(hax(3),time,vd,'DisplayName','$v_d$');
plot(hax(3),time,vq,'DisplayName','$v_q$');

hfig(3) = figure();
figSetting(12,12);
set(hfig(3),'FileName',[pathname 'timeWaveform_mech.fig']);
hax(1) = axes('OuterPosition',[0 1/2 1 1/2]);
hax(2) = axes('OuterPosition',[0 0/2 1 1/2]);
for ii=1:2
    xlabel(hax(ii),'$t$ [s]');
    set(hax(ii),'XLim',[0 time(end)+dt]);
    plot(hax(ii),[0 time(end)+dt],[0 0],'-k','lineWidth',1,'HandleVisibility','off');
    set(hax(ii),'ColorOrderIndex',1)
end
ylabel(hax(1),'$T$ [Nm]')
ylabel(hax(2),'$\theta_m$ [$^\circ$]')

plot(hax(1),time,T);

plot(hax(2),time,wf.th_m);

for ii=1:length(hfig)
    savePrintFigure(hfig(ii));
end

save([pathname 'wafevormVStime.mat'],'wf');


