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

function plot_singt_FFT(filename,pathname)

%close all
if nargin()<2
    [filename, pathname, ~] = uigetfile([cd '\.mat'], 'LOAD DATA');
end

load([pathname filename]);
if exist('per','var')
    delta_sim_singt = per.delta_sim_singt;
else
    delta_sim_singt = 360;
end

if isstruct(out.SOL)
    t60 = abs(out.SOL.T);
else
    t60 = abs(out.SOL(:,6))';
end

t = [repmat(t60,1,360/delta_sim_singt) t60(1)];

f_d = out.SOL.fd;
fd60=f_d;

fd = [repmat(fd60,1,360/delta_sim_singt) fd60(1)];
f_q = out.SOL.fq;
fq60=f_q;

fq = [repmat(fq60,1,360/delta_sim_singt) fq60(1)];
gamma = mean(atan2(out.SOL.iq,out.SOL.id)) * 180/pi;

delta = atan2(fq,fd) * 180/pi;
IPF = cosd(delta-gamma);
th = linspace(0,360,length(t));

% creating folder
folder = 'fft';
[~,message,~] = mkdir(pathname,folder);
if not(isempty(message))
    disp('Warning : existing folder')
end

pathfolder=[pathname folder '\'];

% torque plot
[Tcont,Tharm] = plot_figure_FFT_singt(th,abs(t),50,'$\theta_e$ [$^\circ$]','$T$ [Nm]');
saveas(gcf,[pathfolder 'Torque.fig']);
save([pathfolder 'fft_out'],'out','Tcont','Tharm','th','t');
if exist('geo','var')
    save([pathfolder 'fft_out'],'geo','per','-append');
end

% flux plot
plot_figure_FFT_singt(th,abs(fd+j*fq),50,'$\theta_e$ [$^\circ$]','$|\lambda_{dq}|$ [Vs]')
saveas(gcf,[pathfolder 'Flux.fig']);

plot_figure_FFT_singt(th,fd,50,'$\theta_e$ [$^\circ$]','$\lambda_{d}$ [Vs]')
saveas(gcf,[pathfolder 'FluxD.fig']);

plot_figure_FFT_singt(th,fq,50,'$\theta_e$ [$^\circ$]','$\lambda_{q}$ [Vs]')
saveas(gcf,[pathfolder 'FluxQ.fig']);



