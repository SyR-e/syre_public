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

function plot_singm(F_map,resFolder)

% Interp the flux linkage maps over a very dense grid (256 x 256)

n_interp = 256;    % number of points in [0 Imax] for data interpolation
        
Id = linspace(min(F_map.Id,[],'all'),max(F_map.Id,[],'all'),n_interp);
Iq = linspace(min(F_map.Iq,[],'all'),max(F_map.Iq,[],'all'),n_interp);

[Id,Iq] = meshgrid(Id,Iq);

% refine maps to the 256 x 256 standard resolution
Fd = interp2(F_map.Id,F_map.Iq,F_map.Fd,Id,Iq,'spline');
Fq = interp2(F_map.Id,F_map.Iq,F_map.Fq,Id,Iq,'spline');
T  = interp2(F_map.Id,F_map.Iq,F_map.T,Id,Iq,'spline');

if isfield(F_map,'dT')
    dT = interp2(F_map.Id,F_map.Iq,F_map.dT,Id,Iq,'spline');
end

if isfield(F_map,'dTpp')
    dTpp = interp2(F_map.Id,F_map.Iq,F_map.dTpp,Id,Iq,'spline');
end

if isfield(F_map,'Pfes_h')
    Pfes_h = interp2(F_map.Id,F_map.Iq,F_map.Pfes_h,Id,Iq,'spline');
    Pfes_c = interp2(F_map.Id,F_map.Iq,F_map.Pfes_c,Id,Iq,'spline');
    Pfer_h = interp2(F_map.Id,F_map.Iq,F_map.Pfer_h,Id,Iq,'spline');
    Pfer_c = interp2(F_map.Id,F_map.Iq,F_map.Pfer_c,Id,Iq,'spline');
    velDim = F_map.velDim;
    Pfe    = interp2(F_map.Id,F_map.Iq,F_map.Pfe,Id,Iq,'spline');
    Ppm    = interp2(F_map.Id,F_map.Iq,F_map.Ppm,Id,Iq,'spline');
end

if isfield(F_map,'IM')
    F_IM     = F_map.IM;
    IM.Ir    = interp2(F_map.Id,F_map.Iq,F_IM.Ir,Id,Iq,'spline');
    IM.kr    = interp2(F_map.Id,F_map.Iq,F_IM.kr,Id,Iq,'spline');
    IM.ks    = interp2(F_map.Id,F_map.Iq,F_IM.ks,Id,Iq,'spline');
    IM.Rs    = interp2(F_map.Id,F_map.Iq,F_IM.Rs,Id,Iq,'spline');
    IM.Rr    = interp2(F_map.Id,F_map.Iq,F_IM.Rr,Id,Iq,'spline');
    IM.Ls    = interp2(F_map.Id,F_map.Iq,F_IM.Ls,Id,Iq,'spline');
    IM.Lr    = interp2(F_map.Id,F_map.Iq,F_IM.Lr,Id,Iq,'spline');
    IM.Lm    = interp2(F_map.Id,F_map.Iq,F_IM.Lm,Id,Iq,'spline');
    IM.Lt    = interp2(F_map.Id,F_map.Iq,F_IM.Lt,Id,Iq,'spline');
    IM.Lls   = interp2(F_map.Id,F_map.Iq,F_IM.Lls,Id,Iq,'spline');
    IM.LLr   = interp2(F_map.Id,F_map.Iq,F_IM.Llr,Id,Iq,'spline');
    IM.sigma = interp2(F_map.Id,F_map.Iq,F_IM.sigma,Id,Iq,'spline');
    IM.tr    = interp2(F_map.Id,F_map.Iq,F_IM.tr,Id,Iq,'spline');
    IM.wslip = interp2(F_map.Id,F_map.Iq,F_IM.wslip,Id,Iq,'spline');
end


% save data
if isoctave()  %OCT
    name_file = strcat(resFolder, 'fdfq_idiq_n',num2str(n_interp),'.mat');
    save ('-mat7-binary', name_file,'Fd','Fq','Id','Iq','T');
    if isfield(F_map,'dT')
        save ('-mat7-binary', name_file,'T','-append');
    end
    if isfield(F_map,'dTpp')
        save ('-mat7-binary', name_file,'dTpp','-append');
    end
    if isfield(F_map,'Pfe')
        save ('-mat7-binary', name_file,'Pfe','-append');
    end
    if isfield(F_map,'Pfes_h')
        save ('-mat7-binary', name_file,'Pfes_h','Pfes_c','Pfer_h','Pfer_c','-append');
    end
    if isfield(F_map,'Ppm')
        save ('-mat7-binary', name_file,'Ppm','-append');
    end
    clear name_file
else
    nameFile = [resFolder 'fdfq_idiq_n' int2str(n_interp) '.mat'];
    save (nameFile,'Fd','Fq','Id','Iq','T');
    if isfield(F_map,'dT')
        save (nameFile,'dT','-append');
    end
    if isfield(F_map,'dTpp')
        save (nameFile,'dTpp','-append');
    end
    if isfield(F_map,'Pfe')
        save (nameFile,'Pfe','-append');
    end
    if isfield(F_map,'Pfes_h')
        save (nameFile,'Pfes_h','Pfes_c','Pfer_h','Pfer_c','velDim','-append');
    end
    if isfield(F_map,'Ppm')
        save (nameFile,'Ppm','-append');
    end
    if isfield(F_map,'speed')
        velDim = F_map.speed;
        save (nameFile,'velDim','-append');
    end
    if exist('IM','var')
        save (nameFile,'IM','-append');
    end
end

%% Figures
FigDir=[resFolder,'fig - flux maps'];
mkdir(FigDir);
FigDir = [FigDir '\'];

% flux maps
figure
plot(Id(1,:),Fd([1 end],:),F_map.Id(1,:),F_map.Fd([1 end],:),'kx'), grid on, hold on
plot(Iq(:,1),Fq(:,[1 end]),F_map.Iq(:,1),F_map.Fq(:,[1 end]),'kx'),
xlabel('id,iq [A]'), ylabel('\lambda_d, \lambda_q [Vs]'), %zlabel('\lambda_d')
h=gcf(); %OCT
if isoctave()
    fig_name=strcat(FigDir, 'Curves.fig');
    hgsave(h,[fig_name]);
    clear fig_name
else
    saveas(gcf,[FigDir 'Curves.fig'])
end

figure
surfc(Id,Iq,Fd), grid on, xlabel('id'), ylabel('iq'), zlabel('\lambda_d')
h=gcf(); %OCT
if isoctave()
    fig_name=strcat(FigDir, '\Fdsurf.fig');
    hgsave(h,[fig_name]);
    clear fig_name
else
    saveas(gcf,[FigDir 'Fdsurf.fig'])
end

figure
surfc(Id,Iq,Fq), grid on, xlabel('id'), ylabel('iq'), zlabel('\lambda_q')
h=gcf(); %OCT
if isoctave()
    fig_name=strcat(FigDir, 'Fqsurf.fig');
    hgsave(h,[fig_name]);
    clear fig_name
else
    saveas(gcf,[FigDir 'Fqsurf.fig'])
end


% TORQUE MAP
figure
surf(Id,Iq,abs(T)), grid on, xlabel('id [A]'), ylabel('iq [A]'), zlabel('Torque [Nm]')
h=gcf(); %OCT
if isoctave()
    fig_name=strcat(FigDir, 'Tsurf.fig');
    hgsave(h,[fig_name]);
    clear fig_name
else
    saveas(gcf,[FigDir 'Tsurf.fig'])
end

if isfield(F_map,'dT')
    figure
    surf(Id,Iq,dT), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Torque ripple (std) [Nm]')
    h=gcf(); %OCT
    if isoctave()
        fig_name=strcat(FigDir, 'dTsurf.fig');
        hgsave(h,[fig_name]);
        clear fig_name
    else
        saveas(gcf,[FigDir 'dTsurf.fig'])
    end
end
if isfield(F_map,'dTpp')
    figure
    surf(Id,Iq,dTpp), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Torque ripple (pk-pk) [Nm]')
    h=gcf(); %OCT
    if isoctave()
        fig_name=strcat(FigDir, 'dTppsurf.fig');
        hgsave(h,[fig_name]);
        clear fig_name
    else
        saveas(gcf,[FigDir 'dTppsurf.fig'])
    end
    
end

if isfield(F_map,'Pfe')
    figure
    surf(Id,Iq,Pfe), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Iron Loss [W]')
    h=gcf(); %OCT
    if isoctave()
        fig_name=strcat(FigDir, 'Pfesurf.fig');
        hgsave(h,[fig_name]);
        clear fig_name
    else
        saveas(gcf,[FigDir 'Pfesurf.fig'])
    end
end

if isfield(F_map,'IM')
    figure
    surf(Id,Iq,IM.Ir), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('I_r [A]')
    h=gcf(); %OCT
    if isoctave()
        fig_name=strcat(FigDir, 'Irsurf.fig');
        hgsave(h,[fig_name]);
        clear fig_name
    else
        saveas(gcf,[FigDir 'Irsurf.fig'])
    end
end
