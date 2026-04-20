% Copyright 2025
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

function plot_singm3D(F_map,resFolder,flagPlot)

if nargin()==2
    flagPlot=1;
end

% Interp the flux linkage maps over a very dense grid (41 x 41)

n_interp = 41;    % number of points in [0 Imax] for data interpolation
        
Id = linspace(min(F_map.Id,[],'all'),max(F_map.Id,[],'all'),n_interp);
Iq = linspace(min(F_map.Iq,[],'all'),max(F_map.Iq,[],'all'),n_interp);
Ir = linspace(min(F_map.If,[],'all'),max(F_map.If,[],'all'),n_interp);

[Id,Iq,Ir] = meshgrid(Id,Iq,Ir);

Fd = interp3(F_map.Id,F_map.Iq,F_map.If,F_map.Fd,Id,Iq,Ir,'spline');
Fq = interp3(F_map.Id,F_map.Iq,F_map.If,F_map.Fq,Id,Iq,Ir,'spline');
T  = interp3(F_map.Id,F_map.Iq,F_map.If,F_map.T,Id,Iq,Ir,'spline');

if isfield(F_map,'dT')
    dT  = interp3(F_map.Id,F_map.Iq,F_map.If,F_map.dT,Id,Iq,Ir,'spline');
end

if isfield(F_map,'dTpp')
    dTpp = interp3(F_map.Id,F_map.Iq,F_map.If,F_map.dTpp,Id,Iq,Ir,'spline');
end

if isfield(F_map,'We')
    We  = interp3(F_map.Id,F_map.Iq,F_map.If,F_map.We,Id,Iq,Ir,'spline');
    Wc  = interp3(F_map.Id,F_map.Iq,F_map.If,F_map.Wc,Id,Iq,Ir,'spline');
end

if isfield(F_map,'Pfes_h')
    Pfes_h = interp3(F_map.Id,F_map.Iq,F_map.If,F_map.Pfes_h,Id,Iq,Ir,'spline');
    Pfes_c = interp3(F_map.Id,F_map.Iq,F_map.If,F_map.Pfes_c,Id,Iq,Ir,'spline');
    Pfer_h = interp3(F_map.Id,F_map.Iq,F_map.If,F_map.Pfer_h,Id,Iq,Ir,'spline');
    Pfer_c = interp3(F_map.Id,F_map.Iq,F_map.If,F_map.Pfer_c,Id,Iq,Ir,'spline');
    velDim = F_map.velDim;
    Pfe    = interp3(F_map.Id,F_map.Iq,F_map.If,F_map.Pfe,Id,Iq,Ir,'spline');
    Ppm    = interp3(F_map.Id,F_map.Iq,F_map.If,F_map.Ppm,Id,Iq,Ir,'spline');
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
    if isfield(F_map,'We')
        save ('-mat7-binary', name_file,'We','Wc','-append');
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
    if isfield(F_map,'If')
        save ('-mat7-binary', name_file,'Ir','-append');
    end
    if isfield(F_map,'Fr')
        save ('-mat7-binary', name_file,'Fr','-append');
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
    if isfield(F_map,'We')
        save (nameFile,'We','Wc','-append');
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
    if isfield(F_map,'If')
        save (nameFile,'Ir','-append');
    end
    if isfield(F_map,'Fr')
        save (nameFile,'Fr','-append');
    end
    if exist('IM','var')
        save (nameFile,'IM','-append');
    end
end

%% Figures
if flagPlot
    FigDir = [resFolder,'fig - flux maps'];
    mkdir(FigDir);
    FigDir = checkPathSyntax([FigDir '\']);
    
    % flux maps
    figure
    plot(Id(1,:,end),Fd([1 end],:,end),F_map.Id(1,:,end),F_map.Fd([1 end],:,end),'kx'), grid on, hold on
    plot(Iq(:,1,end),Fq(:,[1 end],end),F_map.Iq(:,1,end),F_map.Fq(:,[1 end],end),'kx'),
    xlabel('id,iq [A]'), ylabel('\lambda_d, \lambda_q [Vs]'), title('Curves at Ir_{max}')%zlabel('\lambda_d')
    h=gcf(); %OCT
    if isoctave()
        fig_name=strcat(FigDir, 'Curves.fig');
        hgsave(h,[fig_name]);
        clear fig_name
    else
        saveas(gcf,[FigDir 'Curves.fig'])
    end
    
    figure
    surfc(Id(:,:,end),Iq(:,:,end),Fd(:,:,end)), grid on, xlabel('id'), ylabel('iq'), zlabel('\lambda_d (Ir_{max})')
    h=gcf(); %OCT
    if isoctave()
        fig_name = checkPathSyntax(strcat(FigDir, '\Fdsurf.fig'));
        hgsave(h,[fig_name]);
        clear fig_name
    else
        saveas(gcf,[FigDir 'Fdsurf.fig'])
    end
    
    figure
    surfc(Id(:,:,end),Iq(:,:,end),Fq(:,:,end)), grid on, xlabel('id'), ylabel('iq'), zlabel('\lambda_q (Ir_{max})')
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
    surf(Id(:,:,end),Iq(:,:,end),abs(T(:,:,end))), grid on, xlabel('id [A]'), ylabel('iq [A]'), zlabel('Torque [Nm]')
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
        surf(Id(:,:,end),Iq(:,:,end),dT(:,:,end)), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Torque ripple (std) [Nm]')
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
        surf(Id(:,:,end),Iq(:,:,end),dTpp(:,:,end)), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Torque ripple (pk-pk) [Nm]')
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
        surf(Id(:,:,end),Iq(:,:,end),Pfe(:,:,end)), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Iron Loss [W]')
        h=gcf(); %OCT
        if isoctave()
            fig_name=strcat(FigDir, 'Pfesurf.fig');
            hgsave(h,[fig_name]);
            clear fig_name
        else
            saveas(gcf,[FigDir 'Pfesurf.fig'])
        end
    end
    
    
    figure
    surf(Id(:,:,end),Iq(:,:,end),Ir(:,:,end)), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Ir [A]')
    h=gcf(); %OCT
    if isoctave()
        fig_name=strcat(FigDir, 'Irsurf.fig');
        hgsave(h,[fig_name]);
        clear fig_name
    else
        saveas(gcf,[FigDir 'Irsurf.fig'])
    end
end

