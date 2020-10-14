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

function plot_singm(F_map,NewDir,filename)

% Interp the flux linkage maps over a very dense grid (256 x 256)

n_interp = 256;    % number of points in [0 Imax] for data interpolation
klength = 1; kturns = 1; n2=n_interp;
        
id = F_map.Id; iq = F_map.Iq;
Fd = F_map.Fd; Fq = F_map.Fq;

T = F_map.T;
if isfield(F_map,'dT')
    dT = F_map.dT;
end
if isfield(F_map,'dTpp')
    dTpp = F_map.dTpp;
end
if isfield(F_map,'Pfe')
    Pfe = F_map.Pfe;
end
if isfield(F_map,'Pfes_h')
    Pfes_h = F_map.Pfes_h;
    Pfes_c = F_map.Pfes_c;
    Pfer_h = F_map.Pfer_h;
    Pfer_c = F_map.Pfer_c;
end
if isfield(F_map,'Ppm')
    Ppm = F_map.Ppm;
end

i_d=linspace(id(1),id(end),n2);
i_q=linspace(iq(1),iq(end),n2);

[Id,Iq]=meshgrid(i_d,i_q);

% refine maps to the 256 x 256 standard resolution
Fd = interp2(F_map.Id(1,:,1),F_map.Iq(:,1,1)',mean(Fd,3),Id,Iq,'spline')*klength;
Fq = interp2(F_map.Id(1,:,1),F_map.Iq(:,1,1)',mean(Fq,3),Id,Iq,'spline')*klength;
T = interp2(F_map.Id(1,:,1),F_map.Iq(:,1,1)',mean(T,3),Id,Iq,'spline')*klength;
if isfield(F_map,'dT')
    dT = interp2(F_map.Id(1,:),F_map.Iq(:,1)',mean(dT,3),Id,Iq,'spline');
    dT = dT*klength;
end
if isfield(F_map,'dTpp')
    dTpp = interp2(F_map.Id(1,:,1),F_map.Iq(:,1,1)',mean(dTpp,3),Id,Iq,'spline');
    dTpp = dTpp*klength;
end
if isfield(F_map,'Pfes_h')
    Pfes_h = interp2(F_map.Id(1,:,1),F_map.Iq(:,1,1)',mean(Pfes_h,3),Id,Iq,'spline')*klength;
    Pfes_c = interp2(F_map.Id(1,:,1),F_map.Iq(:,1,1)',mean(Pfes_c,3),Id,Iq,'spline')*klength;
    Pfer_h = interp2(F_map.Id(1,:,1),F_map.Iq(:,1,1)',mean(Pfer_h,3),Id,Iq,'spline')*klength;
    Pfer_c = interp2(F_map.Id(1,:,1),F_map.Iq(:,1,1)',mean(Pfer_c,3),Id,Iq,'spline')*klength;
    velDim = F_map.velDim;
    Pfe    = interp2(F_map.Id(1,:),F_map.Iq(:,1)',mean(Pfe,3),Id,Iq,'spline');
    Ppm    = interp2(F_map.Id(1,:,1),F_map.Iq(:,1,1)',mean(Ppm,3),Id,Iq,'spline')*klength;
end

%% rewind
Id=Id/kturns;
Iq=Iq/kturns;
Fd=Fd*kturns;
Fq=Fq*kturns;

% %% add end-connections term
% Fd = Fd + Lld * Id;
% Fq = Fq + Llq * Iq;

if isoctave()  %OCT
    name_file = strcat(NewDir, 'fdfq_idiq_n',num2str(n2),'.mat');
    save ('-mat7-binary', name_file,'Fd','Fq','Id','Iq');
    save ('-mat7-binary', name_file,'T','-append');
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
    save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'Fd','Fq','Id','Iq');
    save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'T','-append');
    if isfield(F_map,'dT')
        save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'dT','-append');
    end
    if isfield(F_map,'dTpp')
        save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'dTpp','-append');
    end
    if isfield(F_map,'Pfe')
        save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'Pfe','-append');
    end
    if isfield(F_map,'Pfes_h')
        save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'Pfes_h','Pfes_c','Pfer_h','Pfer_c','velDim','-append');
    end
    if isfield(F_map,'Ppm')
        save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'Ppm','-append');
    end
    if isfield(F_map,'speed')
        velDim = F_map.speed;
        save ([NewDir 'fdfq_idiq_n' num2str(n2) '.mat'],'velDim','-append');
    end
end

% Figures
FigDir=[NewDir,'fig - flux maps'];
[success,message,messageid] = mkdir(FigDir);
FigDir = [FigDir '\'];

% flux maps
figure
plot(Id(1,:),Fd([1 end],:),F_map.Id(1,:),F_map.Fd([1 end],:),'kx'), grid on, hold on
plot(Iq(:,1),Fq(:,[1 end]),F_map.Iq(:,1),F_map.Fq(:,[1 end]),'kx'),
xlabel('id,iq [A]'), ylabel('\lambda_d, \lambda_q [Vs]'), %zlabel('\lambda_d')
h=gcf(); %OCT
if isoctave()
    fig_name=strcat(FigDir, 'Curves_', strrep(filename,'.mat',''));
    hgsave(h,[fig_name]);
    clear fig_name
else
    saveas(gcf,[FigDir 'Curves_' strrep(filename,'.mat','.fig')])
end

figure
surfc(Id,Iq,Fd), grid on, xlabel('id'), ylabel('iq'), zlabel('\lambda_d')
if not(kturns == 1)
    title(['Rewind factor = ' num2str(kturns)])
end
h=gcf(); %OCT
if isoctave()
    fig_name=strcat(FigDir, '\Fdsurf_', strrep(filename,'.mat',''));
    hgsave(h,[fig_name]);
    clear fig_name
else
    saveas(gcf,[FigDir 'Fdsurf_' strrep(filename,'.mat','.fig')])
end

figure
surfc(Id,Iq,Fq), grid on, xlabel('id'), ylabel('iq'), zlabel('\lambda_q')
if not(kturns == 1)
    title(['Rewind factor = ' num2str(kturns)])
end
h=gcf(); %OCT
if isoctave()
    fig_name=strcat(FigDir, 'Fqsurf_', strrep(filename,'.mat',''));
    hgsave(h,[fig_name]);
    clear fig_name
else
    saveas(gcf,[FigDir 'Fqsurf_' strrep(filename,'.mat','.fig')])
end


% TORQUE MAP
figure
surf(Id,Iq,abs(T)), grid on, xlabel('id [A]'), ylabel('iq [A]'), zlabel('Torque [Nm]')
if not(kturns == 1)
    title(['Rewind factor = ' num2str(kturns)])
end
h=gcf(); %OCT
if isoctave()
    fig_name=strcat(FigDir, 'Tsurf_', strrep(filename,'.mat',''));
    hgsave(h,[fig_name]);
    clear fig_name
else
    saveas(gcf,[FigDir 'Tsurf_' strrep(filename,'.mat','.fig')])
end

if isfield(F_map,'dT')
    figure
    surf(Id,Iq,dT), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Torque ripple (std) [Nm]')
    if not(kturns == 1)
        title(['Rewind factor = ' num2str(kturns)])
    end
    h=gcf(); %OCT
    if isoctave()
        fig_name=strcat(FigDir, 'dTsurf_', strrep(filename,'.mat',''));
        hgsave(h,[fig_name]);
        clear fig_name
    else
        saveas(gcf,[FigDir 'dTsurf_' strrep(filename,'.mat','.fig')])
    end
end
if isfield(F_map,'dTpp')
    figure
    surf(Id,Iq,dTpp), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Torque ripple (pk-pk) [Nm]')
    if not(kturns == 1)
        title(['Rewind factor = ' num2str(kturns)])
    end
    h=gcf(); %OCT
    if isoctave()
        fig_name=strcat(FigDir, 'dTppsurf_', strrep(filename,'.mat',''));
        hgsave(h,[fig_name]);
        clear fig_name
    else
        saveas(gcf,[FigDir 'dTppsurf_' strrep(filename,'.mat','.fig')])
    end
    
end

if isfield(F_map,'Pfe')
    figure
    surf(Id,Iq,Pfe), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('Iron Loss [W]')
    if not(kturns == 1)
        title(['Rewind factor = ' num2str(kturns)])
    end
    h=gcf(); %OCT
    if isoctave()
        fig_name=strcat(FigDir, 'Pfesurf_', strrep(filename,'.mat',''));
        hgsave(h,[fig_name]);
        clear fig_name
    else
        saveas(gcf,[FigDir 'Pfesurf' strrep(filename,'.mat','.fig')])
    end
    
end
