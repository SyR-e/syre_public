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

function eval_vonMisesStress(dataIn)

% Run Structural Simulation in Matlab, at the considered rotor speed.

custom = dataIn.custom;

pathname = dataIn.currentpathname;
filename = dataIn.currentfilename;
load([pathname filename])

dataSet.EvalSpeed = dataIn.EvalSpeed;

materialCodes;

evalSpeed = dataIn.EvalSpeed;
if evalSpeed==0
    error('Please set a speed higher than zero!!!')
end

clc

resFolder = [filename(1:end-4) '_results\FEA results\'];
if ~exist([pathname resFolder],'dir')
    mkdir([pathname resFolder]);
end

newDir = ['structural_' int2str(evalSpeed) 'rpm\'];

newDir = [pathname resFolder newDir];
mkdir(newDir);

simSetup.filename  = filename;
simSetup.pathname  = pathname;
simSetup.evalSpeed = evalSpeed;
simSetup.meshSize  = 'fine';    % fine or coarse
simSetup.flagFull  = 0;         % 0-->Qs simulation / 1-->full motor simulation
simSetup.shaftBC   = 1;         % 1-->locked shaft / 0-->free shaft

warning('off')
disp(['Creation of the PDE model...'])
tic

% if (custom)
    [structModel,data4GeoMat] = femm2pde(geo,mat,simSetup);
% else
%     [structModel] = syre2pde(geo,mat,simSetup);
% end

tEnd = toc();
warning('on')
disp(['PDE model created in ' num2str(tEnd) ' s'])
save([newDir filename(1:end-4) '_structModel.mat'],'structModel','dataSet','geo','per','mat','data4GeoMat');


hfig = figure();
figSetting();
pdeplot(structModel,'Mesh','on');
saveas(gcf,[newDir 'StructMesh.fig']);

warning('off')
disp('Solving the PDE model...')
tic
[sVonMises,R,structModel] = calcVonMisesStress(structModel);
tEnd = toc();
warning('on')
disp(['PDE model solved in ' num2str(tEnd) ' s'])
save([newDir filename(1:end-4) '_structModel.mat'],'structModel','sVonMises','R','dataSet','geo','per','mat');

[out] = eval_maxStress(structModel,sVonMises,geo,mat);

save([newDir filename(1:end-4) '_structModel.mat'],'out','-append');

figure();
figSetting();
set(gca,'DataAspectRatio',[1 1 1]);
pdeplot(structModel,'XYData',abs(R.NodalSolution(:,1)+j*R.NodalSolution(:,2)),'ZData',abs(R.NodalSolution(:,1)+j*R.NodalSolution(:,2)))
colormap jet
xlabel('[m]')
ylabel('[m]')
view(2)
title('Displacement [m]')
saveas(gcf,[newDir 'Displacement.fig'])


figure();
figSetting();
set(gca,'DataAspectRatio',[1 1 1]);
tmp.ux = R.NodalSolution(:,1);
tmp.uy = R.NodalSolution(:,2);
pdeplot(structModel,'XYData',sVonMises/1e6,'ZData',sVonMises/1e6);
colormap jet
xlabel('[m]')
ylabel('[m]')
view(2)
title('Von Mises Stress [MPa]')
%legend ('Location','northwest')
set(gca,'CLim',[0 max(sVonMises)]/1e6)
saveas(gcf,[newDir 'VonMisesStress.fig'])

figure();
figSetting();
set(gca,'DataAspectRatio',[1 1 1]);
tmp.ux = R.NodalSolution(:,1);
tmp.uy = R.NodalSolution(:,2);
pdeplot(structModel,'XYData',sVonMises/1e6,'ZData',sVonMises/1e6,'Deformation',tmp,'DeformationScaleFactor',100);
colormap jet
xlabel('[m]')
ylabel('[m]')
view(2)
title('Von Mises Stress [MPa] - deformation scale=100')
%legend ('Location','northwest')
set(gca,'CLim',[0 max(sVonMises)]/1e6)
saveas(gcf,[newDir 'VonMisesStressDeformation.fig'])

if ~isempty(out.x_over)
    figure()
    figSetting()
    hax=axes('OuterPosition',[0 0 1 1]);
    set(gca,'DataAspectRatio',[1 1 1]);
    set(gca,'XLim',[geo.r-geo.r*sin(pi/geo.p)-1 geo.r+1]/1e3,'YLim',[-1 geo.r*sin(pi/geo.p)+1]/1e3);
    xlabel('[m]')
    ylabel('[m]')
    %set(gca,'XTick',[],'YTick',[]);
    title(['Von Mises Stress Limit Exceeded'])
    
    x_over = out.x_over;
    y_over = out.y_over;
    x_max = out.x_max;
    y_max = out.y_max;
    
    rotorplot = geo.rotor;
    rotorplot(:,1:6) = rotorplot(:,1:6)/1e3;
    tmp = rotorplot(:,8);
    index = 1:1:numel(tmp);
    index = index(tmp~=codMatShaft);
    rotorplot = rotorplot(index,:);
    GUI_Plot_Machine(hax,rotorplot);
    
    plot(x_over/1e3,y_over/1e3,'r.','MarkerSize',5,'DisplayName', 'Stress Exceeded');
    plot(x_max/1e3,y_max/1e3,'bo','MarkerSize',5.5,'DisplayName', 'Max Stress')
    grid off
    %legend ('Location','northwest')
    saveas(gcf,[newDir 'VonMisesStressLimit.fig'])
end


disp([num2str(out.nodesOver),' nodes out of ',num2str(length(structModel.Mesh.Nodes)), ' exceed the maximum material stress.'])









