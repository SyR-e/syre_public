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

function [sVonMises,R,structModel] = calcVonMisesStress(structModel,flagPlot)

% sigma_max = mat.Rotor.sigma_max;
% SF = 0.8;

if nargin==1
    flagPlot=0;
end

R = solvepde(structModel);

[cgradx,cgrady,~] = evaluateCGradient(R);
sxx = cgradx(:,1);
sxy = cgradx(:,2);
syx = cgrady(:,1);
syy = cgrady(:,2);
sxz = 0;
syz = 0;
szx = 0;
szy = 0;
szz = 0;
sVonMises = sqrt( 0.5*( (sxx-syy).^2 + (syy -szz).^2 +(szz-sxx).^2) + 3*(sxy.^2 + syz.^2 + szx.^2));

sVM_max = max(sVonMises);

if flagPlot
    figure();
    figSetting();
    set(gca,'DataAspectRatio',[1 1 1]);
    pdeplot(structModel,'XYData',abs(R.NodalSolution(:,1)+j*R.NodalSolution(:,2)),'ZData',abs(R.NodalSolution(:,1)+j*R.NodalSolution(:,2)))
    colormap turbo
    view(2)
    title('Displacement [m]')


    figure();
    figSetting();
    set(gca,'DataAspectRatio',[1 1 1]);
    pdeplot(structModel,'XYData',sVonMises/1e6,'ZData',sVonMises/1e6)
    colormap turbo
    view(2)
    title('Von Mises Stress [MPa]')
    set(gca,'CLim',[0 sVM_max]/1e6)
end

% figure();
% figSetting();
% set(gca,'DataAspectRatio',[1 1 1]);
% pdeplot(structModel,'XYData',sVonMises/1e6,'ZData',sVonMises/1e6)
% colormap turbo
% view(2)
% title('Von Mises Stress [MPa]')
% % set(gca,'CLim',[0 sVM_max]/1e6)
% 
% 
% if sVM_max/1e6>sigma_max
%     set(gca,'CLim',[sigma_max sVM_max]/1e6)
% elseif sVM_max/1e6>sigma_max*SF
%     set(gca,'CLim',[sigma_max*SF sigma_max]/1e6);
% else
%     set(gca,'CLim',[0 sigma_max*SF]/1e6);
% end
    
    