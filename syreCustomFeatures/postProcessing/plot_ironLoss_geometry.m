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

function plot_ironLoss_geometry(geo,per,mat,out,pathname)

SOL = evalIronLossFEMM(geo,per,mat,out.SOL,2);

gifName = 'ironFluxDensity';
figName = 'localLoss';

hfig(1) = figure();
figSetting();
hax(1) = axes('OuterPosition',[0 0 1 1]);
axis equal
colormap turbo
title('Iron loss (W)')
set(hfig(1),'FileName',[pathname figName '.fig'])

hfig(2) = figure();
figSetting();
hax(2) = axes('OuterPosition',[0 0 1 1]);
axis equal
colormap turbo
zLim = max(abs(out.SOL.bs+out.SOL.br),[],'all');
set(gca,'CLim',[0 zLim]);
set(hfig(2),'FileName',[pathname gifName '.fig'])


for ii=1:length(hfig)
    GUI_Plot_Machine(hax(ii),[geo.stator;geo.rotor]);
    hchild = get(hax(ii),'Children');
    for jj=1:length(hchild)
        set(hchild(jj),'Color','k','LineWidth',1,'HandleVisibility','off');
    end
end

scatter(hax(1),real(SOL.pos),imag(SOL.pos),10,(SOL.psh+SOL.psc+SOL.prc+SOL.prh),'filled')
hgif = scatter(hax(2),real(SOL.pos),imag(SOL.pos),10,abs(SOL.bs(1,:)+SOL.br(1,:)),'filled');

colorbar(hax(1));
colorbar(hax(2));

savePrintFigure(hfig(1));

v = VideoWriter([pathname gifName '.avi']);
v.Quality = 95;
open(v);

for ii=1:length(SOL.th)
    set(hgif,'CData',abs(SOL.bs(ii,:)+SOL.br(ii,:)));
    title(['$\theta = ' int2str(SOL.th(ii)-SOL.th(1)) '^\circ$'])
    drawnow
    
    frame = getframe(hfig(2));
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if ii==1
        imwrite(imind,cm,[pathname gifName '.gif'],'gif','LoopCount',inf,'DelayTime',0.01);
    else
        imwrite(imind,cm,[pathname gifName '.gif'],'gif','WriteMode','append','DelayTime',0.01);
    end
    writeVideo(v,im);
end

close(v);










