% Copyright 2026
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

function [out,Ppm] = recompute_PMloss(geo0,per0,mat0,out0,nL,nT,debug)

if ((nargin()==0)||(nargin()==1))
    if nargin()==1
        pathname = geo0;
    else
        pathname = [cd '\'];
    end
    [fileName,pathname] = uigetfile([pathname '*.mat'], 'Select source .mat file');
    load(fullfile(pathname, fileName));
    geo0 = geo;
    per0 = per;
    mat0 = mat;
    out0 = out;
    clear geo per mat out;
    prompt = {...
        'Axial segments',...
        'Tangential segments'};
    dlgTitle = 'PM segmentation inputs';
    defaultAns = {'1','1'};
    answer = inputdlg(prompt, dlgTitle, 1, defaultAns);
    nL = round(str2double(answer{1}));
    nT = round(str2double(answer{1}));
    debug = 1;
end


if ~exist('debug','var')
    debug = 0;
end

% compute geometry with updated PMs segmentation
geo = geo0;
mat = mat0;
per = per0;
geo.PMNa      = nL;
geo.PMNc(2,1) = nT;
eval_type = 'singt';    % default 
fem = dimMesh(geo,eval_type);


[rotor,~,geo,mat] = ROTmatr(geo,fem,mat);

SOL = out0.SOL;

groNo = SOL.groNo;
groNo(groNo>200) = 200;
pos = SOL.pos;

rotorTmp = rotor(rotor(:,8)==6,:);
ids = unique(rotorTmp(:,9));
numPM = numel(ids);

for ii=1:numPM
    rotorTmp = rotor(rotor(:,9)==ids(ii),:);
    X = rotorTmp(:,1);
    Y = rotorTmp(:,2);
    ps = polyshape(X,Y);
    index = isinterior(ps,real(pos),imag(pos));
    groNo(index==1) = 200+ii;
end

out = out0;
out.SOL.groNo = groNo;

SOL = evalIronLossFEMM(geo,per,mat,out.SOL,2);

out.SOL           = SOL;
out.Ppm           = sum(sum(SOL.ppm))*(2*geo.p/geo.ps);
out.ppm_no3D      = sum(sum(SOL.ppm_no3D))*(2*geo.p/geo.ps);
out.ppm_noRFno3D  = sum(sum(SOL.ppm_noRFno3D))*(2*geo.p/geo.ps);
out.Ppm_breakdown = SOL.ppm_PM*(2*geo.p/geo.ps);

Ppm = out.Ppm;

flagPlot = ones(size(SOL.groNo));
flagPlot(SOL.groNo==22) = NaN;
flagPlot(SOL.groNo==12) = NaN;


if debug
    tmp = sum(out0.SOL.ppm);
    hfig(1) = figure();
    figSetting()
    GUI_Plot_Machine(gca,geo.rotor);
    GUI_Plot_Machine(gca,geo.stator);
    axis equal
    title('monolitic')
    scatter(real(SOL.pos(flagPlot==1)),imag(SOL.pos(flagPlot==1)),10,tmp(flagPlot==1),'filled');
    colormap turbo
    colorbar
    hax(1) = gca;
    cLim1 = get(gca,'CLim');


    tmp = sum(SOL.ppm);
    hfig(2) = figure();
    figSetting()
    GUI_Plot_Machine(gca,geo.rotor);
    GUI_Plot_Machine(gca,geo.stator);
    axis equal
    title('segmented')
    scatter(real(SOL.pos(flagPlot==1)),imag(SOL.pos(flagPlot==1)),10,tmp(flagPlot==1),'filled');
    colormap turbo
    colorbar
    hax(2) = gca;
    cLim2 = get(gca,'CLim');

    cLim = [0 max([cLim1(2) cLim2(2)])];
    set(hax(1),'CLim',cLim);
    set(hax(2),'CLim',cLim);

    hfig(3) = figure();
    figSetting()
    bar([1 2],[out0.Ppm out.Ppm],'BarWidth',0.8);
    set(gca,'XLim',[0.5 2.5],'XTick',[1 2],'XTickLabel',{'monolitic','segmented'})
    ylabel('$P_{PM}$ (W)')
end

if nargout()==0
    % answer = questdlg('Load results on WorkSpace','Save','Yes','No','Yes');
    % if strcmp(answer,'Yes')
    %     assignin('base','out',out);
    %     assignin('base','Ppm',Ppm);
    % end
    answer = questdlg('Save results','Save','Yes','No','Yes');
    if strcmp(answer,'Yes')
        resFolder = checkPathSyntax([pathname 'PMseg_' int2str(nL) 'axial_' int2str(nT) 'tangential\']);
        mkdir(resFolder)
        set(hfig(1),'FileName',[resFolder 'PMloss_monolitic.fig']);
        set(hfig(2),'FileName',[resFolder 'PMloss_segmented.fig']);
        set(hfig(3),'FileName',[resFolder 'PMloss_comparison.fig']);
        for ii=1:3
            savePrintFigure(hfig(ii));
        end
        save([resFolder 'out.mat'],'out','geo0','per0','mat0','out0');
    end
    clear out Ppm
end






