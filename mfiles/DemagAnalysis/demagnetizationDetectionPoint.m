% Copyright 2018
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
%    See the License for the spcific language governing permissions and
%    limitations under the License.

function [DemagOut]=demagnetizationDetectionPoint(setup)
%
% [DemagOut]=demagnetizationDetectionPoint(setup)
%
% Version 2 of the demag detector. This function run multiple simulations
% at different current level, but at the same temperature. The results are
% the percentage of demagnetized point on each PMs.

if nargin<1
    
    [filename,pathname,~]=uigetfile([cd '\.mat'],'Select a PM machine');
    load([pathname filename]);
    
    i0 = per.i0;
    
    prompt={'Maximum current [A]',...
        'Number of current levels',...
        'PM temperature',...
        'plot figures? (yes=1/no=0)'};
    name='Demagnetization input';
    defaultanswer={ num2str(i0),...
        num2str(dataSet.NumGrid),...
        num2str(20),...
        '1'};
    answer=inputdlg(prompt,name,1,defaultanswer);
    
    setup.filename = filename;
    setup.pathname = pathname;
    setup.geo      = geo;
    setup.per      = per;
    setup.mat      = mat;
    setup.per.overload = linspace(0,eval(answer{1}),eval(answer{2}))/i0;
    setup.per.tempPP   = eval(answer{3});
%     setup.Imax     = eval(answer{1});
%     setup.ncur     = eval(answer{2});
%    setup.temp     = eval(answer{3});
    setup.figFlag  = eval(answer{4});
    
end

%load([setup.pathname setup.filename]);

geo     = setup.geo;
per     = setup.per;
mat     = setup.mat;
dataSet = setup.dataSet;

if strcmp(mat.LayerMag.MatName,'Air')
    error('Select a PM machine!');
elseif ~isfield(mat.LayerMag,'temp') 
    error('Please select a PMs with thermal properties or update the PM properties');
end

filename = setup.filename;
pathname = setup.pathname;
i0 = per.i0;
% xdeg     = 360/(6*geo.p*geo.q)*geo.p;
% nsim     = 2;
Imax     = max(per.overload*i0);
ncur     = length(per.overload);
Ivect    = per.overload*i0;
figFlag  = setup.figFlag;

if per.tempPP>max(mat.LayerMag.temp.temp)
    error('Input temperature higher than the maximum allowed!')
end

Br = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,per.tempPP);
Bd = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Bd,per.tempPP);

resFolder=[filename(1:end-4) '_DemagPoint_' int2str(per.tempPP) 'degC\'];
if ~isfolder([pathname resFolder])
    mkdir(pathname,resFolder);
end
resFolder=[pathname resFolder];
save([resFolder 'setupData.mat'],'setup');

motFolder=[resFolder 'critical machines\'];
if ~isfolder(motFolder)
    mkdir(motFolder);
end

%Ivect=linspace(Imax/2,Imax,ncur);

if figFlag
    
    hfig{1}=figure();
    figSetting(16,10);
    hax{1}=axes('OuterPosition',[0 0 1 1]);
    xlabel('$I_{demag}$ [A]')
    ylabel('$B_{PMs}$ [T]')
    title(['$\theta_{PM}=' num2str(per.tempPP) '^\circ C$'])
    set(hax{1},'XLim',[0 Imax]);
    set(hax{1},'YLim',[0 Br*1.1]);
    plot(hax{1},[0 Imax],Br*[1,1],'g','DisplayName',['$B_r=' num2str(Br) '$ [T]']);
    plot(hax{1},[0 Imax],Bd*[1,1],'r','DisplayName',['$B_d=' num2str(Bd) '$ [T]']);
    set(hax{1},'ColorOrderIndex',1);
    hleg=legend('show','Location','northeastoutside');
    
    hfig{2}=figure();
    figSetting(16,10);
    hax{2}=axes('OuterPosition',[0 0 1 1]);
    xlabel('$I_{demag}$ [A]')
    ylabel('$PM_{demag}$ [pu]')
    title(['$\theta_{PM}=' num2str(per.tempPP) '^\circ C$'])
    set(hax{2},'XLim',[0 Imax]);
    set(hax{2},'YLim',[-0.1 1.1]);
    plot(hax{2},[0 Imax],[1,1],'r','HandleVisibility','off');
    plot(hax{2},[0 Imax],[0 0],'g','HandleVisibility','off');
    set(hax{2},'ColorOrderIndex',1);
    hleg=legend('show','Location','northeastoutside');
    
    for ii=1:2
        hplot{1,ii}=plot(hax{1},0,0,'-o','DisplayName',['pos ' int2str(ii)]);
        hplot{2,ii}=plot(hax{2},0,0,'-o','DisplayName',['pos ' int2str(ii)]);
    end
    drawnow
end

disp('Starting FEMM simulations...')

openfemm(1)
opendocument([pathname filename(1:end-4) '.fem']);

% if (strcmp(eval_type,'idemag') || strcmp(eval_type,'idemagmap'))
    % the PMs group must be changed, from 201 on, based on the number of PMs area.
    tmp=geo.BLKLABELS.rotore.xy;
    index=201;
    BrDir=[];
    for ii=1:length(tmp(:,1))
        if tmp(ii,3)==6
            mi_selectlabel(tmp(ii,1),tmp(ii,2));
            mi_setgroup(index);
            BrDir=[BrDir atan2(tmp(ii,7),tmp(ii,6))];
            mi_clearselected;
            index=index+1;
        end
    end
    clear index
    BrGro=200+[1:1:numel(BrDir)];
% end

% change the PMs group: group 200
tmp=geo.BLKLABELS.rotore.xy;
for ii=1:length(tmp(:,1))
    if tmp(ii,3)==6
        mi_selectlabel(tmp(ii,1),tmp(ii,2));
        mi_setgroup(200);
        mi_clearselected;
    end
end
mi_saveas([resFolder filename(1:end-4) '.fem']);
mi_close;

mat.LayerMag.Hc = Br/(4e-7*pi*mat.LayerMag.mu);
%per.tempPP = temp;

Bmin     = Br*ones(ncur,2); % the dimensions are current on the rows and position on the lines
dPM      = ones(ncur,2);
VolDem   = zeros(ncur,2);
VolTot   = zeros(ncur,2);
xyDemagC = cell(ncur,2);
xyDemagV = cell(ncur,2);
xyDemagB = cell(ncur,2);

%[i0,~]=calc_io(geo,per);

for ii=1:ncur
    
    disp(['Current level ' int2str(ii) ' of ' int2str(ncur)]);
    
    if (strcmp(geo.RotType,'SPM')||strcmp(geo.RotType,'Vtype'))
        per.gamma = 180;
    else
        per.gamma = 90;
    end
    per.overload = Ivect(ii)/i0;
    nsim = 2;
    
    for jj = 1:nsim
        
        [~,pathname]=createTempDir();
        
        copyfile([resFolder filename(1:end-4) '.fem'],[pathname filename(1:end-4) '.fem'])
        
        per.nsim_singt = 1;
        per.delta_sim_singt = (0.5*360/(6*geo.q))*(jj-1);    % 1) theta = 0, 2) theta = 1/2 slot pitch
        
        SOL = simulate_xdeg(geo,per,mat,'idemagmap',pathname,[filename(1:end-4) '.fem'],1);
        
        if isfolder([pathname 'critical machines'])
            copyfile([pathname 'critical machines'],motFolder)
        end
        
        Bmin(ii,jj)      = SOL.Bmin;
        dPM(ii,jj)       = SOL.dPM;
        VolDem(ii,jj)    = SOL.VolDem;
        VolTot(ii,jj)    = SOL.VolTot;
        xyDemagC{ii,jj}  = SOL.xyDemagTmpC;
        xyDemagV{ii,jj}  = SOL.xyDemagTmpV;
        xyDemagB{ii,jj}  = SOL.xyDemagTmpB;
        
    end
    
    %update figure
    if figFlag
        for hh=1:2
            set(hplot{1,hh},'XData',Ivect(1:ii),'YData',Bmin(1:ii,hh),'MarkerIndices',ii);
            set(hplot{2,hh},'XData',Ivect(1:ii),'YData',dPM(1:ii,hh),'MarkerIndices',ii);
        end
        drawnow
    end
    
end



DemagOut.Ivect     = Ivect;
DemagOut.Bmin      = Bmin;
DemagOut.dPM       = dPM;
DemagOut.VolDem    = VolDem;
DemagOut.VolTot    = VolTot;
DemagOut.xyDemag.C = xyDemagC;
DemagOut.xyDemag.V = xyDemagV;
DemagOut.xyDemag.B = xyDemagB;
DemagOut.setup     = setup;

save([resFolder 'DemagAnalysisResult.mat'],'DemagOut','dataSet','geo','per','mat');

if figFlag
    saveas(hfig{1},[resFolder 'DemagLimit.fig']);
    print(hfig{1},[resFolder 'DemagLimit.png'],'-dpng');
    
    saveas(hfig{2},[resFolder 'DemagLimitPU.fig']);
    print(hfig{2},[resFolder 'DemagLimitPU.png'],'-dpng');
end

% [status,~,~] = rmdir(tmpFolder,'s');
% if status
%     disp('Temporary folder removed')
% else
%     disp('Temporary folder not removed')
% end


