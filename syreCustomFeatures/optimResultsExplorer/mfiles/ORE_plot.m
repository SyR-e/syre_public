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

function ORE_plot(optRes,OREsetup,hax)

% optRes = app.optRes;
% OREsetup = app.OREsetup;
% hax = app.UIAxes;

if ischar(hax)
    flagUI = 0;
    if strcmp(hax,'2D')
        flag3D=0;
    else
        flag3D = 1;
    end
else
    flagUI = 1;
    flag3D = 0;
end


if isempty(optRes)
    cla(hax)
else
    optRes = ORE_removeUnfeasible(optRes,OREsetup.removeUnfeasible);

    vartot = horzcat('none','motorID',optRes.varNames,optRes.objNames);
    tmp = [ones(size(optRes.motNum))' optRes.motNum' optRes.varParetoFront optRes.objParetoFront];
    limits = [nan nan nan(size(optRes.varNames)) optRes.objLimits'];
    
    tf = strcmp(vartot,OREsetup.xaxis);
    idx = find(tf==1);
    xData = tmp(:,idx)'.*optRes.filt;
    xSet = limits(idx);

    tf = strcmp(vartot,OREsetup.yaxis);
    idx = find(tf==1);
    yData = tmp(:,idx)'.*optRes.filt;
    ySet = limits(idx);

    tf = strcmp(vartot,OREsetup.zaxis);
    idx = find(tf==1);
    zData = tmp(:,idx)'.*optRes.filt;
    zSet = limits(idx);

    tf = strcmp(vartot,OREsetup.caxis);
    idx = find(tf==1);
    cData = tmp(:,idx)'.*optRes.filt;
    cSet = limits(idx);

end

if flagUI
    cla(hax)
    set(hax,'ClimMode','auto','ZLimMode','auto')
    scatter3(hax,xData,yData,zData,[],cData,'filled','MarkerEdgeColor','k')
    xLim = get(hax,'XLim');
    yLim = get(hax,'YLim');
    zLim = get(hax,'ZLim');
    cLim = get(hax,'CLim');
    if ~isnan(xSet)
        plot3(hax,xSet*[1 1],yLim,zLim(1)*[1 1],'-r')
    end
    if ~isnan(ySet)
        plot3(hax,xLim,ySet*[1 1],zLim(1)*[1 1],'-r')
    end
    if ~isnan(zSet)
        set(hax,'ZLim',[zLim(1) zSet])
    end
    if ~isnan(cSet)
        set(hax,'CLim',[cLim(1) cSet])
    end
else
    figure()
    figSetting();
    scatter3(xData,yData,zData,[],cData,'filled','MarkerEdgeColor','k')
    colormap turbo
    xlabel(OREsetup.xaxis)
    ylabel(OREsetup.yaxis)
    zlabel(OREsetup.zaxis)
    c = colorbar(gca);
    c.Label.String = OREsetup.caxis;
    xLim = get(gca,'XLim');
    yLim = get(gca,'YLim');
    zLim = get(gca,'ZLim');
    cLim = get(gca,'CLim');
    if flag3D
        view(3)
        if ~isnan(xSet)
            surf(gca,xSet*[1 1;1 1],[yLim;yLim],[zLim' zLim'],'EdgeColor','r','FaceColor','r','FaceAlpha',0.5)
        end
        if ~isnan(ySet)
            surf(gca,[xLim;xLim],ySet*[1 1;1 1],[zLim' zLim'],'EdgeColor','r','FaceColor','r','FaceAlpha',0.5)
        end
        if ~isnan(zSet)
            surf(gca,[xLim;xLim],[yLim' yLim'],zSet*[1 1;1 1],'EdgeColor','r','FaceColor','r','FaceAlpha',0.5)
        end
        if ~isnan(cSet)
            set(gca,'CLim',[cLim(1) cSet])
        end
    else
        view(2)
        if ~isnan(xSet)
            plot3(gca,xSet*[1 1],yLim,zLim(1)*[1 1],'-r')
        end
        if ~isnan(ySet)
            plot3(gca,xLim,ySet*[1 1],zLim(1)*[1 1],'-r')
        end
        if ~isnan(zSet)
            set(gca,'ZLim',[zLim(1) zSet])
        end
        if ~isnan(cSet)
            set(gca,'CLim',[cLim(1) cSet])
        end
    end
end