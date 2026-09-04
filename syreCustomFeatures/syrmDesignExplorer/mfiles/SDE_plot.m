% Copyright 2021
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

function [hfig] = SDE_plot(map,flag3D)


figure()
figSetting();
xlabel('$x$');
ylabel('$b$');
set(gca,'XLim',[min(map.xx,[],'all') max(map.xx,[],'all')])
set(gca,'YLim',[min(map.bb,[],'all') max(map.bb,[],'all')])
if flag3D
    view(3)
    legend('show','Location','north');
else
    plot(map.xRaw,map.bRaw,'o','Color',[0 0.5 0],'MarkerFaceColor',[0 0.5 0],'DisplayName','FEAfix')
    plot(map.xx(isnan(map.T)),map.bb(isnan(map.T)),'rx','DisplayName','unfeasible','MarkerSize',8)
    legend('show','Location','northeastoutside');
end

colors = get(gca,'ColorOrder');
jj=0;
[~,gui2map] = SDE_namesMapping();
for ii=1:length(map.dataSelect)
    dataSelectGUI = map.dataSelect{ii};
    if isKey(gui2map,map.dataSelect{ii})
        map.dataSelect{ii} = gui2map(map.dataSelect{ii});
    else
        map.dataSelect{ii} = map.dataSelect{ii}(7:end);
    end
    command = ['tmp = map.' map.dataSelect{ii} ';'];
    eval(command);

    %if sum(size(tmp{1,1}))>1
    if iscell(tmp)
        for jj=1:numel(tmp{1,1})
            command = ['tmp = map.' map.dataSelect{ii} ';'];
            eval(command);
            tmp = cellfun(@(v)v(jj),tmp);

            tmp = tmp.*ones(size(map.xx));
            if flag3D
                surf(map.xx,map.bb,tmp,'DisplayName',[map.dataSelect{ii} num2str(jj)]);
            else
                indexColor = ii+jj;
                if indexColor>size(colors,1)
                    indexColor=indexColor-floor(indexColor/size(colors,1))*size(colors,1);
                end
                contour(map.xx,map.bb,tmp,'LineWidth',1,'LineColor',colors(indexColor,:),'DisplayName',[map.dataSelect{ii} num2str(jj)],'ShowText','on');
            end
        end
    else
        tmp = tmp.*ones(size(map.xx));
        if flag3D
            surf(map.xx,map.bb,tmp,'DisplayName',dataSelectGUI);
        else
            indexColor = ii+jj;
            if indexColor>size(colors,1)
                indexColor=indexColor-floor(indexColor/size(colors,1))*size(colors,1);
            end
            contour(map.xx,map.bb,tmp,'LineWidth',1,'LineColor',colors(indexColor,:),'DisplayName',dataSelectGUI,'ShowText','on');
        end
    end

end


map = rmfield(map,'dataAvailable');
map = rmfield(map,'dataSelect');
set(gcf,'UserData',map);

if nargout==1
    hfig = gcf;
end


