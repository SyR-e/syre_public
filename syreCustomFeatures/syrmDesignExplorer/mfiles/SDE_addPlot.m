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

function SDE_addPlot(app)

map = app.map;

hax = app.UIAxes;

colors{1} = [0.0 0.0 1.0];
colors{2} = [1.0 0.0 0.0];
colors{3} = [0.0 0.8 0.0];
colors{4} = [1.0 0.5 0.0];
colors{5} = [0.0 0.8 0.8];
colors{6} = [0.8 0.0 0.8];
colors{7} = [1.0 0.8 0.0];

for ii=1:length(map.dataSelect)
    command = ['tmp = map.' map.dataSelect{ii} ';'];
    eval(command);
    tmp = tmp.*ones(size(map.xx));
    hchild = get(hax,'Children');

    indexColor = length(hchild);
    if indexColor>size(colors)
        indexColor=indexColor-floor(indexColor/size(colors))*size(colors);
    end
    contour(app.UIAxes,map.xx,map.bb,tmp,'LineWidth',1,'LineColor',colors{indexColor},'DisplayName',map.dataSelect{ii},'ShowText','on');
end


