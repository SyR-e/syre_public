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

function draw_lines_arcs(Mat,gruppo,res)
% Importa_Mat(Mat,gruppo)
[nrig,ncol] = size(Mat);
% stop=130;
for ii=1:nrig
    if Mat(ii,7)==0
        mi_drawline(Mat(ii,3),Mat(ii,4),Mat(ii,1),Mat(ii,2));
        mi_selectsegment(mean([Mat(ii,1) Mat(ii,3)]),mean([Mat(ii,2) Mat(ii,4)]));
        
        mi_setsegmentprop('None', res, 0, 0, gruppo);
        mi_clearselected
        mi_selectnode(Mat(ii,1),Mat(ii,2));
        mi_selectnode(Mat(ii,3),Mat(ii,4));
        mi_setnodeprop('None',gruppo);
        mi_clearselected
        %         keyboard
    elseif Mat(ii,7)== eps || Mat(ii,7)== -eps  %mod walter
        mi_drawline(Mat(ii,1),Mat(ii,2),Mat(ii,3),Mat(ii,4));
        mi_selectsegment(mean([Mat(ii,1) Mat(ii,3)]),mean([Mat(ii,2) Mat(ii,4)]));
        
        mi_setsegmentprop('None', res, 0, 0, gruppo);
        mi_clearselected
        mi_selectnode(Mat(ii,1),Mat(ii,2));
        mi_selectnode(Mat(ii,3),Mat(ii,4));
        mi_setnodeprop('None',gruppo);
        mi_clearselected
        
    else
        
        [maxsegdeg,raggio,ang1,ang]=Disegna_Arco(Mat(ii,:),res);
        [x_temp,y_temp]=pol2cart((ang1+0.5*ang)*pi/180,raggio);
        x=x_temp+Mat(ii,1);
        y=y_temp+Mat(ii,2);
        % mi_selectarcsegment( Mat(ii,5)*0.01+Mat(ii,3)*0.99 , Mat(ii,6)*0.01+Mat(ii,4)*0.99 );
        % mi_setarcsegmentprop(res, 'None', 0, gruppo);
        mi_selectarcsegment(x,y);
        mi_setarcsegmentprop(maxsegdeg, 'None', 0, gruppo);
        mi_clearselected
        mi_selectnode(Mat(ii,5),Mat(ii,6));
        mi_selectnode(Mat(ii,3),Mat(ii,4));
        mi_setnodeprop('None',gruppo);
        mi_clearselected
        %        keyboard
    end
end
% mi_saveas('testZero.fem')
% system(['move testZero.fem ..'])
% cd('..')
% error('OK')



