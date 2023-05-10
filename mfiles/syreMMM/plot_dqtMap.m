
% Copyright 2018
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_dqtMap(hax,data,filename)

[xS,yS,zS]=size(data.x);
th=data.th;

hfig=get(hax,'Parent');

v = VideoWriter([filename(1:end-4) '.avi']);
v.Quality = 95;
open(v);

for ii=1:zS
    cla(hax)
    surf(hax,data.x(:,:,ii),data.y(:,:,ii),data.z(:,:,ii));
    title(hax,['$\theta_e = ' num2str(th(ii)) ' ^\circ$']);
    drawnow
    
    frame=getframe(hfig);
    im=frame2im(frame);
    [imind,cm]=rgb2ind(im,256);
    if ii==1
        imwrite(imind,cm,filename,'gif','LoopCount',inf,'DelayTime',0.01);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.01);
    end
    writeVideo(v,im);
end

close(v);