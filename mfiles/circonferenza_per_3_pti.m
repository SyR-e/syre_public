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

% Determinazione del centro e raggio di una circonferenza per 3 pti noti

function [xc,yc,r,tA,tB]=circonferenza_per_3_pti(xA,yA,xB,yB,xC,yC)

delta=det([xA,yA,1; xB,yB,1;xC,yC,1]);
deltaA=det([-(xA^2+yA^2),yA,1;-(xB^2+yB^2),yB,1;-(xC^2+yC^2),yC,1]);
deltaB=det([xA,-(xA^2+yA^2),1;xB,-(xB^2+yB^2),1;xC,-(xC^2+yC^2),1]);
deltaC=det([xA,yA,-(xA^2+yA^2);xB,yB,-(xB^2+yB^2);xC,yC,-(xC^2+yC^2)]);
a=deltaA/delta; b=deltaB/delta; c=deltaC/delta;
xc=-a/2;
yc=-b/2;
r=sqrt(xc^2+yc^2-c);
tA=atan2((yA-yc),(xA-xc));
tB=atan2((yB-yc),(xB-xc));
% angle=(tA-tB)*180/pi;

end