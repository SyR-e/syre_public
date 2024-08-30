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

function [geo,mat,temp] = nodes_rotor_SPM_Halbach(geo,mat)

% parameters

r = geo.r;                      % rotor radius
p = geo.p;                      % pole pairs
% phi = geo.phi/p;                % angle range of permanent magnet
phi = geo.dalpha_pu*180/p;
% lm = geo.lm;                    % the thickness of permant magnet
lm  = geo.hc_pu*geo.g;
% hc = 0;
seg = geo.dx;                   % the number of segments of magnet
hs  = geo.hs;                   % rotor sleeve thickness 

r = r-hs;   % sleeve thickness inside the rotor space and not in the airgap

if seg ~= 2
    seg = 2;
end

if seg == 1
    PMdir = 'p';
else
    PMdir = 'r';
end

PMregular = geo.betaPMshape;
if PMregular > 1
    PMregular =1;                      % limit to per unit
end

xPMco = r;
yPMco = 0;

xPMregular = r-lm + PMregular*lm;
yPMregular = 0;

[xPMo,yPMo] = rot_point(xPMregular,yPMregular,phi/2*pi/180);    % PM edge point
[xAiro, yAiro] = rot_point(xPMco,yPMco,90/p*pi/180);                  % Air zone point
[x6, y6] = rot_point(xPMco,yPMco,phi/2*pi/180);                 % Air zone point

xPMci = r-lm; yPMci = 0;
[xPMi,yPMi] = rot_point(xPMci,yPMci,phi/2*pi/180);              % PM edge point on the steel

[xAiri, yAiri] = rot_point(xPMci,yPMci,90/p*pi/180);                  % Rotor edge point
%% chao 2017.01.09 use an arc to assume sinusoidal
xArccenter = (xPMco + xPMo - (yPMo^2/(xPMco-xPMo)))/2;          % find arc center location
yArccenter = 0;

% mean line in the PM (for block labels)
xLabCenter = (mean([xPMco xPMci]) + mean([xPMo xPMi]) - (mean([yPMo yPMi])^2/(mean([xPMco xPMci])-mean([xPMo xPMi]))))/2;          % find arc center location
yLabCenter = 0;


%% calculate PM area
area2 = phi/2/360*pi*r^2 - 0.5*xArccenter*xPMo;
area3 = atan(yPMo/(xPMo-xArccenter))/2*(xPMco-xArccenter)^2 - area2;
PMarea = 2*(area3 + phi/2/360*pi*xPMregular^2 - phi/2/360*pi*xPMci^2);
geo.PMarea = PMarea * 2*p*1e-6;
geo.PMvol = geo.PMarea*geo.l*1e-3;

geo.PMmass = geo.PMarea*geo.l*mat.LayerMag.kgm3*1e-3;
% end

%% segmentation point build by outer circular wave
% if seg~=1
% %     % a quarter pole is considered
% %     NoSeg = 1:floor(seg/2);
% %     SegAngle = phi/seg;                         % angle span of each segment
% %     % coodinate displacement to Arccenter to calculate y on outside shape of PM
% %     [xPMso,yPMso] = rot_point(xPMo-xArccenter,yPMo,-atan(yPMo/(xPMo-xArccenter))*2/seg*NoSeg);
% %     % coodinate recovery
% %     xPMso = xPMso + xArccenter;
% %     % get positions on inside shape of PM
% %     [xPMsi,yPMsi] = rot_point(xPMi,yPMi,-(NoSeg*SegAngle)*pi/180);
% %     
% %     
%     SegAngle = phi/seg;
%     SegAngle = cumsum(SegAngle*ones(1,floor(seg)));
%     %SegAngle = cumsum(SegAngle*ones(1,floor(seg/2)));
%     SegAngle = [0 SegAngle]*pi/180;
%     [xPMso,yPMso] = intersezione_retta_circonferenza(xArccenter,yArccenter,xPMco-xArccenter,tan(SegAngle),0);
%     [xPMsi,yPMsi] = intersezione_retta_circonferenza(0,0,xPMci,tan(SegAngle),0);
%     
% else
%     xPMsi = [];
%     yPMsi = [];
%     xPMso = [];
%     yPMso = [];
% end

SegAngle = phi/seg/2;
SegAngle = cumsum(SegAngle*ones(1,floor(seg)));
%SegAngle = cumsum(SegAngle*ones(1,floor(seg/2)));
SegAngle = [0 SegAngle]*pi/180;                         
[xPMso,yPMso] = intersezione_retta_circonferenza(xArccenter,yArccenter,xPMco-xArccenter,tan(SegAngle),0);
[xPMsi,yPMsi] = intersezione_retta_circonferenza(0,0,xPMci,tan(SegAngle),0);

LabAngle = SegAngle(1:end-1)+diff(SegAngle)/2;        
[xc,yc] = intersezione_retta_circonferenza(xLabCenter,yLabCenter,mean([xPMci xPMco])-xLabCenter,tan(LabAngle),0);

xc = (xPMsi(1:end-1)+diff(xPMsi)/2+xPMso(1:end-1)+diff(xPMso)/2)/2;
yc = (yPMsi(1:end-1)+diff(yPMsi)/2+yPMso(1:end-1)+diff(yPMso)/2)/2;

xmag = cos(atan2(yc,xc)+pi);
ymag = sin(atan2(yc,xc)+pi);

xmag_rot = xmag(1,2);
ymag_rot = ymag(1,2);

xmag(1,2) = ymag_rot;              % 90° rotation of the magn direction
ymag(1,2) = -xmag_rot;

zmag = zeros(size(xc));

%%  SAVE THE FINAL DATA:
temp.xPMco      = xPMco;
temp.yPMco      = yPMco;
temp.xPMo       = xPMo;
temp.yPMo       = yPMo;
temp.xPMci      = xPMci;
temp.yPMci      = yPMci;
temp.xPMi       = xPMi;
temp.yPMi       = yPMi;
temp.xAiri      = xAiri;
temp.yAiri      = yAiri;
temp.xAiro      = xAiro;
temp.yAiro      = yAiro;
temp.xArccenter = xArccenter;
temp.yArccenter = yArccenter;
temp.x6         = x6;
temp.y6         = y6;
temp.xPMso      = xPMso(2:end-1);
temp.yPMso      = yPMso(2:end-1);
temp.xPMsi      = xPMsi(2:end-1);
temp.yPMsi      = yPMsi(2:end-1);

temp.xc   = xc;        
temp.yc   = yc;
temp.xmag = xmag;
temp.ymag = ymag;
temp.zmag = zmag;

temp.mirrorFlag    = ones(1,length(xc));
temp.mirrorFlagAir = ones(1,length(xc));

geo.hc = lm;