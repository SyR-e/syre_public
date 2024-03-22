% Copyright 2023
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

function [geo,mat,temp] = nodes_rotor_Spoke(geo,mat)

% Thesis Automotive Engineering
% Matrix created by Franco Porras


%% Inizializzazione dei dati in ingresso
% angle=geo.VanglePM;    % Inclinazione barriera [radianti]

r               = geo.r;                % Inner airgap radius
rshaft          = geo.Ar;               % Shaft Radio
hc              = geo.hc;               % Thickness of permanent magnet
g               = geo.g;                % Airgap thickness
wm              = geo.PMdim(1,1);       % PM radial thickness
% geo.PMdim(2,1)  = 0;
pontT           = geo.pontT;            % rib thickness
p               = geo.p;                % Poles pairs
hfe_min         = geo.hfe_min;          % Minimum distance between shaft and magnet

% nlay = geo.nlay;    % N° layers

%% Miss typos inputs
if pontT <= 0
    pontT = 1;
end

if hfe_min <= 0
    hfe_min = 0.2;
end

if hc <= 0
    hc = 1;
end
%% Radial axis points of permanent magnet

xxB1k = r - pontT; 
yyB1k = 0;   
xxB2k = xxB1k - wm;
yyB2k = 0;

%% Corner points of barrier

xxD1k = xxB1k;              % Upper left point of permanent magnet - x
yyD1k = yyB1k + hc/2;       % Upper left point of permanent magnet - y
xxD2k = xxB2k;              % Bottom left point of permanent magnet - x
yyD2k = yyB2k + hc/2;       % Bottom left point of permanent magnet - y

% Check rib direction & tangential
hipo = sqrt(xxB1k^2+(hc/2)^2);          % radio from origin to a outer corner of PM
y_limit = abs(sqrt((r-hfe_min)^2-xxB1k^2));  % limit width of PM
tan_out = r - hipo;                     % difference between corner and r-line



if tan_out < hfe_min
    yyD1k = y_limit;
    yyD2k = yyD1k;
    hc = abs(y_limit*2);
    hipo = sqrt(xxB1k^2+(hc/2)^2);              % just to check
    tan_out = r - hipo;                         % just to check
end


% Check inner radio and ribs

rad_in  = xxB2k - rshaft;

if rad_in < hfe_min 
    xxB2k = rshaft + hfe_min;
    xxD2k = xxB2k;
    wm = xxB1k - xxB2k;
end
    % compute the side line of the draw machine
x1 = 0; y1=0;                                               % origin
[x2,y2] = rot_point(r,0,pi/(2*p));                          % intersection vertice sideline/circunference
[m,q] = calc_line2points(x1,y1,x2,y2);
[a,b,c] = calc_retta_offset(m,-1,q,-hfe_min);               % a*x+b*y+c=0, offset line limit
tan_rib = calc_distanza_retta_punto(a,b,c,xxD2k,yyD2k);     % distance from one corner to limit offset
flag_ang = atan(yyD2k/xxD2k);                               % slope angle

if tan_rib < hfe_min || flag_ang > pi/(2*p)                 % if distance is in range, or slope is bigger
    yyD2k   = (-c-a*xxD2k)/b;                               % limit y point of the offset line and magnet corner
    yyD1k   = yyD2k;
    hc      = 2*yyD2k;
end

%% Magnet Geometry
% Center of the permanent magnet

% Comment SF: start do draw the geometry from the right side (towards
% airgap)

xc      = xxB1k - wm/2;           % Center point in magnet axis - x
yc      = yyB1k;                  % Center point in magnet axis - Y

%% Magnet Direction

xmag = cos(atan2(yc,xc)-pi/2);
ymag = sin(atan2(yc,xc)-pi/2);
zmag = 0;

%% MirrorFlag

mirrorFlag = zeros(size(xc));
mirrorFlagAir = zeros(size(xc));


% 1) PM radial dimension --> if the PM stay inside the rotor and not in the shaft
% 2) PM tangential dimension --> if the PM stay inside the pole

% If it helps
% Check retta_per_2pti.m function
% Check calc_retta_offset.m function 


%% Assignment variable - output

geo.hc          = abs(hc);
geo.PMdim(1,1)  = wm;
geo.pontT       = pontT;
geo.hfe_min     = hfe_min;


temp.xxB1k      = xxB1k;
temp.yyB1k      = yyB1k;
temp.xxB2k      = xxB2k;
temp.yyB2k      = yyB2k;
temp.xxD1k      = xxD1k;
temp.yyD1k      = yyD1k;
temp.xxD2k      = xxD2k;
temp.yyD2k      = yyD2k;
temp.xc         = xc;
temp.yc         = yc;
temp.xmag       = xmag;
temp.ymag       = ymag;
temp.zmag       = zmag;
temp.mirrorFlag = mirrorFlag;
temp.mirrorFlagAir = mirrorFlagAir;





