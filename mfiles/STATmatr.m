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

function [geo,stator,BLKLABELSstat] = STATmatr(geo,fem)
%%
%% Stator
%%
% described by matrixes contatining
% the master nodes
% the description of straight lines
% the description of arches
% the nodes where block labels must be placed

Qs           = geo.Qs;
ps           = geo.ps;
avv          = geo.win.avv;
lt           = geo.lt;
p            = geo.p;
q            = geo.q;
R            = geo.R;
r            = geo.r;
g            = geo.g;
n3ph         = geo.win.n3phase;
pont0        = geo.pont0;
flagDivision = geo.statorYokeDivision;  % division line for end-winding extrusion (GalFer Challenge 2024)

% mesh resolution
res     = fem.res;
% res_traf = fem.res_traf;

alpha_slot = 2*pi/(6*p*q*n3ph);   % angolo di passo cava
RSI=r+g;                     % r traferro statore
RSE=R;                        % r esterno statore
RS5=r+g+lt;                  % r esterno cave
PassoCava=(RSI)*2*pi/(Qs*ps);  % passo cava di statore


% draw one slot
[geo,temp] = drawSlot(geo);

CavaMat = temp.CavaMat;
xAir = temp.xAir;
yAir = temp.yAir;
xCu  = temp.xCu;
yCu  = temp.yCu;
xFe  = temp.xFe;
yFe  = temp.yFe;


materialCodes;
% codMatAir = 2;
% codMatFeSta  = 4;
% codMatCu  = 3;

nCu  = size(xCu,2);
nAir = size(xAir,2);
nFe  = size(xFe,2);

if ~isnan(xAir)
    xy0 = [
        xCu'  yCu'  codMatCu*ones(nCu,1)      res*ones(nCu,1)  ones(nCu,1)
        xAir' yAir' codMatAirSta*ones(nAir,1) res*ones(nAir,1) ones(nAir,1)
        xFe'  yFe'  codMatFeSta*ones(nFe,1)   res*ones(nFe,1)  ones(nFe,1)
        ];
else
    xy0 = [
        xCu'  yCu'  codMatCu*ones(nCu,1)    res*ones(nCu,1)  ones(nCu,1)
        xFe'  yFe'  codMatFeSta*ones(nFe,1) res*ones(nFe,1)  ones(nFe,1)
        ];
end

[nRow_mat,nCol_mat] = size(CavaMat);
[nRow_xy,nCol_xy]   = size(xy0);

stator  = zeros(Qs*nRow_mat,nCol_mat);
xyLabel = zeros(Qs*nRow_xy,nCol_xy);

% replicate/rotate one slot to obtain all slots
for ii=1:Qs
    % geometry
     CavaTmp = rotateMatrix(CavaMat,(ii-1)*alpha_slot);
     indexEle = CavaTmp(:,9);
     indexEle(indexEle==0)=-(ii-1);
     CavaTmp(:,9) = indexEle+(ii-1);
     stator(1+nRow_mat*(ii-1):nRow_mat*ii,:) = CavaTmp;
    % labels
    [xtemp,ytemp] = rot_point(xy0(:,1),xy0(:,2),(ii-1)*alpha_slot);
    xyLabel(1+nRow_xy*(ii-1):nRow_xy*ii,:) = [xtemp ytemp xy0(:,3:end)];
end



% draw the tooth divisions

thFe = 0:alpha_slot:2*pi/(6*p*q*n3ph)*Qs;
nFe = length(thFe);

stator = [stator
    (r+g)*cos(thFe') (r+g)*sin(thFe') R*cos(thFe') R*sin(thFe') nan(nFe,1) nan(nFe,1) zeros(nFe,1) codMatFeSta*ones(nFe,1) zeros(nFe,1)
    ];


% complete the stator with the outer lines and the slot divisions
indexEle = max(stator(:,9));

if (ps<2*p)
    % partial machine
    xse1 = R;
    yse1 = 0;
    xse2 = R*cos(pi/p*ps);
    yse2 = R*sin(pi/p*ps);
    
    xsi1 = r+g;
    ysi1 = 0;
    xsi2 = (r+g)*cos(pi/p*ps);
    ysi2 = (r+g)*sin(pi/p*ps);

    xsm1 = RS5+pont0;
    ysm1 = 0;
    xsm2 = (RS5+pont0)*cos(pi/p*ps);
    ysm2 = (RS5+pont0)*sin(pi/p*ps);
    
    stator = [stator
        xsi1 ysi1 xse1 yse1 NaN  NaN   0 codMatFeSta    indexEle+1
        0    0    xse1 yse1 xse2 yse2 +1 codMatFeSta    indexEle+1
        xse2 yse2 xsi2 ysi2 NaN  NaN   0 codMatFeSta    indexEle+1
        0    0    xsi2 ysi2 xsi1 ysi1 -1 codMatFeSta    indexEle+1
        ];
    if flagDivision
        stator = [stator
            0    0    xsm1 ysm1 xsm2 ysm2 +1 codMatFeSta    0
            ];
    end
else
    % full machine
    xse1 = R;
    yse1 = 0;
    xse2 = -R;
    yse2 = 0;
    
    xsi1 = r+g;
    ysi1 = 0;
    xsi2 = -(r+g);
    ysi2 = 0;

    xsm1 = RS5+pont0;
    ysm1 = 0;
    xsm2 = -(RS5+pont0);
    ysm2 = 0;
    
    stator = [stator
        0 0 xse1 yse1 xse2 yse2 1 codMatFeSta indexEle+1
        0 0 xse2 yse2 xse1 yse1 1 codMatFeSta indexEle+1
        0 0 xsi1 ysi1 xsi2 ysi2 1 codMatFeSta indexEle+1
        0 0 xsi2 ysi2 xsi1 ysi1 1 codMatFeSta indexEle+1
        ];
    if flagDivision
        stator = [stator
            0 0 xsm1 ysm1 xsm2 ysm2 1 codMatFeSta 0
            0 0 xsm2 ysm2 xsm1 ysm1 1 codMatFeSta 0
            ];
    end
end

% reorder the stator matrix (before the boundary conditions)

stator = [stator(end-3:end,:); stator(1:end-4,:)];

% block labels names (not used for FEMM)
nome_Cu_slot  = cell(size(avv,1)*Qs,1);
nome_FeYoke   = cell(Qs,1);
nome_air_slot = cell(Qs,1);
index = 1;
for num_cave=1:Qs
    for num_str=1:size(avv,1)
        nome_Cu_slot{index} = {['slot_',num2str(num_str),'_',num2str(num_cave)]};
        index=index+1;
    end
    nome_FeYoke{num_cave}   = {['statore_',num2str(num_cave)]};
    nome_air_slot{num_cave} = {['slot_air_',num2str(num_cave)]};
end

% boundary conditions (FEMM-style codes)
codBound_FluxTan  = 0;      % 0 flux tangential
if (Qs<geo.p*geo.q*6*geo.win.n3phase)
    codBound_periodic = 10;     % 10 odd or even periodicity
else
    codBound_periodic = -10;    % -10 no periodicity, simulate full machine
end

xyTmp = xyLabel(xyLabel(:,3)==codMatFeSta,:);
thBoundRSE = atan2(xyTmp(:,2),xyTmp(:,1));
xBoundRSE = RSE*cos(thBoundRSE);
yBoundRSE = RSE*sin(thBoundRSE);


%[xBoundRSE,yBoundRSE]=rot_point(RSE,0,2*pi/(6*p*q)*Qs/2);  % rotazione di 1/4 di passo cava

BoundRSE = [xBoundRSE,yBoundRSE,codBound_FluxTan*ones(size(thBoundRSE))];

[xBoundLAT1,yBoundLAT1] = rot_point(mean([RSE,(RS5+pont0)]),0,0);
[xBoundLAT2,yBoundLAT2] = rot_point(mean([RSE,(RS5+pont0)]),0,Qs*alpha_slot);

[xBoundLAT3,yBoundLAT3] = rot_point(mean([RSI,(RS5+pont0)]),0,0);
[xBoundLAT4,yBoundLAT4] = rot_point(mean([RSI,(RS5+pont0)]),0,Qs*alpha_slot);


BoundLAT=[
    xBoundLAT1 yBoundLAT1 codBound_periodic
    xBoundLAT2 yBoundLAT2 codBound_periodic
    xBoundLAT3 yBoundLAT3 codBound_periodic
    xBoundLAT4 yBoundLAT4 codBound_periodic
    ];

%% OUTPUT
geo.stator = stator;
geo.ly     = geo.R-geo.r-geo.g-geo.lt;

BLKLABELSstat.xy             = xyLabel;
BLKLABELSstat.names.air_slot = nome_air_slot;
BLKLABELSstat.names.Cu_slot  = nome_Cu_slot;
BLKLABELSstat.names.FeYoke   = nome_FeYoke;
BLKLABELSstat.names.legend   = {'air_slot','Cu_slot','FeYoke'};
BLKLABELSstat.boundary       = [BoundLAT;BoundRSE];

