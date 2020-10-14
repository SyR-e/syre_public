% Copyright 2020
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

function [mFeS,mFeR] = calcMassFe(geo,mat)

%% stator

stator  = geo.stator;
kgm3    = mat.Stator.kgm3;
p       = geo.p;
ps      = geo.ps;
Ro      = geo.R;
Ri      = (geo.r+geo.g);
r       = geo.r;
l       = geo.l;
Qs      = geo.Qs;
Ar      = geo.Ar;

area = pi*(Ro^2-Ri^2)*ps/(2*p);

% Compure the area of just one slot (+slot air) and then subctract to the
% annulus portion area (iron)

areaSlot = geo.Aslot;

nAir = sum(stator(:,8)==2)/Qs;
tmp = stator(:,9);
index = 1:1:size(stator,1);
tmp = tmp(stator(:,8)==2);
index = index(stator(:,8)==2);
val = tmp(1);

indTmp = find(tmp==val,nAir,'first');
index = index(indTmp);



X = [];
Y = [];

for ii=1:length(index)
    X = [X stator(index(ii),1) stator(index(ii),3)];
    Y = [Y stator(index(ii),2) stator(index(ii),4)];
end

X = [X X(1)];
Y = [Y Y(1)];

areaAir = polyarea(X,Y);

area = (area-Qs*(areaSlot+areaAir))*2*p/ps;

mFeS = area/1e6*l/1e3*kgm3;


%% rotor

rotor = geo.rotor;
kgm3  = mat.Rotor.kgm3;

nEle = max(rotor(:,9));

mFeR = pi*((r/1e3)^2-(Ar/1e3)^2)*ps/(2*p)*l/1e3*kgm3;

for ii=1:nEle
    rotTmp = rotor(rotor(:,9)==ii,:);
    if rotTmp(1,8)==1 % air
        areaTmp = calcAreaShape(rotTmp);
        mFeR = mFeR-(areaTmp/1e6*l/1e3*kgm3);
    end
end

mFeR = mFeR*2*p/ps;




