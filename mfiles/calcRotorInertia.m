% Copyright 2020
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

function [J] = calcRotorInertia(geo,mat)

materialCodes;

rotor = geo.rotor;
l     = geo.l;
p     = geo.p;
ps    = geo.ps;

rhoFe = mat.Rotor.kgm3;
rhoPM = mat.LayerMag.kgm3;
rhoAir = -rhoFe;
if strcmp(mat.Shaft.MatName,'Air')
    rhoShaft = 0;
else
    rhoShaft = mat.Shaft.kgm3;
end

nEle = max(rotor(:,9));

r = nan(1,nEle);
M = nan(1,nEle);

for ii=1:nEle
    rotTmp = rotor(rotor(:,9)==ii,:);
    [M(ii),r(ii)] = calcAreaShape(rotTmp);
    switch rotTmp(1,8)
        case codMatAirRot
            M(ii) = (M(ii)*l)/1e9*rhoAir;
        case codMatBar
            M(ii) = (M(ii)*l)/1e9*rhoPM;
        case codMatFeRot
            M(ii) = (M(ii)*l)/1e9*rhoFe;
        case codMatShaft
            M(ii) = (M(ii)*l)/1e9*rhoShaft;
    end
end

r = r/1e3;

J = sum(M.*r.^2)*2*p/ps;


