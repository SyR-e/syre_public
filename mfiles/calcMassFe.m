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

areaS = pi*(Ro^2-Ri^2)*ps/(2*p);

% Compure the area of just one slot (+slot air) and then subctract to the
% annulus portion area (iron)

areaSlot = geo.Aslot;
if geo.acs~=0
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
else
    areaAir = 0;
end

areaS = (areaS-Qs*(areaSlot+areaAir))*2*p/ps;

mFeS = areaS/1e6*l/1e3*kgm3;

if isfield(geo,'pShape')
    if geo.pShape.flag
        areaS = area(geo.pShape.stator);
        mFeS = areaS/1e6*l/1e3*kgm3*2*p/ps;
    end
end




%% rotor

rotor = geo.rotor;
kgm3  = mat.Rotor.kgm3;

if strcmp(geo.RotType,'EESM')
    ry = Ar + geo.lyr;
    xR(1) = Ar*cos(pi/p);
    yR(1) = Ar*sin(pi/p);
    alpha = atan(rotor(1,6)/rotor(1,5));
    theta = linspace(pi/p,alpha,10);
    for ii = 2:11
        xR(ii) = ry*cos(theta(ii-1));
        yR(ii) = ry*sin(theta(ii-1));
    end
    pos = 12;
    for ii = 2:length(rotor)
        if rotor(ii,9) == 1
            if rotor(ii,7) == 0
                xR(pos) = rotor(ii,1);
                yR(pos) = rotor(ii,2);
            elseif rotor(ii,7) == -1
                xR(pos) = rotor(ii,3);
                yR(pos) = rotor(ii,4);
            else % ( == 1)
                %warning('Error in rotor iron mass calculation');
            end
            pos = pos + 1;
        end
    end
    flag = 0;
    for ii = length(rotor):-1:1
        if rotor(ii,9) == 3
            if rotor(ii,7) == 0
                xR(pos) = rotor(ii,1);
                yR(pos) = rotor(ii,2);
            elseif rotor(ii,7) == 1
                xR(pos) = rotor(ii,3);
                yR(pos) = rotor(ii,4);
            else % ( == 1)
                %warning('Error in rotor iron mass calculation');
            end
            pos = pos + 1;
        elseif rotor(ii,9) == 2 && flag == 0
            flag = 1;
            pos = pos-1;
            alpha = atan(rotor(1,6)/-rotor(1,5))+pi/p;
            theta = linspace(alpha,0,10);
            for ii = 1:10
                xR(pos+ii-1) = ry*cos(theta(ii));
                yR(pos+ii-1) = ry*sin(theta(ii));
            end
            pos = pos + ii - 1;
        end
    end
    theta = linspace(0,pi/p,20);
    for ii = 1:20
        xR(pos+ii) = Ar*cos(theta(ii));
        yR(pos+ii) = Ar*sin(theta(ii));
    end
    
    % % Check Geometry
    % figure
    % figSetting
    % plot(xR,yR,'b-')

    AR = polyarea(xR,yR);

    mFeR = round(AR/1e6*2*p*l/1e3*kgm3,2);

    %warning('Rotor iron mass for EESM not yet computed')
else
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
    
    if isfield(geo,'pShape')
        if geo.pShape.flag
            areaR = area(geo.pShape.rotor);
            mFeR = areaR/1e6*l/1e3*kgm3*2*p/ps;
        end
    end
end







