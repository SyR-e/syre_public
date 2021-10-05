% Copyright 2021
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

function [out]=eval_maxStress(structModel,sVonMises,geo,mat)
x         = structModel.Mesh.Nodes(1,:);
y         = structModel.Mesh.Nodes(2,:);
sigma_max = mat.Rotor.sigma_max*10^6;
MaxStress = max(sVonMises);
pos_max = find(sVonMises==MaxStress);

x_max = x(pos_max)*10^3;
y_max = y(pos_max)*10^3;

nodesOver = sum (sVonMises>sigma_max);
pos_over = find(sVonMises>sigma_max);
x_over = x(pos_over)*10^3;
y_over = y(pos_over)*10^3;

stress_R = nan(1,geo.nlay);
stress_T = nan(1,geo.nlay);

if any(geo.pontR & strcmp(geo.RotType,'Seg'))
    ii = (geo.pontR>0 & geo.radial_ribs_split==0);
    xR_xo(ii) = (geo.B2k(ii)-geo.hc(ii)/2);
    yR_xo(ii) = 0;
    
    jj = (geo.pontR>0 & geo.radial_ribs_split==1);
    xR_xo(jj) = (geo.B2k(jj)-geo.hc(jj)/2);
    yR_xo(jj) = geo.YpBar1(jj)-geo.pontR(jj)/4;
    
    [xR,yR] = rot_point(xR_xo,yR_xo,pi/2/geo.p);
    
    for kk=1:geo.nlay
        if xR(kk)>0
            [~,pos_mid]  = min(abs((x-xR(kk)*10^-3))+abs((y-yR(kk)*10^-3)));
            stress_R(kk) = sVonMises(pos_mid);
        end
    end
end

if strcmp(geo.RotType,'Seg')
    xT_xo = geo.xpont;
    yT_xo = geo.ypont;
    
    [xT,yT] = rot_point(xT_xo,yT_xo,pi/2/geo.p);
    
    for kk=1:geo.nlay
        if xT(kk)>0
            [~,pos_mid]  = min(abs((x-xT(kk)*10^-3))+abs((y-yT(kk)*10^-3)));
            stress_T(kk) = sVonMises(pos_mid);
        end
    end
end

out.MaxStress = MaxStress;
out.x_max     = x_max;
out.y_max     = y_max;
out.x_over    = x_over;
out.y_over    = y_over;
out.nodesOver = nodesOver;
out.stress_T  = stress_T;
out.stress_R  = stress_R;

