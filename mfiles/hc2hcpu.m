% Copyright 2022
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

function [hc_pu,dataSet] = hc2hcpu(dataSet)


[~,~,geo,~,~] = data0(dataSet);

if strcmp(dataSet.TypeOfRotor,'SPM')
    geo.hc = dataSet.HCmm;
    hc_pu = geo.hc/geo.g;
else

    hc_mm   = dataSet.HCmm;
    r       = geo.r;
    R       = geo.R;
    p       = geo.p;
    pont0   = geo.pont0;
    g       = geo.g;
    nlay    = geo.nlay;
    lt      = geo.lt;
    hfe_min = geo.hfe_min;
    x0      = geo.x0;
    alpha   = cumsum(geo.dalpha);
    hfe_min = geo.hfe_min;
    
    
    cos_x0 = cos(pi/2/p);
    sin_x0 = sin(pi/2/p);
    ArLim = r*(1/cos_x0-sqrt(((1-cos_x0^2)/cos_x0)^2+sin_x0^2));    % (max) shaft radius [mm]
    
    
    if nlay==1
        hcMax1 = 2*(alpha*pi/180/(1+alpha*pi/180)*(r-pont0));
        temp_alpha_hfemin = hfe_min/r; % rad
        temp_alpha_hc_2 = pi/(2*p) - alpha*pi/180 - temp_alpha_hfemin;
        hcMax2 = 2*(temp_alpha_hc_2/(1+temp_alpha_hc_2)*(r-pont0));
        hcMax = min([hcMax1 hcMax2]);
        hc_pu = hc_mm/hcMax;
    else
        beta = calc_apertura_cerchio(alpha*pi/180,r,x0);
        rbeta = (x0-r*cosd(alpha))./(cos(beta));
        [xpont,ypont] = calc_intersezione_cerchi(r-pont0,rbeta,x0);
        rpont_x0 = sqrt(ypont.^2+(x0-xpont).^2);
        Bx0 = x0-(rpont_x0);
        hc_min = (r-ArLim-(R-r-g-lt))/nlay/4;
        delta = (0.5*hc_mm(1)+sum(hc_mm(2:end-1))+0.5*hc_mm(end)-hc_min*(nlay-1))/(Bx0(1)-Bx0(end)-hfe_min*(nlay-1)-hc_min*(nlay-1));
        hc_pu = hc_mm*(delta*nlay)/sum(hc_mm);
    end
end

dataSet.HCpu = hc_pu;




