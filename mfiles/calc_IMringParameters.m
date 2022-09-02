% Copyright 2022
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

function [Rring,Lring] = calc_IMringParameters(geo,per,mat)

% The parameters refer to one portion of the ring between 2 bars
% Ref for Lring: Boldea, The Induction Machine Handbook

wt    = geo.IM.wt;
lt    = geo.IM.lt;
Nbars = geo.IM.Nbars;
r     = geo.r;
ws    = 2*pi*r/Nbars-wt;
l     = geo.l;
p     = geo.p;

mu0 = 4*pi*1e-7;
% rho = 1/mat.SlotCond.sigma;
% warning('bar conductor not considered!!!')
% rho = 1/mat.BarCond.sigma;

tempRot = per.tempPP;

if exist('mat','var')
    ro0 = 1/mat.BarCond.sigma;
    alphaCond = mat.BarCond.alpha;
    rho = ro0*(1+alphaCond*(tempRot-20));
else
    rho = (1.7241e-08)*(234.5+tempRot)/(234.5+20);
    warning('Copper windings computation');
end


Rring = rho*2*pi*(r-lt/2)/1000/(Nbars*lt/1e3*ws/1e3);

Lring = mu0*(2*pi*(r-lt/2)/1e3)/Nbars*(2.3*2*(r-lt/2)/1e3)/(4*Nbars*l/1e3*sin(pi*p/Nbars)^2)*log(4.7*2*(r-lt/2)/1e3/(ws/1e3+2*lt/1e3));



