% Copyright 2025
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


function [per] = calc_if(geo,per,mat)

Jf  = per.Jf;
if0 = per.if0;


l      = geo.l/1e3;          % stack length [m]
lendf  = geo.lendf/1e3;      % end-winding length [m]
Nf     = geo.win.Nf;         % number of turns in series per pole
if isfield(geo,'Acoilf')
    Acoilf = geo.Acoilf/1e6;      % slot area [m^2]
else
    Acoilf = NaN;
end

kcuf   = geo.win.kcuf;        % slot filling factor
p      = geo.p;              % pole pairs number

tempCu = per.tempcu;        % target copper temperature [°C]


if exist('mat','var')
    ro0 = 1/mat.BarCond.sigma;
    alphaCond = mat.BarCond.alpha;
    rocu = ro0*(1+alphaCond*(tempCu-20));
else
    rocu = (1.7241e-08)*(234.5+tempCu)/(234.5+20);
    warning('Copper rotor winding computation');
end

flag = 1;
if ~isnan(Jf)
    if0 = Jf*(Acoilf*1e6*kcuf)/Nf;
elseif ~isnan(if0)
    Jf = if0/(Acoilf*1e6*kcuf)*Nf;
    flag = 0;
else
    warning('Wrong field currrent input!')
end


Rf = rocu*2*(l+lendf)/(Acoilf*kcuf)*Nf^2*(2*p);

per.Rf = Rf;
per.if0 = if0;
per.Jf = Jf;











