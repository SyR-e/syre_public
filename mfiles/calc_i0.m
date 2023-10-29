% Copyright 2021
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

function [per] = calc_i0(geo,per,mat)
% Computation of the rated current and phase resistance from loss, thermal loading or slot
% current density

kj   = per.kj;      % thermal load [W/m^2]
J    = per.J;       % current density [Arms/mm^2]
loss = per.Loss;    % admitted loss [W]
i0   = per.i0;      % current [Apk]
%% rimuovere
% kj   = nan;
% J    = nan;
% loss = nan;
% i0   = per.i0;
%%

lend  = geo.lend/1e3;       % end-winding length [m]
l     = geo.l/1e3;          % stack length [m]
R     = geo.R/1e3;          % outer radius [m]
Ns    = geo.win.Ns;         % number of turns in series per phase
Aslot = geo.Aslot/1e6;      % slot area [m^2]
kcu   = geo.win.kcu;        % slot filling factor
n3ph  = geo.win.n3phase;    % number of three-phase sets
p     = geo.p;              % pole pairs number
q     = geo.q;              % number of slot per pole per phase

tempCu = per.tempcu;        % target copper temperature [°C]

if exist('mat','var')
    ro0 = 1/mat.SlotCond.sigma;
    alphaCond = mat.SlotCond.alpha;
    rocu = ro0*(1+alphaCond*(tempCu-20));
else
    rocu = (1.7241e-08)*(234.5+tempCu)/(234.5+20);
    warning('Copper windings computation');
end
Aslots = Aslot*(6*p*q*n3ph);

flag = 1;
if ~isnan(kj)
    loss = kj.*(2*pi*R*l);
    J = ((kj.*pi*R/3.*l./(l+lend)./n3ph./(p*q*Aslot*kcu)./rocu).^0.5)/1e6;
elseif ~isnan(loss)
    kj = loss./(2*pi*R*l);
    J = ((kj.*pi*R/3.*l./(l+lend)./n3ph./(p*q*Aslot*kcu)./rocu).^0.5)/1e6;
elseif ~isnan(J)
    kj = (J*1e6).^2.*p*q.*Aslot.*kcu.*rocu.*3./(pi*R).*(l+lend)./l.*n3ph;
    loss = kj.*(2*pi*R*l);
elseif ~isnan(i0)
    Rs   = 12*rocu*(l+lend)./(kcu*Aslots)*Ns^2;
    loss = n3ph*3/2*Rs*i0^2;
    kj   = loss./(2*pi*R*l);
    J    = ((kj.*pi*R/3.*l./(l+lend)./n3ph./(p*q*Aslot*kcu)./rocu).^0.5)/1e6;
    flag = 0;
else
    warning('Wrong rated currrent input!')
end

if flag
    i0 = 1/Ns.*(kj.*kcu./rocu.*l./(l+lend).*pi.*R.*Aslots/9).^0.5/n3ph; % [Apk]
    i0 = real(i0);
    Rs = loss./(n3ph.*3/2.*i0.^2);
end

per.Loss = loss;
per.kj   = kj;
per.J    = J;

per.i0 = i0;
per.Rs = Rs;

