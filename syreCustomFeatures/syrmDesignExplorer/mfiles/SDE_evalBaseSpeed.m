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

function [map] = SDE_evalBaseSpeed(map)
p       = map.geo.p;
% kj      = map.kj;
% R       = map.geo.R;
% l       = map.geo.l;
% n3ph    = map.geo.win.n3phase;
% Ns      = map.geo.win.Ns;
% i0      = map.i0;
fd      = map.fd;
fq      = map.fq;
fM      = map.fM;
id      = map.id;
iq      = map.iq;
% q       = map.geo.q;
T       = map.T; 
Vdc     = map.Vdc;
Rs      = map.Rs;

fdq = fd + j*fq;
%Rs = kj.*(2*pi*R/1000*l/1000)./(n3ph.*3/2.*i0.^2);


[w_elet] = calcLimitPulsation(id,iq,fd,fq,Rs,Vdc/sqrt(3));
% f_elet = w_elet/2/pi;
% nbase = f_elet/pi*60;
nbase = w_elet*30/pi/p;
Pbase =  T.*w_elet/p;
% fUGOpu = (abs(fdq))./fM;

% map.fUGOpu = fUGOpu;
map.nbase = nbase;
map.Pbase = Pbase;