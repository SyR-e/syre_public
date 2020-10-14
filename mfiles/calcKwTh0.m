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

function [kw, th0] = calcKwTh0(geo)

% INPUT
% - WinTable     % winding description (one period)

% OUTPUT
% - kw: winding coeff
% - phase1_offset: angular phase of phase 1 respect to the horizon (elt deg)
% refers to phase 1 also for multi-three-phase

WinTable = geo.win.avv;
Q = 6*geo.q*geo.p*geo.win.n3phase; 
p = geo.p;
n3ph = geo.win.n3phase;

% Star of slots
[nLayers,nSlots] = size(WinTable);
Passo_elt = 2*pi*p/Q;
Posiz_cave = zeros(1,nSlots);
for ee = 2:nSlots
    Posiz_cave(ee) = Posiz_cave(ee-1) + Passo_elt;
    while Posiz_cave(ee)>2*pi
        Posiz_cave(ee)=Posiz_cave(ee)-2*pi;
    end
end

% Seleziono i vettori della stella degli avvolgimenti appartenenti alla fase 1 e ne eseguo la somma vettoriale
% La fase del vettore risultante identifica la posizione dell'asse della fase 1

kw = zeros(1,n3ph);
th0 = zeros(1,n3ph);

for ss=1:n3ph
    num = 0;
    den = 0;
    iph = 1+3*(ss-1);
    for tt=1:nLayers
        for yy=1:nSlots
            if WinTable(tt,yy)==iph || WinTable(tt,yy)==-iph
                den=den+1;
                num = num+WinTable(tt,yy)/iph*exp(1i*(Posiz_cave(yy)+0.5*pi));
            end
        end
    end
    
    kw(ss)  = abs(num)/den; % winding factor for each 3-phase set
    th0(ss) = angle(num);   % angular position of the axis of first phase of the 3-phase set respect to the x-axis
end

th0 = th0*180/pi;   % elec rad --> elec deg