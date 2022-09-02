% Copyright 2016
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

function Mass = calcMassCu(geo,mat)
% 
% Mass = calcMassCu(geo)
% 
% Evaluate the copper mass of the stator

% SF - 07/07/2016

% Data
rhoCu=8940; % [kg/m^3]
rhoCu = mat.SlotCond.kgm3;  % [kg/m^3]

% evaluation of lend

avv=geo.win.avv;
[r,c]=size(avv);
ini=0;
fin=0;
for ii=1:c
    if ini==0
        if avv(1,ii)==1 || avv(2,ii)==1
            ini=ii;
        end
    elseif fin==0
        if avv(1,ii)==-1 || avv(2,ii)==-1
            fin=ii;
        end
    end
end

if fin==0
    fin=c+1;
end

yq=fin-ini;

if yq==1  %concentrated winding
    lend=0.5*(geo.wt+pi*(geo.r+geo.g+geo.lt/2)*sin(pi/6/geo.p/geo.q));    % Gamba - A new PMASRM with nonconventional FS pole combination
else
    alpha=yq*2*pi/(6*geo.p*geo.q*geo.win.n3phase);
    lend=2*geo.lt+(geo.r+geo.g+geo.lt/2)*alpha;
    clear alpha
end
if isfield(geo,'pShape')
    if geo.pShape.flag
        geo.Aslot = area(geo.pShape.slot)/c;
    end
end

Mass=rhoCu*geo.Aslot/1e6*(geo.l+lend)/1000*(6*geo.p*geo.q*geo.win.n3phase)*geo.win.kcu; % [kg]


