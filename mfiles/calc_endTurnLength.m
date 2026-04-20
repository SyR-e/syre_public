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

function [lend,geo] = calc_endTurnLength(geo)
% Computation of the end-winding length

wt = geo.wt;
r  = geo.r;
g  = geo.g;
lt = geo.lt;
p  = geo.p;
q  = geo.q;
n3phase = geo.win.n3phase;
avv=geo.win.avv;


if geo.q<1
    %concentrated winding
    if(isnan(n3phase))
        nphases = 5;
        lend=0.5.*(wt+pi*(r+g+lt/2).*sin(pi./(2*p.*q.*nphases)))/1e3;    % [m] from Gamba - A new PMASRM with nonconventional FS pole combination
    else
        lend=0.5.*(wt+pi*(r+g+lt/2).*sin(pi./(6*p.*q.*n3phase)))/1e3;    % [m] from Gamba - A new PMASRM with nonconventional FS pole combination
    end
else
    % distributed winding
    [~,c]=size(avv);
    ini=0;  % initial slot
    fin=0;  % final slot
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

    if(isnan(n3phase))
        nphases = 5;
        alpha=2*pi*yq./(2.*p.*q.*nphases);  % coil pitch (rad)
    else
        alpha=2*pi*yq./(6.*p.*q.*n3phase);  % coil pitch (rad)
    end

    lend=(2*lt+(r+g+lt/2).*alpha)/1e3; %[m]
end

lend = lend*1e3; %m-->mm
geo.lend = lend;

