% Copyright 2019
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


function Lend = calc_Lend(geo)

% calc_Lend.m
% end-winding leakage inductance 

Ns = geo.win.Ns;

if geo.q<1
    % concentrated winding
    warning('Lend for concentrated windings set to 0 .. to be completed!');
    Lend=0;
else
    if geo.win.kracc == 1
        % full pitch (one layer)
        Nm = geo.q * geo.p;     % number of coils per phase
        As = geo.Aslot*1e-6;         % cross-sectional area of the coil (m2)
    else
        % short pitch (two layers)
        Nm = geo.q * geo.p * 2;
        As = 0.5*geo.Aslot*1e-6;
    end
    
    N = Ns/Nm; % turns per coil
  
    avv=geo.win.avv;
    [r,c]=size(avv);
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
    alpha=2*pi*yq/(6*geo.p*geo.q*geo.win.n3phase);  % coil pitch in radians
    taucp = (geo.r+geo.g+geo.lt/2)*alpha*1e-3;   % mean coil pitch (m)
    
    mu0 = 4e-7*pi;
    % Hanselman (from MOTORCAD)
    LendHANS = Nm*mu0*taucp*N^2/2*log(taucp*sqrt(pi)/sqrt(2*As));
    % Rosa and Grover (from MOTORCAD)
    LendROSA = Nm*mu0*taucp*N^2/2*(log(taucp*pi*sqrt(pi)/sqrt(2*As))-0.75);
    % Lipo (Introduction to AC machines Design - 3rd edition)
    a = geo.lt*1e-3;    % overhang length [mm]
    b = taucp;          % coil pitch
    NumCondSlot = round(geo.win.Nbob * 2);  % conductors in one slot
    Acond = geo.Aslot*geo.win.kcu/NumCondSlot*1e-6;  % conductor cross section (m2)
    epsilon = sqrt(Acond/pi);   % radius of the conductor or winding bundle (m)
    Lew1 = mu0/pi*(...
        -2*a*log(2*a/b+sqrt((2*a/b)+1))...
        -b*log(b/(2*a)+sqrt((b/(2*a)+1)))...
        +2*a*log(4*a/epsilon)+b*log(2*b/epsilon)+2*sqrt(4*a^2+b^2)-2*b-4*a); % one turn, both ends
    Linternal = mu0/8/pi*(4*a+2*b); % internal to the conductor or bundle
    LendLIPO = (Lew1+Linternal)*Nm*N^2;
    
    % as LIPO is the most conservative, this is the one in use
    Lend = LendLIPO;
end

