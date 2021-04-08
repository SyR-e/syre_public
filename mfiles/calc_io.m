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


function [i0,Rs,geo] = calc_io(geo,per)

% calc_io.m
% data la geometria e la dissipazione ammessa calcola l'ampiezza del vettore corrente
% input:
% - Kj : W /m2 dissipabili sullo statore
% - g : geometria della macchina (+ car. avvolgimento)
% output:
% - io (A)
% - Rs (Ohm)

% ref. Tutorial Course Notes: Design, Analysis and Control of Interior PM Synchronous Machines
% IEEE IAS Annual Meeting 2004 - Seattle
% cap. 6 - Vagati

% input
loss = per.Loss; % [W]
R = geo.R/1e3; % [m]
l = geo.l/1e3; % [m]
kcu = geo.win.kcu;
N=geo.win.Ns;
% rocu = 17.8*(234.5+per.tempcuest)/(234.5+20)*1e-9;
rocu = 17.8*(234.5+per.tempcu)/(234.5+20)*1e-9;
n3phase=geo.win.n3phase; %AS
Aslots = geo.Aslot*(6*geo.p*geo.q*n3phase)/1e6; % [m^2]

% end turn length
if geo.q<1
    %concentrated winding
    lend=0.5*(geo.wt+pi*(geo.r+geo.g+geo.lt/2)*sin(pi/6/geo.p/geo.q*n3phase))/1e3;    % [m] from Gamba - A new PMASRM with nonconventional FS pole combination
else
    % distributed winding
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
    alpha=2*pi*yq/(6*geo.p*geo.q*n3phase);  % coil pitch (rad)
    
%     if strcmp(geo.win.nCond,'round')
        lend=(2*geo.lt+(geo.r+geo.g+geo.lt/2)*alpha)/1e3; %[m]
%     else
%         lend=(2*geo.wt+(geo.r+geo.g+geo.lt/2)*alpha/cos(atan(1/3)))/1e3; %[m]
%     end
end

kj=loss/(2*pi*R*l);

i0 = (kj*kcu/rocu*l/(l+lend)*2*pi*R*2*Aslots/(4*(n3phase*3)^2)/N^2)^0.5; %AS

Rs = loss/(n3phase*3/2*i0^2); %AS





