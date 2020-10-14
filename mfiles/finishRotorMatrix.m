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

function [rotore] = finishRotorMatrix(rotore,geo)

% input: half pole description of rotor geometry
% output: full rotor description (ps poles)

ps = geo.ps;
p  = geo.p;
r  = geo.r;
Ar = geo.Ar;
lm = geo.lm;
RotType = geo.RotType;


% Check the rotor matrix, to avoid plot error due to arcs
[rotore]=checkPlotMatrix(rotore,1e-9);

if not(isempty(rotore))
    % mirror the rotor points
    rotNeg=rotore;
    rotNeg(:,[2 4 6 7])=-rotore(:,[2 4 6 7]);
    rotore=[rotore;rotNeg];
    
    % replicate the pole ps-1 times (fractional slots have ps>1)
    Temp2=[];
    num_poli=1;
    while num_poli<ps
        Temp1=[];
        for kk=1:2:size(rotore,2)-2
            [xtemp,ytemp]=rot_point(rotore(:,kk),rotore(:,kk+1),num_poli*180/p*pi/180);
            Temp1=[Temp1,xtemp,ytemp];
        end
        Temp2=[Temp2;[Temp1,rotore(:,end)]];
        num_poli=num_poli+1;
    end
    rotore=[rotore;Temp2];
    %     clear Temp1 Temp2 xtemp ytemp;
end
%
% raggio esterno
if strcmp(RotType,'SPM')
    xre1 = r-lm;     yre1 = 0;
else
    xre1 = r;        yre1 = 0;
end
if (ps<2*p)
    [xre3,yre3] = rot_point(xre1,yre1,(ps-1/2)*180/p*pi/180);
else
    [xre3,yre3] = rot_point(xre1,yre1,(ps/2-1/2)*180/p*pi/180);
end
[xre2,yre2] = rot_point(xre1,yre1,-90/p*pi/180);

if (ps<2*p)
    rotore=[rotore;
        0 0 xre2 yre2 xre3 yre3 1];
else
    rotore=[rotore;
        0 0 xre2 yre2 xre3 yre3 1;
        0 0 xre3 yre3 xre2 yre2 1];
end

% Albero
if (ps<2*p)
    [xAl1,yAl1] = rot_point(Ar,0,(ps-1/2)*180/p*pi/180);
else
    [xAl1,yAl1] = rot_point(Ar,0,(ps/2-1/2)*180/p*pi/180);
end
[xAl2,yAl2] = rot_point(Ar,0,-90/p*pi/180);

if (ps<2*p)
    rotore=[rotore;
        0,0,xAl2,yAl2,xAl1,yAl1,1];
else
    rotore=[rotore;
        0,0,xAl2,yAl2,xAl1,yAl1,1;
        0,0,xAl1,yAl1,xAl2,yAl2,1];
end

% Completo il rotore con linee di chiusura laterali
if (ps<2*p)
    rotore=[rotore;
        0,0,xre2,yre2,NaN,NaN,0;
        0,0,xre3,yre3,NaN,NaN,0];
end