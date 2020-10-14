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

%% 02/05/2013 MG
%% Circonferenza(incognita) tangente ad una retta (nota)
% funzione : si determina la posizione del centro della circonferenza il raggio, il pto di tangenza ad una retta nota:
% vincoli si impone il luogo del centro della circonferenza mediante le coordinate di un punto noto (xp,yp) appartenente 
% alla retta su cui si ipotizza giacia il centro della circonferenza; la
% retta richiede in ingresso 2 pti per cui passa
% output: xt2 yt2 pto di tangenza; x2 y2 r2: centro e raggio della
% circonferenza incognita.

function [xt2,yt2,x2,y2,r2]=tg_cir(xr1,yr1,xr2,yr2,xp,yp)
% keyboard
[a,b,c]=retta_per_2pti(xr1,yr1,xr2,yr2);
mp=yp/xp;

A=(- (4*a^2 + 4*b^2)*(b^2 - b^2*(mp^2 + 1) + b^2*mp^2) + (2*b^2 - 2*a*b*mp)^2);
B= (- ((2*xp + 2*mp*yp)*b^2 + 2*c*mp*b)*(4*a^2 + 4*b^2) - 4*a*c*(2*b^2 - 2*a*b*mp));
C=(4*a^2 + 4*b^2)*(b^2*(xp^2 + yp^2) - c^2) + 4*a^2*c^2;

x2_vect=roots([A,B,C]);
x2=min(x2_vect);
y2=mp*x2;
r2=sqrt((xp-x2)^2+(yp-y2)^2);

xt2=-(2*a*c-2*b^2*x2+2*a*b*y2)/(2*(a^2+b^2));
yt2=-(a/b)*xt2-(c/b);

