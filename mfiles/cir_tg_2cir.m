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
%% determinazione pti di tangenza di una circonferenza interna a 2 circonferenze date
%%
% xp yp rappresenta il pto per cui passa la circonferenza incognita e
% tangnente alla C1.
% la retta (0,0)-->(xp,yp) è il luogo su cui giace il centro della
% circonferenza incognita
% x2 y2 r2: centro e raggio della 2° circonferenza considerata
% 0 0 r: 1° circonferenza ha centro nell'origine (il metodo nn è del tutto
% generale)
% output: xt2 yt2 pto di tangenza con C2; xc yc rc= centro e raggio della
% circonferenza incognita

function [xt2,yt2,xc,yc,rc]=cir_tg_2cir(xp,yp,r,x2,y2,r2) 

mp=yp/xp;

A=-(r2^2*(4*mp^2+4)-(2*x2-2*xp+2*mp*y2-2*mp*yp)^2);
B=-(2*(-r^2+r2^2+x2^2+y2^2)*(2*x2-2*xp+2*mp*y2-2*mp*yp)-r2^2*(8*x2+8*mp*y2));
C=(x2^2+y2^2+r2^2-r^2)^2-r2^2*(4*x2^2+4*y2^2);
xc_vect=roots([A,B,C]);
xc=min(xc_vect);
yc=mp*xc;
% [THETA,RHO] = cart2pol(xc,yc);
% THETA=THETA*180/pi;
% RHO;
rc=sqrt((1+mp^2)*xc^2-xc*(2*xp+2*mp*yp)+r^2);

a=2*(xc-x2); b=2*(yc-y2); c=(x2^2-xc^2)+(y2^2-yc^2)+(rc^2-r2^2);

xt2=-(2*a*c-2*b^2*x2+2*b*y2*a)/(2*(a^2+b^2));
yt2=-(a/b)*xt2-(c/b);
