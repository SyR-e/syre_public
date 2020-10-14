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


%% Vengono calcolati i pti di intersezione tra 2 circonferenze, le circonferenze in esame vengono assegnate mediante posizione del centro e raggio
% xc1..2 ,yc1..2 =coordinate del centro circonferenza; rc1..2=raggio.

function [x,y]=intersezione_cerchi(xc1,yc1,rc1,xc2,yc2,rc2)

a1=-2*xc1; b1=-2*yc1;
a2=-2*xc2; b2=-2*yc2;
c1=-rc1^2+xc1^2+yc1^2;
c2=-rc2^2+xc2^2+yc2^2;

%Formula intersezione cerchi:
% ((a1 - a2)^2/(b1 - b2)^2 + 1)*x^2 + (a1 - (b1*(a1 - a2))/(b1 - b2) + (2*(a1 - a2)*(c1 - c2))/(b1 - b2)^2)*x + c1 + (c1 - c2)^2/(b1 - b2)^2 - (b1*(c1 - c2))/(b1 - b2)
A=((a1 - a2)^2/(b1 - b2)^2 + 1);
B=(a1 - (b1*(a1 - a2))/(b1 - b2) + (2*(a1 - a2)*(c1 - c2))/(b1 - b2)^2);
C= c1 + (c1 - c2)^2/(b1 - b2)^2 - (b1*(c1 - c2))/(b1 - b2);

x1=(-B+sqrt(B^2-4*A*C))/(2*A);
x2=(-B-sqrt(B^2-4*A*C))/(2*A);
y1=-(a1-a2)/(b1-b2)*x1-(c1-c2)/(b1-b2);
y2=-(a1-a2)/(b1-b2)*x2-(c1-c2)/(b1-b2);
x=[x1;x2];
y=[y1;y2];