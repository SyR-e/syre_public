% Copyright 2023
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

function [d] = calc_distanza_retta_punto(a,b,c,x0,y0)
% calcola la distanza del punto (x0,y0) dalla retta a*x+b*y+c=0

% calcolo retta perpendicolare
bp = -1;
ap = b/a;
cp = -ap*x0-bp*y0;

[xTmp,yTmp] = intersezione_tra_rette(a,b,c,ap,bp,cp);

d = ((x0-xTmp)^2+(y0-yTmp)^2)^0.5;