% Copyright 2022
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [xf,yf] = proiezione_punto_retta(a,b,c,x0,y0)
% Proiezione del punto (x0,y0) sulla retta a*x+b*y+c=0

if a~=0
    yf = (a^2*y0-b*c-a*b*x0)/(a^2+b^2);
    xf = (-b*yf-c)/a;
else
    xf = (b^2*x0-a*c-a*b*y0)/(a^2+b^2);
    yf = (-a*xf-c)/b;
end


% yf = y0-b/a*x0-(b-b*c)/a^2;
% xf = -(b*yf-c)/a;
