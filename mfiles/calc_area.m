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

%Calcolo area mediante prodotto vettoriale:
function [A]=calc_area(x1,x2,x3,y1,y2,y3)
vx=x2-x1; vy=y2-y1;
wx=x3-x1; wy=y3-y1;
A=abs(1/2*det([vx vy;wx wy]));