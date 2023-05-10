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

% compute a point belonging to a line (y=mx+q) and with a given distance
% (W) from a specific point (x1,y1)

function [x,y] = calc_pointDistanceLine(m,q,W,x1,y1)
A = 1 + m.^2;
B = 2*m.*q-2*y1.*m-2*x1;
C = x1.^2+y1.^2+q.^2-2*y1.*q-W.^2;
x = roots([A,B,C]);
y = m.*x+q;

