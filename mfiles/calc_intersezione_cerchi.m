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

% calc_intersezione_cerchi.m
% calcola l'intersezione tra circonferenza centrata in 0,0 di raggio r e
% circonferenza centrata in x0,0 di raggio r1
% input: raggi r ed r1
% output: x,y del punto

function [x,y] = calc_intersezione_cerchi(r,r1,x0)

x = (r.^2 - r1.^2 + x0^2)/(2*x0);
gamma = acos(x ./ r);
y = r .* sin(gamma);