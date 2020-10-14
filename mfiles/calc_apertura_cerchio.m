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

% calc_apertura_cerchio.m
% calcola l'apertura angolare rispetto al centro (x0,0) di un punto con
% coordinate polari (r, alpha) rispetto al centro (0,0)
% input: raggi r, alpha, x0
% output: angolo rispetto a (x0,0)

function beta = calc_apertura_cerchio(alpha,r,x0)

beta = atan(r * sin(alpha) ./ (x0 - r * cos(alpha)));