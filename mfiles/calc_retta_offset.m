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

function [a,b,c] = calc_retta_offset(a0,b0,c0,d)
% calcola l'equazione della retta parallela, con distanza d.
% Se d<0 --> l'offset è fatto verso l'alto
% Se d>0 --> l'offset è fatto verso il basso

a = a0;
b = b0;
c = c0+d.*(a.^2+b.^2).^0.5;