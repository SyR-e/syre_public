 % Copyright 2025
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

function [lendf,geo] = calc_endTurnFieldLength(geo)
% computation of the field coil end-winding length


wp = geo.wp;
wb = geo.wb;

% lendf = pi*(wp/2+wb/2);
% lendf = wp + 2*wb  % -> wb/2 + (wb/2+wp+wb/2) + wb/2 (rettangolo)
lendf = wp+pi*wb/2; % due quarti di cerchio più ampiezza polo


geo.lendf = lendf;


