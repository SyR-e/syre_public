% Copyright 2021
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

function [map] = SDE_scaleNs(map,kN)

map.id    = map.id./kN;
map.iq    = map.iq./kN;
map.i0    = map.i0./kN;
map.iAmp  = map.iAmp./kN;
map.fd    = map.fd.*kN;
map.fq    = map.fq.*kN;
map.Lbase = map.Lbase.*kN.^2;
map.Lmd   = map.Lbase.*kN.^2;
map.Rs    = map.Rs.*kN.^2;
map.fM    = map.fM.*kN;
map.iHWC  = map.iHWC./kN;
map.ich   = map.ich./kN;

if numel(kN)==1
    map.dataSet.TurnsInSeries = map.dataSet.TurnsInSeries*kN;
end