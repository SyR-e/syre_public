% Copyright 2019
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

function Mass = calcMassPM(geo,mat)

if strcmp(mat.LayerMag.MatName,'Air')
    Mass = 0;
else
    rhoPM = mat.LayerMag.kgm3;
    if strcmp(geo.RotType,'SPM')
        Mass = pi*(geo.r^2-(geo.r-geo.lm)^2)*geo.l*(geo.phi/180)/1e9*rhoPM;
    else
        Mass = (2*geo.p)*rhoPM*(geo.l/1000)*2*sum([geo.AreaC geo.AreaE])/1e6;
    end
end