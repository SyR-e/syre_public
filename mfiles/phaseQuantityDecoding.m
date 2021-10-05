% Copyright 2021
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

function [out] = phaseQuantityDecoding(xa,xb,xc,xdeg)


switch xdeg
    case 360
        out.a = xa;
        out.b = xb;
        out.c = xc;
    case 180
        out.a = [xa -xa];
        out.b = [xb -xb];
        out.c = [xc -xc];
    case 120
        out.a = [xa xc xb];
        out.b = [xb xa xc];
        out.c = [xc xb xa];
    case 60
        out.a = [xa -xb xc -xa xb -xc];
        out.b = [xb -xc xa -xb xc -xa];
        out.c = [xc -xa xb -xc xa -xb];
    otherwise
        out.a = nan(size(xa));
        out.b = nan(size(xb));
        out.c = nan(size(xc));
        out.a = nan(1,round(360/xdeg*size(xa,1)));
end
