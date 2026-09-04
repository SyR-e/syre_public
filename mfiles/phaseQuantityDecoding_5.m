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

function [out] = phaseQuantityDecoding(xa,xb,xc,xd,xe,xdeg)


switch xdeg
    case 360
        out.a = xa;
        out.b = xb;
        out.c = xc;
        out.d = xd;
        out.e = xe;
    case 180
        out.a = [xa -xa];
        out.b = [xb -xb];
        out.c = [xc -xc];
        out.d = [xd -xd];
        out.e = [xe -xe];
    case 72
        % out.a = [xa xb xc xd xe];
        % out.b = [xb xc xd xe xa];
        % out.c = [xc xd xe xa xb];
        % out.d = [xd xe xa xb xc];
        % out.e = [xe xa xb xc xd];
        out.a = [xa xe xd xc xb];
        out.b = [xb xa xe xd xc];
        out.c = [xc xb xa xe xd];
        out.d = [xd xc xb xa xe];
        out.e = [xe xd xc xb xa];
    case 36
        out.a = [xa -xc xe -xb xd -xa xc -xe xb -xd];
        out.b = [xb -xd xa -xc xe -xb xd -xa xc -xe];
        out.c = [xc -xe xb -xd xa -xc xe -xb xd -xa];
        out.d = [xd -xa xc -xe xb -xd xa -xc xe -xb];
        out.e = [xe -xb xd -xa xc -xe xb -xd xa -xc];
    otherwise
        out.a = nan(size(xa));
        out.b = nan(size(xb));
        out.c = nan(size(xc));
        out.d = nan(size(xd));
        out.e = nan(size(xe));
end
