% Copyright 2025
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

function [out] = clearTmpFolder(flagDisp)

if nargin()==0
    flagDisp = 1;
end

try
    rmdir(checkPathSyntax([cd,'\tmp']),'s');
    if flagDisp==1
        msgbox('directory empty');
    elseif flagDisp==2
        disp('tmp directory empty');
    end
catch
    if flagDisp==1
        msgbox('directory empty');
    elseif flagDisp==2
        disp('tmp directory empty');
    end
end
if exist(checkPathSyntax([cd,'\tmp']),'dir') == 0
    mkdir(checkPathSyntax([cd,'\tmp']));
end

if nargout()==0
    clear out
else
    out=1;
end