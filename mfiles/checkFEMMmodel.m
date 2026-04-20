% Copyright 2026
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

function [out] = checkFEMMmodel(femmFile,solver)

femmFile = checkPathSyntax(femmFile);

[~,filename,ext] = fileparts(femmFile);
filename = [filename ext]; % fem file name

[~,pathname]=createTempDir();

testFile = [pathname filename ext];

copyfile(femmFile,testFile);

if strcmp(solver,'FEMM')
    openfemm(1);
    opendocument(testFile);
else
    FemmProblem = loadfemmfile(testFile);
    FemmProblem.ProbInfo.SmartMesh = 0;
    FemmProblem = setcircuitcurrent(FemmProblem,'fase1',1);
    writefemmfile(femmFile,FemmProblem);
end

out = 1;

try
    if strcmp(solver,'FEMM')
        mi_analyze(1);
    else
        testFile = fmesher(testFile);
        ansFile = fsolver(testFile,0,1);
    end
catch
    out = 0;
end

if strcmp(solver,'FEMM')
    closefemm
end




    

