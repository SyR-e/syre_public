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

function [out,dirName] = skinEffect_point(tRef,fRef,geo,per,mat,filename)


[thisfilepath,dirName]=createTempDir();

copyfile(filename,[dirName 'slot_0.fem'])

load([filename(1:end-4) '.mat'],'slot');

nCond = slot.nCond;

iRef = per.i0;

if ~geo.XFEMMsimulation
    openfemm(1)
    opendocument([dirName 'slot_0.fem']);
else
    FemmProblem = loadfemmfile(checkPathSyntax([dirName 'slot_0.fem']));
    femmFileName = [dirName 'slot_0.fem'];
    FemmProblem.ProbInfo.SmartMesh = 0;
end

for ii=1:nCond
    if ~geo.XFEMMsimulation
        mi_modifycircprop(['conductor_' int2str(ii)],1,iRef/nCond);
    else
        FemmProblem = setcircuitcurrent(FemmProblem,['conductor_' int2str(ii)],iRef/nCond);
    end
end

if ~isfield(mat.SlotCond,'alpha')
    mat.SlotCond.alpha = 0.004;
    warning('Resistivity thermal coefficient not set')
end

rho20 = 1/mat.SlotCond.sigma;
rho = rho20*(1+mat.SlotCond.alpha*(tRef-20));

if ~geo.XFEMMsimulation
    mi_modifymaterial(mat.SlotCond.MatName,5,1/rho/1e6);
    mi_probdef(fRef);
else
    FemmProblem = modifymaterial_mfemm(FemmProblem,mat.SlotCond.MatName,'Sigma',1/rho/1e6);
    FemmProblem.ProbInfo.Frequency = fRef;
end

if ~geo.XFEMMsimulation
    mi_analyze(1);
    mi_loadsolution();
else
    writefemmfile(checkPathSyntax(femmFileName), FemmProblem);
        femmFileName = fmesher(checkPathSyntax(femmFileName));
        ansFile = fsolver(femmFileName,0,1);
        myfpproc = fpproc();
        myfpproc.opendocument(ansFile);
end

% post processing


out.T = tRef;
out.f = fRef;
out.I = iRef;

Ptmp = 0;
Qtmp = 0;

for ii=1:nCond
    if ~geo.XFEMMsimulation
        tmp = mo_getcircuitproperties(['conductor_' int2str(ii)]);
    else
        tmp = myfpproc.getcircuitprops(['conductor_' int2str(ii)]);
    end
    Ptmp = Ptmp+real(tmp(1)*tmp(2));
    Qtmp = Qtmp+imag(tmp(1)*tmp(2));
end

out.P = Ptmp;
out.Q = Qtmp;

if ~geo.XFEMMsimulation
    mo_close;
    mi_close;
    closefemm;
end



