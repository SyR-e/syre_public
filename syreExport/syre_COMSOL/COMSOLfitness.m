% Copyright 2024
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

function [geo,mat,out,pathname] = COMSOLfitness(RQ,geo,per,mat,eval_type,pathname,filename)

pathnameIn = pathname;
[~,pathname]=createTempDir();
copyfile([pathnameIn strrep(filename,'.mat','.mph')],[pathname strrep(filename,'.mat','.mph')]); % copy .mph in the temporary folder

[SOL] = simulate_xdeg_COMSOL(geo,per,eval_type,pathname,filename);

out.id = SOL.id;
out.iq = SOL.iq;
out.fd = SOL.fd;
out.fq = SOL.fq;
out.T  = SOL.T;
out.IPF = sin(atan2(out.iq,out.id)-atan2(out.fq,out.fd));
out.th = SOL.th;
out.SOL = SOL;

out.ia = SOL.ia; 
out.ib = SOL.ib;
out.ic = SOL.ic;

out.fa = SOL.fa;
out.fb = SOL.fb;
out.fc = SOL.fc;

out.Ppm    = SOL.Ppm;
out.Pfes   = SOL.Pfes;
out.Pfer   = SOL.Pfer;

out.Pfe    = out.Pfes + out.Pfer;
out.velDim = per.EvalSpeed;

