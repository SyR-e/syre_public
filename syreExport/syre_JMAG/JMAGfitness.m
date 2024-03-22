% Copyright 2024
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

function [geo,mat,out,pathname] = JMAGfitness(RQ,geo,per,mat,pathname,filename)

% currentDir=pwd();
pathnameIn = pathname;
[~,pathname]=createTempDir();
copyfile(fullfile(pathnameIn, strcat(strrep(filename,'.mat','.jmag'),'.jproj')),fullfile (pathname, strcat(strrep(filename,'.mat','.jmag'),'.jproj'))); % copy .jproj in the temporary folder


[SOL] = simulate_xdeg_JMAG(geo,per,pathname,filename);

% Post processing
out.id     = mean(SOL.id);
out.iq     = mean(SOL.iq);
out.fd     = mean(SOL.fd);
out.fq     = mean(SOL.fq);
out.T      = abs(mean(SOL.T));
out.dT     = std(SOL.T);
out.dTpu   = std(SOL.T)/out.T;
out.dTpp   = max(SOL.T)-min(SOL.T);
out.IPF    = sin(atan(out.iq./out.id)-atan(out.fq./out.fd));
out.Ppm    = sum(abs(mean(SOL.Ppm)));
out.Pfes_h = SOL.Pfes_h;
out.Pfer_h = SOL.Pfer_h;
out.Pfes_c = SOL.Pfes_c;
out.Pfer_c = SOL.Pfer_c;
out.Pfe    = out.Pfes_h + out.Pfes_c + out.Pfer_h + out.Pfer_c;
out.velDim = per.EvalSpeed;
out.SOL    = SOL;

end
