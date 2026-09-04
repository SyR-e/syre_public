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

flagCharger = false;

if(~isnan(per.rotorPos))
    flagCharger = true;
end

% currentDir=pwd();
pathnameIn = pathname;
[~,pathname]=createTempDir();
copyfile(fullfile(pathnameIn, strcat(strrep(filename,'.mat','.jmag'),'.jproj')),fullfile (pathname, strcat(strrep(filename,'.mat','.jmag'),'.jproj'))); % copy .jproj in the temporary folder

[SOL] = simulate_xdeg_JMAG(geo,per,pathname,filename);

% Post processing
if(flagCharger)
    if(isnan(geo.win.n3phase))
        out.id1     = mean(SOL.id1);
        out.iq1     = mean(SOL.iq1);
        out.id3     = mean(SOL.id3);
        out.iq3     = mean(SOL.iq3);
        out.fd1     = mean(SOL.fd1);
        out.fq1     = mean(SOL.fq1);
        out.fd3     = mean(SOL.fd3);
        out.fq3     = mean(SOL.fq3);
        out.IPF    = sin(atan(out.iq1./out.id1)-atan(out.fq1./out.fd1));
    else
        for i=1:geo.win.n3phase
            out.id(:,i)     = mean(SOL.id(:,i));
            out.iq(:,i)     = mean(SOL.iq(:,i));
            out.fd(:,i)     = mean(SOL.fd(:,i));
            out.fq(:,i)     = mean(SOL.fq(:,i));
        end

        out.IPF    = sin(atan(out.iq(:,1)./out.id(:,1))-atan(out.fq(:,1)./out.fd(:,1)));
    end
else
    if(isnan(geo.win.n3phase))
        out.id1     = mean(SOL.id1);
        out.iq1     = mean(SOL.iq1);
        out.id3     = mean(SOL.id3);
        out.iq3     = mean(SOL.iq3);
        out.fd1     = mean(SOL.fd1);
        out.fq1     = mean(SOL.fq1);
        out.fd3     = mean(SOL.fd3);
        out.fq3     = mean(SOL.fq3);
        out.IPF    = sin(atan(out.iq1./out.id1)-atan(out.fq1./out.fd1));
    else
        out.id     = mean(SOL.id);
        out.iq     = mean(SOL.iq);
        out.fd     = mean(SOL.fd);
        out.fq     = mean(SOL.fq);
        out.IPF    = sin(atan(out.iq./out.id)-atan(out.fq./out.fd));
    end
end

out.T      = abs(mean(SOL.T));
out.dT     = std(SOL.T);
out.dTpu   = std(SOL.T)/out.T;
out.dTpp   = max(SOL.T)-min(SOL.T);
out.Ppm    = sum(abs(mean(SOL.Ppm)));
out.Pfes_h = SOL.Pfes_h;
out.Pfer_h = SOL.Pfer_h;
out.Pfes_c = SOL.Pfes_c;
out.Pfer_c = SOL.Pfer_c;
out.Pfe    = out.Pfes_h + out.Pfes_c + out.Pfer_h + out.Pfer_c;
out.velDim = per.EvalSpeed;
out.SOL    = SOL;

end
