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

function [OUT] = FEAfixSimulation(RQ,geo,per,mat,eval_type,filemot,gammaFix)

% load simulation (Torque)
if ~gammaFix
    [~,~,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
    OUT.fd = out.fd;
    OUT.fq = out.fq;
    OUT.id = out.id;
    OUT.iq = out.iq;
else
    deltaGamma = 3;
    numSim     = 7;
    
    gVect = (0:deltaGamma:deltaGamma*(numSim-1));
    gVect = gVect-mean(gVect);
    gVect = gVect+RQ(end);
    
    TVect  = nan(size(gVect));
    FdVect = nan(size(gVect));
    FqVect = nan(size(gVect));
    IdVect = nan(size(gVect));
    IqVect = nan(size(gVect));
    for ii=1:length(gVect)
        RQ(end) = gVect(ii);
        [~,~,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
        TVect(ii)  = out.T;
        FdVect(ii) = out.fd;
        FqVect(ii) = out.fq;
        IdVect(ii) = out.id;
        IqVect(ii) = out.iq;
    end
    
    [~,index] = max(TVect);
    OUT.fd = FdVect(index);
    OUT.fq = FqVect(index);
    OUT.id = IdVect(index);
    OUT.iq = IqVect(index);
    
end

% Characteristic current (Vtype and SPM)
if (strcmp(geo.RotType,'SPM')||strcmp(geo.RotType,'Vtype'))
    % characteristic current
    RQ(end) = 180;
    [~,~,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
    OUT.f0 = out.fd;
    % no-load simulation
    per.overload=0;
    RQ(end)=0;
    [~,~,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
    OUT.fM = out.fd;
else
    OUT.fM = 0;
    OUT.f0 = 0;
end

% Flux density in airgap, tooth and stator yoke (debug mode)
if strcmp(eval_type,'flxdn')
    OUT.Bt = max(max(out.SOL.Bt(:,2:end)));
    OUT.By = max(max(out.SOL.By(:,2:end)));
    if rem(geo.ps,2)~=0
        Bg = [out.SOL.Bg(:,2:end);-out.SOL.Bg(:,2:end)];
    else
        Bg = out.SOL.Bg(:,2:end);
    end
    a=fft(Bg,2^nextpow2(length(Bg(:,1))),1);
    harm=2*abs(a(2,:))/r;
    OUT.Bg = mean(harm);
else
    OUT.Bt = 0;
    OUT.By = 0;
    OUT.Bg = 0;
end


