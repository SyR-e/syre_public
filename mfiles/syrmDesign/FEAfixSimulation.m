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

if strcmp(mat.LayerMag.MatName,'Air')
    flagPM=0;
else
    flagPM=1;
end

nFEA = 0;

% load simulation (Torque)
if ~gammaFix
    [~,geo,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
    nFEA = nFEA+1;
    OUT.fd  = out.fd;
    OUT.fq  = out.fq;
    OUT.id  = out.id;
    OUT.iq  = out.iq;
    OUT.mPM = calcMassPM(geo,mat);
    OUT.mCu = calcMassCu(geo,mat);
else
    maxIter   = 20;
    gammaStep = 2;
    % max angle from initial: 36 elt deg
    direction = 0;
    
    gamma0 = RQ(end);

    ii = 1;

    done = 0;

    TVect  = nan(1,maxIter);
    FdVect = nan(1,maxIter);
    FqVect = nan(1,maxIter);
    IdVect = nan(1,maxIter);
    IqVect = nan(1,maxIter);
    gVect  = nan(1,maxIter);

    while ~done
        if ii==1
            gammaSim = gamma0;
        elseif ii==2
            gammaSim = gamma0+gammaStep;
        elseif ii==3
            gammaSim = gamma0-gammaStep;
        else
            gammaSim = gammaSim+direction*gammaStep;
        end
        
        RQ(end) = gammaSim;
        [~,geo,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
        nFEA = nFEA+1;
        TVect(ii)  = out.T;
        FdVect(ii) = out.fd;
        FqVect(ii) = out.fq;
        IdVect(ii) = out.id;
        IqVect(ii) = out.iq;
        gVect(ii)  = gammaSim;
        
        if ii==3
            [~,index] = max(TVect,[],'omitnan');
            if index==1
                done=1;
            elseif index==2
                direction=+1;
            else
                direction=-1;
            end
        elseif ii>3
            if TVect(ii)<TVect(ii-1)
                done=1;
            end
        end

        if ii==maxIter
            done=1;
        end

        ii = ii+1;
    end

    [~,index] = max(TVect,[],'omitnan');
    OUT.fd  = FdVect(index);
    OUT.fq  = FqVect(index);
    OUT.id  = IdVect(index);
    OUT.iq  = IqVect(index);
    OUT.mPM = calcMassPM(geo,mat);
    OUT.mCu = calcMassCu(geo,mat);


%     deltaGamma = 3;
%     numSim     = 7;
%     
%     gVect = (0:deltaGamma:deltaGamma*(numSim-1));
%     gVect = gVect-mean(gVect);
%     gVect = gVect+RQ(end);
%     
%     TVect  = nan(size(gVect));
%     FdVect = nan(size(gVect));
%     FqVect = nan(size(gVect));
%     IdVect = nan(size(gVect));
%     IqVect = nan(size(gVect));
%     for ii=1:length(gVect)
%         RQ(end) = gVect(ii);
%         [~,~,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
%         TVect(ii)  = out.T;
%         FdVect(ii) = out.fd;
%         FqVect(ii) = out.fq;
%         IdVect(ii) = out.id;
%         IqVect(ii) = out.iq;
%     end
%     
%     [~,index] = max(TVect);
%     OUT.fd = FdVect(index);
%     OUT.fq = FqVect(index);
%     OUT.id = IdVect(index);
%     OUT.iq = IqVect(index);
%     
end



% Characteristic current (Vtype and SPM)
if (strcmp(geo.RotType,'SPM')||strcmp(geo.RotType,'Vtype'))
    % characteristic current
    RQ(end) = 180;
    [~,~,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
    nFEA = nFEA+1;
    OUT.f0 = out.fd;
%     % no-load simulation
%     per.overload=0;
%     RQ(end)=0;
%     [~,~,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
%     OUT.fM = out.fd;
else
%     OUT.fM = 0;
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


if flagPM
    % no-load simulation
    per.overload=0;
    RQ(end)=0;
    [~,~,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
    nFEA = nFEA+1;
    OUT.fM = -out.fq;
else
    OUT.fM = 0;
end

OUT.nFEA = nFEA;
