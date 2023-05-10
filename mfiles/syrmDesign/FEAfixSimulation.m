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

function [OUT] = FEAfixSimulation(RQ,geo,per,mat,eval_type,filemot,gammaFix,flagIch,flagSC,flagDemag0,flagDemagHWC)

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
    OUT.T   = out.T;
    gamma0 = RQ(end);
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
    gamma0 = gVect(index);
    OUT.T   = TVect(index);

end

per = calc_i0(geo,per,mat);

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

if flagIch
    MaxIter = 10;
    i0 = per.i0;
    if strcmp(geo.axisType,'PM')
        per.gamma = 180;
        RQ(end) = 180;
    else
        per.gamma = 90;
        RQ(end) = 90;
    end
    ichTest = nan(1,MaxIter);
    FmTest  = nan(1,MaxIter);

    done=0;
    ii=1;
    while ~done
        if ii==1
            ichTest(ii)=0;
            tol=-inf;
        elseif ii==2
            ichTest(ii) = 1;
            tol=FmTest(1)/50;
        else
            ichTest(ii)=ichTest(ii-2)-FmTest(ii-2)/(FmTest(ii-1)-FmTest(ii-2))*(ichTest(ii-1)-ichTest(ii-2));
            tol=FmTest(1)/50;
        end
        if ii==1
            FmTest(ii) = OUT.fM;
        else
            per.overload = ichTest(ii);
            %per.BrPP     = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,per.tempPP);
            [~,~,~,out,~] = FEMMfitness(RQ,geo,per,mat,'singt',filemot);

            if strcmp(geo.axisType,'PM')
                FmTest(ii) = out.fd;
            else
                FmTest(ii) = -out.fq;
            end
        end
        if abs(FmTest(ii))>tol
            ii=ii+1;
            done=0;
        else
            done=1;
        end
        if ii>MaxIter
            done=1;
            ii=ii-1;
        end
    end
        OUT.ich = ichTest(ii)*i0;
    else
        OUT.ich = 0;
end

if flagSC
    setup.RQ               = RQ;
    setup.RQ(end)          = 90;
    setup.flagSave         = 0;
    setup.flagFEAfix       = 1;
    setup.filemot          = filemot;
    setup.geo              = geo;
    setup.mat              = mat;
    setup.per              = per;
    setup.idq0             = OUT.id+j*OUT.iq;
    setup.fdq0             = OUT.fd+j*OUT.fq;
    setup.T0               = OUT.T;

    [pkSCout] = eval_peakShortCircuitCurrent(setup);

    OUT.iHWC = abs(pkSCout.idq);
else
    OUT.iHWC = NaN;
end

if flagDemagHWC
    per.overload = OUT.iHWC/i0;
    if (strcmp(geo.RotType,'SPM')||strcmp(geo.RotType,'Vtype'))
        per.gamma = 180;
    else
        per.gamma = 90;
    end
    RQ(end) = per.gamma;
    [~,~,~,out,~] = FEMMfitness(RQ,geo,per,mat,'demagArea',filemot);
    OUT.dPMHWC  = out.dPM;
    OUT.BminHWC = out.Bmin;
else
    OUT.dPMHWC  = NaN;
    OUT.BminHWC = NaN;
end
if flagDemag0
    per.overload = 1;
    if (strcmp(geo.RotType,'SPM')||strcmp(geo.RotType,'Vtype'))
        per.gamma = 180;
    else
        per.gamma = 90;
    end
    RQ(end) = per.gamma;
    [~,~,~,out,~] = FEMMfitness(RQ,geo,per,mat,'demagArea',filemot);
    OUT.dPM0  = out.dPM;
    OUT.Bmin0 = out.Bmin;
else
    OUT.dPM0  = NaN;
    OUT.Bmin0 = NaN;
end



OUT.nFEA = nFEA;
