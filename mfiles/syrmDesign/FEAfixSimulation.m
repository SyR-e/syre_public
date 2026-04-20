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

function [OUT] = FEAfixSimulation(RQ,geo,per,mat,eval_type,filemot,gammaFix,flagIch,flagSC,flagDemag0,flagDemagHWC,flagDemagUGO,flagMech,flagTherm)

if strcmp(mat.LayerMag.MatName,'Air')
    flagPM=0;
else
    flagPM=1;
end

nFEA = 0;

i0 = per.i0;
% per.flag_OptCurrConst = 0;
% per.nsim_singt = 2;
per.custom_act = 0;

% load simulation (Torque)
stp = false;
if ~gammaFix
    [~,geo,~,out,pathname] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
    nFEA = nFEA+1;
    OUT.fd  = out.fd;
    OUT.fq  = out.fq;
    OUT.id  = out.id;
    OUT.iq  = out.iq;
    OUT.mPM = calcMassPM(geo,mat);
    OUT.mCu = calcMassCu(geo,mat);
    OUT.T   = out.T;

    % gamma0 = RQ(end);
else
    if strcmp(geo.RotType,'EESM') % Metodo della sezione aurea
        

    else
        %Initialization
        maxIter = 5;
        GR       = (sqrt(5)-1)/2; % Golden ratio (0.618034)
        lock    = false;
        gammaVect = nan(1, 3 + maxIter);
        TVect  = nan(1,3 + maxIter);
        FdVect = nan(1,3 + maxIter);
        FqVect = nan(1,3 + maxIter);
        IdVect = nan(1,3 + maxIter);
        IqVect = nan(1,3 + maxIter);
        gVect  = nan(1,3 + maxIter);

        % First 3 fixed tests
        gammaVect(1:3) = [90, 120, 60];
        for ii = 1:3
            gammaSim = gammaVect(ii);
            RQ(end) = gammaSim;
            [~,geo,~,out,pathname] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
            nFEA = nFEA+1;
            TVect(ii)  = out.T;
            FdVect(ii) = out.fd;
            FqVect(ii) = out.fq;
            IdVect(ii) = out.id;
            IqVect(ii) = out.iq;
            gVect(ii)  = gammaSim;
        end

        % Interpolation loop
        for jj = 1:maxIter
            % Sorting gammas
            [g, ord] = sort(gammaVect(1:(3 + jj - 1)));
            T = TVect(ord);
            % Best index in the sorted arrays
            [~, best] = max(TVect(1:(3 + jj - 1)), [], 'omitnan');
            gbest = gammaVect(best);
            [~, index] = min(abs(g - gbest));
            brkt = (index > 1) && (index < numel(g)) && (T(index) >= T(index-1)) && (T(index) >= T(index+1));
            % 1 if max is between two limits

            if lock && index > 1 && index < numel(g)
                brkt = true;
            end

            % Next trial
            if index == 1 || index == numel(g) % maximum at edges - expand
                if index == 1
                    step = g(2) - g(1); if step <= 0, step = 1; end
                    gammaSim = g(1) - step;      % expand left edge
                else
                    step = g(end) - g(end-1); if step <= 0, step = 1; end
                    gammaSim = g(end) + step;    % expand right edge
                end
            elseif brkt % maximum in the middle - refinement inside
                gmin = g(index-1); gmid = g(index); gmax = g(index+1);
                Tmin = T(index-1); Tmid = T(index); Tmax = T(index+1);
                d = ((gmid - gmin)*(Tmid - Tmax) - (gmid - gmax)*(Tmid - Tmin));
                Par = false;
                if d ~= 0 % Parabolic
                    n = ((gmid - gmin)^2)*(Tmid - Tmax) - ((gmid - gmax)^2)*(Tmid - Tmin);
                    u = gmid - 0.5 * (n / d);   % candidate vertex
                    % Safeguards: keep strictly inside, not too close to existing
                    span = gmax - gmin;
                    samp = 0.25;
                    inside = (u > gmin + 0.10*span) && (u < gmax - 0.10*span);  % stay well within bracket
                    farenough = all(abs(u - g) > samp);                        % not too close to any sampled gamma
                    if isfinite(u) && inside && farenough
                        gammaSim = u;
                        Par = true;
                    end
                end
                if ~Par % Golden ratio + "cold-region"
                    dTl  = gmid - gmin; dTu = gmax - gmid; % Distance (in angle) from the middle
                    dgl = Tmid - Tmin; dgU = Tmid - Tmax;  % Distance (in torque) from the middle
                    if dgU >= 1.15*dgl % higher angles 115% further away
                        gammaSim = gmid - GR*dTl;
                    elseif dgl >= 1.15*dgU % lower angles 115% further away
                        gammaSim = gmid + GR*dTu;
                    else % safety
                        if dTu >= dTl
                            gammaSim = gmid + GR*dTu;
                        else
                            gammaSim = gmid - GR*dTl;
                        end
                    end
                end
                lock = true;
            else % midpoint toward the higher neighbor
                if T(index-1) > T(index+1)
                    gammaSim = 0.5*(g(index-1) + g(index));
                else
                    gammaSim = 0.5*(g(index) + g(index+1));
                end
            end
            gammaSim = round(gammaSim);  % force integer angle
            prevg = gammaVect(1:(3 + jj - 1));
            if any(abs(prevg - gammaSim) < 1e-9) % Repeated gamma: loop stops
                [~, final] = max(TVect(1:(3 + jj - 1)), [], 'omitnan');
                OUT.fd  = FdVect(final);
                OUT.fq  = FqVect(final);
                OUT.id  = IdVect(final);
                OUT.iq  = IqVect(final);
                OUT.mPM = calcMassPM(geo,mat);
                OUT.mCu = calcMassCu(geo,mat);
                OUT.T   = TVect(final);
                stp = true;
                break; % exit only the interpolation loop
            end
            % FEMM
            ii = 3 + jj;
            gammaVect(ii) = gammaSim;
            RQ(end) = gammaSim;
            [~,geo,~,out,pathname] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
            nFEA = nFEA+1;
            TVect(ii)  = out.T;
            FdVect(ii) = out.fd;
            FqVect(ii) = out.fq;
            IdVect(ii) = out.id;
            IqVect(ii) = out.iq;
            gVect(ii)  = gammaSim;
        end
        if ~stp
            [~,index] = max(TVect,[],'omitnan');
            OUT.fd  = FdVect(index);
            OUT.fq  = FqVect(index);
            OUT.id  = IdVect(index);
            OUT.iq  = IqVect(index);
            OUT.mPM = calcMassPM(geo,mat);
            OUT.mCu = calcMassCu(geo,mat);
            %gamma0 = gVect(index);
            OUT.T   = TVect(index);
        end
    end
    [~,index] = max(TVect,[],'omitnan');
    OUT.fd  = FdVect(index);
    OUT.fq  = FqVect(index);
    OUT.id  = IdVect(index);
    OUT.iq  = IqVect(index);
    OUT.mPM = calcMassPM(geo,mat);
    OUT.mCu = calcMassCu(geo,mat);
    % gamma0 = gVect(index);
    OUT.T   = TVect(index);

end

per = calc_i0(geo,per,mat);

% Characteristic current (Vtype and SPM)
% if (strcmp(geo.RotType,'SPM')||strcmp(geo.RotType,'Vtype'))
%     % characteristic current
%     RQ(end) = 180;
%     [~,~,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
%     nFEA = nFEA+1;
%     OUT.f0 = out.fd;
%     %     % no-load simulation
%     %     per.overload=0;
%     %     RQ(end)=0;
%     %     [~,~,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
%     %     OUT.fM = out.fd;
% else
%     %     OUT.fM = 0;
OUT.f0 = 0;
% end

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

if flagPM || strcmp(geo.RotType,'EESM')
    
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
            nFEA = nFEA+1;
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
    per.overload     = [1 0];
    setup.RQ         = RQ;
    setup.RQ(end)    = 90;
    setup.flagSave   = 0;
    setup.flagFEAfix = 1;
    setup.filemot    = filemot;
    setup.geo        = geo;
    setup.mat        = mat;
    setup.per        = per;
    setup.idq0       = [OUT.id+1i*OUT.iq 0+j*0]; %#ok<IJCL>
    setup.fdq0       = [OUT.fd+1i*OUT.fq 0-j*OUT.fM]; %#ok<IJCL>
    setup.T0         = [OUT.T 0];
    [pkSCout] = eval_peakShortCircuitCurrent(setup);
    nFEA             = nFEA+pkSCout.count;
    OUT.iHWC         = abs(pkSCout.idq(1));
    OUT.iUGO         = abs(pkSCout.idq(2));
else
    OUT.iHWC = NaN;
    OUT.iUGO = NaN;
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
    nFEA = nFEA+1;
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
    nFEA = nFEA+1;
    OUT.dPM0  = out.dPM;
    OUT.Bmin0 = out.Bmin;
else
    OUT.dPM0  = NaN;
    OUT.Bmin0 = NaN;
end

OUT.dPMUGO  = NaN;
OUT.BminUGO = NaN;


OUT.MaxDef       = NaN;
OUT.agclear      = NaN;
OUT.MaxStress    = NaN;
OUT.TanRibStress = NaN;
OUT.RadRibStress = NaN;
% OUT.MaxStress_prc     = NaN;
OUT.PrcTanStress = NaN;
OUT.PrcRadStress = NaN;


if flagTherm
    mcad=actxserver('MotorCAD.AppAutomation');
    invoke(mcad,'LoadFromFile',[filemot(1:(end-4)) '.mot']);

    %     invoke(mcad,'SetVariable','BPMRotor','0'); %surface radial
    %     invoke(mcad,'SetVariable','Magnet_Thickness',0.1); %to avoid interference

    invoke(mcad,'SetVariable','Stator_Bore',2*(geo.r+geo.g));
    invoke(mcad,'SetVariable','Airgap', geo.g);

    if geo.parallel_slot==0
        invoke(mcad,'SetVariable','Tooth_Width',geo.wt);      %ParallelTooth
    else
        tmp=(geo.r+geo.g+geo.lt/15)*sin(pi/geo.p/geo.Qs)-geo.wt;
        invoke(mcad,'SetVariable','Slot_Width',tmp);       %ParallelSlot
    end

    invoke(mcad,'SetVariable','Slot_Depth', geo.lt);
    if geo.parallel_slot==0
        invoke(mcad,'SetVariable','Slot_Corner_Radius',geo.SFR); %ParallelTooth
    end

    invoke(mcad,'SetVariable','Tooth_Tip_Depth',geo.ttd);

    tmp = (geo.acs)*((geo.r+geo.g)*2*pi/(6*geo.p*geo.q));
    invoke(mcad,'SetVariable','Slot_Opening',tmp);

    invoke(mcad,'SetVariable','Tooth_Tip_Angle', geo.tta);

    if strcmp(geo.win.condType,'Square')

        wCond = geo.st-2*(geo.win.condIns+0.2);
        hCond = geo.win.kcu*geo.Aslot/(geo.win.nCond*wCond);

        invoke(mcad,'SetVariable','Copper_Width',wCond);
        invoke(mcad,'SetVariable','Copper_Height',hCond);
    end

    invoke(mcad,'SetVariable','ShaftSpeed',1);
    invoke(mcad,'SetVariable','Number_Transient_Points',10);
    invoke(mcad,'SetVariable','Transient_Time_Period',30);                            %transient time
    invoke(mcad,'SetVariable','ThermalCalcType', 1);

    invoke(mcad,'SetVariable','Armature_Copper_Loss_@Ref_Speed',per.Loss);            %insert 0rpm loss
    invoke(mcad,'SetVariable','StatorWindingTemperatureAtWhichPcuInput',per.tempcu);  %Cu loss specified @ tempCu

    disp('Thermal Transient Analysis in progress...')
    invoke(mcad,'DoTransientAnalysis');

    % code to iterate until the temperature limit is reached
    %     invoke(mcad,'SetVariable','ThermalCalcType', 0);
    %     TLim_cu = 180;
    %     ii   = 0;
    %     flag = 1;
    %DC_loss = per.Loss;
    %     while flag==1
    %         ii = ii + 1;
    %         flag = 0;
    %     invoke(mcad,'SetVariable','Armature_Copper_Loss_@Ref_Speed', DC_loss);
    %     invoke(mcad,'DoSteadyStateAnalysis');
    %     [~, WindingTemperature_Max] = invoke(mcad,'GetVariable','T_[Winding_Max]');
    %         [~, WindingTemperature_Ave] = invoke(mcad,'GetVariable','T_[Winding_Average]');

    %         if (WindingTemperature_Max<TLim_cu-5 || WindingTemperature_Max>TLim_cu+1)
    %             WindingTemperature_Max(WindingTemperature_Max>2*TLim_cu) = TLim_cu(WindingTemperature_Max>2*TLim_cu)*1.5;
    %             DC_loss = (DC_loss + (TLim_cu-WindingTemperature_Max)/TLim_cu * DC_loss);
    %             flag    = 1;
    %             %             I_th0   = sqrt(DC_loss/3/R*2);
    %             %             disp(['Iteration ' num2str(ii) ': ' num2str(round(I_th0)) ' A - ' num2str(round(WindingTemperature_Max)) '°C copper temperature'])
    %         end
    %     end

    %     tempCu_max = WindingTemperature_Max;
    [~, tempCu_max]    = invoke(mcad,'GetVariable','T_[Winding_Max]');
    [~, tempCuAct_max] = invoke(mcad,'GetVariable','T_[Winding_Active_Max]');
    OUT.outTherm.maxTcu    = tempCu_max;
    OUT.outTherm.maxTcuAct = tempCuAct_max;
    %OUT.outTherm.Loss = DC_loss;
    invoke(mcad,'Quit');
else
    OUT.outTherm.maxTcu    = NaN;
    OUT.outTherm.maxTcuAct = NaN;
end

OUT.nFEA = nFEA;