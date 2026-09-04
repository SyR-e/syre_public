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
per.overload = 1;

tic

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
    % golden ratio search
    switch geo.RotType
        case 'EESM'
            maxIter = 8;
            minGamma = 55;
            maxGamma = 125;
            % max gamma error = 3.16°
        case {'SPM','SPM-Halbach'}
            maxIter = 4;
            minGamma = 90;
            maxGamma = 110;
            % max gamma error = 2.36° (typically not needed)
        otherwise
            if strcmp(mat.LayerMag.MatName,'Air')
                maxIter = 6;
                minGamma = 40;
                maxGamma = 75;
                % max gamma error = 1.56°
            else
                maxIter = 6;
                minGamma = 30;
                maxGamma = 80;
                % max gamma error = 2.25°
            end
    end

    phi  = (1 + sqrt(5))/2;
    invp = 1/phi;
    x = (-1+sqrt(5))/2;

    gammaErrMax = x^(maxIter-3)*(1-x)/2*abs(maxGamma-minGamma);


    %Initialization
    TVect  = nan(1,maxIter);
    FdVect = nan(1,maxIter);
    FqVect = nan(1,maxIter);
    IdVect = nan(1,maxIter);
    IqVect = nan(1,maxIter);
    gVect  = nan(1,maxIter);

    % simulations - golden ratio search
    a = minGamma;
    b = maxGamma;
    c = b - (b-a)*invp;
    d = a + (b-a)*invp;
    for ii=1:maxIter-4
        if ii==1
            gammaSim = a;
            RQ(end) = gammaSim;
            [~,geo,~,out,pathname] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
            nFEA = nFEA+1;
            Ta = out.T;
            TVect(nFEA)  = out.T;
            FdVect(nFEA) = out.fd;
            FqVect(nFEA) = out.fq;
            IdVect(nFEA) = out.id;
            IqVect(nFEA) = out.iq;
            gVect(nFEA)  = gammaSim;
            gammaSim = b;
            RQ(end) = gammaSim;
            [~,geo,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
            nFEA = nFEA+1;
            Tb = out.T;
            TVect(nFEA)  = out.T;
            FdVect(nFEA) = out.fd;
            FqVect(nFEA) = out.fq;
            IdVect(nFEA) = out.id;
            IqVect(nFEA) = out.iq;
            gVect(nFEA)  = gammaSim;
            c = b - (b-a)*invp;
            gammaSim = c;
            RQ(end) = gammaSim;
            [~,geo,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
            nFEA = nFEA+1;
            Tc = out.T;
            TVect(nFEA)  = out.T;
            FdVect(nFEA) = out.fd;
            FqVect(nFEA) = out.fq;
            IdVect(nFEA) = out.id;
            IqVect(nFEA) = out.iq;
            gVect(nFEA)  = gammaSim;
            d = a + (b-a)*invp;
            gammaSim = d;
            RQ(end) = gammaSim;
            [~,geo,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
            nFEA = nFEA+1;
            Td = out.T;
            TVect(nFEA)  = out.T;
            FdVect(nFEA) = out.fd;
            FqVect(nFEA) = out.fq;
            IdVect(nFEA) = out.id;
            IqVect(nFEA) = out.iq;
            gVect(nFEA)  = gammaSim;
        end
        if Tc>Td
            b = d;
            Tb = Td;
            d = c;
            Td = Tc;
            c = b - (b-a)*invp;
            gammaSim = c;
            RQ(end) = gammaSim;
            [~,geo,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
            nFEA = nFEA+1;
            Tc = out.T;
            TVect(nFEA) = out.T;
            FdVect(nFEA) = out.fd;
            FqVect(nFEA) = out.fq;
            IdVect(nFEA) = out.id;
            IqVect(nFEA) = out.iq;
            gVect(nFEA)  = gammaSim;
        else
            a = c;
            Ta = Tc;
            c = d;
            Tc = Td;
            d = a + (b-a)*invp;
            gammaSim = d;
            RQ(end) = gammaSim;
            [~,geo,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
            nFEA = nFEA+1;
            Td = out.T;
            TVect(nFEA) = out.T;
            FdVect(nFEA) = out.fd;
            FqVect(nFEA) = out.fq;
            IdVect(nFEA) = out.id;
            IqVect(nFEA) = out.iq;
            gVect(nFEA)  = gammaSim;
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

timeMTPA = toc;

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
    % no-load simulation
    per.overload=0;
    RQ(end)=0;
    [~,~,~,out,~] = FEMMfitness(RQ,geo,per,mat,eval_type,filemot);
    nFEA = nFEA+1;
    %     OUT.fM = -out.fq;
    OUT.fM = abs(out.fd+j*out.fq);
else
    OUT.fM = 0;
end

timeNoLoad = toc;

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

timeIch = toc;

if flagSC
    per.overload          = [1 0];
    setup.RQ              = RQ;
    setup.RQ(end)         = 90;
    setup.flagSave        = 0;
    setup.flagFEAfix      = 1;
    setup.filemot         = filemot;
    setup.geo             = geo;
    setup.mat             = mat;
    setup.per             = per;
    setup.idq0            = [OUT.id+1i*OUT.iq 0+j*0];
    setup.fdq0            = [OUT.fd+1i*OUT.fq 0-j*OUT.fM];
    setup.T0              = [OUT.T 0];
    % setup.XFEMMsimulation = geo.XFEMMsimulation;
    [pkSCout] = eval_peakShortCircuitCurrent(setup);
    nFEA             = nFEA+pkSCout.count;
    OUT.iHWC         = abs(pkSCout.idq(1));
    OUT.iUGO         = abs(pkSCout.idq(2));
else
    OUT.iHWC = NaN;
    OUT.iUGO = NaN;
end

timeHWC = toc;

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

timeDemagHWC = toc;

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

timeDemag0 = toc;

if flagDemagUGO
    per.overload = OUT.iUGO/i0;
    if (strcmp(geo.RotType,'SPM')||strcmp(geo.RotType,'Vtype'))
        per.gamma = 180;
    else
        per.gamma = 90;
    end
    RQ(end) = per.gamma;
    [~,~,~,out,~] = FEMMfitness(RQ,geo,per,mat,'demagArea',filemot);
    nFEA = nFEA+1;
    OUT.dPMUGO  = out.dPM;
    OUT.BminUGO = out.Bmin;
else
    OUT.dPMUGO  = NaN;
    OUT.BminUGO = NaN;
end

timeDemagUGO = toc;

if flagMech
    [geo,~,mat] = interpretRQ(RQ,geo,mat);

    if ~strcmp(geo.RotType,'EESM')
        % Set up info for solving structural model
        simSetup.evalSpeed  = geo.nmax;
        simSetup.meshSize   = 'PDE fine';
        simSetup.meshShaft  = 0;
        simSetup.flagFull   = 0;
        simSetup.shaftBC    = 1;
        simSetup.flagSoftPM = 0;
        simSetup.flagResin  = 0;
        geo.custom = 0;
        warning('off')
        [~, filename, ~] = fileparts(filemot);

        % Geometry/points generation
        [geo,mat] = draw_motor_in_FEMM(geo,mat,pathname,filename);
        simSetup.pathname  = pathname;
        simSetup.filename  = [filename '.mat'];

        % Convertion FEMM - PDE
        [out.structModel] = femm2pde(geo,mat,simSetup);
        %Solving the PDE model
        [out.sVonMises,R,out.structModel] = calcVonMisesStress(out.structModel);
        % Deformation, clerance and stress calculations
        [out.Stress] = eval_maxStress(out.structModel,out.sVonMises,R,geo,mat);

        if (strcmp(geo.RotType,'Circular')||strcmp(geo.RotType,'Seg'))
            % Deformation magnitude
            OUT.MaxDef       = out.Stress.MaxDef;
            % Air-gap clearance
            OUT.agclear      = out.Stress.agclear;
            % Maximum stress values
            OUT.MaxStress    = out.Stress.MaxStress/1e6;
            OUT.TanRibStress = (max(out.Stress.sigmaTanMax))/1e6;
            OUT.RadRibStress = (max(out.Stress.sigmaRadMax))/1e6;
            % Percentile stress values
            %OUT.MaxStress_prc     = out.Stress.sigmaTotPrc/1e6;
            OUT.PrcTanStress = out.Stress.Tan_sigmaTotPrc/1e6;
            OUT.PrcRadStress = out.Stress.Rad_sigmaTotPrc/1e6;
        else
            OUT.MaxDef            = out.Stress.MaxDef;
            OUT.agclear           = out.Stress.agclear;
            OUT.MaxStress         = out.Stress.MaxStress/1e6;
            % OUT.MaxStress_prc     = out.Stress.sigmaTotPrc/1e6;
        end
    else    
        dataSet.EvalSpeed  = geo.nmax;
        dataSet.Mesh       = 4;
        nameIn             = 'Structural';

        [model,geo] = draw_motor_in_COMSOL(geo,mat,pathname,nameIn);
        dataSet.infoComsol = geo.infoComsol;
        [dataSet,model] = eval_vonMisesStress_COMSOL(dataSet,model);
        OUT.MaxDef       = dataSet.infoComsol.structural.maxDisplacementMM;
        OUT.agclear      = NaN;
        OUT.MaxStress    = dataSet.infoComsol.structural.maxElementVonMisesMPa;
    end




    nFEA = nFEA+1;
else
    OUT.MaxDef       = NaN;
    OUT.agclear      = NaN;
    OUT.MaxStress    = NaN;
    OUT.TanRibStress = NaN;
    OUT.RadRibStress = NaN;
    % OUT.MaxStress_prc     = NaN;
    OUT.PrcTanStress = NaN;
    OUT.PrcRadStress = NaN;
end

timeStruct = toc;


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