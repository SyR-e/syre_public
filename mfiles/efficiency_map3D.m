% Copyright 2025
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


% Special thanks goes to Sandro Rubino, who made the biggest contribution to this script.

function [out,outD] = efficiency_map3D(motorModel,TwData)

path_i   = motorModel.data.pathname;
filename = motorModel.dataSet.currentfilename;

%map_type = 2; % Map Type (1 -> Motoring, 2 -> Motoring/Generator)

if ~exist([path_i filename(1:end-4) '_results\MMM results\TwMap_\IsoTorque'], 'dir')
    mkdir([path_i filename(1:end-4) '_results\MMM results\TwMap_\IsoTorque']);
end
outputFolder = [path_i filename(1:end-4) '_results\MMM results\TwMap_\IsoTorque'];

Vs_lim     = motorModel.data.Vdc/sqrt(3);       % Peak Stator Voltage Limit, (Vpk)
Is_lim     = motorModel.data.i0;                % Peak Stator Current Limit, (Apk)
If_lim     = motorModel.data.if0;  % Rotor Current Limit, (A)
nm_min     = TwData.nmin;
nm_lim     = TwData.nmax;                             % Speed Limit, (rpm)
Tm_min     = TwData.Tmin;
Tm_lim     = TwData.Tmax;                               % Torque Limit, (A)

% Torque-speed map setting
nm_pts = TwData.nstep;              % Speed Step
Tm_pts = TwData.Tstep;              % Torque Step  

% % Temperature Settings
% Ts_0    = 23.3;                             % Reference Stator Temperature (°C)
% Tr_0    = 22.1;                             % Reference Field Temperature (°C)
% Ts_ref  = 140;                              % Operative Stator Temperature (°C)
% Tr_ref  = 140;                              % Operative Field Temperature (°C)
% Rs      = motorModel.data.Rs;%17.55e-3; 
% %Rs      = 8.19*2;
% %Rr     = 15.51; %@20°C
% Rr      = 148/142*2.75/1.68*1.1436*8; % @22,1 °C (mockup con 142 spire in Cu -> avv. finale con 148 spire in Al)
% % Starting Maps Resolution (minumum 128, 256 is the best compromise)
% Np_res  = 256;    % Maps Resolution (128 (~2 mins) - 256 (~5 mins)  - 512 (~45 mins))

% Temperature Settings
Ts_0    = motorModel.data.tempCu;     % Reference Stator Temperature (°C)
Tr_0    = motorModel.data.tempCoRo;   % Reference Field Temperature (°C)
Ts_ref  = TwData.temperature;         % Operative Stator Temperature (°C)
Tr_ref  = TwData.temperature;         % Operative Field Temperature (°C)

alpha_s = motorModel.mat.SlotCond.alpha;
alpha_r = motorModel.mat.BarCond.alpha;

[Rs,~] = calcRsTempFreq(motorModel.data.Rs, Ts_0, alpha_s, motorModel.data.l, motorModel.data.lend , [], '0', Ts_ref, 0);
[Rr,~] = calcRsTempFreq(motorModel.data.Rf, Tr_0, alpha_r, motorModel.data.l, motorModel.data.lendf, [], '0', Tr_ref, 0);
  
% Starting Maps Resolution (minumum 128, 256 is the best compromise)
Np_res  = 256;      % Maps Resolution (128 (~2 mins) - 256 (~5 mins)  - 512 (~45 mins))

% Optimization Case

TnSetup = motorModel.TnSetup;
IronLossFlag = TnSetup.IronLossFlag;

tit = TnSetup.Control;

if     strcmp(tit,'Min Overall Losses')
    opt_case = 1;

elseif strcmp(tit,'Min Stator Losses')
    opt_case = 2;

elseif strcmp(tit,'Min Rotor Losses')
    opt_case = 3;
    
elseif strcmp(tit,'Min Overall Joule Losses')
    opt_case = 4;

elseif strcmp(tit,'Min Stator Joule Losses')
    opt_case = 5;

elseif strcmp(tit,'Min Rotor Joule Losses')
    opt_case = 6;

elseif strcmp(tit,'Max Power Factor')
    opt_case = 7;
    
end

%--------------------------------------------------------------------------
% Magnetic Model Loading & Map Resolution Extension
%--------------------------------------------------------------------------
pp = motorModel.geo.p;
Fd = motorModel.FluxMap_dq.Fd;
Fq = motorModel.FluxMap_dq.Fq;
Fr = motorModel.FluxMap_dq.Fr;
Ir = motorModel.FluxMap_dq.Ir;
Id = motorModel.FluxMap_dq.Id;
Iq = motorModel.FluxMap_dq.Iq;

Id_vct_ref = unique(Id);
Iq_vct_ref = unique(Iq);
Ir_vct_ref = unique(Ir);

Id_new = linspace(min(Id_vct_ref),max(Id_vct_ref),Np_res/2);
Iq_new = linspace(min(Iq_vct_ref),max(Iq_vct_ref),Np_res);
Ir_new = linspace(min(Ir_vct_ref),max(Ir_vct_ref),Np_res/3);

[Iq,Id,Ir]  = ndgrid(Iq_new, Id_new, Ir_new);
Is          = hypot(Id,Iq);
As          = atan2d(Iq,Id);

FD_GI     = griddedInterpolant({Iq_vct_ref,Id_vct_ref,Ir_vct_ref},Fd,"spline");
FQ_GI     = griddedInterpolant({Iq_vct_ref,Id_vct_ref,Ir_vct_ref},Fq,"spline");
FF_GI     = griddedInterpolant({Iq_vct_ref,Id_vct_ref,Ir_vct_ref},Fr,"spline");

Fd = FD_GI({Iq_new,Id_new,Ir_new});
Fq = FQ_GI({Iq_new,Id_new,Ir_new});
Fr = FF_GI({Iq_new,Id_new,Ir_new});

Te = 1.5 * pp * (Fd .* Iq - Fq .* Id);

%--------------------------------------------------------------------------
% Loss Model Loading & Map Resolution Extension
%--------------------------------------------------------------------------

    f0 = 0;
    ksh = motorModel.mat.Stator.alpha;
    krh = motorModel.mat.Rotor.alpha;
    kse = 2.0;
    kre = 2.0;

if strcmp(IronLossFlag,'Yes')
    f0  = motorModel.IronPMLossMap_dq.f0;
    kil = TnSetup.IronLossFactor;

    Psh = kil*motorModel.IronPMLossMap_dq.Pfes_h;
    Pse = kil*motorModel.IronPMLossMap_dq.Pfes_c;
    Prh = kil*motorModel.IronPMLossMap_dq.Pfer_h;
    Pre = kil*motorModel.IronPMLossMap_dq.Pfer_c;

    Id_vct_ref_losses = unique(motorModel.IronPMLossMap_dq.Id);
    Iq_vct_ref_losses = unique(motorModel.IronPMLossMap_dq.Iq);
    If_vct_ref_losses = unique(motorModel.IronPMLossMap_dq.Ir);

    PSH_GI = griddedInterpolant({Iq_vct_ref_losses,Id_vct_ref_losses,If_vct_ref_losses},Psh,"spline");
    PSE_GI = griddedInterpolant({Iq_vct_ref_losses,Id_vct_ref_losses,If_vct_ref_losses},Pse,"spline");
    PRH_GI = griddedInterpolant({Iq_vct_ref_losses,Id_vct_ref_losses,If_vct_ref_losses},Prh,"spline");
    PRE_GI = griddedInterpolant({Iq_vct_ref_losses,Id_vct_ref_losses,If_vct_ref_losses},Pre,"spline");
end

%--------------------------------------------------------------------------
% Reference Torque-Speed Point
%--------------------------------------------------------------------------

Tm_vct = linspace(Tm_min,Tm_lim,Tm_pts);
nm_vct = linspace(nm_min,nm_lim,nm_pts);
[nm_map,Tm_map] = meshgrid(nm_vct,Tm_vct); 
[~,Tm_map_idx]  = meshgrid(nm_vct,1:1:length(Tm_vct));
N_pnts = numel(nm_map);

%--------------------------------------------------------------------------
% Iso Torque
%--------------------------------------------------------------------------
disp('Computation of Iso-Torque...')
Tm_vct = linspace(0,Tm_lim,Tm_pts);
for i = 1 : length(Tm_vct)

    % Compute IsoTorque
    [~,T_iso] = isosurface(Id,Iq,Ir,Te,Tm_vct(1,i));

    % Load IsoTorque Currents
    IsoT.Id_vct = T_iso(:,1);
    IsoT.Iq_vct = T_iso(:,2);
    IsoT.If_vct = T_iso(:,3);

    % Load IsoTorque Fluxes
    IsoT.Is_vct = hypot(IsoT.Id_vct,IsoT.Iq_vct);
    IsoT.As_vct = atan2d(IsoT.Iq_vct,IsoT.Id_vct);
    IsoT.Fd_vct = FD_GI(IsoT.Iq_vct,IsoT.Id_vct,IsoT.If_vct);
    IsoT.Fq_vct = FQ_GI(IsoT.Iq_vct,IsoT.Id_vct,IsoT.If_vct);
    IsoT.Ff_vct = FF_GI(IsoT.Iq_vct,IsoT.Id_vct,IsoT.If_vct);

    % Load IsoTorque Torque
    IsoT.Te_vct = 1.5 * pp * (IsoT.Fd_vct .* IsoT.Iq_vct - IsoT.Fq_vct .* IsoT.Id_vct);

    % Load IsoTorque Ref Losses
    if strcmp(IronLossFlag,'Yes')
        IsoT.Psh_vct = PSH_GI(IsoT.Iq_vct,IsoT.Id_vct,IsoT.If_vct);
        IsoT.Pse_vct = PSE_GI(IsoT.Iq_vct,IsoT.Id_vct,IsoT.If_vct);
        IsoT.Prh_vct = PRH_GI(IsoT.Iq_vct,IsoT.Id_vct,IsoT.If_vct);
        IsoT.Pre_vct = PRE_GI(IsoT.Iq_vct,IsoT.Id_vct,IsoT.If_vct);
    end

    % Save Results
    torquefilename = ['IsoT_' num2str(round(Tm_vct(i))) '_Nm.mat'];
    save([outputFolder '\' torquefilename],'IsoT','-v7.3')

    % Print Progress Status
    %fprintf(['Computation of IsoT ' num2str(i) ' of ' num2str(length(Tm_vct)) ' - Done \n'])

end

%--------------------------------------------------------------------------
% Output Matrices
%--------------------------------------------------------------------------
disp('Efficiency Mapping...')
% dqf Currents
Id_map      = zeros(size(nm_map));
Iq_map      = zeros(size(nm_map));
Ir_map      = zeros(size(nm_map));
% magnetic dqf Currents
Idm_map     = zeros(size(nm_map));
Iqm_map     = zeros(size(nm_map));
% dqf Fluxes
Fd_map      = zeros(size(nm_map));
Fq_map      = zeros(size(nm_map));
Fr_map      = zeros(size(nm_map));
% dqf Voltages
Vd_map      = zeros(size(nm_map));
Vq_map      = zeros(size(nm_map));
Vr_map      = zeros(size(nm_map));
% Losses
Plt_map     = zeros(size(nm_map));
Pjs_map     = zeros(size(nm_map));
Pjr_map     = zeros(size(nm_map));
Pfes_map    = zeros(size(nm_map));
Pfer_map    = zeros(size(nm_map));
Psh_map     = zeros(size(nm_map));
Pse_map     = zeros(size(nm_map));
Prh_map     = zeros(size(nm_map));
Pre_map     = zeros(size(nm_map));

% Ghost Matrix
ghs_map = zeros(size(nm_map));

% Define Temp Vectors
% dqf Currents
Id_vct  = zeros(1,size(nm_map,2));
Iq_vct  = zeros(1,size(nm_map,2));
If_vct  = zeros(1,size(nm_map,2));
% magnetic dqf Currents
Idm_vct = zeros(1,size(nm_map,2));
Iqm_vct = zeros(1,size(nm_map,2));
% dqf Fluxes
Fd_vct  = zeros(1,size(nm_map,2));
Fq_vct  = zeros(1,size(nm_map,2));
Ff_vct  = zeros(1,size(nm_map,2));
% dqf Voltages
Vd_vct  = zeros(1,size(nm_map,2));
Vq_vct  = zeros(1,size(nm_map,2));
Vf_vct  = zeros(1,size(nm_map,2));
% Losses
Plt_vct = zeros(1,size(nm_map,2));
Pjs_vct = zeros(1,size(nm_map,2));
Pjf_vct = zeros(1,size(nm_map,2));
Pis_vct = zeros(1,size(nm_map,2));
Pif_vct = zeros(1,size(nm_map,2));
Psh_vct = zeros(1,size(nm_map,2));
Pse_vct = zeros(1,size(nm_map,2));
Prh_vct = zeros(1,size(nm_map,2));
Pre_vct = zeros(1,size(nm_map,2));
% Ghost
ghs_vct = zeros(1,size(nm_map,2));

% Efficency Mapping

% EM Computation
%addpath iso_torque\
for j = 1 : size(Tm_map,1)

    % Select Torque Point
    Tm_ref = Tm_map(j,1);

    % Load IsoTorque Data
    torquefilename = [outputFolder '\IsoT_' num2str(round(abs(Tm_map(j,1)))) '_Nm.mat'];

    out      = load(torquefilename);
    
    % Set Speeed
    parfor i = 1 : size(Tm_map,2)
        
        idx_min = 0;
        % Select Speed Point
        nm_ref = nm_vct(i);

        % Read IsoTorque
        % IsoT  = out(i).IsoT;
        IsoT = out.IsoT;

        % Negative Torque Correction
        if (Tm_ref < 0)
            IsoT.Iq_vct = - IsoT.Iq_vct;
            IsoT.Fq_vct = - IsoT.Fq_vct;
            IsoT.Te_vct = - IsoT.Te_vct;
        end

        % Compute Frequency
        fe_ref = (pp * nm_ref) / 60;

        % Evaluate Iron Losses
        if strcmp(IronLossFlag,'Yes')
            Psh_iso = IsoT.Psh_vct * (abs(fe_ref / f0) ^ ksh);
            Pse_iso = IsoT.Pse_vct * (abs(fe_ref / f0) ^ kse);
            Prh_iso = IsoT.Prh_vct * (abs(fe_ref / f0) ^ krh);
            Pre_iso = IsoT.Pre_vct * (abs(fe_ref / f0) ^ kre);


            % Overall Stator and Rotor Iron Losses
            Pis_iso = Psh_iso + Pse_iso;
            Pir_iso = Prh_iso + Pre_iso;
        else
            Pis_iso = 0;
            Pir_iso = 0;
        end

        % Compute Synchronous Speed (rad/s)
        we_iso = 2 * pi * fe_ref;

        % Compute Back-Emf
        Ed_iso   = - we_iso .* IsoT.Fq_vct;
        Eq_iso   = + we_iso .* IsoT.Fd_vct;
        Es_iso   = hypot(Ed_iso,Eq_iso);

        % Compute Back-Emf Powers
        Ied_iso = IsoT.Id_vct;
        Ieq_iso = IsoT.Iq_vct;
        Pe_iso  = (3 / 2) * (Ied_iso .* Ed_iso + Ieq_iso .* Eq_iso) + (Pis_iso + Pir_iso);
        Qe_iso  = (3 / 2) * (Ied_iso .* Eq_iso - Ieq_iso .* Ed_iso);

        % Compute Stator dq Currents
        if (nm_ref == 0)
            Isd_iso  = Ied_iso;
            Isq_iso  = Ieq_iso;
        else
            Isd_iso  = (2 / 3) * (1 ./ (Es_iso .* Es_iso)) .* (Ed_iso .* Pe_iso + Eq_iso .* Qe_iso);
            Isq_iso  = (2 / 3) * (1 ./ (Es_iso .* Es_iso)) .* (Eq_iso .* Pe_iso - Ed_iso .* Qe_iso);
        end
        Is_iso   = hypot(Isd_iso,Isq_iso);

        % Compute Stator Joule Losses
        Pjs_iso = 1.5 * Rs * (Is_iso .* Is_iso);

        % Compute Field Joule Losses
        Pjf_iso = Rr * (IsoT.If_vct .* IsoT.If_vct);

        % Compute Overall Stator Powers
        Ps_iso  = Pe_iso + Pjs_iso;
        Qs_iso  = Qe_iso;

        % Compute Stator dq Voltages
        Vsd_iso = (2 / 3) * (1 ./ (Is_iso .* Is_iso)) .* (Isd_iso .* Ps_iso - Isq_iso .* Qs_iso);
        Vsq_iso = (2 / 3) * (1 ./ (Is_iso .* Is_iso)) .* (Isq_iso .* Ps_iso + Isd_iso .* Qs_iso);

        % Apply Low Stator Current Range Correction
        is_th = 0.1;
        Vsd_iso(Is_iso < is_th) = Ed_iso(Is_iso < is_th);
        Vsq_iso(Is_iso < is_th) = Eq_iso(Is_iso < is_th);

        % Stator Voltage Amplitude
        Vs_iso  = hypot(Vsd_iso,Vsq_iso);

        % Field Voltage
        Vf_iso = Rr * IsoT.If_vct;
        
        

        % Overall Losses
        Plt_iso = Pjs_iso + Pjf_iso + Pis_iso + Pir_iso;

        % Power Factor
        pf_iso = cos(atan2(Vsq_iso,Vsd_iso) - atan2(Isq_iso,Isd_iso));

        % Compute Acceptable Points
        range = find((Vs_iso <= Vs_lim) & (Is_iso <= Is_lim) & (abs(IsoT.If_vct) <= If_lim));
        if (isempty(range))
            % NaN on Matrix
            % dqf Currents
            Id_vct(i)  = NaN;
            Iq_vct(i)  = NaN;
            If_vct(i)  = NaN;
            % magnetic dqf Currents
            Idm_vct(i) = NaN;
            Iqm_vct(i) = NaN;
            % dqf Fluxes
            Fd_vct(i)  = NaN;
            Fq_vct(i)  = NaN;
            Ff_vct(i)  = NaN;
            % dqf Voltages
            Vd_vct(i)  = NaN;
            Vq_vct(i)  = NaN;
            Vf_vct(i)  = NaN;
            % Losses
            Plt_vct(i) = NaN;
            Pjs_vct(i) = NaN;
            Pjf_vct(i) = NaN;
            Pis_vct(i) = NaN;
            Pif_vct(i) = NaN;
            Psh_vct(i) = NaN;
            Pse_vct(i) = NaN;
            Prh_vct(i) = NaN;
            Pre_vct(i) = NaN;
            % Ghost
            ghs_vct(i) = 1;
        else 
            % Optimization Case
            % Case 1 = Min Overall Losses
            % Case 2 = Min Stator Losses
            % Case 3 = Min Rotor Losses
            % Case 5 = Max Power Factor
            switch opt_case
                case 1
                    % Minimum Overall Losses Point
                    [~,idx]   = min(Plt_iso(range));
                    % Select Optimal Global Index
                    idx_min = range(idx);
                case 2
                    % Minimum Stator Losses Point
                    if strcmp(IronLossFlag,'Yes')
                        [~,idx]   = min(Pjs_iso(range) + Pis_iso(range));
                    else
                        [~,idx]   = min(Pjs_iso(range));
                    end
                    % Select Optimal Global Index
                    idx_min = range(idx);
                case 3
                    % Minimum Rotor Losses Point
                    if strcmp(IronLossFlag,'Yes')
                        [~,idx]   = min(Pjf_iso(range) + Pir_iso(range));
                    else
                        [~,idx]   = min(Pjf_iso(range));
                    end
                    
                    % Select Optimal Global Index
                    idx_min = range(idx);
                case 4
                    % Minimum Joule Losses Point
                    [~,idx]   = min(Pjs_iso(range) + Pjf_iso(range));
                    % Select Optimal Global Index
                    idx_min = range(idx);
                case 5
                    % Minimum Stator Joule Losses Point
                    [~,idx]   = min(Pjs_iso(range));
                    % Select Optimal Global Index
                    idx_min = range(idx);
                case 6
                    % Minimum Rotor Joule Losses Point
                    [~,idx]   = min(Pjf_iso(range));
                    % Select Optimal Global Index
                    idx_min = range(idx);
                case 7
                    % Maximum Power Factor Point
                    tmp_pf_iso     = round(abs(pf_iso(range)),4);
                    tmp_pl_iso     = Plt_iso(range);
                    range_pf       = find(tmp_pf_iso == max(tmp_pf_iso));
                    [~,idx_min_pl] = min(tmp_pl_iso(range_pf));
                    idx            = range_pf(idx_min_pl);
                    % Select Optimal Global Index
                    idx_min = range(idx); 
            end  
            % Fill Matrices
            % dqf Currents
            Id_vct(i)  = Isd_iso(idx_min);
            Iq_vct(i)  = Isq_iso(idx_min);
            If_vct(i)  = IsoT.If_vct(idx_min);
            % magnetic dqf Currents
            Idm_vct(i) = Ied_iso(idx_min);
            Iqm_vct(i) = Ieq_iso(idx_min);
            % dqf Fluxes
            Fd_vct(i)  = IsoT.Fd_vct(idx_min);
            Fq_vct(i)  = IsoT.Fq_vct(idx_min);
            Ff_vct(i)  = IsoT.Ff_vct(idx_min);
            % dqf Voltages
            Vd_vct(i)  = Vsd_iso(idx_min);
            Vq_vct(i)  = Vsq_iso(idx_min);
            Vf_vct(i)  = Vf_iso(idx_min);
            % Losses
            Plt_vct(i) = Plt_iso(idx_min);
            Pjs_vct(i) = Pjs_iso(idx_min);
            Pjf_vct(i) = Pjf_iso(idx_min);
            if strcmp(IronLossFlag,'Yes')
                Pis_vct(i) = Pis_iso(idx_min);
                Pif_vct(i) = Pir_iso(idx_min);
                Psh_vct(i) = Psh_iso(idx_min);
                Pse_vct(i) = Pse_iso(idx_min);
                Prh_vct(i) = Prh_iso(idx_min);
                Pre_vct(i) = Pre_iso(idx_min);
            end
            % Ghost
            ghs_vct(i) = 0;
        end

        % Plot Executed Point
        %fprintf(['Speed = ' num2str(nm_ref) ' rpm, Torque = ' num2str(Tm_ref) ' Nm - Done!\n'])

    end

    % Assign Vectors to Maps
    % dqf Currents
    Id_map(j,:) = Id_vct;
    Iq_map(j,:) = Iq_vct;
    Ir_map(j,:) = If_vct;
    % magnetic dqf Currents
    Idm_map(j,:) = Idm_vct;
    Iqm_map(j,:) = Iqm_vct;
    % dqf Fluxes
    Fd_map(j,:) = Fd_vct;
    Fq_map(j,:) = Fq_vct;
    Fr_map(j,:) = Ff_vct;
    % dqf Voltages
    Vd_map(j,:) = Vd_vct;
    Vq_map(j,:) = Vq_vct;
    Vr_map(j,:) = Vf_vct;
    % Losses
    Plt_map(j,:)    = Plt_vct;
    Pjs_map(j,:)    = Pjs_vct;
    Pjr_map(j,:)    = Pjf_vct;
    Pfes_map(j,:)   = Pis_vct;
    Pfer_map(j,:)   = Pif_vct;
    Psh_map(j,:)    = Psh_vct;
    Pse_map(j,:)    = Pse_vct;
    Prh_map(j,:)    = Prh_vct;
    Pre_map(j,:)    = Pre_vct;
    % Ghost Matrix
    ghs_map(j,:) = ghs_vct;
end

% Compute Powers
Pm_map = Tm_map .* ((pi / 30) * nm_map);
Pm_map(ghs_map > 0.5) = NaN;
Pe_map = 1.5 * (Id_map .* Vd_map + Iq_map .* Vq_map);
Qe_map = 1.5 * (Id_map .* Vq_map - Iq_map .* Vd_map);
pf_map = cosd(atan2d(Qe_map,Pe_map));

% Efficency Map
Eff_map              = Pm_map ./ (Pe_map + Pjr_map);
Eff_map(Tm_map < 0)  = abs(Pe_map(Tm_map < 0)) ./ (abs(Pm_map(Tm_map < 0)) + abs(Pjr_map(Tm_map < 0)));
Eff_map(Eff_map > 1) = 1;
Eff_map(Eff_map < 0) = 0;
Eff_map(Tm_map == 0) = 0;
Eff_map(nm_map == 0) = 0;

outD.n      = nm_map;
outD.T      = Tm_map;
outD.Id     = Id_map;
outD.Iq     = Iq_map;
outD.Idm    = Idm_map;
outD.Iqm    = Iqm_map;
outD.Ir     = Ir_map;
outD.Io     = hypot(Id_map,Iq_map);
outD.Vd     = Vd_map;
outD.Vq     = Vq_map;
outD.Vr     = Vr_map;
outD.Vo     = hypot(Vd_map,Vq_map);
outD.Fd     = Fd_map;
outD.Fq     = Fq_map;
outD.Fr     = Fr_map;
outD.P      = Pm_map;
outD.Pe     = Pe_map;
outD.Qe     = Qe_map;
outD.eff    = Eff_map;
outD.PF     = abs(pf_map);
outD.Ploss  = Plt_map;
outD.Pfes   = Pfes_map;
outD.Pfer   = Pfer_map;
outD.Pjs    = Pjs_map;
outD.Pjr    = Pjr_map;
if strcmp(IronLossFlag,'Yes')
    outD.Psh    = Psh_map;
    outD.Pse    = Pse_map;
    outD.Prh    = Prh_map;
    outD.Pre    = Pre_map;
end
outD.Teta_s = Ts_ref;
outD.Teta_r = Tr_ref;


Fnan      = ~isnan(Fd_map);
outD.Tout = Tm_map.*Fnan;
outD.Pfe  = Pfes_map+Pfer_map;

% % Output data
% emName = 'EffyMap';
% motorModel.(emName).strategy    = opt_case; 
% motorModel.(emName).nm          = nm_map;
% motorModel.(emName).Tm          = Tm_map;
% motorModel.(emName).Id          = Id_map;
% motorModel.(emName).Iq          = Iq_map;
% motorModel.(emName).Idm         = Idm_map;
% motorModel.(emName).Iqm         = Iqm_map;
% motorModel.(emName).Ir          = Ir_map;
% motorModel.(emName).Is          = hypot(Id_map,Iq_map);
% motorModel.(emName).Vd          = Vd_map;
% motorModel.(emName).Vq          = Vq_map;
% motorModel.(emName).Vr          = Vr_map;
% motorModel.(emName).Vs          = hypot(Vd_map,Vq_map);
% motorModel.(emName).Fd          = Fd_map;
% motorModel.(emName).Fq          = Fq_map;
% motorModel.(emName).Fr          = Fr_map;
% motorModel.(emName).Pm          = Pm_map;
% motorModel.(emName).Pe          = Pe_map;
% motorModel.(emName).Qe          = Qe_map;
% motorModel.(emName).eta         = Eff_map;
% motorModel.(emName).PF          = abs(pf_map);
% motorModel.(emName).Ploss       = Plt_map;
% motorModel.(emName).Pfes        = Pfes_map;
% motorModel.(emName).Pfer        = Pfer_map;
% motorModel.(emName).Pjs         = Pjs_map;
% motorModel.(emName).Pjr         = Pjr_map;
% motorModel.(emName).Psh         = Psh_map;
% motorModel.(emName).Pse         = Pse_map;
% motorModel.(emName).Prh         = Prh_map;
% motorModel.(emName).Pre         = Pre_map;
% motorModel.(emName).Teta_s      = Ts_ref;
% motorModel.(emName).Teta_r      = Tr_ref;
% 
% save([path_i '\' filename], 'motorModel','-append');
% %save([path_i '\' filename '_' emName '_AllVariables'])



% hfig(1)=figure();
% grid on
% figSetting(15,8)
% plot(unique(nm_map), max(Tm_map.*Fnan), 'k', LineWidth=3)
% if min(Tm_map.*Fnan) < 0
%     plot(unique(nm_map), min(Tm_map.*Fnan), 'k', LineWidth=3)
% end
% contourf(nm_map, Tm_map, Eff_map, linspace(0, 1, 1000), 'LineColor','none');
% contour(nm_map, Tm_map, Eff_map,[64:4:80 82:2:88 89:1:100]/100,'k','ShowText','on');
% title(['Efficiency (p.u.) @ Ts= ' num2str(Ts_ref) '°C, Tr= ' num2str(Tr_ref) '°C'], 'Interpreter', 'none')
% ylabel('Torque, (Nm)')
% xlabel('Speed, (rpm)')
% ylim([1.05*min(Tm_map(:)) 1.05*max(Tm_map(:))])
% colorbar
% set(hfig(1),'FileName',[filename '.fig'],'Name','Efficiency Map');
% %saveas(gcf, [path_i emName '_eta.jpg']);
% %saveas(gcf, [path_i emName '_eta.fig']);
% 
% hfig(2)=figure();
% grid on
% figSetting(15,8)
% plot(unique(nm_map), max(Tm_map.*Fnan), 'k', LineWidth=3)
% if min(Tm_map.*Fnan) < 0
%     plot(unique(nm_map), min(Tm_map.*Fnan), 'k', LineWidth=3)
% end
% contourf(nm_map, Tm_map, Id_map, 'LineColor','none');
% contour(nm_map, Tm_map, Id_map,'k','ShowText','on');
% title(['Id (A) @ Ts= ' num2str(Ts_ref) '°C, Tr= ' num2str(Tr_ref) '°C'], 'Interpreter', 'none')
% ylabel('Torque, (Nm)')
% xlabel('Speed, (rpm)')
% ylim([1.05*min(Tm_map(:)) 1.05*max(Tm_map(:))])
% colorbar
% set(hfig(2),'FileName',[filename '.fig'],'Name','d-axis Current Map');
% 
% hfig(3)=figure();
% grid on
% figSetting(15,8)
% plot(unique(nm_map), max(Tm_map.*Fnan), 'k', LineWidth=3)
% if min(Tm_map.*Fnan) < 0
%     plot(unique(nm_map), min(Tm_map.*Fnan), 'k', LineWidth=3)
% end
% contourf(nm_map, Tm_map, Iq_map, 'LineColor','none');
% contour(nm_map, Tm_map, Iq_map,'k','ShowText','on');
% title(['Iq (A) @ Ts= ' num2str(Ts_ref) '°C, Tr= ' num2str(Tr_ref) '°C'], 'Interpreter', 'none')
% ylabel('Torque, (Nm)')
% xlabel('Speed, (rpm)')
% ylim([1.05*min(Tm_map(:)) 1.05*max(Tm_map(:))])
% colorbar
% set(hfig(3),'FileName',[filename '.fig'],'Name','q-axis Current Map')
% 
% hfig(4)=figure();
% grid on
% figSetting(15,8)
% plot(unique(nm_map), max(Tm_map.*Fnan), 'k', LineWidth=3)
% if min(Tm_map.*Fnan) < 0
%     plot(unique(nm_map), min(Tm_map.*Fnan), 'k', LineWidth=3)
% end
% contourf(nm_map, Tm_map, Ir_map, 'LineColor','none');
% contour(nm_map, Tm_map, Ir_map,'k','ShowText','on');
% title(['Ir (A) @ Ts= ' num2str(Ts_ref) '°C, Tr= ' num2str(Tr_ref) '°C'], 'Interpreter', 'none')
% ylabel('Torque, (Nm)')
% xlabel('Speed, (rpm)')
% ylim([1.05*min(Tm_map(:)) 1.05*max(Tm_map(:))])
% colorbar
% set(hfig(4),'FileName',[filename '.fig'],'Name','Rotor Current Map')
% 
% hfig(5)=figure();
% grid on
% figSetting(15,8)
% plot(unique(nm_map), max(Tm_map.*Fnan), 'k', LineWidth=3)
% if min(Tm_map.*Fnan) < 0
%     plot(unique(nm_map), min(Tm_map.*Fnan), 'k', LineWidth=3)
% end
% contourf(nm_map, Tm_map, Vr_map, 'LineColor','none');
% contour(nm_map, Tm_map, Vr_map,'k','ShowText','on');
% title(['Vr (V) @ Ts= ' num2str(Ts_ref) '°C, Tr= ' num2str(Tr_ref) '°C'], 'Interpreter', 'none')
% ylabel('Torque, (Nm)')
% xlabel('Speed, (rpm)')
% ylim([1.05*min(Tm_map(:)) 1.05*max(Tm_map(:))])
% colorbar
% set(hfig(5),'FileName',[filename '.fig'],'Name','Rotor Voltage Map')
% 
% hfig(6)=figure();
% grid on
% figSetting(15,8)
% plot(unique(nm_map), max(Tm_map.*Fnan), 'k', LineWidth=3)
% if min(Tm_map.*Fnan) < 0
%     plot(unique(nm_map), min(Tm_map.*Fnan), 'k', LineWidth=3)
% end
% contourf(nm_map, Tm_map, Plt_map/1000, 'LineColor','none');
% contour(nm_map, Tm_map, Plt_map/1000,'k','ShowText','on');
% title(['Ploss (kW) @ Ts= ' num2str(Ts_ref) '°C, Tr= ' num2str(Tr_ref) '°C'], 'Interpreter', 'none')
% ylabel('Torque, (Nm)')
% xlabel('Speed, (rpm)')
% ylim([1.05*min(Tm_map(:)) 1.05*max(Tm_map(:))])
% colorbar
% set(hfig(6),'FileName',[filename '.fig'],'Name','Total Losses Map')
% 
% % hfig(7)=figure();
% % grid on
% % figSetting(15,8)
% % plot(unique(nm_map), max(Tm_map.*Fnan), 'k', LineWidth=3)
% % if min(Tm_map.*Fnan) < 0
% %     plot(unique(nm_map), min(Tm_map.*Fnan), 'k', LineWidth=3)
% % end
% % contourf(nm_map, Tm_map, pf_map, 'LineColor','none');
% % contour(nm_map, Tm_map, pf_map,'k','ShowText','on');
% % title(['cos$\phi$ @ Ts= ' num2str(Ts_ref) '°C, Tr= ' num2str(Tr_ref) '°C'], 'Interpreter', 'none')
% % ylabel('Torque, (Nm)')
% % xlabel('Speed, (rpm)')
% % ylim([1.05*min(Tm_map(:)) 1.05*max(Tm_map(:))])
% % colorbar
% % set(hfig(7),'FileName',[filename '.fig'],'Name','Power Factor Map')
% 
% 
% hfig(7)=figure();
% grid on
% figSetting(15,8)
% plot(unique(nm_map), max(Tm_map.*Fnan), 'k', LineWidth=2)
% % if min(Tm_map.*Fnan) < 0
% %     plot(unique(nm_map), min(Tm_map.*Fnan), 'k', LineWidth=3)
% % end
% % contourf(nm_map, Tm_map, Eff_map, linspace(0, 1, 1000), 'LineColor','none');
% % contour(nm_map, Tm_map, Eff_map,[64:4:80 82:2:88 89:1:100]/100,'k','ShowText','on');
% title(['Torque Limit @ $I_{max}$ = ' num2str(round(motorModel.data.i0)) 'A, $V_{DC}$ = ' num2str(round(motorModel.data.Vdc)) 'V'], 'Interpreter','latex')
% ylabel('Torque, (Nm)')
% xlabel('Speed, (rpm)')
% ylim([0 1.05*max(Tm_map(:))])
% set(hfig(7),'FileName',[filename '.fig'],'Name','Torque Limit');
% 
% hfig(8)=figure();
% grid on
% figSetting(15,8)
% plot(unique(nm_map), max(Tm_map.*Fnan)'.*unique(nm_map)*pi/30, 'k', LineWidth=2)
% % if min(Tm_map.*Fnan) < 0
% %     plot(unique(nm_map), min(Tm_map.*Fnan), 'k', LineWidth=3)
% % end
% % contourf(nm_map, Tm_map, Eff_map, linspace(0, 1, 1000), 'LineColor','none');
% % contour(nm_map, Tm_map, Eff_map,[64:4:80 82:2:88 89:1:100]/100,'k','ShowText','on');
% title(['Power Limit @ $I_{max}$ = ' num2str(round(motorModel.data.i0)) 'A, $V_{DC}$ = ' num2str(round(motorModel.data.Vdc)) 'V'], 'Interpreter','latex')
% ylabel('Torque, (Nm)')
% xlabel('Speed, (rpm)')
% ylim([0 1.05*max(max(Tm_map.*Fnan)'.*unique(nm_map)*pi/30)])
% set(hfig(8),'FileName',[filename '.fig'],'Name','Power Limit');