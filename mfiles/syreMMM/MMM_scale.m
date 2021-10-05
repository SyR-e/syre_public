% Copyright 2020
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

function [motorModel] = MMM_scale(motorModel,scaleFactors)
% Change number of turns, stack length and dq leakage inductances of a
% motor model.

kN  = scaleFactors.Ns/motorModel.data.Ns;
kL  = scaleFactors.l/motorModel.data.l;
kD  = 1;
Lld = scaleFactors.Lld;
Llq = scaleFactors.Llq;

% 1)  change fdfq
%  a) change number of turns
motorModel.FluxMap_dq.Id = motorModel.FluxMap_dq.Id/kN;
motorModel.FluxMap_dq.Iq = motorModel.FluxMap_dq.Iq/kN;
motorModel.FluxMap_dq.Fd = motorModel.FluxMap_dq.Fd*kN;
motorModel.FluxMap_dq.Fq = motorModel.FluxMap_dq.Fq*kN;
%  b) change stack length
motorModel.FluxMap_dq.Fd = motorModel.FluxMap_dq.Fd*kL;
motorModel.FluxMap_dq.Fq = motorModel.FluxMap_dq.Fq*kL;
motorModel.FluxMap_dq.T  = motorModel.FluxMap_dq.T*kL;
%  c) change stator radius
motorModel.FluxMap_dq.Id = motorModel.FluxMap_dq.Id*kD;
motorModel.FluxMap_dq.Iq = motorModel.FluxMap_dq.Iq*kD;
motorModel.FluxMap_dq.Fd = motorModel.FluxMap_dq.Fd*kD;
motorModel.FluxMap_dq.Fq = motorModel.FluxMap_dq.Fq*kD;
motorModel.FluxMap_dq.T  = motorModel.FluxMap_dq.T*kD^2;
%  d) add leakage inductances
if Lld || Llq
    motorModel.FluxMap_dq.Fd = motorModel.FluxMap_dq.Fd+Lld*motorModel.FluxMap_dq.Id;
    motorModel.FluxMap_dq.Fq = motorModel.FluxMap_dq.Fq+Llq*motorModel.FluxMap_dq.Iq;
    Tnew = 3/2*motorModel.data.p*motorModel.data.n3phase*(motorModel.FluxMap_dq.Fd.*motorModel.FluxMap_dq.Iq-motorModel.FluxMap_dq.Fq.*motorModel.FluxMap_dq.Id);
    motorModel.FluxMap_dq.T = Tnew;
end

% 2)  check for iron loss update
if ~isempty(motorModel.IronPMLossMap_dq)
    if strcmp(motorModel.IronPMLossMap_dq.type,'map')
        ironLoss = motorModel.IronPMLossMap_dq;
        %  a) change number of turns
        ironLoss.Id = ironLoss.Id/kN;
        ironLoss.Iq = ironLoss.Iq/kN;
        %  b) change stack length
        ironLoss.Pfes_h = ironLoss.Pfes_h*kL;
        ironLoss.Pfes_c = ironLoss.Pfes_c*kL;
        ironLoss.Pfer_h = ironLoss.Pfer_h*kL;
        ironLoss.Pfer_c = ironLoss.Pfer_c*kL;
        ironLoss.Ppm    = ironLoss.Ppm*kL;
        %  c) change stator radius
        ironLoss.Id = ironLoss.Id*kD;
        ironLoss.Iq = ironLoss.Iq*kD;
        ironLoss.Pfes_h = ironLoss.Pfes_h*kD^2;
        ironLoss.Pfes_c = ironLoss.Pfes_c*kD^2;
        ironLoss.Pfer_h = ironLoss.Pfer_h*kD^2;
        ironLoss.Pfer_c = ironLoss.Pfer_c*kD^2;
        ironLoss.Ppm    = ironLoss.Ppm*kD^4;
        
        motorModel.IronPMLossMap_dq = ironLoss;
    else
        motorModel.IronPMLossMap_dq = [];
    end
end

% 3) check for dqtMap model update
if ~isempty(motorModel.FluxMap_dqt)
    dqtMap = motorModel.FluxMap_dqt;
    %  a) change number of turns
    dqtMap.Id = dqtMap.Id/kN;
    dqtMap.Iq = dqtMap.Iq/kN;
    dqtMap.data.Id = dqtMap.data.Id/kN;
    dqtMap.data.Iq = dqtMap.data.Iq/kN;
    if isfield(dqtMap.data,'Ia')
        dqtMap.data.Ia = dqtMap.data.Ia/kN;
        dqtMap.data.Ib = dqtMap.data.Ib/kN;
        dqtMap.data.Ic = dqtMap.data.Ic/kN;
    end
    dqtMap.data.Fd = dqtMap.data.Fd*kN;
    dqtMap.data.Fq = dqtMap.data.Fq*kN;
    if isfield(dqtMap.data,'Fa')
        dqtMap.data.Fa = dqtMap.data.Fa*kN;
        dqtMap.data.Fb = dqtMap.data.Fb*kN;
        dqtMap.data.Fc = dqtMap.data.Fc*kN;
    end
    %  b) change stack length
    dqtMap.data.Fd = dqtMap.data.Fd*kL;
    dqtMap.data.Fq = dqtMap.data.Fq*kL;
    if isfield(dqtMap.data,'Fa')
        dqtMap.data.Fa = dqtMap.data.Fa*kL;
        dqtMap.data.Fb = dqtMap.data.Fb*kL;
        dqtMap.data.Fc = dqtMap.data.Fc*kL;
    end
    dqtMap.data.T  = dqtMap.data.T*kL;
    %  c) change stator radius
    dqtMap.Id = dqtMap.Id*kD;
    dqtMap.Iq = dqtMap.Iq*kD;
    dqtMap.data.Id = dqtMap.data.Id*kD;
    dqtMap.data.Iq = dqtMap.data.Iq*kD;
    if isfield(dqtMap.data,'Ia')
        dqtMap.data.Ia = dqtMap.data.Ia*kD;
        dqtMap.data.Ib = dqtMap.data.Ib*kD;
        dqtMap.data.Ic = dqtMap.data.Ic*kD;
    end
    dqtMap.data.Fd = dqtMap.data.Fd*kD;
    dqtMap.data.Fq = dqtMap.data.Fq*kD;
    if isfield(dqtMap.data,'Fa')
        dqtMap.data.Fa = dqtMap.data.Fa*kD;
        dqtMap.data.Fb = dqtMap.data.Fb*kD;
        dqtMap.data.Fc = dqtMap.data.Fc*kD;
    end
    %  d) add leakage inductance
    if Lld || Llq
        dqtMap.data.Fd = dqtMap.data.Fd+Lld*dqtMap.data.Id;
        dqtMap.data.Fq = dqtMap.data.Fq+Llq*dqtMap.data.Iq;
        Told = mean(dqtMap.data.T,3);
        Tnew = 3/2*motorModel.data.p*motorModel.data.n3phase*(mean(dqtMap.data.Fd,3).*mean(dqtMap.data.Iq,3)-mean(dqtMap.data.Fq,3).*mean(dqtMap.data.Id,3));
        dqtMap.data.T = dqtMap.data.T.*repmat(Tnew,[1,1,size(dqtMap.data.T,3)])./repmat(Told,[1,1,size(dqtMap.data.T,3)]);
        dqtMap.data.T(isnan(dqtMap.data.T)) = 0;
    end
    % re-compute fInt
    dqtMap.fInt.Id = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Id,'spline');
    dqtMap.fInt.Iq = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Iq,'spline');
    dqtMap.fInt.th = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.th,'spline');
    dqtMap.fInt.Fd = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fd,'spline');
    dqtMap.fInt.Fq = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fq,'spline');
    dqtMap.fInt.T  = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.T,'spline');
    if isfield(dqtMap.data,'Fa')
        dqtMap.fInt.Fa = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fa,'spline');
        dqtMap.fInt.Fb = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fb,'spline');
        dqtMap.fInt.Fc = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fc,'spline');
    end
    motorModel.FluxMap_dqt = dqtMap;
end

% 4) update motor data
motorModel.tmpScale = scaleFactors;
% 
motorModel.data.Ns = motorModel.data.Ns*kN;
motorModel.data.l  = motorModel.data.l*kL;
motorModel.data.R  = motorModel.data.R*kD;


[motorModel.dataSet] = back_compatibility(motorModel.dataSet,motorModel.geo,motorModel.per,0);

%Radial scaling
motorModel.dataSet.StatorOuterRadius    =   motorModel.data.R;
motorModel.dataSet.ShaftRadius          =   motorModel.dataSet.ShaftRadius*kD;
motorModel.dataSet.AirGapRadius         =   motorModel.dataSet.AirGapRadius*kD;
motorModel.dataSet.AirGapThickness      =   motorModel.dataSet.AirGapThickness*kD;
motorModel.dataSet.ToothLength          =   motorModel.dataSet.ToothLength*kD;
motorModel.dataSet.ToothWidth           =   motorModel.dataSet.ToothWidth*kD;
motorModel.dataSet.ToothTangDepth       =   motorModel.dataSet.ToothTangDepth*kD;
motorModel.dataSet.FilletCorner         =   motorModel.dataSet.FilletCorner*kD;
motorModel.dataSet.RadShiftInner        =   motorModel.dataSet.RadShiftInner*kD;
motorModel.dataSet.TanRibEdit           =   motorModel.dataSet.TanRibEdit*kD;
motorModel.dataSet.RotorFilletTan1      =   motorModel.dataSet.RotorFilletTan1*kD;
motorModel.dataSet.RotorFilletTan2      =   motorModel.dataSet.RotorFilletTan2*kD;
motorModel.dataSet.RadRibEdit           =   motorModel.dataSet.RadRibEdit*kD;
motorModel.dataSet.pontRoffsetEdit      =   motorModel.dataSet.pontRoffsetEdit*kD;
motorModel.dataSet.RotorFilletIn        =   motorModel.dataSet.RotorFilletIn*kD;
motorModel.dataSet.RotorFilletOut       =   motorModel.dataSet.RotorFilletOut*kD;
motorModel.dataSet.Mesh                 =   motorModel.dataSet.Mesh*kD;
motorModel.dataSet.Mesh_MOOA            =   motorModel.dataSet.Mesh_MOOA*kD;
motorModel.dataSet.MinMechTol           =   motorModel.dataSet.MinMechTol*kD ;
motorModel.dataSet.PMdim                =   motorModel.dataSet.PMdim*kD;

motorModel.geo.R                        =   motorModel.dataSet.StatorOuterRadius;
motorModel.geo.Ar                       =   motorModel.dataSet.ShaftRadius;
motorModel.geo.r                        =   motorModel.dataSet.AirGapRadius;
motorModel.geo.g                        =   motorModel.dataSet.AirGapThickness;
motorModel.geo.lt                       =   motorModel.dataSet.ToothLength;
motorModel.geo.wt                       =   motorModel.dataSet.ToothWidth;
motorModel.geo.ttd                      =   motorModel.dataSet.ToothTangDepth;
motorModel.geo.SFR                      =   motorModel.dataSet.FilletCorner;
motorModel.geo.dxIB                     =   motorModel.dataSet.RadShiftInner; 
motorModel.geo.pontT                    =   motorModel.dataSet.TanRibEdit;
motorModel.geo.RotorFilletTan1          =   motorModel.dataSet.RotorFilletTan1; 
motorModel.geo.RotorFilletTan2          =   motorModel.dataSet.RotorFilletTan2; 
motorModel.geo.pontR                    =   motorModel.dataSet.RadRibEdit; 
motorModel.geo.pontRoffset              =   motorModel.dataSet.pontRoffsetEdit; 
motorModel.geo.RotorFillet1             =   motorModel.dataSet.RotorFilletIn; 
motorModel.geo.RotorFillet2             =   motorModel.dataSet.RotorFilletOut; 
motorModel.geo.mesh_K                   =   motorModel.dataSet.Mesh; 
motorModel.geo.mesh_K_MOOA              =   motorModel.dataSet.Mesh_MOOA; 
motorModel.geo.pont0                    =   motorModel.dataSet.MinMechTol; 
motorModel.geo.PMdim                    =   motorModel.dataSet.PMdim;

motorModel.geo.Aslot                    =   motorModel.geo.Aslot*kD^2;
motorModel.geo.win.Ns                   =   motorModel.geo.win.Ns*kN;
motorModel.geo.l                        =   motorModel.geo.l*kL;

motorModel.per.tempcu                   =   motorModel.data.tempCu;

geo = motorModel.geo;
per = motorModel.per;
per.Loss = per.Loss*kD*kL;
% geo.l = geo.l*kL;
% geo.win.Ns = geo.win.Ns*kN;
% geo.R = geo.R*kD;
% geo.Aslot = geo.Aslot*kD^2;

[i0,Rs] = calc_io(geo,per);
motorModel.data.i0 = i0;
motorModel.data.Rs = Rs;

motorModel.dataSet.StackLength   = motorModel.data.l;
motorModel.dataSet.TurnsInSeries = motorModel.data.Ns;
motorModel.dataSet.RatedCurrent  = i0;
motorModel.dataSet.Rs            = Rs;
% motorModel.geo.l                 = motorModel.data.l;
% motorModel.geo.win.Ns            = motorModel.data.Ns;
motorModel.per.Rs                = Rs;
motorModel.per.i0                = i0;


motorModel.data.Imax = motorModel.data.Imax/kN*kD;
motorModel.data.Vdc  = motorModel.data.Vdc*kN*kL*kD;



% Mass and Inertia computation
fem.res = 0;
fem.res_traf = 0;

% nodes
[motorModel.geo.rotor,~,motorModel.geo] = ROTmatr(motorModel.geo,fem,motorModel.mat);
[motorModel.geo,motorModel.geo.stator,~] = STATmatr(motorModel.geo,fem);

motorModel.geo.mCu  = calcMassCu(motorModel.geo,motorModel.mat);
motorModel.geo.mPM  = calcMassPM(motorModel.geo,motorModel.mat);
[motorModel.geo.mFeS,motorModel.geo.mFeR] = calcMassFe(motorModel.geo,motorModel.mat);
motorModel.geo.J    = calcRotorInertia(motorModel.geo,motorModel.mat);

motorModel.dataSet.MassWinding    = motorModel.geo.mCu;
motorModel.dataSet.MassMagnet     = motorModel.geo.mPM;
motorModel.dataSet.MassStatorIron = motorModel.geo.mFeS;
motorModel.dataSet.MassRotorIron  = motorModel.geo.mFeR;
motorModel.dataSet.RotorInertia   = motorModel.geo.J;
motorModel.data.J                 = motorModel.geo.J;



% 5) delete or recompute other models
motorModel.acLossFactor = [];

if ~isempty(motorModel.controlTrajectories)
    motorModel.controlTrajectories = MMM_eval_AOA(motorModel,motorModel.controlTrajectories.method);
end

if ~isempty(motorModel.IncInductanceMap_dq)
    motorModel.IncInductanceMap_dq = MMM_eval_inductanceMap(motorModel);
end

if ~isempty(motorModel.FluxMapInv_dq)
    motorModel.FluxMapInv_dq = MMM_eval_inverseModel_dq(motorModel);
end

if ~isempty(motorModel.FluxMapInv_dqt)
    motorModel.FluxMapInv_dqt = MMM_eval_inverse_dqtMap(motorModel);
end

