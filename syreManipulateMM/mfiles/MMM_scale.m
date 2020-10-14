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

kr  = scaleFactors.Ns/motorModel.data.Ns;
kl  = scaleFactors.l/motorModel.data.l;
Lld = scaleFactors.Lld;
Llq = scaleFactors.Llq;

% 1)  change fdfq
%  a) change number of turns
motorModel.fdfq.Id = motorModel.fdfq.Id/kr;
motorModel.fdfq.Iq = motorModel.fdfq.Iq/kr;
motorModel.fdfq.Fd = motorModel.fdfq.Fd*kr;
motorModel.fdfq.Fq = motorModel.fdfq.Fq*kr;
%  b) change stack length
motorModel.fdfq.Fd = motorModel.fdfq.Fd*kl;
motorModel.fdfq.Fq = motorModel.fdfq.Fq*kl;
motorModel.fdfq.T  = motorModel.fdfq.T*kl;
%  c) add leakage inductances
motorModel.fdfq.Fd = motorModel.fdfq.Fd+Lld*motorModel.fdfq.Id;
motorModel.fdfq.Fq = motorModel.fdfq.Fq+Llq*motorModel.fdfq.Iq;
Tnew = 3/2*motorModel.data.p*motorModel.data.n3phase*(motorModel.fdfq.Fd.*motorModel.fdfq.Iq-motorModel.fdfq.Fq.*motorModel.fdfq.Id);
motorModel.fdfq.T = Tnew;

% 2)  check for iron loss update
if ~isempty(motorModel.ironLoss)
    if strcmp(motorModel.ironLoss.type,'map')
        ironLoss = motorModel.ironLoss;
        %  a) change number of turns
        ironLoss.Id = ironLoss.Id/kr;
        ironLoss.Iq = ironLoss.Iq/kr;
        %  b) change stack length
        ironLoss.Pfes_h = ironLoss.Pfes_h*kl;
        ironLoss.Pfes_c = ironLoss.Pfes_c*kl;
        ironLoss.Pfer_h = ironLoss.Pfer_h*kl;
        ironLoss.Pfer_c = ironLoss.Pfer_c*kl;
        ironLoss.Ppm    = ironLoss.Ppm*kl;

        motorModel.ironLoss = ironLoss;
    else
        motorModel.ironLoss = [];
    end
end

% 3) check for dqtMap model update
if ~isempty(motorModel.dqtMap)
    dqtMap = motorModel.dqtMap;
    %  a) change number of turns
    dqtMap.Id = dqtMap.Id/kr;
    dqtMap.Iq = dqtMap.Iq/kr;
    dqtMap.data.Id = dqtMap.data.Id/kr;
    dqtMap.data.Iq = dqtMap.data.Iq/kr;
    dqtMap.data.Ia = dqtMap.data.Ia/kr;
    dqtMap.data.Ib = dqtMap.data.Ib/kr;
    dqtMap.data.Ic = dqtMap.data.Ic/kr;
    dqtMap.data.Fd = dqtMap.data.Fd*kr;
    dqtMap.data.Fq = dqtMap.data.Fq*kr;
    dqtMap.data.Fa = dqtMap.data.Fa*kr;
    dqtMap.data.Fb = dqtMap.data.Fb*kr;
    dqtMap.data.Fc = dqtMap.data.Fc*kr;
    %  b) change stack length
    dqtMap.data.Fd = dqtMap.data.Fd*kl;
    dqtMap.data.Fq = dqtMap.data.Fq*kl;
    dqtMap.data.Fa = dqtMap.data.Fa*kl;
    dqtMap.data.Fb = dqtMap.data.Fb*kl;
    dqtMap.data.Fc = dqtMap.data.Fc*kl;
    dqtMap.data.T  = dqtMap.data.T*kl;
    %  c) add leakage inductance
    dqtMap.data.Fd = dqtMap.data.Fd+Lld*dqtMap.data.Id;
    dqtMap.data.Fq = dqtMap.data.Fq+Llq*dqtMap.data.Iq;
    Told = mean(dqtMap.data.T,3);
    Tnew = 3/2*motorModel.data.p*motorModel.data.n3phase*(mean(dqtMap.data.Fd,3).*mean(dqtMap.data.Iq,3)-mean(dqtMap.data.Fq,3).*mean(dqtMap.data.Id,3));
    dqtMap.data.T = dqtMap.data.T.*repmat(Tnew,[1,1,size(dqtMap.data.T,3)])./repmat(Told,[1,1,size(dqtMap.data.T,3)]);
    % re-compute fInt
    dqtMap.fInt.Id = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Id,'spline');
    dqtMap.fInt.Iq = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Iq,'spline');
    dqtMap.fInt.th = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.th,'spline');
    dqtMap.fInt.Fd = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fd,'spline');
    dqtMap.fInt.Fq = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fq,'spline');
    dqtMap.fInt.T  = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.T,'spline');
    dqtMap.fInt.Fa = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fa,'spline');
    dqtMap.fInt.Fb = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fb,'spline');
    dqtMap.fInt.Fc = griddedInterpolant(dqtMap.data.Id,dqtMap.data.Iq,dqtMap.data.th,dqtMap.data.Fc,'spline');
    
    motorModel.dqtMap = dqtMap;
end

% 4) update motor data
motorModel.scale = scaleFactors;
% 
motorModel.data.Ns = motorModel.data.Ns*kr;
motorModel.data.l  = motorModel.data.l*kl;




% motorModel.data.pathname = 'unvalidPath';

geo = motorModel.geo;
per = motorModel.per;
geo.l = geo.l*kl;
geo.win.Ns = geo.win.Ns*kr;
[i0,Rs] = calc_io(geo,per);
motorModel.data.i0 = i0;
motorModel.data.Rs = Rs;

motorModel.dataSet.StackLength   = motorModel.data.l;
motorModel.dataSet.TurnsInSeries = motorModel.data.Ns;
motorModel.dataSet.RatedCurrent  = i0;
motorModel.dataSet.Rs            = Rs;
motorModel.geo.l                 = motorModel.data.l;
motorModel.geo.win.Ns            = motorModel.data.Ns;
motorModel.per.Rs                = Rs;
motorModel.per.i0                = i0;

% motorModel.data.Imax = motorModel.data.Imax/kr;
motorModel.data.Vdc = motorModel.data.Vdc*kr*kl;

[motorModel] = MMM_eval_Ratings(motorModel,0);

% 5) delete other models
motorModel.AOA        = [];
motorModel.Inductance = [];
motorModel.idiq       = [];
motorModel.dqtMapF    = [];
motorModel.skinEffect = [];







