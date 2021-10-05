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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [AOA] = MMM_eval_AOA(motorModel,method)

%% Load data
fdfq = motorModel.FluxMap_dq;
axisType = motorModel.data.axisType;

Id   = fdfq.Id;
Iq   = fdfq.Iq;
Fd   = fdfq.Fd;
Fq   = fdfq.Fq;
T    = fdfq.T;
dTpp = fdfq.dTpp;
PF   = sin(atan2(Iq,Id)-atan2(Fq,Fd));


%% Extract curves

%MTPA
[id,iq] = calcOptCtrl(Id,Iq,T,abs(Id+j*Iq),axisType);
if ((id(1)~=0)&&(iq(1)~=0))
    id = [0 id];
    iq = [0 iq];
end

MTPA.id   = id;
MTPA.iq   = iq;

% MTPA.fd   = interp2(Id,Iq,Fd,id,iq);
% MTPA.fq   = interp2(Id,Iq,Fq,id,iq);
% MTPA.T    = interp2(Id,Iq,T,id,iq);
% MTPA.dTpp = interp2(Id,Iq,dTpp,id,iq);

% MTPV
[id,iq] = calcOptCtrl(Id,Iq,T,abs(Fd+j*Fq),axisType);

MTPV.id = id;
MTPV.iq = iq;

% MPFPA
[id,iq] = calcOptCtrl(Id,Iq,PF,abs(Id+j*Iq),axisType);

MPFPA.id = id;
MPFPA.iq = iq;

%% interpolate (if needed)

nPoints = 101;


if strcmp(method,'Fit')
    ft = fittype('poly8');
    if strcmp(axisType,'SR')
        opts = fitoptions(...
            'Method','LinearLeastSquares',...
            'Lower',[zeros(1,8) 0],...
            'Upper',[inf*ones(1,8),0]);
        [xData,yData] = prepareCurveData(MTPA.id,MTPA.iq);
        fitFun = fit(xData,yData,ft,opts);
        MTPA.id = linspace(0,max(MTPA.id),nPoints)';
        MTPA.iq = fitFun(MTPA.id);
        if ~isempty(MTPV.id)
            opts = fitoptions(...
                'Method','LinearLeastSquares',...
                'Lower',[zeros(1,8)],...
                'Upper',[inf*ones(1,8)]);
            [xData,yData] = prepareCurveData(MTPV.id,MTPV.iq);
            fitFun = fit(xData,yData,ft,opts);
            MTPV.id = linspace(0,max(MTPV.id),nPoints)';
            MTPV.iq = fitFun(MTPV.id);
        end
    else
        opts = fitoptions(...
            'Method','LinearLeastSquares',...
            'Upper',[zeros(1,8) 0],...
            'Lower',[-inf*ones(1,8),0]);
        [xData,yData] = prepareCurveData(MTPA.iq,MTPA.id);
        fitFun = fit(xData,yData,ft,opts);
        MTPA.iq = linspace(0,max(MTPA.iq),nPoints)';
        MTPA.id = fitFun(MTPA.iq);
        if ~isempty(MTPV.iq)
            opts = fitoptions(...
                'Method','LinearLeastSquares',...
                'Upper',[zeros(1,8)],...
                'Lower',[-inf*ones(1,8)]);
            [xData,yData] = prepareCurveData(MTPV.iq,MTPV.id);
            fitFun = fit(xData,yData,ft,opts);
            MTPV.iq = linspace(0,max(MTPV.iq),nPoints)';
            MTPV.id = fitFun(MTPV.iq);
        end
    end
    AOA.method = 'Fit';
else
    AOA.method = 'LUT';
end


%% Interp flux linkages and torque

MTPA.fd    = interp2(Id,Iq,Fd,MTPA.id,MTPA.iq);
MTPA.fq    = interp2(Id,Iq,Fq,MTPA.id,MTPA.iq);
MTPA.T     = interp2(Id,Iq,T,MTPA.id,MTPA.iq);
MTPA.dTpp  = interp2(Id,Iq,dTpp,MTPA.id,MTPA.iq);

MTPV.fd    = interp2(Id,Iq,Fd,MTPV.id,MTPV.iq);
MTPV.fq    = interp2(Id,Iq,Fq,MTPV.id,MTPV.iq);
MTPV.T     = interp2(Id,Iq,T,MTPV.id,MTPV.iq);
MTPV.dTpp  = interp2(Id,Iq,dTpp,MTPV.id,MTPV.iq);

MPFPA.fd   = interp2(Id,Iq,Fd,MPFPA.id,MPFPA.iq);
MPFPA.fq   = interp2(Id,Iq,Fq,MPFPA.id,MPFPA.iq);
MPFPA.T    = interp2(Id,Iq,T,MPFPA.id,MPFPA.iq);
MPFPA.dTpp = interp2(Id,Iq,dTpp,MPFPA.id,MPFPA.iq);

AOA.MTPA  = MTPA;
AOA.MTPV  = MTPV;
AOA.MPFPA = MPFPA;




















