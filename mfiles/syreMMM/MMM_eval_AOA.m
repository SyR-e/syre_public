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
fdfq = motorModel.fdfq;
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
[id,iq] = calcOptCtrl(Id,Iq,T,abs(Id+j*Iq));
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
[id,iq] = calcOptCtrl(Id,Iq,T,abs(Fd+j*Fq));

MTPV.id = id;
MTPV.iq = iq;

% MPFPA
[id,iq] = calcOptCtrl(Id,Iq,PF,abs(Id+j*Iq));

MPFPA.id = id;
MPFPA.iq = id;

%% interpolate (if needed)

nPoints = 101;
nPoly = 7;

if strcmp(method,'Fit')
    if strcmp(axisType,'SR')
        p = polyfit(MTPA.id,MTPA.iq,nPoly);
        MTPA.id = linspace(0,max(MTPA.id),nPoints);
        MTPA.iq = polyval(p,MTPA.id);
        
        if ~isempty(MTPV.id)
            p = polyfit(MTPV.id,MTPV.iq,nPoly);
            MTPV.id = linspace(0,max(MTPV.id),nPoints);
            MTPV.iq = polyval(p,MTPV.id);
        end
        
        p = polyfit(MPFPA.id,MPFPA.iq,nPoly);
        MPFPA.id = linspace(0,max(MPFPA.id),nPoints);
        MPFPA.iq = polyval(p,MPFPA.id);
    else
        p = polyfit(MTPA.iq,MTPA.id,nPoly);
        MTPA.iq = linspace(0,max(MTPA.iq),nPoints);
        MTPA.id = polyval(p,MTPA.iq);
        
        if ~isempty(MTPV.id)
            p = polyfit(MTPV.iq,MTPV.id,nPoly);
            MTPV.iq = linspace(0,max(MTPV.iq),nPoints);
            MTPV.id = polyval(p,MTPV.iq);
        end
        
        p = polyfit(MPFPA.iq,MPFPA.id,nPoly);
        MPFPA.iq = linspace(0,max(MPFPA.iq),nPoints);
        MPFPA.id = polyval(p,MPFPA.iq);
        
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




















