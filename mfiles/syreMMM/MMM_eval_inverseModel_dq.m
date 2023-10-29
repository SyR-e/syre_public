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

function [idiq] = MMM_eval_inverseModel_dq(motorModel,setup)

if nargin==1
    setup.flagE = 0;
    setup.extrapolationMethod = 'none';
    % Extrapolation of the flux map controlled from flagE:
    % - flagE==0 --> no extrapolation, flux map is reduced to have a complete rectangular domain
    % - flagE==1 --> flux map domain is not reduced and flux maps are extrapolated according to extrapolationMethod
    %
    % extrapolationMethod='none' --> flux maps are not extrapolated and NaN are added at the corners
end

flagE = setup.flagE;
extrapolationMethod = setup.extrapolationMethod;

Id = motorModel.FluxMap_dq.Id;
Iq = motorModel.FluxMap_dq.Iq;
Fd = motorModel.FluxMap_dq.Fd;
Fq = motorModel.FluxMap_dq.Fq;
T  = motorModel.FluxMap_dq.T;

% flagE = 0; % Extrapolation flag. If 1, the flux linkage limits are not restricted

% compute the regular grid limits
if flagE
    tmp   = max(Fd,[],2,'omitnan');     % max on d-axis
    FdMax = max(tmp,[],1,'omitnan');    % min on q-axis

    tmp   = min(Fd,[],2,'omitnan');     % min on d-axis
    FdMin = min(tmp,[],1,'omitnan');    % max on q-axis

    tmp   = max(Fq,[],1,'omitnan');     % max on q-axis
    FqMax = max(tmp,[],2,'omitnan');    % min on d-axis

    tmp   = min(Fq,[],1,'omitnan');     % min on q-axis
    FqMin = min(tmp,[],2,'omitnan');    % max on d-axis
    
%     extrapolationMethod = 'linear';
%     extrapolationMethod = 'none';
else
    tmp   = max(Fd,[],2,'omitnan');     % max on d-axis
    FdMax = min(tmp,[],1,'omitnan');    % min on q-axis

    tmp   = min(Fd,[],2,'omitnan');     % min on d-axis
    FdMin = max(tmp,[],1,'omitnan');    % max on q-axis

    tmp   = max(Fq,[],1,'omitnan');     % max on q-axis
    FqMax = min(tmp,[],2,'omitnan');    % min on d-axis

    tmp   = min(Fq,[],1,'omitnan');     % min on q-axis
    FqMin = max(tmp,[],2,'omitnan');    % max on d-axis

%     FdMax = max(Fd(:),[],'omitnan');
%     FdMin = min(Fd(:),[],'omitnan');
%     FqMax = max(Fq(:),[],'omitnan');
%     FqMin = min(Fq(:),[],'omitnan');
    
%     extrapolationMethod = 'none';
end

fD = linspace(FdMin,FdMax,256);
fQ = linspace(FqMin,FqMax,256);
[fD,fQ]=meshgrid(fD,fQ);

% matrix-->vector for scatteredInterpolant
index=1:1:numel(Id);

IdVect = Id(index)';
IqVect = Iq(index)';
FdVect = Fd(index)';
FqVect = Fq(index)';
TVect  = T(index)';


% filt NaN
IdVect = IdVect(~isnan(FdVect));
IqVect = IqVect(~isnan(FdVect));
FqVect = FqVect(~isnan(FdVect));
FdVect = FdVect(~isnan(FdVect));
TVect  = TVect(~isnan(FdVect));

% interpolant functions
intD = scatteredInterpolant(FdVect,FqVect,IdVect,'linear',extrapolationMethod);
intQ = scatteredInterpolant(FdVect,FqVect,IqVect,'linear',extrapolationMethod);
intT = scatteredInterpolant(FdVect,FqVect,TVect,'linear',extrapolationMethod);

iD = intD(fD,fQ);
iQ = intQ(fD,fQ);
Tf = intT(fD,fQ);

% output data
idiq.Fd = fD;
idiq.Fq = fQ;
idiq.Id = iD;
idiq.Iq = iQ;
idiq.T  = Tf;









