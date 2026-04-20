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
if strcmp(motorModel.dataSet.TypeOfRotor,'EESM')
    Ir = motorModel.FluxMap_dq.Ir;
    Fr = motorModel.FluxMap_dq.Fr;
    Ir_length = length(Ir(1,1,:));
    intDIM=41;
else
    Ir_length = 1;
    intDIM=256;
end

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
if strcmp(motorModel.dataSet.TypeOfRotor,'EESM')
    warning('Fr not included yet')
end
%     FdMax = max(Fd(:),[],'omitnan');
%     FdMin = min(Fd(:),[],'omitnan');
%     FqMax = max(Fq(:),[],'omitnan');
%     FqMin = min(Fq(:),[],'omitnan');
    
%     extrapolationMethod = 'none';
end

fD = zeros([intDIM intDIM Ir_length]);
fQ = zeros([intDIM intDIM Ir_length]);
iD = zeros([intDIM intDIM Ir_length]);
iQ = zeros([intDIM intDIM Ir_length]);
Tf = zeros([intDIM intDIM Ir_length]);
% if strcmp(motorModel.dataSet.TypeOfRotor,'EESM')
    for ii = 1:Ir_length
        tmpD = linspace(FdMin(ii),FdMax(ii),intDIM);
        tmpQ = linspace(FqMin(ii),FqMax(ii),intDIM);
        [fD(:,:,ii),fQ(:,:,ii)]=meshgrid(tmpD,tmpQ);
    end
% else
%     fD = linspace(FdMin,FdMax,intDIM);
%     fQ = linspace(FqMin,FqMax,intDIM);
%     [fD,fQ]=meshgrid(fD,fQ);
% end


for ii = 1:Ir_length
    % matrix-->vector for scatteredInterpolant
    index=1:1:numel(Id(:,:,ii));

    IdVect = reshape(Id(:,:,ii),1,[]);
    IqVect = reshape(Iq(:,:,ii),1,[]);
    FdVect = reshape(Fd(:,:,ii),1,[]);
    FqVect = reshape(Fq(:,:,ii),1,[]);
    TVect  = reshape(T(:,:,ii),1,[]);

    IdVect = IdVect';
    IqVect = IqVect';
    FdVect = FdVect';
    FqVect = FqVect';
    TVect  = TVect';


    % filt NaN
    IdVect = IdVect(~isnan(FdVect(:)));
    IqVect = IqVect(~isnan(FdVect));
    FqVect = FqVect(~isnan(FdVect));
    FdVect = FdVect(~isnan(FdVect));
    TVect  = TVect(~isnan(FdVect));

    % interpolant functions
    intD = scatteredInterpolant(FdVect,FqVect,IdVect,'linear',extrapolationMethod);
    intQ = scatteredInterpolant(FdVect,FqVect,IqVect,'linear',extrapolationMethod);
    intT = scatteredInterpolant(FdVect,FqVect,TVect,'linear',extrapolationMethod);

    iD(:,:,ii) = intD(fD(:,:,ii),fQ(:,:,ii));
    iQ(:,:,ii) = intQ(fD(:,:,ii),fQ(:,:,ii));
    Tf(:,:,ii) = intT(fD(:,:,ii),fQ(:,:,ii));
end

% output data
idiq.Fd = fD;
idiq.Fq = fQ;
idiq.Id = iD;
idiq.Iq = iQ;
idiq.T  = Tf;









