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

function [dqtMapF] = MMM_eval_inverse_dqtMap(motorModel)

dqtMap = motorModel.FluxMap_dqt;

[xS,yS,zS]=size(dqtMap.data.Id);

data=dqtMap.data;
dataF.th = data.th;
dataF.Id = zeros([xS,yS,zS]);
dataF.Iq = zeros([xS,yS,zS]);
dataF.Fd = zeros([xS,yS,zS]);
dataF.Fq = zeros([xS,yS,zS]);
dataF.T  = zeros([xS,yS,zS]);

% compute the regular grid limits

tmp   = min(data.Fd,[],3,'omitnan');    % min on th
tmp   = max(tmp,[],1,'omitnan');        % max on d-axis
FdMax = min(tmp,[],2,'omitnan');        % min on q-axis

tmp   = max(data.Fd,[],3,'omitnan');    % max on th
tmp   = min(tmp,[],1,'omitnan');        % min on d-axis
FdMin = max(tmp,[],2,'omitnan');        % max on q-axis

tmp   = min(data.Fq,[],3,'omitnan');    % min on th
tmp   = max(tmp,[],2,'omitnan');        % max on q-axis
FqMax = min(tmp,[],1,'omitnan');        % min on d-axis

tmp   = max(data.Fq,[],3,'omitnan');    % max on th
tmp   = min(tmp,[],2,'omitnan');        % min on q-axis
FqMin = max(tmp,[],1,'omitnan');        % max on d-axis


dqtMapF.Fd=linspace(FdMin,FdMax,xS);
dqtMapF.Fq=linspace(FqMin,FqMax,yS);
dqtMapF.th=dqtMap.th;
[dataF.Fd,dataF.Fq,dataF.th]=ndgrid(dqtMapF.Fd,dqtMapF.Fq,dqtMapF.th);

index = 1:1:numel(data.Id(:,:,1));

% disp('Inversion process running...')
for zz=1:zS    
    % extract the layer th=th(zz) of the 3D matrix
    IdVect = data.Id(:,:,zz);
    IqVect = data.Iq(:,:,zz);
    FdVect = data.Fd(:,:,zz);
    FqVect = data.Fq(:,:,zz);
    TVect  = data.T(:,:,zz);
    FdGrid = dataF.Fd(:,:,zz);
    FqGrid = dataF.Fq(:,:,zz);
    
    % matrix-->vector (for scatteredInterpolant)
    IdVect = IdVect(index)';
    IqVect = IqVect(index)';
    FdVect = FdVect(index)';
    FqVect = FqVect(index)';
    TVect  = TVect(index)';
    
    % filt NaN
    IdVect = IdVect(~isnan(FdVect));
    IqVect = IqVect(~isnan(FdVect));
    TVect  = TVect(~isnan(FdVect));
    FqVect = FqVect(~isnan(FdVect));
    FdVect = FdVect(~isnan(FdVect));
    
    % inversion with scatteredInterpolant
    intD = scatteredInterpolant(FdVect,FqVect,IdVect,'linear','none');
    intQ = scatteredInterpolant(FdVect,FqVect,IqVect,'linear','none');
    intT = scatteredInterpolant(FdVect,FqVect,TVect,'linear','none');
    
    % fill the layer th=th(zz) of the inverted 3D matrix
    dataF.Id(:,:,zz) = intD(dataF.Fd(:,:,zz),dataF.Fq(:,:,zz));
    dataF.Iq(:,:,zz) = intQ(dataF.Fd(:,:,zz),dataF.Fq(:,:,zz));
    dataF.T(:,:,zz)  = intT(dataF.Fd(:,:,zz),dataF.Fq(:,:,zz));

%     disp([num2str(zz,'%03.0f') ' of ' num2str(zS,'%03.0f') ' done'])

end

% disp('Inversion done!')

%% Interpolant
% fInt.Id = griddedInterpolant(dataF.Fd,dataF.Fq,dataF.th,dataF.Id,'linear','none');
% fInt.Iq = griddedInterpolant(dataF.Fd,dataF.Fq,dataF.th,dataF.Iq,'linear','none');
% fInt.th = griddedInterpolant(dataF.Fd,dataF.Fq,dataF.th,dataF.th,'linear','none');
% fInt.Fd = griddedInterpolant(dataF.Fd,dataF.Fq,dataF.th,dataF.Fd,'linear','none');
% fInt.Fq = griddedInterpolant(dataF.Fd,dataF.Fq,dataF.th,dataF.Fq,'linear','none');
% fInt.T  = griddedInterpolant(dataF.Fd,dataF.Fq,dataF.th,dataF.T,'linear','none');

%% Save data
dqtMapF.th   = dqtMap.th;
dqtMapF.Fd   = unique(dataF.Fd);
dqtMapF.Fq   = unique(dataF.Fq);
% dqtMapF.fInt = fInt;
dqtMapF.dataF = dataF;















