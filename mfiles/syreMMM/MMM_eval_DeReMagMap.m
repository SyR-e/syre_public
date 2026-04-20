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

function [Control] = MMM_eval_DeReMagMap(motorModel,n)

if nargin == 1
    n = 51;
end
%% Load data

DemId = motorModel.VFMdata.data.Demag.id;
DemIq = motorModel.VFMdata.data.Demag.iq;
DemMS = motorModel.VFMdata.data.Demag.MS;
RemId = motorModel.VFMdata.data.Remag.id;
RemIq = motorModel.VFMdata.data.Remag.iq;
RemMS = motorModel.VFMdata.data.Remag.MS;
FlxMap = motorModel.FluxMap_dq;

if any(RemMS>2)
    RemMS = RemMS/100;
end 

if any(DemMS>2)
    DemMS = DemMS/100;
end 

if any(FlxMap.MS>2)
    FlxMap.MS = FlxMap.MS/100;
end 
 

% Flux interpolant function geneartion
index = 1:1:numel(FlxMap.Fd);
FlxId_vect = reshape(FlxMap.Id(index), [],1);
FlxIq_vect = reshape(FlxMap.Iq(index), [],1);
FlxMS_vect = reshape(FlxMap.MS(index), [],1);
FlxFd_vect = reshape(FlxMap.Fd(index), [],1);
FlxFq_vect = reshape(FlxMap.Fq(index), [],1);

FlxFuncFd = scatteredInterpolant(FlxId_vect, FlxIq_vect, FlxMS_vect, FlxFd_vect,'linear','linear');
FlxFuncFq = scatteredInterpolant(FlxId_vect, FlxIq_vect, FlxMS_vect, FlxFq_vect, 'linear', 'linear');


%% Fluxes for the demag and remag currents

for ii=1:numel(DemId)
    DemFd(ii) = FlxFuncFd(DemId(ii), DemIq(ii), DemMS(ii));
    DemFq(ii) = FlxFuncFq(DemId(ii), DemIq(ii), DemMS(ii));
end

MS_threshold = min(FlxMap.MS, [] ,'all');                      % Avoid NaN after the interpolation by limiting the minimum MS  
valid_indices = find(RemMS >= MS_threshold);
RemMS_valid = RemMS(valid_indices);     % Extract only valid data points
RemId_valid = RemId(valid_indices);
RemIq_valid = RemIq(valid_indices);
RemFd = zeros(size(RemMS_valid));       % Pre-allocation to have the same size of RemId_valid for the torque computation
RemFq = zeros(size(RemMS_valid));
for jj=1:numel(RemIq_valid)
    RemFd(jj) = FlxFuncFd(RemId_valid(jj), RemIq_valid(jj), RemMS_valid(jj));
    RemFq(jj) = FlxFuncFq(RemId_valid(jj), RemIq_valid(jj), RemMS_valid(jj));
end

DemFd = reshape(DemFd, size(DemId));  % Reshape to the origianl matrix size for torque computation
DemFq = reshape(DemFq, size(DemId));

F_Dem = hypot(DemFd, DemFq);                                                % Flux for all the id iq demag couples
T_Dem = 1.5*motorModel.data.p*(DemIq.*DemFd - DemId.*DemFq);                % Torque for all the id iq demag couples

F_Rem = hypot(RemFd, RemFq);                                                % Flux for all the id iq remag valid couples
T_Rem = 1.5*motorModel.data.p*(RemIq_valid.*RemFd-RemId_valid.*RemFq);      % Torque for all the id iq remag valid couples


%% DFVC Controller
% DFVC DEMAG MAP
iindex = find(DemFd >= 0);                     % DFVC cannot reach negative flux d value while keeping the desired torque 
F_DemDFVC = F_Dem(iindex);       
T_DemDFVC = T_Dem(iindex);
MS_DemDFVC = DemMS(iindex);

T_DemDomainDFVC = linspace(max(min(T_DemDFVC,[],'all'), 0), max(T_DemDFVC,[],'all'), n);
MS_DemDomainDFVC = linspace(min(MS_DemDFVC,[],'all'), max(MS_DemDFVC,[],'all'), n);
index = 1:1:numel(MS_DemDFVC);
DFVC_DemMS_vect = reshape(MS_DemDFVC(index), [], 1);                        % reshape for having column vector 
DFVC_DemT_vect = reshape(T_DemDFVC(index), [], 1);
DFVC_DemF_vect = reshape(F_DemDFVC(index), [], 1);
DFVC_DemFunc = scatteredInterpolant(DFVC_DemMS_vect, DFVC_DemT_vect,DFVC_DemF_vect, 'linear', 'none');

[Control.DFVC.Demag.MS, Control.DFVC.Demag.T] = meshgrid(MS_DemDomainDFVC, T_DemDomainDFVC);
Control.DFVC.Demag.F = DFVC_DemFunc(Control.DFVC.Demag.MS, Control.DFVC.Demag.T);

% 1D Torque limit LUT, DFVC DEMAG
Control.DFVC.Demag.Limit.MSlim = MS_DemDomainDFVC;
[~, cols] = size(Control.DFVC.Demag.F);
last_valid_idx = arrayfun(@(col) find(~isnan(Control.DFVC.Demag.F(:, col)), 1, 'last'), 1:cols, 'UniformOutput', false);
% Handle columns with all NaN by checking if result is empty
last_valid_idx(cellfun(@isempty, last_valid_idx)) = {1};
last_valid_idx = cell2mat(last_valid_idx);
Control.DFVC.Demag.Limit.Tlim = Control.DFVC.Demag.T(last_valid_idx,1);

% DFVC REMAG MAP
T_RemDomainDFVC = linspace(max(min(T_Rem,[],'all'), 0), max(T_Rem,[],'all'), n);
MS_RemDomainDFVC = linspace(min(FlxMap.MS,[],'all'), max(RemMS,[],'all'), n);
DFVC_RemFunc = scatteredInterpolant(RemMS_valid, T_Rem, F_Rem, 'linear', 'none');

[Control.DFVC.Remag.MS, Control.DFVC.Remag.T] = meshgrid(MS_RemDomainDFVC, T_RemDomainDFVC);
Control.DFVC.Remag.F = DFVC_RemFunc(Control.DFVC.Remag.MS, Control.DFVC.Remag.T);

% 1D Torque limit LUT, DFVC DEMAG
Control.DFVC.Remag.Limit.MSlim = MS_RemDomainDFVC;
[~, cols] = size(Control.DFVC.Remag.F);
last_valid_idx = arrayfun(@(col) find(~isnan(Control.DFVC.Remag.F(:, col)), 1, 'last'), 1:cols, 'UniformOutput', false);
% Handle columns with all NaN by checking if result is empty
last_valid_idx(cellfun(@isempty, last_valid_idx)) = {1};
last_valid_idx = cell2mat(last_valid_idx);
Control.DFVC.Remag.Limit.Tlim = Control.DFVC.Remag.T(last_valid_idx,1);


%% FOC Controller
% FOC DEMAG MAP
T_DemDomainIdT = linspace(max(min(T_Dem, [], 'all'), 0),max(T_Dem, [], 'all'),n);
MS_DemDomainIdT = linspace(min(DemMS, [], 'all'),max(DemMS, [], 'all'),n);

index = 1:1:numel(DemMS);
FOC_DemMS_vect = reshape(DemMS(index), [], 1);
FOC_DemT_vect = reshape(T_Dem(index), [], 1);
FOC_DemId_vect = reshape(DemId(index), [], 1);
FOC_DemIq_vect = reshape(DemIq(index), [], 1);
FOC_DemFuncId = scatteredInterpolant(FOC_DemMS_vect, FOC_DemT_vect, FOC_DemId_vect, 'linear', 'none');
FOC_DemFuncIq = scatteredInterpolant(FOC_DemMS_vect, FOC_DemT_vect, FOC_DemIq_vect, 'linear', 'none');


[Control.FOC.Demag.MS, Control.FOC.Demag.T] = meshgrid(MS_DemDomainIdT, T_DemDomainIdT);
Control.FOC.Demag.id = FOC_DemFuncId(Control.FOC.Demag.MS, Control.FOC.Demag.T);
Control.FOC.Demag.iq = FOC_DemFuncIq(Control.FOC.Demag.MS, Control.FOC.Demag.T);

% 1D Torque limit LUT, FOC DEMAG
Control.FOC.Demag.Limit.MSlim = MS_DemDomainIdT;
[~, cols] = size(Control.FOC.Demag.id);
last_valid_idx = arrayfun(@(col) find(~isnan(Control.FOC.Demag.id(:, col)), 1, 'last'), 1:cols, 'UniformOutput', false);
% Handle columns with all NaN by checking if result is empty
last_valid_idx(cellfun(@isempty, last_valid_idx)) = {1};
last_valid_idx = cell2mat(last_valid_idx);
Control.FOC.Demag.Limit.Tlim = Control.FOC.Demag.T(last_valid_idx,1);

% FOC REMAG MAP 
T_RemDomainIdT = linspace(max(min(T_Rem, [], 1),0),max(T_Rem, [], 1),n);
MS_RemDomainIdT = linspace(min(FlxMap.MS,[],'all'), max(RemMS,[],'all'), n);

FOC_RemFuncId = scatteredInterpolant(RemMS_valid, T_Rem, RemId_valid, 'linear', 'none');
FOC_RemFuncIq = scatteredInterpolant(RemMS_valid, T_Rem, RemIq_valid, 'linear', 'none');

[Control.FOC.Remag.MS, Control.FOC.Remag.T] = meshgrid(MS_RemDomainIdT, T_RemDomainIdT);
Control.FOC.Remag.id = FOC_RemFuncId(Control.FOC.Remag.MS, Control.FOC.Remag.T);
Control.FOC.Remag.iq = FOC_RemFuncIq(Control.FOC.Remag.MS, Control.FOC.Remag.T);

% 1D Torque limit LUT, FOC REMAG
Control.FOC.Remag.Limit.MSlim = MS_RemDomainIdT;
[~, cols] = size(Control.FOC.Remag.id);
last_valid_idx = arrayfun(@(col) find(~isnan(Control.FOC.Remag.id(:, col)), 1, 'last'), 1:cols, 'UniformOutput', false);  % store the last not NaN value for each column 
last_valid_idx(cellfun(@isempty, last_valid_idx)) = {1};                    % Handle columns with all NaN and set the index to 1
last_valid_idx = cell2mat(last_valid_idx);
Control.FOC.Remag.Limit.Tlim = Control.FOC.Remag.T(last_valid_idx,1);

% % Output Variables
% 
% motorModel.VFMdata.Control = Control;