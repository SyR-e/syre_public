% Copyright 2014
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

function BarName = defineBlockNames(BarCenter,geo)

% define names of rotor barriers (not used in FEMM, needed for
% compatibility with other finite-element codes)

%% Main revision
% 01/08/2018 (Simone Ferrari): now BarName is write starting from
% BarCenter. This method is more robust respect the previous one and
% doesn't need code changes if the geometry is changed.

indAir=1;
indPMs=1;
indBar=1;
ind=1;

rotType=geo.RotType;

% Material codes: refer to BLKLABELS.materials positions
materialCodes;
% codMatFe    = 5;
% codMatBar   = 6;
% codMatShaft = 7;
% codMatPM    = 6;
% codMatAir   = 1;

for ii=1:length(BarCenter(:,3))
    if (BarCenter(ii,3)==codMatAirRot)
        BarName{ind}=['Air_Bar_' int2str(indAir)];
        ind=ind+1;
        indAir=indAir+1;
    elseif BarCenter(ii,3)==codMatBar
        BarName{ind}=['Magnet_Bar_' int2str(indPMs)];
        ind=ind+1;
        indPMs=indPMs+1;
    elseif BarCenter(ii,3)==codMatCuRot
        BarName{ind}=['Bar_' int2str(indBar)];
        ind=ind+1;
        indBar=indBar+1;
    end
end

