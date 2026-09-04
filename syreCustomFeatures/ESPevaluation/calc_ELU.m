% Copyright 2026
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
%    limitations under the Licens

function [total,partial,ELUkg] = calc_ELU(dataSet)

% Computation of the Environmental Load Units to assess the environmental
% impact of the active parts of the motor
% ref : https://www.ivl.se/download/18.7342a03f17582337c283720/1605276804178/EPS2015d-Including_climate_impacts_from_secondary_particles_Damage_costs.xlsx

% Stator and rotor cores
% FeSi iron laminations are composed of about
% - 96% Fe (0.810039525197803 ELU/kg)
% - 2.8-3.8% Si (0 ELU/kg)
% --> 0.7776 ELU/kg
% FeCo iron lamination are composed of about
% - 50% Fe (0.810039525197803 ELU/kg)
% - 50% Co (178.776470588235 ELU/kg)
% --> 89.7933 ELU/kg
%
% Winding
% copper: 90.948 ELU/kg (insulation not considered)
%
% PMs
% - ferrite: iron oxide (rust) and Barium or Strontium oxide. No data available on the %. Hypothesis: full iron --> (0.810039525197803% ELU/kg)
% - NdFeB: typical:
%   - 32% Nd (139.923076923077 ELU/kg) / Pr (512.394366197183 ELU/kg) (facciamo 16% + 16%)
%   -  2% Dy (1039.42857142857 ELU/kg)
%   -  1% B  (0.05 ELU/kg)
%   -  1% Al (0.35 ELU/kg)
%   -  1% Co (178.7764706 ELU/kg)
%   - 63% Fe (0.810039525197803 ELU/kg)
%   --> total: 170.4615 ELU/kg
% - SmCo (series 2:17)
%   - 25% Sm (808.4444444 ELU/kg)
%   - 75% Co (178.776470588235 ELU/kg)
%   --> total 336.18 ELU/kg

ELUkg.FeSi         =   0.7776;
ELUkg.FeCo         =  89.7933;
ELUkg.Cu           =  90.948;
ELUkg.Al           =   0.3473;
ELUkg.NdFeB        = 170.4615;
ELUkg.NdFeB_DyFree = 104.9216;
ELUkg.SmCo         = 336.1800;
ELUkg.ferrite      =   0.81;


xNames{1} = 'FeSi';
xNames{2} = 'FeCo';
xNames{3} = 'Cu';
xNames{4} = 'Al';
xNames{5} = 'NdFeB';
xNames{6} = 'Dy-free NdFeB';
xNames{7} = 'SmCo';
xNames{8} = 'Ferrite';

if nargin~=0

    if strcmp(dataSet.FluxBarrierMaterial,'Air')
        ELUkg.PM = 0;
    elseif ~isempty(strfind(dataSet.FluxBarrierMaterial,'N'))
        ELUkg.PM = ELUkg.NdFeB;
    elseif ~isempty(strfind(dataSet.FluxBarrierMaterial,'Recoma'))
        ELUkg.PM = ELUkg.SmCo;
    elseif ~isempty(strfind(dataSet.FluxBarrierMaterial,'SmCo'))
        ELUkg.PM = ELUkg.SmCo;
    else
        ELUkg.PM = ELUkg.ferrite;
    end
    
    partial.PM = dataSet.MassMagnet*ELUkg.PM;
    partial.Cu = dataSet.MassWinding*ELUkg.Cu;
    partial.Fe = (dataSet.MassStatorIron+dataSet.MassRotorIron)*ELUkg.FeSi;
    
    total = partial.PM+partial.Cu+partial.Fe;
else
    partial = [];
    total = [];
end





