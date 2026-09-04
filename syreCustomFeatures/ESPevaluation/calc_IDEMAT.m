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
%    limitations under the License

function [CL,OL,ECOCOSTkg] = calc_IDEMAT(dataSet)

% Computation of the Eco Cost according to IDEMAT DataBase
% ref : https://www.ecocostsvalue.com/data-tools-books

%% Materials list OpenLoop
% Laminations
ECOCOSTkg.OL.FeSi               = 0.633808513043164; % 0.97Fe+0.027Si+0.003Mn
ECOCOSTkg.OL.FeCo               = 28.0903819696835;  % 0.49Fe+0.49Co+0.02V
% Conductors
ECOCOSTkg.OL.Cu                 = 5.219829219;
ECOCOSTkg.OL.Al                 = 1.91350744506325;
% Permanent Magnets
ECOCOSTkg.OL.PM_Fer             = 1.4670649623861;
ECOCOSTkg.OL.PM_NdFeB           = 83.3702266280105;
ECOCOSTkg.OL.PM_NdFeB_DyFree    = 60.6570706652032;
ECOCOSTkg.OL.PM_SmCo            = 89.338597567697;
ECOCOSTkg.OL.PM_FeN             = 1.58876414388147;

%% Materials list ClosedLoop (recycling)
% Laminations
ECOCOSTkg.CL.FeSi               = 0.163587840434164; % 0.97Fe+0.027Si+0.003Mn
ECOCOSTkg.CL.FeCo               = 8.03500769848352;  % 0.49Fe+0.49Co+0.02V
% Conductors        
ECOCOSTkg.CL.Cu                 = 0.0998292190000001;
ECOCOSTkg.CL.Al                 = 0.13350744506325;
% Permanent Magnets     
ECOCOSTkg.CL.PM_Fer             = 0.4003094623861;
ECOCOSTkg.CL.PM_NdFeB           = 3.925839061742;
ECOCOSTkg.CL.PM_NdFeB_DyFree    = 2.43815176446715;
ECOCOSTkg.CL.PM_SmCo            = 21.067195465529;
ECOCOSTkg.CL.PM_FeN             = 1.36525952188147;


%% computation from dataSet

if nargin~=0
    
    if strcmp(dataSet.FluxBarrierMaterial,'Air')
        OL.PM = 0;
        CL.PM = 0;
    elseif ~isempty(strfind(dataSet.FluxBarrierMaterial,'N'))
        OL.PM = ECOCOSTkg.OL.PM_NdFeB;
        CL.PM = ECOCOSTkg.CL.PM_NdFeB;
    elseif ~isempty(strfind(dataSet.FluxBarrierMaterial,'Recoma'))
        OL.PM = ECOCOSTkg.OL.PM_SmCo;
        CL.PM = ECOCOSTkg.CL.PM_SmCo;
    elseif ~isempty(strfind(dataSet.FluxBarrierMaterial,'SmCo'))
        OL.PM = ECOCOSTkg.OL.PM_SmCo;
        CL.PM = ECOCOSTkg.CL.PM_SmCo;
    else
        OL.PM = ECOCOSTkg.OL.PM_Fer;
        CL.PM = ECOCOSTkg.CL.PM_Fer;
    end
    
    OL.PM = OL.PM*dataSet.MassMagnet;
    OL.Cu = ECOCOSTkg.OL.Cu*dataSet.MassWinding;
    OL.Fe = ECOCOSTkg.OL.FeSi*(dataSet.MassStatorIron+dataSet.MassRotorIron);

    CL.PM = CL.PM*dataSet.MassMagnet;
    CL.Cu = ECOCOSTkg.CL.Cu*dataSet.MassWinding;
    CL.Fe = ECOCOSTkg.CL.FeSi*(dataSet.MassStatorIron+dataSet.MassRotorIron);
    
    OL.tot = OL.PM+OL.Cu+OL.Fe;
    CL.tot = CL.PM+CL.Cu+CL.Fe;
else
    OL = [];
    CL = [];
end






