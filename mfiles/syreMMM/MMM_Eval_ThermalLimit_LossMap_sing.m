% Copyright 2021
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

function [out] = MMM_Eval_ThermalLimit_LossMap_sing(mcad,loss,TCu_lim,Tpm_lim,speed)

invoke(mcad,'ShowThermalContext');
MagnetTemperature = nan;
invoke(mcad,'SetVariable','ShaftSpeed',speed);

Pfe_s = loss(:,2);
Pfe_s = Pfe_s(isfinite(Pfe_s));
index = floor(length(Pfe_s)/2);
Pfe_r = loss(:,3);

[~,wt]  = invoke(mcad,'GetVariable','Tooth_Width');
[~,sd]  = invoke(mcad,'GetVariable','Slot_Depth');
[~,rse] = invoke(mcad,'GetVariable','Stator_Lam_Dia');
[~,rsi] = invoke(mcad,'GetVariable','Stator_Bore');

wy = (rse/2-rsi/2-sd);
ratio_pfe = wt/wy;
%%%%
% ratio_pfe =0.5; 
%%%%
flag   = 1;
flagPM = 0;

while flag==1
    
    flag = 0;
    Pfe_init  = loss(index,2);
    Pfet_init = ratio_pfe*Pfe_init;
    Pfey_init = (1-ratio_pfe)*Pfe_init;
    Pfer_init = loss(index,3);
    Ppm_init = loss(index,4);
    Pjs_init = loss(index,5);
    invoke(mcad,'SetVariable','Stator_Iron_Loss_@Ref_Speed_[Tooth]', Pfet_init);
    invoke(mcad,'SetVariable','Stator_Iron_Loss_@Ref_Speed_[Back_Iron]', Pfey_init);
    invoke(mcad,'SetVariable','Magnet_Iron_Loss_@Ref_Speed', Ppm_init);
    invoke(mcad,'SetVariable','Armature_Copper_Loss_@Ref_Speed', Pjs_init);
    invoke(mcad,'SetVariable','Rotor_Iron_Loss_@Ref_Speed_[Embedded_Magnet_Pole]', Pfer_init/2);
    invoke(mcad,'SetVariable','Rotor_Iron_Loss_@Ref_Speed_[Back_Iron]', Pfer_init/2);
    
    invoke(mcad,'DoSteadyStateAnalysis');
    invoke(mcad,'SetVariable','MagneticThermalCoupling',2);
    invoke(mcad,'DoMagneticCalculation');
    [~, WindingTemperature_Max] = invoke(mcad,'GetVariable','T_[Winding_Max]');
    [~, MagnetTemperature]      = invoke(mcad,'GetVariable','Magnet_Temperature');
    
    if WindingTemperature_Max<TCu_lim-5 && flagPM == 0
        index = index+1;
        flag = 1;
    elseif WindingTemperature_Max>TCu_lim+1 && flagPM == 0
        index = index-1;
        flag = 1;
    end
    
    if MagnetTemperature>Tpm_lim && flag == 0
        index = index-1;
        flag = 1;
        flagPM = 1;
    end
    
    if index==1
        Loss_Tot               = nan;
        deltaWJ                = nan;
        WindingTemperature_Max = nan;
        
        disp(['No continous point feasible at ' num2str(speed) ' rpm']);
        break
    end
end

if flag==0
    [~, WJ_out]  = invoke(mcad,'GetVariable','WJ_Fluid_Outlet_Temp_[Active]');
    [~, WJ_in]   = invoke(mcad,'GetVariable','WJ_Fluid_Inlet_Temperature');
    deltaWJ      = WJ_out - WJ_in;
    Loss_Tot     = loss(index,end);
end

%% Output
out.Ploss   = Loss_Tot;
out.deltaWJ = deltaWJ;
out.TempCu  = WindingTemperature_Max;
out.TempPM  = MagnetTemperature;
