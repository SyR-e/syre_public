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

function dataSet = Eval_ThermalLimit0rpm(dataSet,SimType)

mcad = actxserver('MotorCAD.AppAutomation');

path      = dataSet.currentpathname;
file      = dataSet.currentfilename;
TLim_cu   = dataSet.TempCuLimit;
I_th_calc = dataSet.RatedCurrent;
Time      = dataSet.TransTime;
R_0       = dataSet.Rs;
T_ini     = dataSet.InitTemp;
T_amb     = dataSet.AmbTemp;

R = R_0 * (234.5+dataSet.TargetCopperTemp)/(234.5+TLim_cu);

tmp = load([path file]);
geo = tmp.geo;

% Load Motor-CAD model
filenameMot = [file(1:(end-3)) 'mot'];
[fail] = invoke(mcad,'LoadFromFile',[path filenameMot]);
if fail == 1
    error('The requested MotorCAD does not exist. Export the model first!')
end

tic

%% Simulation settings
invoke(mcad,'SetVariable','MagneticThermalCoupling',1);
invoke(mcad,'SetVariable','ProximityLossModel',0);
invoke(mcad,'SetVariable','MagneticSolver',1);

invoke(mcad,'SetVariable','TorqueCalculation','False');
invoke(mcad,'SetVariable','BackEMFCalculation','False');
invoke(mcad,'SetVariable','CoggingTorqueCalculation','False');
invoke(mcad,'SetVariable','TorqueSpeedCalculation','False');
invoke(mcad,'SetVariable','DemagnetizationCalc','False');
invoke(mcad,'SetVariable','ElectromagneticForcesCalc_Load','False');

invoke(mcad,'SetVariable','StatorCopperLossesVaryWithTemp', 'False');
invoke(mcad,'SetVariable','Ambient_Temperature',T_amb);

%% Steady
if strcmp(SimType,'Steady')
    invoke(mcad,'SetVariable','PeakCurrent',I_th_calc);
    speed_sim = 1;
    invoke(mcad,'SetVariable','ShaftSpeed',speed_sim);
    invoke(mcad,'DoMagneticCalculation');
    invoke(mcad,'ClearMessages');
    DC_loss = 3/2*R*I_th_calc^2;
    
    invoke(mcad,'SetVariable','ThermalCalcType', 0);
    
    ii   = 0;
    flag = 1;
    
    while flag==1
        ii = ii + 1;
        flag = 0;
        invoke(mcad,'SetVariable','Armature_Copper_Loss_@Ref_Speed', DC_loss);
        invoke(mcad,'DoSteadyStateAnalysis');
        [~, WindingTemperature_Max] = invoke(mcad,'GetVariable','T_[Winding_Max]');
        [~, WindingTemperature_Ave] = invoke(mcad,'GetVariable','T_[Winding_Average]');
        
        if (WindingTemperature_Max<TLim_cu-5 || WindingTemperature_Max>TLim_cu+1)
            WindingTemperature_Max(WindingTemperature_Max>2*TLim_cu) = TLim_cu(WindingTemperature_Max>2*TLim_cu)*1.5;
            DC_loss = (DC_loss + (TLim_cu-WindingTemperature_Max)/TLim_cu * DC_loss);
            flag    = 1;
            I_th0   = sqrt(DC_loss/3/R*2);
            disp(['Iteration ' num2str(ii) ': ' num2str(round(I_th0)) ' A - ' num2str(round(WindingTemperature_Max)) '°C copper temperature'])
        end
        
    end
    
    if flag == 0
        R = R_0 * (234.5+dataSet.TargetCopperTemp)/(234.5+WindingTemperature_Ave);
        I_th0       = sqrt(DC_loss/3/R*2);
        [~,cu_cond] = invoke(mcad,'GetVariable','ArmatureConductorCSA');
        [~,n_cond]  = invoke(mcad,'GetVariable','ConductorsPerSlot');
        [~,n_par]   = invoke(mcad,'GetVariable','ParallelPaths');
        [~,turn]    = invoke(mcad,'GetVariable','MagTurnsConductor');
        [~,lay]     = invoke(mcad,'GetVariable','WindingLayers');
        cu_area     = cu_cond*double(n_cond);
        I_den       = I_th0*(double(n_par)/double(turn)*double(lay))/(cu_area*sqrt(2));
        
        [~, WJ_out]  = invoke(mcad,'GetVariable','WJ_Fluid_Outlet_Temp_[Active]');
        [~, WJ_in]   = invoke(mcad,'GetVariable','WJ_Fluid_Inlet_Temperature');
        deltaWJ      = WJ_out - WJ_in;
        
        [~, deltaHousingTemperature_Max] = invoke(mcad,'GetVariable','dT_[Housing_-_Ambient]');
        HousingTemperature_Max = 80 + deltaHousingTemperature_Max;
        
    end
    
    kj_eval  = DC_loss/(2*pi*geo.R*geo.l*10^-6);
    
    disp(['Continuous limit found at ' num2str(round(I_th0)) ' A - ' num2str(round(WindingTemperature_Max)) '°C copper temperature after ' num2str(round(ii)) ' iterations'])
    
    %Save output
    clock_1 = fix(clock);
    time = join((string(clock_1(1:end-1))),"");
    resFolder = [path file(1:(end-4)) '_results\FEA results\Thermal\Steady_Cu' num2str(TLim_cu) 'C_' num2str(time)];
    mkdir(resFolder);
    
    dataSet.SimIth0 = round(I_th0);
end

%% Transient
if strcmp(SimType,'Transient')
    
    invoke(mcad,'SetVariable','InitialTransientTemperatureOption',3);
    invoke(mcad,'SetVariable','Initial_Machine_Temperature',T_ini);
    
    I_trans_calc = 2*I_th_calc;
    invoke(mcad,'SetVariable','PeakCurrent',I_trans_calc);
    speed_sim = 1;
    invoke(mcad,'SetVariable','ShaftSpeed',speed_sim);
    invoke(mcad,'DoMagneticCalculation');
    invoke(mcad,'ClearMessages');
    DC_loss = 3/2*R*I_th_calc^2;
    
    invoke(mcad,'SetVariable','ThermalCalcType', 1);
    
    invoke(mcad,'SetVariable','Transient_Time_Period',Time);
    invoke(mcad,'SetVariable','Number_Transient_Points',4);
    ii   = 1;
    flag = 1;
    
    while flag == 1
        flag = 0;
        invoke(mcad,'SetVariable','Armature_Copper_Loss_@Ref_Speed', DC_loss);
        invoke(mcad,'DoTransientAnalysis');
        [~, WindingTemperature_Max] = invoke(mcad,'GetVariable','T_[Winding_Max]');
        [~, WindingTemperature_Ave] = invoke(mcad,'GetVariable','T_[Winding_Average]');
        if (WindingTemperature_Max<TLim_cu-5 || WindingTemperature_Max>TLim_cu+1)
            WindingTemperature_Max(WindingTemperature_Max>2*TLim_cu) = TLim_cu(WindingTemperature_Max>2*TLim_cu)*1.5;
            DC_loss = (DC_loss + (TLim_cu-WindingTemperature_Max)/TLim_cu * DC_loss);
            I_th0   = sqrt(DC_loss/3/R*2);
            disp(['Iteration ' num2str(ii) ': ' num2str(round(I_th0)) ' A - ' num2str(round(WindingTemperature_Max)) '°C copper temperature'])
            
            flag = 1;
        end
        ii = ii + 1;
    end
    
    if flag == 0
        R = R_0 * (234.5+dataSet.TargetCopperTemp)/(234.5+WindingTemperature_Ave);
        I_th0       = sqrt(DC_loss/3/R*2);
        [~,cu_cond] = invoke(mcad,'GetVariable','ArmatureConductorCSA');
        [~,n_cond]  = invoke(mcad,'GetVariable','ConductorsPerSlot');
        [~,n_par]   = invoke(mcad,'GetVariable','ParallelPaths');
        [~,turn]    = invoke(mcad,'GetVariable','MagTurnsConductor');
        [~,lay]     = invoke(mcad,'GetVariable','WindingLayers');
        cu_area     = cu_cond*double(n_cond);
        I_den       = I_th0*(double(n_par)/double(turn)*double(lay))/(cu_area*sqrt(2));
        
        [~, WJ_out]  = invoke(mcad,'GetVariable','WJ_Fluid_Outlet_Temp_[Active]');
        [~, WJ_in]   = invoke(mcad,'GetVariable','WJ_Fluid_Inlet_Temperature');
        deltaWJ      = WJ_out - WJ_in;
    end
    
    
    kj_eval  = DC_loss/(2*pi*geo.R*geo.l*10^-6);
    
    disp(['Transient limit found at ' num2str(round(I_th0)) ' A - ' num2str(round(WindingTemperature_Max)) '°C copper temperature after ' num2str(ii) ' iterations'])
    
    %Save output
    clock_1 = fix(clock);
    time = join((string(clock_1(1:end-1))),"");
    resFolder = [path file(1:(end-4)) '_results\FEA results\Thermal\Transient_Cu' num2str(TLim_cu) 'C_' num2str(time)];
    mkdir(resFolder);
    
    dataSet.SimIthpk = round(I_th0);
end
toc
%% Output
out.Loss     = DC_loss;
out.kj       = kj_eval;
out.I_th     = I_th0;
out.I_den    = I_den;
out.deltaTWJ = deltaWJ;
out.TempWind = WindingTemperature_Max;

%% Update the Thermal parameter panel
if strcmp(SimType,'Steady')
    
    dataSet.EstimatedCopperTemp = round(WindingTemperature_Ave);
    dataSet.AdmiJouleLosses     = round(DC_loss);
    dataSet.ThermalLoadKj       = round(kj_eval);
    dataSet.TargetCopperTemp    = round(WindingTemperature_Ave);
    dataSet.HousingTemp         = round(HousingTemperature_Max);
    dataSet.RatedCurrent        = round(I_th0,2);

end

save([resFolder '\ThLim0rpm' SimType '.mat'], 'out');