% Copyright 2020
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

function dataSet = ExportThermalParameters_MCAD(dataSet)

% == NEW SCRIPT FOR SAVING A MACHINE FROM GUI MANUAL ENTRIES ==
% ==== TAKE THE DIR AND THE STRUCT DATASET ===============================

filemot = strrep(dataSet.currentfilename,'.mat','.mot');
tmp = exist([dataSet.currentpathname filemot],'file');

if tmp == 2
    mcad=actxserver('MotorCAD.AppAutomation');
    
    if nargin < 2
        %load Syr-e and MCAD model
        file_mot=[dataSet.currentfilename(1:(end-4)) '.mot'];
        invoke(mcad,'LoadFromFile',[dataSet.currentpathname file_mot]);
    end
    
    %Setting Cooling System
    switch  dataSet.HousingType
        case 'Axial fins (Servo)'
            invoke(mcad,'SetVariable','HousingType', 9 );
            tmp=dataSet.InletTemperature;
            tmp=num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','Ambient_Temperature',tmp);
        case 'Water Jacket (Axial)'
            invoke(mcad,'SetVariable','HousingType', 11 );
            invoke(mcad,'SetVariable','Housing_Water_Jacket', 1 );
            tmp=dataSet.InletTemperature;
            tmp=num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','WJ_Fluid_Inlet_Temperature',tmp);
            
            invoke(mcad,'SetVariable','WJ_Channel-Lam',2);
            invoke(mcad,'SetVariable','WJ_Channel_Spacing',2);
            invoke(mcad,'SetVariable','WJ_Channel_Height',5);
            
            
            tmp = dataSet.FlowRate/(5.988*10^4);
            tmp = num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','WJ_Fluid_Volume_Flow_Rate', tmp);   %l/min
            
            %%Water-Glycol properties 50-50
            tmp = 0.399;
            tmp = num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','WJ_Fluid_Thermal_Conductivity', tmp);
            tmp =1065;
            tmp = num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','WJ_Fluid_Density', tmp);
            tmp = 3364;
            tmp = num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','WJ_Fluid_Specific_Heat', tmp);
            tmp = 2.25*10^-6;
            tmp = num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','WJ_Fluid_Kinematic_Viscosity', tmp);
            tmp = 0.002396;
            tmp = num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','WJ_Fluid_Dynamic_Viscosity', tmp);
            tmp = 20.2;
            tmp = num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','WJ_Fluid_Prandlt_Number', tmp);
            
        case 'Water Jacket (Spiral)'
            invoke(mcad,'SetVariable','HousingType', 12 );
            invoke(mcad,'SetVariable','Housing_Water_Jacket', 1 );
            tmp=dataSet.InletTemperature;
            tmp=num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','WJ_Fluid_Inlet_Temperature',tmp);
            
            invoke(mcad,'SetVariable','WJ_Channel-Lam',5);
            invoke(mcad,'SetVariable','WJ_Channel_Height',5);
            
            
            
            tmp = dataSet.FlowRate/(5.988*10^4);
            tmp = num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','WJ_Fluid_Volume_Flow_Rate', tmp);   %l/min
            
            if strcmp(dataSet.Fluid,'W/G 50/50')
                %%Water-Glycol 50 50 properties
                tmp = 0.399;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Thermal_Conductivity', tmp);
                tmp =1065;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Density', tmp);
                tmp = 3364;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Specific_Heat', tmp);
                tmp = 2.25*10^-6;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Kinematic_Viscosity', tmp);
                tmp = 0.002396;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Dynamic_Viscosity', tmp);
                tmp = 20.2;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Prandlt_Number', tmp);
            end
            
            if strcmp(dataSet.Fluid,'W/G 60/40')
                %%Water-Glycol 60 40 properties
                tmp = 0.361;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Thermal_Conductivity', tmp);
                tmp =1079;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Density', tmp);
                tmp = 3262;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Specific_Heat', tmp);
                tmp = 2.79*10^-6;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Kinematic_Viscosity', tmp);
                tmp = 0.00301;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Dynamic_Viscosity', tmp);
                tmp = 27.2;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Prandlt_Number', tmp);
            end
            
            if strcmp(dataSet.Fluid,'Water')
                %%Water-Glycol 60 40 properties
                tmp = 0.6233;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Thermal_Conductivity', tmp);
                tmp =992.3;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Density', tmp);
                tmp = 4178;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Specific_Heat', tmp);
                tmp = 6.59*10^-7;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Kinematic_Viscosity', tmp);
                tmp = 0.0006539;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Dynamic_Viscosity', tmp);
                tmp = 4.383;
                tmp = num2str(tmp);
                tmp(tmp=='.')=',';
                invoke(mcad,'SetVariable','WJ_Fluid_Prandlt_Number', tmp);
            end
            
            
        case 'None'
            invoke(mcad,'SetVariable','HousingType', 13 );
            tmp=dataSet.InletTemperature;
            tmp=num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','Ambient_Temperature',tmp);
    end
    
    %Setting calculation type
    invoke(mcad,'SetVariable','Transient_Calculation_Type', 0); 

    %Show thermal context
    invoke(mcad,'ShowThermalContext');
    
    %save MCAD model
    %invoke(mcad,'SaveToFile',[pathname filename]);
    invoke(mcad,'SaveToFile',[dataSet.currentpathname file_mot]);
    
    %Display save
    disp('Motor-CAD Thermal Model file saved in:')
    disp([dataSet.currentpathname file_mot])
    disp(' ')
    
    %Close Motor-CAD
    %     invoke(mcad,'Quit');
else
    error('Error: File .mot not found...')
end

end


