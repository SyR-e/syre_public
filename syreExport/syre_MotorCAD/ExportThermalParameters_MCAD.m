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
    
    %Setting type Housing
    switch  dataSet.HousingType
        case 'Axial fins (Servo)'
            invoke(mcad,'SetVariable','HousingType', 9 );
        case 'Water Jacket (Axial)'
            invoke(mcad,'SetVariable','HousingType', 11);
        case 'Water Jacket (Spiral)'
            invoke(mcad,'SetVariable','HousingType', 12);
        case 'None'
            invoke(mcad,'SetVariable','HousingType', 13);
    end
    
    % Setting inlet temperature
    switch  dataSet.HousingType
        case 'Axial fins (Servo)'
            tmp=dataSet.InletTemperature;
            tmp=num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','Ambient_Temperature',tmp);
        case 'Water Jacket (Axial)'
            tmp=dataSet.InletTemperature;
            tmp=num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','WJ_Fluid_Inlet_Temperature',tmp);
        case 'Water Jacket (Spiral)'
            tmp=dataSet.InletTemperature;
            tmp=num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','WJ_Fluid_Inlet_Temperature',tmp);
        case 'None'
            tmp=dataSet.InletTemperature;
            tmp=num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','Ambient_Temperature',tmp);
    end
    
    %Setting calculation type
    invoke(mcad,'SetVariable','Transient_Calculation_Type', 0);
    
    %Setting transient period
    if dataSet.TransientPeriod ~= inf
        tmp=dataSet.TransientPeriod;
        tmp=num2str(tmp);
        tmp(tmp=='.')=',';
        invoke(mcad,'SetVariable','Transient_Time_Period',tmp);
    %Setting number of Point    
        tmp=dataSet.TransientTimeStep;
        tmp=num2str(tmp);
        tmp(tmp=='.')=',';
        invoke(mcad,'SetVariable','Number_Transient_Points',tmp);
    end
    
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
    invoke(mcad,'Quit');
else
    error('Error: File .mot not found...')
end

end


