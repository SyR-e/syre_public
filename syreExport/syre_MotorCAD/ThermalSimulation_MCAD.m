%    Copyright 2020
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

function dataSet = ThermalSimulation_MCAD(dataSet)

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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Setting Simulation Parameters %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % save output into individual folders
    FILENAME = [filemot(1:end-4),'_ThermalSimulation_',num2str(dataSet.TransientPeriod) 's_',num2str(dataSet.InletTemperature) 'C','_MCAD'];
    mkdir(dataSet.currentpathname,FILENAME);
    newDir=[dataSet.currentpathname,FILENAME,'\'];
    
    %Show thermal context
    invoke(mcad,'ShowThermalContext');
    
    %Check type of calculation
    if dataSet.TransientPeriod==inf
        invoke(mcad,'SetVariable','ThermalCalcType', 0);
        disp('Thermal Steady State Analysis in progress...')
        success=invoke(mcad,'DoSteadyStateAnalysis');
        if success==0
            disp('Thermal calculation successfully completed')
        else
            disp('Thermal calculation failed')
        end
        
        % Plot steadystate temperature
        [tmp, WindingTemperature_Min]=invoke(mcad,'GetVariable','T_[Winding_Min]');
        [tmp, WindingTemperature_Max]=invoke(mcad,'GetVariable','T_[Winding_Max]');
        [tmp, WindingTemperature_Average]=invoke(mcad,'GetVariable','T_[Winding_Average]');
        T=[WindingTemperature_Min,WindingTemperature_Max,WindingTemperature_Average];
        c=categorical({'Winding Min','Winding Max','Winding Average'});
        figure()
        figSetting()
        bar(c,T)
        title('Temperature in Steady-State Analysis')
        ylabel('T [Celsius]')
        savefig([newDir, 'SteadyState_temperature.fig']);
        outTherm.WindingTempMeanSteadyState=WindingTemperature_Average;
    else
        invoke(mcad,'SetVariable','ThermalCalcType', 1);
        disp('Thermal Transient Analysis in progress...')
        success=invoke(mcad,'DoTransientAnalysis');
        if success==0
            disp('Thermal calculation successfully completed')
        else
            disp('Thermal calculation failed')
        end
        
        %Save Transient Solution .csv
        tmp=[newDir 'TransientTemperatures.csv'];
        invoke(mcad,'SaveTransientTemperatures',tmp);
        
        %Set Number Points
        tmp=dataSet.TransientTimeStep;
        tmp=num2str(tmp);
        tmp(tmp=='.')=',';
        invoke(mcad,'SetVariable','Number_Transient_Points',tmp);
        
        %Save winding average
        nPoints = dataSet.TransientTimeStep;
        WindingTemp_Average_Transient = zeros(nPoints,1);
        Time = zeros(nPoints,1);
        for timestep=1: 1 : nPoints
            [success,x,y]=invoke(mcad,'GetTemperatureGraphPoint','Winding (Average)',timestep);
            if success ==0
                Time(timestep)=x;
                WindingTemp_Average_Transient(timestep)=y;
            end 
        end
        
        outTherm.timeTransient=Time;
        outTherm.WindingTempTransient=WindingTemp_Average_Transient;
        
        % Plot winding Average
        figure
        figSetting
        figSetting(16,10)
        plot(Time,WindingTemp_Average_Transient,'b*','DisplayName','Winding Average')
        xlabel('Time [s]')
        ylabel('Temperature [Celsius]')
        legend show
        grid on
        saveas(gcf,[newDir 'TransientTemperature'])
        print(gcf,[newDir 'TransientTemperature.png'],'-dpng','-r300')
    end
    
    %save MCAD model
    %invoke(mcad,'SaveToFile',[pathname filename]);
    invoke(mcad,'SaveToFile',[dataSet.currentpathname file_mot]);
    
    %Display save
    disp('Motor-CAD Thermal Simulation file saved in:')
    disp([dataSet.currentpathname file_mot])
    disp(' ')
    
    %Save variable
    save([newDir,FILENAME,'.mat'],'outTherm','dataSet');
    
    %Close Motor-CAD
    invoke(mcad,'Quit');
    
else
    error('Error: File .mot not found...')
end

end




