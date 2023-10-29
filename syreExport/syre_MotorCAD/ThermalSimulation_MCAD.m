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
%dataSet.th_eval_type = 'Steady State';
th_eval_type = dataSet.th_eval_type;

pathname = dataSet.currentpathname;
filename = dataSet.currentfilename;

filemot = strrep(filename,'.mat','.mot');
tmp = exist([pathname filemot],'file');


if tmp == 2

    mcad = actxserver('MotorCAD.AppAutomation');

    if nargin < 2
        %load SyR-e and MCAD model
        file_mot=[filename(1:(end-4)) '.mot'];
        invoke(mcad,'LoadFromFile',[pathname file_mot]);
    end

    %% Input Loss
    if dataSet.CustomLossMCADCheck
        prompt = {'Copper loss [W]:','Stator loss [W]:','Rotor loss [W]:','Magnet loss [W]:'};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'2000','1000','500','30'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);

        Pjs  = str2double(answer{1});
        Pfes = str2double(answer{2});
        Pfer = str2double(answer{3});
        Ppm  = str2double(answer{4});
    else
        FEAfolder = [pathname filename(1:(end-4)) '_results\FEA results\'];
        if ~exist(FEAfolder,'dir')
            FEAfolder = pathname;
        end

        [filename1,pathname1] = uigetfile([FEAfolder '\*.mat'],'Load loss data');
        data = load([pathname1 filename1]);
        I = abs(data.out.id+1i*data.out.iq);

        Pfes = data.out.Pfes_h + data.out.Pfes_c;
        Pfer = data.out.Pfer_h + data.out.Pfer_c;
        Ppm  = data.out.Ppm;

        clear data

        if dataSet.CustomCurrentEnable
            [outPcu] = skinEffect_evalCustomCurrent(dataSet);
            Pjs = outPcu.Ploss;
        else
            Rs   = dataSet.Rs;
            l    = dataSet.StackLength;
            lend = dataSet.EndWindingsLength;
            dataSet.SlotConductorFrequency   = dataSet.NumOfPolePairs*dataSet.EvalSpeed/60;
            dataSet.SlotConductorTemperature = dataSet.TargetCopperTemp;

            [outPcu] = skinEffect_eval(dataSet);
            Pjs = 3/2*Rs.*(outPcu.k*l/(lend+l)+lend/(lend+l))*I.^2;
        end
    end

    [~,wt]  = invoke(mcad,'GetVariable','Tooth_Width');
    [~,sd]  = invoke(mcad,'GetVariable','Slot_Depth');
    [~,rse] = invoke(mcad,'GetVariable','Stator_Lam_Dia');
    [~,rsi] = invoke(mcad,'GetVariable','Stator_Bore');
    wy = (rse/2-rsi/2-sd);
    ratio_pfet = wt/(wy+wt);
    Pfest = ratio_pfet*Pfes;
    Pfesy = (1-ratio_pfet)*Pfes;

    invoke(mcad,'SetVariable','Stator_Iron_Loss_@Ref_Speed_[Tooth]', Pfest);
    invoke(mcad,'SetVariable','Stator_Iron_Loss_@Ref_Speed_[Back_Iron]', Pfesy);
    invoke(mcad,'SetVariable','Magnet_Iron_Loss_@Ref_Speed', Ppm);
    invoke(mcad,'SetVariable','Armature_Copper_Loss_@Ref_Speed', Pjs);
    invoke(mcad,'SetVariable','Rotor_Iron_Loss_@Ref_Speed_[Embedded_Magnet_Pole]', Pfer/2);
    invoke(mcad,'SetVariable','Rotor_Iron_Loss_@Ref_Speed_[Back_Iron]', Pfer/2);
    invoke(mcad,'SetVariable','StatorCopperLossesVaryWithTemp', 'True');
    invoke(mcad,'SetVariable','StatorWindingTemperatureAtWhichPcuInput', dataSet.TargetCopperTemp);

    %% Setting inlet temperature
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

    invoke(mcad,'SetVariable','Ambient_Temperature',dataSet.AmbTemp);

    %     %Setting calculation type
    %     invoke(mcad,'SetVariable','Transient_Calculation_Type', 0);

    % save output into individual folders
    FILENAME = [filemot(1:end-4),'_ThermalSimulation_',num2str(dataSet.TransientPeriod) 's_',num2str(dataSet.InletTemperature) 'C','_MCAD'];
    mkdir(dataSet.currentpathname,FILENAME);
    newDir=[dataSet.currentpathname,FILENAME,'\'];

    switch th_eval_type
        case char('Transient')

            %Whole machine at specified temperature
            invoke(mcad,'SetVariable','InitialTransientTemperatureOption',3);
            tmp = dataSet.InitTemp;
            tmp = num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','Initial_Machine_Temperature',tmp);

            tmp=dataSet.TransientPeriod;
            tmp=num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','Transient_Time_Period',tmp);


            %Setting number of Point
            tmp=dataSet.TransientTimeStep;
            tmp=num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','Number_Transient_Points',tmp);


            invoke(mcad,'SetVariable','ThermalCalcType', 1);
            disp('Thermal Transient Analysis in progress...')
            success=invoke(mcad,'DoTransientAnalysis');
            if success==0
                disp('Thermal calculation successfully completed')
            else
                disp('Thermal calculation failed')
            end

            %Save Transient Solution .csv
            tmp = [newDir 'TransientTemperatures.csv'];
            invoke(mcad,'SaveTransientTemperatures',tmp);

            %Set Number Points
            tmp=dataSet.TransientTimeStep;
            tmp=num2str(tmp);
            tmp(tmp=='.')=',';
            invoke(mcad,'SetVariable','Number_Transient_Points',tmp);

            %Save winding average
            nPoints = dataSet.TransientTimeStep;
            WindingTemp_Average_Transient = zeros(nPoints,1);
            WindingTemp_Coolspot_Transient = zeros(nPoints,1);
            WindingTemp_Hotspot_Transient = zeros(nPoints,1);
            MagnetTemp_Transient = zeros(nPoints,1);
            Time = zeros(nPoints,1);
            for timestep = 1: 1 : nPoints
                [success,x1,y1] = invoke(mcad,'GetTemperatureGraphPoint','Winding (Avg)',timestep);
                [~,~,y2] = invoke(mcad,'GetTemperatureGraphPoint','Winding (Coolspot)',timestep);
                [~,~,y3] = invoke(mcad,'GetTemperatureGraphPoint','Winding (Hotspot)',timestep);
                [~,~,y4] = invoke(mcad,'GetTemperatureGraphPoint','Magnet [R]',timestep);
                if success == 0
                    Time(timestep) = x1;
                    WindingTemp_Average_Transient(timestep)  = y1;
                    WindingTemp_Coolspot_Transient(timestep) = y2;
                    WindingTemp_Hotspot_Transient(timestep)  = y3;
                    MagnetTemp_Transient(timestep)           = y4;
                end
            end

            Time = [0 ;Time];
            WindingTemp_Average_Transient  = [ dataSet.InitTemp; WindingTemp_Average_Transient];
            WindingTemp_Coolspot_Transient = [ dataSet.InitTemp; WindingTemp_Coolspot_Transient];
            WindingTemp_Hotspot_Transient  = [ dataSet.InitTemp; WindingTemp_Hotspot_Transient];
            MagnetTemp_Transient           =  [ dataSet.InitTemp; MagnetTemp_Transient];

            outTherm.timeTransient                 = Time;
            outTherm.WindingAvgTempTransient       = WindingTemp_Average_Transient;
            outTherm.WindingTemp_Hotspot_Transient = WindingTemp_Hotspot_Transient;
            outTherm.MagnetTemp_Transient          = MagnetTemp_Transient;

            % Plot winding Average
            figure
            figSetting
            figSetting(16,10,10)
            T_limit = dataSet.TargetCopperTemp*ones(length(Time),1);
            plot(Time,T_limit,':r','LineWidth',0.5,'DisplayName','Target Copper Temperature')
            plot(Time,WindingTemp_Average_Transient,'b','DisplayName','Winding Average')
            plot(Time,WindingTemp_Coolspot_Transient,'Color',[0 0.8 0],'DisplayName','Winding Coolspot')
            plot(Time,WindingTemp_Hotspot_Transient,'r','DisplayName','Winding Hotspot')
            plot(Time,MagnetTemp_Transient,'k','DisplayName','Magnet')
            xlabel('Time [s]')
            xlim([0 Time(end)])
            ylabel('$\Theta$ [$^\circ$C]')
            title('Temperature - Transient')
            legend('Location','southeast')
            saveas(gcf,[newDir 'TransientTemperature'])
            print(gcf,[newDir 'TransientTemperature.png'],'-dpng','-r300')

        case char('Steady State')

            invoke(mcad,'SetVariable','ThermalCalcType', 0);
            disp('Thermal Steady State Analysis in progress...')

            invoke(mcad,'SetVariable','MagneticThermalCoupling',2);
            success = invoke(mcad,'DoSteadyStateAnalysis');
            if success==0
                disp('Thermal calculation successfully completed')
            else
                disp('Thermal calculation failed')
            end

            % Plot steadystate temperature
            [~, WindingTemperature_Min] = invoke(mcad,'GetVariable','T_[Winding_Min]');
            [~, WindingTemperature_Max] = invoke(mcad,'GetVariable','T_[Winding_Max]');
            [~, WindingTemperature_Average] = invoke(mcad,'GetVariable','T_[Winding_Average]');

            %             invoke(mcad,'DoMagneticCalculation');
            %             [~, MagnetTemperature] = invoke(mcad,'GetVariable','Magnet_Temperature');


            T = [WindingTemperature_Min,WindingTemperature_Max,WindingTemperature_Average];
            c = categorical({'Min','Max','Average'});
            figure()
            figSetting(8,10,10)
            bar(c(1),T(1),0.4,'FaceColor',[0 0.8 0]);
            bar(c(2),T(2),0.4,'r');
            bar(c(3),T(3),0.4,'b');
            title('Winding Temperature - Steady-State')
            ylabel('$\Theta_{Cu}$ [$^\circ$C]')
            savefig([newDir, 'SteadyState_temperature.fig']);

            %             c1 = categorical({'Magnet'});
            %             figure()
            %             figSetting(8,10,10)
            %             bar(c1,MagnetTemperature,0.4,'r');
            %             title('Magnet Temperature - Steady-State')
            %             ylabel('$\Theta_{PM}$ [$^\circ$C]')


            outTherm.WindingTempMeanSteadyState = WindingTemperature_Average;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Show thermal context
    invoke(mcad,'ShowThermalContext');

    %save MCAD model
    invoke(mcad,'SaveToFile',[dataSet.currentpathname file_mot]);

    %Display save
    disp('Motor-CAD Thermal Simulation file saved in:')
    disp([dataSet.currentpathname file_mot])
    disp(' ')

    %Save variable
    save([newDir,FILENAME,'.mat'],'outTherm','dataSet');

    %Close Motor-CAD
    %invoke(mcad,'Quit');

else
    error('Error: File Motor-CAD (.mot) not found')
end

end