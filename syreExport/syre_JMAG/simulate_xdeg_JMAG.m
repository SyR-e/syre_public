% Copyright 2024
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
% '------------------------------------------------------------------------
function [SOL] = simulate_xdeg_JMAG(geo,per,pathname,filename)

n3ph = geo.win.n3phase;

flagCharger = false;

if(~isnan(per.rotorPos))
    flagCharger = true;
end

custom_act = per.custom_act;

if(flagCharger)
    curr_pos = per.rotorPos;
    custom_Amp = per.custom_Amp(:,curr_pos);
    custom_Ph = per.custom_Ph(:,curr_pos);
    rotorPos = per.custom_rotorPos(curr_pos);
    if(isnan(n3ph))
        win_delta = NaN;
    else
        win_delta  = per.custom_win_delta;
    end
end

%'----------------------------Open JMAG Designer---------------------------
 JMAGversion='250'; 
% '# Create an "app" application object to launch JMAG-Designer
JDesigner = actxserver(strcat('designer.Application.',JMAGversion));
JDesigner.Show(); % Show the JMAG interface
% '------------------------------------------------------------------------
% 'Motor Winding diagram definitions :
Ncoilseries = geo.win.Ns;%number of series turns
if(isnan(n3ph))
    coilTurn = 5*Ncoilseries/((geo.Qs)*(geo.p*2));%stator coil number set in external circuit
else
    coilTurn = Ncoilseries/((geo.Qs/(3*n3ph))*(geo.p*2));%stator coil number set in external circuit
end

coilRes = per.Rs;% 'Stator coil resistance in ohms
coilLind = geo.win.Lend;% 'Stator coil LeakageInductance in Henry
phase_order = 1;% 'Phase order: 1--> UWV and 0 --> UVW
% '------------------------------------------------------------------------
% 'Motion Setting definitions :
th0 = (-geo.th0)/(geo.p);% 'Initial position of rotor in mechanical degrees
% '------------------------------------------------------------------------
% 'Study property settings: Step Control and Resolution
% 'Resolution or number of Step division per cycle of calculations
if(flagCharger)
    if (length(per.custom_rotorPos) == 1)
        CalcDivision = 9;
    else
        CalcDivision = size(per.custom_rotorPos,2)-1; 
    end
else
    CalcDivision = per.nsim_singt*geo.p*2;
end
%'#Number of electric periods to be calculed
% ncylce=per.delta_sim_singt*geo.p*2/360;
ncylce = 1;
% '------------------------------------------------------------------------
if strcmp(geo.RotType,'Hybrid')
    MagnetCG=geo.BLKLABELS.rotore.xy(geo.BLKLABELS.rotore.xy(:,3)==5,1:2);%center of gravity(G) of Magnet
else
    MagnetCG=geo.BLKLABELS.rotore.xy(geo.BLKLABELS.rotore.xy(:,3)==6,1:2);%center of gravity(G) of Magnet
end
Magnetxy=geo.rotor(geo.rotor(:,end-1)==6,:);
% '________________________________________________________________________
% '# Open jproj (Motor Drawing File)

JDesigner.NewProject('Untitled');
JDesigner.Load (fullfile(pathname, strcat(strrep(filename,'.mat','.jmag'),'.jproj')));

pwm_mod = false;
num_studies = 1;

if(JDesigner.NumStudies() == 2)
    if strcmp(geo.RotType,'Hybrid')
        num_studies = 2;
    else
        if(flagCharger)
            % pwm_mod = true;
            % num_studies = 2;
        end
    end
end

SOL_old = [];

for ii = 1:num_studies
    
    if strcmp(geo.RotType,'Hybrid')
        if(ii == 1)
            study_name = 'Demag';
        else
            study_name = 'Load';
        end
    else
        if(ii == 1)
            study_name = 'Load';
        else
            study_name = 'PWM';
        end
    end
    
    JDesigner.SetCurrentStudy(study_name);
    
    study = JDesigner.GetCurrentStudy();
    % '------------------------------------------------------------------------
    % 'Create the Design Parameters:
    Design_parameters = study.GetDesignTable();
    

    if(flagCharger)
        % 'initial position(motion) of the rotor in degree
        createDesignParameters_JMAG(Design_parameters,'InitPosition',0,th0(1)+rotorPos-1*180/geo.p,'Initial position(motion) of the rotor (degrees)')
        % 'rotation speed in rpm
        createDesignParameters_JMAG(Design_parameters,'rspeed',0,0,'rotation speed (rpm)')
    else
        % 'initial position(motion) of the rotor in degree
        createDesignParameters_JMAG(Design_parameters,'InitPosition',0,th0,'Initial position(motion) of the rotor (degrees)')
        % 'rotation speed in rpm
        createDesignParameters_JMAG(Design_parameters,'rspeed',0,per.EvalSpeed,'rotation speed (rpm)')
    end
    % 'Stator coil turn number
    createDesignParameters_JMAG(Design_parameters,'coilTurn',0,coilTurn,'Coil turn number')
    % 'Stator coil resistance in ohms
    createDesignParameters_JMAG(Design_parameters,'coilRes',0,coilRes,'Coil resistance (ohms)')
    % 'Stator coil LeakageInductance in Henry
    createDesignParameters_JMAG(Design_parameters,'coilLind',0,coilLind,'coil LeakageInductance (Henry)')
    
    if(flagCharger)
    
        if(strcmp(study_name,'PWM'))
            
            grid_amp = 230*sqrt(2)*1;
    
            % % 'Voltage_Grid in volt for grid voltage souces
            createDesignParameters_JMAG(Design_parameters,'Voltage_Grid',0,grid_amp,'Amplitude of grid voltage (V)')
            % % 'Freq_Grid in Hz for grid voltage souces
            createDesignParameters_JMAG(Design_parameters,'Freq_Grid',0,50,'Frequency of grid voltage (Hz)')
    
            for i=1:3
                % % 'Phase_Grid in degrees for grid voltage sources
                createDesignParameters_JMAG(Design_parameters,strcat('Phase_Grid','_',num2str(i)),0,90+120*(i-1),'Phase of grid voltages (deg)')
            end
        end
    
        if(isnan(n3ph))
            nph = 5;
        else
            nph = 3*n3ph;
        end

        for i=1:nph
            % % 'Amplitude in amps for stator coil current sources
            createDesignParameters_JMAG(Design_parameters,strcat('Aumont_Source','_',num2str(i)),0,custom_Amp(i),'Amplitude of motor drive (A or V)')
            % 'Phase in degrees for stator coil current sources
            createDesignParameters_JMAG(Design_parameters,strcat('Phase_Source','_',num2str(i)),0,90-custom_Ph(i)*180/pi,'Phase angle of motor drive (degrees)')
        end
    else
        
        src_amp = per.i0*per.overload;
    
        if strcmp(geo.RotType,'Hybrid')

            if(strcmp(study_name,'Demag'))
                src_amp = 0.1;
            end

            for i=1:3
                % % 'Amplitude in amps for stator coil current sources
                createDesignParameters_JMAG(Design_parameters,strcat('Aumont_Source','_',num2str(i)),0,src_amp,'Amplitude of motor drive (A or V)')
                % 'Phase in degrees for stator coil current sources
                createDesignParameters_JMAG(Design_parameters,strcat('Phase_Source','_',num2str(i)),0,per.gamma+90-120*(i-1),'Phase angle of motor drive (degrees)')
            end
        else
    
            if(isnan(n3ph))
                % % 'Aumont of drive in volt for voltage type and in ampere for current type for the stator coils
                createDesignParameters_JMAG(Design_parameters,'Aumont_Source',0,src_amp,'Amplitude of motor drive (A or V)')
                % 'Phase angle of drive for voltage type or current type for the stator coils
                createDesignParameters_JMAG(Design_parameters,'Phase_Source',0,per.gamma+90,'Phase angle of motor drive (degrees)')
            else
                % % 'Aumont of drive in volt for voltage type and in ampere for current type for the stator coils
                createDesignParameters_JMAG(Design_parameters,'Aumont_Source',0,src_amp,'Amplitude of motor drive (A or V)')
                % 'Phase angle of drive for voltage type or current type for the stator coils
                createDesignParameters_JMAG(Design_parameters,'Phase_Source',0,per.gamma+90,'Phase angle of motor drive (degrees)')
            end
        end
    end

    if(flagCharger)
        % 'Frequency of drive for voltage type or current type for the stator coils
        createDesignParameters_JMAG(Design_parameters,'Freq_Source',0,50,'Frequncy of motor drive (Hz)')
    else    
        % 'Frequency of drive for voltage type or current type for the stator coils
        createDesignParameters_JMAG(Design_parameters,'Freq_Source',0,per.EvalSpeed*geo.p/60,'Frequncy of motor drive (Hz)')
    end
    
    % 'Resolution or number of Step division per cycle of calculations
    createDesignParameters_JMAG(Design_parameters,'CalcDivision',0,CalcDivision,'Resolution or Step division per cycle of calculations')
    %'#Number of electric periods to be calculed
    createDesignParameters_JMAG(Design_parameters,'ncylce',0,ncylce,'Number of electric periods to be calculed')
    
    %% 'Study property settings
    % '------------------------------------------------------------------------
    if(flagCharger)
        %'Calculation time for one electric period / poles:
        CalcTime='ncylce/50';   % in second
    else
        %'Calculation time for one electric period / poles:
        CalcTime='ncylce/(rspeed*PolePair/60)';   % in second
    end
    
    % '# 'Step Control settings
    StepControl = study.GetStep();
    
    if(strcmp(study_name,'PWM'))
        JDesigner.SetCurrentStudy(study_name);
        JDesigner.GetDataManager().CreatePointArray('point_array/timevsnonlinear', 'SwitchingTable');
    
        [sarray_data, parray_data, pwm_pts] = getPWMSwitchingTable(n3ph, CalcDivision, SOL_old);
    
        CalcStep=num2str(length(sarray_data));
        StepControl.SetValue('Step', strcat(CalcStep));
    
        JDesigner.GetDataManager().GetDataSet('SwitchingTable').SetTable(sarray_data);
        StepControl.SetValue("StepType", 2);
        StepControl.SetTableProperty("Nonlinear", JDesigner.GetDataManager().GetDataSet('SwitchingTable'));
    else
        % ''#Number of calculation steps for Step Control
        CalcStep='CalcDivision+1';% in order to have ncycle of electric period
        StepControl.SetValue('Step', strcat(CalcStep));
        StepControl.SetValue('StepType', 1); % 0 = 'Simple', 1 = 'Uniform', 2 = 'Time'
        StepControl.SetValue('StepDivision ', 'CalcDivision');
        StepControl.SetValue('EndPoint', strcat(CalcTime)); %for 'StepType' = 1 
    end
    
    JDesigner.View().SetCurrentCase(study.GetCurrentCase());
    % '------------------------------------------------------------------------
    % %% 'Parallel Computing Setting
    % study.GetStudyProperties().SetValue('UseMultiCPU', 1);  
    % study.GetStudyProperties().SetValue('MultiCPU',8);
    % '#'______________________________________________________________________
    %% Electric Circuit of windings
    circuit=study.CreateCircuit();%Create the circuit
    % #shift between the components and FEM coils on X axis and Y axis
    dx=2; dy=3;
    % #Shift between the phase connections U - V - W
    duvw=dy;
    % #reference position of the first FEM coil
    X0=0; Y0=0;
    % #create the FEM Coil components
    % #FEM Coils U, V, W
    
    if(flagCharger)
        if(isnan(n3ph))
            Gu(1,1)=X0; Gu(1,2)=Y0;
            Gv(1,1)=X0; Gv(1,2)=Y0-duvw;
            Gw(1,1)=X0; Gw(1,2)=Y0-2*duvw;
            Gx(1,1)=X0; Gx(1,2)=Y0-3*duvw;
            Gy(1,1)=X0; Gy(1,2)=Y0-4*duvw;
            phaseLabel = {'U', 'V', 'W', 'X', 'Y'}; % or any other label
        else
            Gx = NaN;
            Gy = NaN;

            for i=1:1:n3ph
                % if(n3ph == 3)
                %     Gu(i,1)=X0; Gu(i,2)=Y0-3*(3*i-3)*duvw + 8*(i-1)*duvw;
                %     Gv(i,1)=X0; Gv(i,2)=Y0-3*(3*i-2)*duvw + 8*(i-1)*duvw;
                %     Gw(i,1)=X0; Gw(i,2)=Y0-3*(3*i-1)*duvw + 8*(i-1)*duvw;
                % else
                    Gu(i,1)=X0; Gu(i,2)=Y0-(3*i-3)*duvw;
                    Gv(i,1)=X0; Gv(i,2)=Y0-(3*i-2)*duvw;
                    Gw(i,1)=X0; Gw(i,2)=Y0-(3*i-1)*duvw;
                % end
                phaseLabel(3*i-2) = {strcat('U',num2str(i))}; %
                phaseLabel(3*i-1) = {strcat('V',num2str(i))}; %
                phaseLabel(3*i) = {strcat('W',num2str(i))}; %
            end

        end
    else
        if(isnan(n3ph))
            Gu(1,1)=X0; Gu(1,2)=Y0;
            Gv(1,1)=X0; Gv(1,2)=Y0-duvw;
            Gw(1,1)=X0; Gw(1,2)=Y0-2*duvw;
            Gx(1,1)=X0; Gx(1,2)=Y0-3*duvw;
            Gy(1,1)=X0; Gy(1,2)=Y0-4*duvw;
            phaseLabel = {'U', 'V', 'W', 'X', 'Y'}; % or any other label
        else
            Gu(1,1)=X0; Gu(1,2)=Y0;
            Gv(1,1)=X0; Gv(1,2)=Y0-duvw;
            Gw(1,1)=X0; Gw(1,2)=Y0-2*duvw;
            phaseLabel = {'U', 'V', 'W'}; % or any other label
        end
    end
    
    if strcmp(geo.RotType,'Hybrid')
        Grc(1,1)=X0; Grc(1,2)=Y0-6*duvw;
        phaseLabel = [phaseLabel, 'RC']; % or any other label
    end
   
    if(flagCharger)
        if(isnan(n3ph))
            positionU = [Gu(1,1), Gu(1,2)];
            positionV = [Gv(1,1), Gv(1,2)];
            positionW = [Gw(1,1), Gw(1,2)];
            positionX = [Gx(1,1), Gx(1,2)];
            positionY = [Gy(1,1), Gy(1,2)];
        else
            for i=1:1:n3ph
                positionU(i,:) = [Gu(i,1), Gu(i,2)];
                positionV(i,:) = [Gv(i,1), Gv(i,2)];
                positionW(i,:) = [Gw(i,1), Gw(i,2)];
            end
        end
    else
        positionU = [Gu(1,1), Gu(1,2)];
        positionV = [Gv(1,1), Gv(1,2)];
        positionW = [Gw(1,1), Gw(1,2)];
        if(isnan(n3ph))
            positionX = [Gx(1,1), Gx(1,2)];
            positionY = [Gy(1,1), Gy(1,2)];
        end
    end
    
    if strcmp(geo.RotType,'Hybrid')
        positionRC = [Grc(1,1), Grc(1,2)];
    end
    
    if(flagCharger)
        if(isnan(n3ph))
            FEMCoil_Creation_JMAG(circuit, phaseLabel{1}, positionU);
            FEMCoil_Creation_JMAG(circuit, phaseLabel{2}, positionV);
            FEMCoil_Creation_JMAG(circuit, phaseLabel{3}, positionW);
            FEMCoil_Creation_JMAG(circuit, phaseLabel{4}, positionX);
            FEMCoil_Creation_JMAG(circuit, phaseLabel{5}, positionY);
        else
            for i=1:1:n3ph
                FEMCoil_Creation_JMAG(circuit, phaseLabel{3*i-2}, positionU(i,:));
                FEMCoil_Creation_JMAG(circuit, phaseLabel{3*i-1}, positionV(i,:));
                FEMCoil_Creation_JMAG(circuit, phaseLabel{3*i}, positionW(i,:));
            end
        end
    else
        FEMCoil_Creation_JMAG(circuit, phaseLabel{1}, positionU);
        FEMCoil_Creation_JMAG(circuit, phaseLabel{2}, positionV);
        FEMCoil_Creation_JMAG(circuit, phaseLabel{3}, positionW);
        if(isnan(n3ph))
            FEMCoil_Creation_JMAG(circuit, phaseLabel{4}, positionX);
            FEMCoil_Creation_JMAG(circuit, phaseLabel{5}, positionY);
        end
    end
    
    if strcmp(geo.RotType,'Hybrid')
        FEMCoil_Creation_JMAG(circuit, phaseLabel{end}, positionRC);
    end
    
    if(flagCharger)
        if(strcmp(study_name,'PWM'))
            drawPWMcircuit(n3ph, Gu, Gv, Gw, dx, parray_data, circuit, JDesigner);
        else
            if(isnan(n3ph))
                nc = 5;
            else
                nc = 3*n3ph;
            end

            dx1=Gw(1,1)-4*dx; dy1=Gw(1,1);
            dx2=Gv(1,1)-4*dx; dy2=Gv(1,2);
            dx3=Gw(1,1)+4*dx;

            for i=1:nc
                circuit.CreateComponent('CurrentSource', strcat('CS',num2str(i)));
                circuit.CreateInstance(strcat('CS',num2str(i)), dx1, dy1-(i-1)*3);
                circuit.GetComponent(strcat('CS',num2str(i))).SetName(strcat(num2str(nc),'Phase-','Current','_',num2str(i)));
                circuit.GetComponent(strcat(num2str(nc),'Phase-','Current','_',num2str(i))).SetValue("FunctionType", 1);
                func = JDesigner.FunctionFactory().Sin(strcat('Aumont_Source','_',num2str(i)), 'Freq_Source', strcat('Phase_Source','_',num2str(i)), false);
                circuit.GetComponent(strcat(num2str(nc),'Phase-','Current','_',num2str(i))).SetFunction(func);        
            end
        end
    else
        if(isnan(n3ph))
            % # Add Electric source components into the circuit
            dx2=Gw(1,1)-4*dx; dy2=Gw(1,2);
            circuit.CreateComponent('MultiPhaseCurrentSource', 'CS1');
            circuit.CreateInstance('CS1', dx2, dy2-1);
            circuit.GetComponent('CS1').SetName(strcat('5Phase-','Current'));
            circuit.GetComponent(strcat('5Phase-','Current')).SetValue('Num Phases ', 5);
            circuit.GetComponent(strcat('5Phase-','Current')).SetValue('PhaseWindingType ', 0);
            %circuit.GetComponent(strcat('5Phase-','Current')).SetValue('Num Phases ', 5);
            circuit.GetComponent(strcat('5Phase-','Current')).SetValue('XType ', 'Time');
            % circuit.GetComponent(strcat('5Phase-','Current')).SetValue('CommutatingSequence', phase_order);
            circuit.GetComponent(strcat('5Phase-','Current')).SetValue('Amplitude', 'Aumont_Source');
            circuit.GetComponent(strcat('5Phase-','Current')).SetValue('Frequency', 'Freq_Source');
            circuit.GetComponent(strcat('5Phase-','Current')).SetValue('PhaseU', 'Phase_Source');
        else
            % # Add Electric source components into the circuit
            dx2=Gv(1,1)-4*dx; dy2=Gv(1,2);
            
            if strcmp(geo.RotType,'Hybrid')
                for i=1:3
                    circuit.CreateComponent('CurrentSource', strcat('CS',num2str(i)));
                    circuit.CreateInstance(strcat('CS',num2str(i)), dx2, dy2-(i-2)*3);
                    circuit.GetComponent(strcat('CS',num2str(i))).SetName(strcat(num2str(3),'Phase-','Current','_',num2str(i)));
                    circuit.GetComponent(strcat(num2str(3),'Phase-','Current','_',num2str(i))).SetValue("FunctionType", 1);
                    func = JDesigner.FunctionFactory().Sin(strcat('Aumont_Source','_',num2str(i)), 'Freq_Source', strcat('Phase_Source','_',num2str(i)), false);
                    func.SetValue("UseRangeEnd", 1);
                    func.SetValue("End", num2str(0.7*ncylce/(per.EvalSpeed*geo.p/60)));
                    circuit.GetComponent(strcat(num2str(3),'Phase-','Current','_',num2str(i))).SetFunction(func);       
                end
            else
                circuit.CreateComponent('3PhaseCurrentSource', 'CS1');
                circuit.CreateInstance('CS1', dx2, dy2);
                circuit.GetComponent('CS1').SetName(strcat('3Phase-','Current'));
                circuit.GetComponent(strcat('3Phase-','Current')).SetValue('XType ', 'Time');
                circuit.GetComponent(strcat('3Phase-','Current')).SetValue('CommutatingSequence', phase_order);
                circuit.GetComponent(strcat('3Phase-','Current')).SetValue('Amplitude', 'Aumont_Source');
                circuit.GetComponent(strcat('3Phase-','Current')).SetValue('Frequency', 'Freq_Source');
                circuit.GetComponent(strcat('3Phase-','Current')).SetValue('PhaseU', 'Phase_Source');
            end
        end
    end
    
    if strcmp(geo.RotType,'Hybrid')
        [dxrc, dyrc] = drawHybridCircuit(Grc, dx, ii, circuit, JDesigner);
    end

    if(flagCharger)
        if(isnan(n3ph))
            %Add Voltage probe for U-Phase
            AddVoltageProbe_JMAG (circuit, 1, positionU+[-dx, 2], phaseLabel{1})
            %Add Voltage probe for V-Phase
            AddVoltageProbe_JMAG (circuit, 2, positionV+[-dx, 2], phaseLabel{2})
            %Add Voltage probe for W-Phase
            AddVoltageProbe_JMAG (circuit, 3, positionW+[-dx, 2], phaseLabel{3})
            %Add Voltage probe for X-Phase
            AddVoltageProbe_JMAG (circuit, 4, positionX+[-dx, 2], phaseLabel{4})
            %Add Voltage probe for Y-Phase
            AddVoltageProbe_JMAG (circuit, 5, positionY+[-dx, 2], phaseLabel{5})
        else
            for i=1:1:n3ph
                %Add Voltage probe for X-Phase
                AddVoltageProbe_JMAG (circuit, 3*i-2, positionU(i,:)+[-dx, 2], phaseLabel{3*i-2})
                %Add Voltage probe for Y-Phase
                AddVoltageProbe_JMAG (circuit, 3*i-1, positionV(i,:)+[-dx, 2], phaseLabel{3*i-1})
                %Add Voltage probe for Z-Phase
                AddVoltageProbe_JMAG (circuit, 3*i, positionW(i,:)+[-dx, 2], phaseLabel{3*i})     
            end
        end
    else
        %Add Voltage probe for U-Phase
        AddVoltageProbe_JMAG (circuit, 1, positionU+[-dx, 2], phaseLabel{1})
        %Add Voltage probe for V-Phase
        AddVoltageProbe_JMAG (circuit, 2, positionV+[-dx, 2], phaseLabel{2})
        %Add Voltage probe for W-Phase
        AddVoltageProbe_JMAG (circuit, 3, positionW+[-dx, 2], phaseLabel{3})
        if(isnan(n3ph))
            %Add Voltage probe for X-Phase
            AddVoltageProbe_JMAG (circuit, 4, positionX+[-dx, 2], phaseLabel{4})
            %Add Voltage probe for Y-Phase
            AddVoltageProbe_JMAG (circuit, 5, positionY+[-dx, 2], phaseLabel{5})
        end
    end
    
    if(strcmp(study_name,'PWM'))
        if(isnan(n3ph))
            %Add Voltage probe for U-Phase
            AddVoltageProbe_JMAG (circuit, 1, positionU+[dx, 2], [phaseLabel{1} '1'])
            %Add Voltage probe for V-Phase
            AddVoltageProbe_JMAG (circuit, 2, positionV+[dx, 2], [phaseLabel{2} '1'])
            %Add Voltage probe for W-Phase
            AddVoltageProbe_JMAG (circuit, 3, positionW+[dx, 2], [phaseLabel{3} '1'])
            %Add Voltage probe for X-Phase
            AddVoltageProbe_JMAG (circuit, 4, positionX+[dx, 2], [phaseLabel{4} '1'])
            %Add Voltage probe for Y-Phase
            AddVoltageProbe_JMAG (circuit, 5, positionY+[dx, 2], [phaseLabel{5} '1'])
        else
            for i=1:n3ph
                %Add Voltage probe for U-Phase
                AddVoltageProbe_JMAG (circuit, 3*i-2, positionU(i,:)+[dx, 2], [phaseLabel{3*i-2} '1'])
                %Add Voltage probe for V-Phase
                AddVoltageProbe_JMAG (circuit, 3*i-1, positionV(i,:)+[dx, 2], [phaseLabel{3*i-1} '1'])
                %Add Voltage probe for W-Phase
                AddVoltageProbe_JMAG (circuit, 3*i, positionW(i,:)+[dx, 2], [phaseLabel{3*i} '1'])
            end
        end
    end
    
    if strcmp(geo.RotType,'Hybrid')
        AddVoltageProbe_JMAG (circuit, 5, positionRC+[-dx, 2], phaseLabel{end})
    end
    
    %Add a Ground component
    circuit.CreateComponent('Ground', 'Ground');
    if(flagCharger)
        if(isnan(n3ph))
            if(strcmp(study_name,'PWM'))
                circuit.CreateInstance('Ground', Gw(1,1)+7*dx, Gw(1,2)-dy);
            else
                circuit.CreateInstance('Ground', Gw(1,1)+3*dx, Gw(1,2)-dy);
            end
        else
            if(strcmp(study_name,'PWM'))
                if(n3ph == 3)
                    circuit.CreateInstance('Ground', Gw(1,1)+7*dx, Gu(3,2)-dy);
                else
                    circuit.CreateInstance('Ground', Gw(1,1)+7*dx, Gw(1,2)-dy);
                end
            else
                circuit.CreateInstance('Ground', Gw(1,1)+3*dx, Gw(1,2)-dy);
            end        
        end
    else
        if(isnan(n3ph))
            circuit.CreateInstance('Ground', Gw(1,1)+3*dx, Gw(1,2)-dy);
        else
            circuit.CreateInstance('Ground', Gv(1,1)+3*dx, Gv(1,2)-dy);
        end
    end

    if strcmp(geo.RotType,'Hybrid')
        circuit.CreateInstance('Ground', Grc(1,1)+3*dx, Grc(1,2)-2);
    end
    
    if(strcmp(study_name,'PWM'))
        if(isnan(n3ph))
            % #create star connection of FEM Coils
            circuit.CreateWire(Gu(1,1)+2+4*dx, Gu(1,2), Gv(1,1)+2+4*dx, Gv(1,2));
            circuit.CreateWire(Gv(1,1)+2+4*dx, Gv(1,2), Gw(1,1)+2+4*dx, Gw(1,2));
            circuit.CreateWire(Gw(1,1)+2+4*dx, Gw(1,2), Gx(1,1)+2+4*dx, Gx(1,2));
            circuit.CreateWire(Gx(1,1)+2+4*dx, Gx(1,2), Gy(1,1)+2+4*dx, Gy(1,2));
        else
            Gs = [Gu; Gv; Gw];
            for i=1:1:3*n3ph-1
                ind = floor((i-1)/3)+n3ph*mod(i-1,3)+1;
                ind1 = floor((i)/3)+n3ph*mod(i,3)+1;
                % #create star connection of FEM Coils
                circuit.CreateWire(Gs(ind,1)+2+4*dx, Gs(ind,2), Gs(ind1,1)+2+4*dx, Gs(ind1,2));
            end
        end
    else

        if(flagCharger)
            if(isnan(n3ph))
                % #create star connection of FEM Coils
                circuit.CreateWire(Gu(1,1)+2, Gu(1,2), Gv(1,1)+2, Gv(1,2));
                circuit.CreateWire(Gv(1,1)+2, Gv(1,2), Gw(1,1)+2, Gw(1,2));
                circuit.CreateWire(Gw(1,1)+2, Gw(1,2), Gx(1,1)+2, Gx(1,2));
                circuit.CreateWire(Gx(1,1)+2, Gx(1,2), Gy(1,1)+2, Gy(1,2));
            else
                Gs = [Gu; Gv; Gw];
                for i=1:1:3*n3ph-1
                    ind = floor((i-1)/3)+n3ph*mod(i-1,3)+1;
                    ind1 = floor((i)/3)+n3ph*mod(i,3)+1;

                    % #create star connection of FEM Coils
                    circuit.CreateWire(Gs(ind,1)+2, Gs(ind,2), Gs(ind1,1)+2, Gs(ind1,2));
                end
            end
        else
            % #create star connection of FEM Coils
            circuit.CreateWire(Gu(1,1)+2, Gu(1,2), Gv(1,1)+2, Gv(1,2));
            circuit.CreateWire(Gv(1,1)+2, Gv(1,2), Gw(1,1)+2, Gw(1,2));

            if(isnan(n3ph))
                % #create star connection of FEM Coils
                circuit.CreateWire(Gw(1,1)+2, Gw(1,2), Gx(1,1)+2, Gx(1,2));
                circuit.CreateWire(Gx(1,1)+2, Gx(1,2), Gy(1,1)+2, Gy(1,2));
            end
        end
    end
    
    if(flagCharger)
        %Connect Electric source to the Coils.
        if(strcmp(study_name,'PWM'))
            for i=1:1:n3ph
                if(n3ph == 3)
                    Gu(i,1)=X0; Gu(i,2)=Y0-3*(3*i-3)*duvw + 8*(i-1)*duvw;
                    Gv(i,1)=X0; Gv(i,2)=Y0-3*(3*i-2)*duvw + 8*(i-1)*duvw;
                    Gw(i,1)=X0; Gw(i,2)=Y0-3*(3*i-1)*duvw + 8*(i-1)*duvw;
                end
            end

            drawPWMConnections(n3ph, Gu, Gv, Gw, Gx, Gy, dx, dx2, circuit);
        else
            if(isnan(n3ph))
                circuit.CreateWire(dx2+2, Gu(1,2), Gu(1,1)-2, Gu(1,2));
                circuit.CreateWire(dx2+2, Gv(1,2), Gv(1,1)-2, Gv(1,2));
                circuit.CreateWire(dx2+2, Gw(1,2), Gw(1,1)-2, Gw(1,2));
                circuit.CreateWire(dx2+2, Gx(1,2), Gx(1,1)-2, Gx(1,2));
                circuit.CreateWire(dx2+2, Gy(1,2), Gy(1,1)-2, Gy(1,2));
            else
                for i=1:1:n3ph
                    circuit.CreateWire(dx2+2, Gu(i,2), Gu(i,1)-2, Gu(i,2));
                    circuit.CreateWire(dx2+2, Gv(i,2), Gv(i,1)-2, Gv(i,2));
                    circuit.CreateWire(dx2+2, Gw(i,2), Gw(i,1)-2, Gw(i,2));
                end
            end
        end
    else
        if(isnan(n3ph))
            circuit.CreateWire(dx2+2, dy2+4, Gu(1,1)-2, Gu(1,2));
            circuit.CreateWire(dx2+2, dy2+2, Gv(1,1)-2, Gv(1,2));
            circuit.CreateWire(dx2+2, dy2, Gw(1,1)-2, Gw(1,2));
            circuit.CreateWire(dx2+2, dy2-2, Gx(1,1)-2, Gx(1,2));
            circuit.CreateWire(dx2+2, dy2-4, Gy(1,1)-2, Gy(1,2));
        else
            if strcmp(geo.RotType,'Hybrid')
                %Connect Electric source to the Coils.
                circuit.CreateWire(dx2+2, dy2+3, Gu(1,1)-2, Gu(1,2));
                circuit.CreateWire(dx2+2, dy2, Gv(1,1)-2, Gv(1,2));
                circuit.CreateWire(dx2+2, dy2-3, Gw(1,1)-2, Gw(1,2));
            else
                %Connect Electric source to the Coils.
                circuit.CreateWire(dx2+2, dy2+2, Gu(1,1)-2, Gu(1,2));
                circuit.CreateWire(dx2+2, dy2, Gv(1,1)-2, Gv(1,2));
                circuit.CreateWire(dx2+2, dy2-2, Gw(1,1)-2, Gw(1,2));
            end
        end
    end
    
    if strcmp(geo.RotType,'Hybrid')
        %Connect Electric source to the Coils.
        circuit.CreateWire(dxrc+2, dyrc, Grc(1,1)-2, Grc(1,2));   
    end
    
    if(flagCharger)
        if(isnan(n3ph))
            %Connect Ground
            if(strcmp(study_name,'PWM'))
                circuit.CreateWire(Gw(1,1)+2+4*dx, Gw(1,2), Gw(1,1)+3*dx+4*dx, Gw(1,2)-dy+2);
            else
                circuit.CreateWire(Gw(1,1)+2, Gw(1,2), Gw(1,1)+3*dx, Gw(1,2)-dy+2);
            end
        else
            %Connect Ground
            if(strcmp(study_name,'PWM'))
                circuit.CreateWire(Gw(1,1)+2+4*dx, Gw(1,2), Gw(1,1)+3*dx+4*dx, Gw(1,2)-dy+2);
            else
                circuit.CreateWire(Gw(1,1)+2, Gw(1,2), Gw(1,1)+3*dx, Gw(1,2)-dy+2);
            end
        end
    else
        if(isnan(n3ph))
            %Connect Ground
            circuit.CreateWire(Gw(1,1)+2, Gw(1,2), Gw(1,1)+3*dx, Gw(1,2)-dy+2);
        else        
            %Connect Ground
            circuit.CreateWire(Gv(1,1)+2, Gv(1,2), Gv(1,1)+3*dx, Gv(1,2)-dy+2);
        end
    end
    
    if strcmp(geo.RotType,'Hybrid')
        %Connect Ground
        circuit.CreateWire(Grc(1,1)+2, Grc(1,2), Grc(1,1)+3*dx, Grc(1,2));
    end
    
    % '------------------------------------------------------------------------
    %% Apply Condictions
    % '#'______________________________________________________________________
    % '# Apply Motion Condition
    study.CreateCondition('RotationMotion', 'Speed');
    Motion_Condition=study.GetCondition('Speed');
    Motion_Condition.SetValue('AngularVelocity', 'rspeed');
    Motion_Condition.SetValue('InitialRotationAngle', 'InitPosition');
    Motion_Condition.ClearParts();
    sel = Motion_Condition.GetSelection();
    sel.SelectPart('Rotor Core');
    if any(unique(Magnetxy(:,9)) ~= 0)
    for ID = 1:1:length(MagnetCG)
    sel.SelectPart(strcat('Magnet',num2str(ID)));
    end
    end
    
    if strcmp(geo.RotType,'Hybrid')
    
        RCoilCG=geo.BLKLABELS.rotore.xy(geo.BLKLABELS.rotore.xy(:,3)==8,1:2);%center of gravity(G) of Copper Coils
    
        sel.SelectPart('RCoil');
    end
    
    Motion_Condition.AddSelected(sel);
    JDesigner.SetCurrentStudy(study_name);
    
    % '#'_______________________________________________________________
    % '# Apply Torque Condition
    study.CreateCondition('Torque', 'Electromagnetic_Torque');
    Torque_Condition=study.GetCondition('Electromagnetic_Torque');
    Torque_Condition.SetValue('TargetType', 0);
    Torque_Condition.ClearParts();
    Torque_Condition.SetValue('TargetType', 1);
    Torque_Condition.SetLinkWithType('LinkedMotion', 'Speed');
    
    JDesigner.SetCurrentStudy(study_name);
    % '#'_______________________________________________________________
    % '# Apply Iron loss condition on Stator Core
    Iron_losscondition_JMAG (study,'Stator Core','Stator_Core_Loss',1,2,1, flagCharger)
    % '#'_______________________________________________________________
    % '# Apply Iron loss condition on Rotor Core
    Iron_losscondition_JMAG (study,'Rotor Core','Rotor_Core_Loss',1,2,1, flagCharger)
    % '#'_______________________________________________________________
    % '##FEM Coil conditions for Winding Connections
    % '_______________________________________________________________
    % coil direction Definition: Upward or Downward
    % coildirectionp = 0; % upward coil direction, coildirectionn = 1; % downward coil direction
    
    if(flagCharger)
        if(isnan(n3ph))
            % 'Create FEM Coil conditions for U-phase
            Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{1},+1,-1)
            % '_______________________________________________________________
            % 'Create FEM Coil conditions for X-phase
            Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{4},+4,-4)
            % '_______________________________________________________________
            % 'Create FEM Coil conditions for V-phase
            Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{2},+2,-2)
            % '_______________________________________________________________
            % 'Create FEM Coil conditions for Y-phase
            Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{5},+5,-5)
            % '_______________________________________________________________
            % 'Create FEM Coil conditions for W-phase
            Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{3},+3,-3)
        else
            for i=1:n3ph
                % 'Create FEM Coil conditions for U-phase
                Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{3*i-2},+(3*i-2),-(3*i-2))
                % '_______________________________________________________________
                % 'Create FEM Coil conditions for V-phase
                Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{3*i-1},+(3*i-1),-(3*i-1))
                % '_______________________________________________________________
                % 'Create FEM Coil conditions for W-phase
                Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{3*i},+3*i,-3*i)
            end
        end
    else
        if(isnan(n3ph))
            % 'Create FEM Coil conditions for U-phase
            Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{1},+1,-1)
            % '_______________________________________________________________
            % 'Create FEM Coil conditions for X-phase
            Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{4},+4,-4)
            % '_______________________________________________________________
            % 'Create FEM Coil conditions for V-phase
            Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{2},+2,-2)
            % '_______________________________________________________________
            % 'Create FEM Coil conditions for Y-phase
            Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{5},+5,-5)
            % '_______________________________________________________________
            % 'Create FEM Coil conditions for W-phase
            Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{3},+3,-3)
        else
            % 'Create FEM Coil conditions for U-phase
            Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{1},+1,-1)
            % '_______________________________________________________________
            % 'Create FEM Coil conditions for W-phase
            Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{3},+3,-3)
            % '_______________________________________________________________
            % 'Create FEM Coil conditions for V-phase
            Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{2},+2,-2)
        end
    end
    
    if strcmp(geo.RotType,'Hybrid')
        RCoil_DirectionDefinition_JMAG (geo,study,phaseLabel{end})
    end
    
    JDesigner.SetCurrentStudy(study_name);
    
    if strcmp(geo.RotType,'Hybrid')
        setMagParams(ii, study_name, study, JDesigner);
    end
    
    %% -------------------------Generate Mesh----------------------------------
    % 'Air region Mesh Scale factor:
    air_mesh_scale=1.2;
    
    % '#'_______________________________________________________________
    % '#MESH
    Mesh = study.GetMeshControl();
    Mesh.SetValue('MeshType', 1);
    Mesh.SetValue('AutoGapDivision', 0);
    Mesh.SetValue('AutoDivision', 0);
    Mesh.SetValue('RadialDivision', 'RDivision');
    Mesh.SetValue('CircumferentialDivision', 'CDivision');
    Mesh.SetValue('AirRegionScale', air_mesh_scale);
    
    % '#'Apply the rotation periodic mesh
    Mesh.CreateCondition('RotationPeriodicMeshAutomatic', "symmetric mesh");
    
    % '#Stator core mesh size
    Create_Mesh_JMAG(Mesh,'Stator','Stator Core','stator_mesh_size')
    
    % '#Rotor core mesh size
    Create_Mesh_JMAG(Mesh,'Rotor','Rotor Core','rotor_mesh_size')
    
    % '#Magnet mesh size
    if any(unique(Magnetxy(:,9)) ~= 0)
    Part_mesh = Mesh.CreateCondition('Part', 'Magnet');
    Part_mesh.SetName('Magnet');
    Part_mesh.SetValue('Size', 'magnet_mesh_size');
    Part_mesh.ClearParts();
    sel = Part_mesh.GetSelection();
    for ID = 1:1:length(MagnetCG)
    sel.SelectPart(strcat('Magnet',num2str(ID)));
    end
    Part_mesh.AddSelected(sel);
    end
    
    
    % '#Coil mesh size
    Coil_mesh = Mesh.CreateCondition('Part', 'Coil');
    Coil_mesh.SetName('Coil');
    Coil_mesh.SetValue('Size', 'coil_mesh_size');
    Coil_mesh.ClearParts();
    sel = Coil_mesh.GetSelection();

    if(flagCharger)
        if(isnan(n3ph))
            sel.SelectPart(strcat('Coil','_U'));
            sel.SelectPart(strcat('Coil','_V'));
            sel.SelectPart(strcat('Coil','_W'));
            sel.SelectPart(strcat('Coil','_X'));
            sel.SelectPart(strcat('Coil','_Y'));
        else
            for i=1:1:n3ph
                sel.SelectPart(strcat('Coil','_U',num2str(i)));
                sel.SelectPart(strcat('Coil','_V',num2str(i)));
                sel.SelectPart(strcat('Coil','_W',num2str(i)));  
            end
        end
    else
        sel.SelectPart(strcat('Coil','_U'));
        sel.SelectPart(strcat('Coil','_V'));
        sel.SelectPart(strcat('Coil','_W'));

        if(isnan(n3ph))
            sel.SelectPart(strcat('Coil','_X'));
            sel.SelectPart(strcat('Coil','_Y'));
        end
    end

    Coil_mesh.AddSelected(sel);
    
    JDesigner.SetCurrentStudy(study_name);
    
    %Generate the mesh
    study.CreateMesh();
    JDesigner.View().ShowMesh();
    
    
    %% Run SMULATION
    study.Run();
    JDesigner.Save();
    
    
    %% Post processing
    % '#'______________________________________________________________________
    % Results Graphs Data
    Data_Manager = JDesigner.GetDataManager();
    
    if(flagCharger)
        if(isnan(n3ph))
            Uphase_name = 'U-Phase Coil';Vphase_name = 'V-Phase Coil';Wphase_name = 'W-Phase Coil';
            Xphase_name = 'X-Phase Coil'; Yphase_name = 'Y-Phase Coil';
        else
            for i=1:1:n3ph
                Uphase_name(i,:) = strcat('U',num2str(i),'-Phase Coil');
                Vphase_name(i,:) = strcat('V',num2str(i),'-Phase Coil');
                Wphase_name(i,:) = strcat('W',num2str(i),'-Phase Coil');
            end
        end
    else
        Uphase_name = 'U-Phase Coil';Vphase_name = 'V-Phase Coil';Wphase_name = 'W-Phase Coil';

        if(isnan(n3ph))
            Xphase_name = 'X-Phase Coil'; Yphase_name = 'Y-Phase Coil';
        end
    end
    
    JDesigner.SetCurrentStudy(study_name);
    
    % Results Graphs of Coil Flux-Linkage
    CoilFluxLinkage_result = study.GetDataSet('Coil Flux-Linkage',1);
    Data_Manager.CreateAllCasesGraphModel(CoilFluxLinkage_result,'Coil Flux-Linkage');
    % '------------------------------------------------------------------------
    %% Results Graphs of CoilFluxLinkages
    NTL = CoilFluxLinkage_result.GetRows();
    
    Fluxa = zeros(NTL,1);
    Fluxb = zeros(NTL,1);
    Fluxc = zeros(NTL,1);
    if(flagCharger)
        if(isnan(n3ph))
            Fluxd = zeros(NTL,1);
            Fluxe = zeros(NTL,1);         
        else
            Fluxd = zeros(NTL,n3ph);
            Fluxq = zeros(NTL,n3ph);
            Flux0 = zeros(NTL,n3ph);
        end
    else
        if(isnan(n3ph))
            Fluxd = zeros(NTL,1);
            Fluxe = zeros(NTL,1);
        end
    end

    if strcmp(geo.RotType,'Hybrid')
        cir_pos = 2;
    else
        cir_pos = 1;
    end
    
    for row =1:1:NTL
    
        if(flagCharger)
            if(isnan(n3ph))
                if strcmp(CoilFluxLinkage_result.GetColumnName(cir_pos), Uphase_name) == 1
                    Fluxa(row,1) = CoilFluxLinkage_result.GetValue(row-1, cir_pos);
                end
                if strcmp(CoilFluxLinkage_result.GetColumnName(cir_pos+1), Vphase_name) == 1
                    Fluxb(row,1) = CoilFluxLinkage_result.GetValue(row-1, cir_pos+1);
                end
                if strcmp(CoilFluxLinkage_result.GetColumnName(cir_pos+2), Wphase_name) == 1
                    Fluxc(row,1) = CoilFluxLinkage_result.GetValue(row-1, cir_pos+2);
                end
                if strcmp(CoilFluxLinkage_result.GetColumnName(cir_pos+3), Xphase_name) == 1
                    Fluxd(row,1) = CoilFluxLinkage_result.GetValue(row-1, cir_pos+3);
                end
                if strcmp(CoilFluxLinkage_result.GetColumnName(cir_pos+4), Yphase_name) == 1
                    Fluxe(row,1) = CoilFluxLinkage_result.GetValue(row-1, cir_pos+4);
                end
            else
                for i=1:1:n3ph
                    if strcmp(CoilFluxLinkage_result.GetColumnName(cir_pos+i-1), Uphase_name(i,:)) == 1
                        Fluxa(row,i) = CoilFluxLinkage_result.GetValue(row-1, cir_pos + i-1);
                    end
                    if strcmp(CoilFluxLinkage_result.GetColumnName(cir_pos+n3ph+i-1), Vphase_name(i,:)) == 1
                        Fluxb(row,i) = CoilFluxLinkage_result.GetValue(row-1, cir_pos + n3ph+i-1);
                    end   
                    if strcmp(CoilFluxLinkage_result.GetColumnName(cir_pos+2*n3ph+i-1), Wphase_name(i,:)) == 1
                        Fluxc(row,i) = CoilFluxLinkage_result.GetValue(row-1, cir_pos + 2*n3ph+i-1);
                    end    
                end
            end
        else
            if strcmp(CoilFluxLinkage_result.GetColumnName(cir_pos), Uphase_name) == 1
                Fluxa(row,1) = CoilFluxLinkage_result.GetValue(row-1, cir_pos);
            end
            if strcmp(CoilFluxLinkage_result.GetColumnName(cir_pos+1), Vphase_name) == 1
                Fluxb(row,1) = CoilFluxLinkage_result.GetValue(row-1, cir_pos+1);
            end
            if strcmp(CoilFluxLinkage_result.GetColumnName(cir_pos+2), Wphase_name) == 1
                Fluxc(row,1) = CoilFluxLinkage_result.GetValue(row-1, cir_pos+2);
            end
            if(isnan(n3ph))
                if strcmp(CoilFluxLinkage_result.GetColumnName(cir_pos+3), Xphase_name) == 1
                    Fluxd(row,1) = CoilFluxLinkage_result.GetValue(row-1, cir_pos+3);
                end
                if strcmp(CoilFluxLinkage_result.GetColumnName(cir_pos+4), Yphase_name) == 1
                    Fluxe(row,1) = CoilFluxLinkage_result.GetValue(row-1, cir_pos+4);
                end
            end
        end
    end
    % '------------------------------------------------------------------------
    %% Results Graphs of Circuit Current
  
    if(flagCharger)
        if(isnan(n3ph))
            ia = zeros(NTL,1);
            ib = zeros(NTL,1);
            ic = zeros(NTL,1);
            ua = zeros(NTL,1);
            ub = zeros(NTL,1);
            uc = zeros(NTL,1);
            id = zeros(NTL,1);
            ie = zeros(NTL,1);
            ud = zeros(NTL,1);
            ue = zeros(NTL,1);
        else
            ia = zeros(NTL,n3ph);
            ib = zeros(NTL,n3ph);
            ic = zeros(NTL,n3ph);
            ua = zeros(NTL,n3ph);
            ub = zeros(NTL,n3ph);
            uc = zeros(NTL,n3ph);
            id = zeros(NTL,n3ph);
            iq = zeros(NTL,n3ph);
            i0 = zeros(NTL,n3ph);
        end
    else
        ia = zeros(NTL,1);
        ib = zeros(NTL,1);
        ic = zeros(NTL,1);
        ua = zeros(NTL,1);
        ub = zeros(NTL,1);
        uc = zeros(NTL,1);

        if(isnan(n3ph))
            id = zeros(NTL,1);
            ie = zeros(NTL,1);
            ud = zeros(NTL,1);
            ue = zeros(NTL,1);
        end
    end

    if(flagCharger)
        if(isnan(n3ph))
            if(strcmp(study_name,'PWM'))
                cir_pos = 6;
            else
                cir_pos = 1;
            end
        else
            if(strcmp(study_name,'PWM'))
                cir_pos = 3*n3ph+1;
            else
                cir_pos = 1;
            end
        end
    else
        cir_pos = 2;
    end
    
    % if strcmp(geo.RotType,'Hybrid')
    %     cir_pos = cir_pos + 1;
    % end
    
    CircuitCurrent_result = study.GetDataSet('Circuit Current',1);
    Data_Manager.CreateAllCasesGraphModel(CircuitCurrent_result,'Circuit Current Graph');
    for row =1:1:NTL
        if(flagCharger)
            if(isnan(n3ph))
                if strcmp(CircuitCurrent_result.GetColumnName(cir_pos), Uphase_name) == 1
                    ia(row,1) = CircuitCurrent_result.GetValue(row-1, cir_pos);
                end
                if strcmp(CircuitCurrent_result.GetColumnName(cir_pos+1), Vphase_name) == 1
                    ib(row,1) = CircuitCurrent_result.GetValue(row-1, cir_pos+1);
                end
                if strcmp(CircuitCurrent_result.GetColumnName(cir_pos+2), Wphase_name) == 1
                    ic(row,1) = CircuitCurrent_result.GetValue(row-1, cir_pos+2);
                end
                if strcmp(CircuitCurrent_result.GetColumnName(cir_pos+3), Xphase_name) == 1
                    id(row,1) = CircuitCurrent_result.GetValue(row-1, cir_pos+3);
                end
                if strcmp(CircuitCurrent_result.GetColumnName(cir_pos+4), Yphase_name) == 1
                    ie(row,1) = CircuitCurrent_result.GetValue(row-1, cir_pos+4);
                end
            else
                for i=1:1:n3ph
                    if strcmp(CircuitCurrent_result.GetColumnName(cir_pos+i-1), Uphase_name(i,:)) == 1
                        ia(row,i) = CircuitCurrent_result.GetValue(row-1, cir_pos+i-1);
                    end
                    if strcmp(CircuitCurrent_result.GetColumnName(cir_pos+n3ph+i-1), Vphase_name(i,:)) == 1
                        ib(row,i) = CircuitCurrent_result.GetValue(row-1, cir_pos+n3ph+i-1);
                    end
                    if strcmp(CircuitCurrent_result.GetColumnName(cir_pos+2*n3ph+i-1), Wphase_name(i,:)) == 1
                        ic(row,i) = CircuitCurrent_result.GetValue(row-1, cir_pos+2*n3ph+i-1);
                    end
                end
            end
        else
            if(isnan(n3ph))
                if strcmp(CircuitCurrent_result.GetColumnName(cir_pos), Uphase_name) == 1
                    ia(row,1) = CircuitCurrent_result.GetValue(row-1, cir_pos);
                end
                if strcmp(CircuitCurrent_result.GetColumnName(cir_pos+1), Vphase_name) == 1
                    ib(row,1) = CircuitCurrent_result.GetValue(row-1, cir_pos+1);
                end
                if strcmp(CircuitCurrent_result.GetColumnName(cir_pos+2), Wphase_name) == 1
                    ic(row,1) = CircuitCurrent_result.GetValue(row-1, cir_pos+2);
                end
                if strcmp(CircuitCurrent_result.GetColumnName(cir_pos+3), Xphase_name) == 1
                    id(row,1) = CircuitCurrent_result.GetValue(row-1, cir_pos+3);
                end
                if strcmp(CircuitCurrent_result.GetColumnName(cir_pos+4), Yphase_name) == 1
                    ie(row,1) = CircuitCurrent_result.GetValue(row-1, cir_pos+4);
                end
            else
                if strcmp(CircuitCurrent_result.GetColumnName(cir_pos), Uphase_name) == 1
                    ia(row,1) = CircuitCurrent_result.GetValue(row-1, cir_pos);
                end
                if strcmp(CircuitCurrent_result.GetColumnName(cir_pos+1), Vphase_name) == 1
                    ib(row,1) = CircuitCurrent_result.GetValue(row-1, cir_pos+1);
                end
                if strcmp(CircuitCurrent_result.GetColumnName(cir_pos+2), Wphase_name) == 1
                    ic(row,1) = CircuitCurrent_result.GetValue(row-1, cir_pos+2);
                end
            end
        end
    end
    
    cir_pos = 1;
    
    if strcmp(geo.RotType,'Hybrid')
        cir_pos = cir_pos + 1;
    end
    
    CircuitVoltage_result = study.GetDataSet('Circuit Voltage',1);
    Data_Manager.CreateAllCasesGraphModel(CircuitVoltage_result,'Circuit Voltage Graph');
    if(flagCharger)
        for row =1:1:NTL
            if(isnan(n3ph))
                if strcmp(CircuitVoltage_result.GetColumnName(cir_pos), Uphase_name(1)) == 1
                    ua(row,1) = CircuitVoltage_result.GetValue(row-1, cir_pos);
                end
                if strcmp(CircuitVoltage_result.GetColumnName(cir_pos+1), Vphase_name(1)) == 1
                    ub(row,1) = CircuitVoltage_result.GetValue(row-1, cir_pos+1);
                end
                if strcmp(CircuitVoltage_result.GetColumnName(cir_pos+2), Wphase_name(1)) == 1
                    uc(row,1) = CircuitVoltage_result.GetValue(row-1, cir_pos+2);
                end
                if strcmp(CircuitVoltage_result.GetColumnName(cir_pos+3), Xphase_name(1)) == 1
                    ud(row,1) = CircuitVoltage_result.GetValue(row-1, cir_pos+3);
                end
                if strcmp(CircuitVoltage_result.GetColumnName(cir_pos+4), Yphase_name(1)) == 1
                    ue(row,1) = CircuitVoltage_result.GetValue(row-1, cir_pos+4);
                end
            else
                for i=1:1:n3ph
                    if strcmp(CircuitVoltage_result.GetColumnName(cir_pos+i-1), Uphase_name(i,:)) == 1
                        ua(row,i) = CircuitVoltage_result.GetValue(row-1, cir_pos+i-1);
                    end
                    if strcmp(CircuitVoltage_result.GetColumnName(cir_pos+n3ph+i-1), Vphase_name(i,:)) == 1
                        ub(row,i) = CircuitVoltage_result.GetValue(row-1, cir_pos+n3ph+i-1);
                    end
                    if strcmp(CircuitVoltage_result.GetColumnName(cir_pos+2*n3ph+i-1), Wphase_name(i,:)) == 1
                        uc(row,i) = CircuitVoltage_result.GetValue(row-1, cir_pos+2*n3ph+i-1);
                    end
                end
            end
        end
    else
        for row =1:1:NTL
            if(isnan(n3ph))
                if strcmp(CircuitVoltage_result.GetColumnName(cir_pos), Uphase_name(1)) == 1
                    ua(row,1) = CircuitVoltage_result.GetValue(row-1, cir_pos);
                end
                if strcmp(CircuitVoltage_result.GetColumnName(cir_pos+1), Vphase_name(1)) == 1
                    ub(row,1) = CircuitVoltage_result.GetValue(row-1, cir_pos+1);
                end
                if strcmp(CircuitVoltage_result.GetColumnName(cir_pos+2), Wphase_name(1)) == 1
                    uc(row,1) = CircuitVoltage_result.GetValue(row-1, cir_pos+2);
                end
                if strcmp(CircuitVoltage_result.GetColumnName(cir_pos+3), Xphase_name(1)) == 1
                    ud(row,1) = CircuitVoltage_result.GetValue(row-1, cir_pos+3);
                end
                if strcmp(CircuitVoltage_result.GetColumnName(cir_pos+4), Yphase_name(1)) == 1
                    ue(row,1) = CircuitVoltage_result.GetValue(row-1, cir_pos+4);
                end
            else
                if strcmp(CircuitVoltage_result.GetColumnName(cir_pos), Uphase_name(1)) == 1
                    ua(row,1) = CircuitVoltage_result.GetValue(row-1, cir_pos);
                end
                if strcmp(CircuitVoltage_result.GetColumnName(cir_pos+1), Vphase_name(1)) == 1
                    ub(row,1) = CircuitVoltage_result.GetValue(row-1, cir_pos+1);
                end
                if strcmp(CircuitVoltage_result.GetColumnName(cir_pos+2), Wphase_name(1)) == 1
                    uc(row,1) = CircuitVoltage_result.GetValue(row-1, cir_pos+2);
                end
            end
        end
    end

    % '------------------------------------------------------------------------
    %% Results Graphs of Torque
    electromagnetic_torque = zeros(NTL,1);
    torque_result = study.GetDataSet('Torque',1);
    Data_Manager.CreateAllCasesGraphModel(torque_result,'Torque Graph');
    % #Returns the number of rows in the data set as an integer.
    tempo = zeros(NTL,1);
    for row = 1:NTL
        tempo(row) = torque_result.GetValue(row-1, 0);  %time in s
        electromagnetic_torque(row) = torque_result.GetValue(row-1, 1);
    end
    
    % '------------------------------------------------------------------------
    %% Iron Loss
    % % --> Standard Calculation Method
    IronLoss_result = study.GetDataSet('Iron Loss (Iron loss)',1);
    Data_Manager.CreateGraphModel(IronLoss_result, 'Iron Core Losses (W)');
    if strcmp(IronLoss_result.GetColumnName(2), 'Stator Core') == 1
        IronLoss_stator(1,1) = IronLoss_result.GetValue(0, 2);
    end
    if strcmp(IronLoss_result.GetColumnName(1), 'Rotor Core') == 1
        IronLoss_rotor(1,1) = IronLoss_result.GetValue(0, 1);
    end
    % '------------------------------------------------------------------------
    %% Hystersis Loss
    % % --> Standard Calculation Method
    HystersisLoss_result = study.GetDataSet('Hysteresis Loss (Iron loss)',1);
    Data_Manager.CreateGraphModel(HystersisLoss_result, 'Hysteresis Loss (W)');
    if strcmp(HystersisLoss_result.GetColumnName(2), 'Stator Core') == 1
        HystersisLoss_stator(1,1) = HystersisLoss_result.GetValue(0, 2);
    end
    if strcmp(HystersisLoss_result.GetColumnName(1), 'Rotor Core') == 1
        HystersisLoss_rotor(1,1) = HystersisLoss_result.GetValue(0, 1);
    end
    % '------------------------------------------------------------------------
    %% Eddy Current Loss
    % % --> Standard Calculation Method
    EddyLoss_result = study.GetDataSet('Joule Loss (Iron loss)',1);
    Data_Manager.CreateGraphModel(EddyLoss_result, 'Eddy Losses (W)');
    if strcmp(EddyLoss_result.GetColumnName(2), 'Stator Core') == 1
        EddyLoss_stator(1,1) = EddyLoss_result.GetValue(0, 2);
    end
    if strcmp(EddyLoss_result.GetColumnName(1), 'Rotor Core') == 1
        EddyLoss_rotor(1,1) = EddyLoss_result.GetValue(0, 1);
    end
    % '------------------------------------------------------------------------
    %% PM Losses
    Loss_PM = zeros(NTL,length(MagnetCG));
    
    PMLoss_result = study.GetDataSet('Joule Loss',1);
    Data_Manager.CreateGraphModel(PMLoss_result, 'PM Losses (W)');
    if any(unique(Magnetxy(:,9)) ~= 0)
        for ID = 1:1:length(MagnetCG)
            for row = 1:NTL
                if strcmp(PMLoss_result.GetColumnName(ID), strcat('Magnet',num2str(ID))) == 1
                    Loss_PM(row,ID) = PMLoss_result.GetValue(row-1, ID);
                end
            end
        end
    end
    
    if(flagCharger)
        theta = rotorPos*geo.p*(pi/180);
    else
        % '------------------------------------------------------------------------
        %% Results Graphs of rotor electrical position in deg
        pos_Mov = tempo.*per.EvalSpeed*2*pi/60; % rotor position (rad mec)
        theta = ( pos_Mov) .* geo.p; 
    end
    
    %% Results Graphs of d-q Flux Linkages
    % Fluxq = 2/3*(+(Fluxa-0.5*Fluxb-0.5*Fluxc).*sin(theta)-sqrt(3)/2*(Fluxb-Fluxc).*cos(theta));
    % Fluxd = -2/3*((Fluxa-0.5*Fluxb-0.5*Fluxc).*cos(theta)+sqrt(3)/2*(Fluxb-Fluxc).*sin(theta));
    
    if(flagCharger)
        if(isnan(n3ph))
            for i=1:NTL
                Fluxdq(:,i) = abcde2dqd3q30([Fluxa(i);Fluxb(i);Fluxc(i);Fluxd(i);Fluxe(i)],theta);
            end
            Fluxd1 = Fluxdq(1,:);
            Fluxq1 = Fluxdq(2,:);
            Fluxd3 = Fluxdq(3,:);
            Fluxq3 = Fluxdq(4,:);
            Flux0 = Fluxdq(5,:);
        else    
            for i=1:NTL
                Fluxdq(:,i) = abc2dq0_g(n3ph, Fluxa(i,:),Fluxb(i,:),Fluxc(i,:),theta,win_delta);
            end
            for i=1:n3ph
                Fluxd(:,i) = Fluxdq(3*i-2,:);
                Fluxq(:,i) = Fluxdq(3*i-1,:);
                Flux0(:,i) = Fluxdq(3*i,:);
            end
        end
    else
        if(isnan(n3ph))
            for i=1:NTL
                Fluxdq(:,i) = abcde2dqd3q30([Fluxa(i);Fluxb(i);Fluxc(i);Fluxd(i);Fluxe(i)],theta);
            end
            Fluxd1 = Fluxdq(1,:);
            Fluxq1 = Fluxdq(2,:);
            Fluxd3 = Fluxdq(3,:);
            Fluxq3 = Fluxdq(4,:);
        else    
            [Fluxdq] = abc2dq(Fluxa',Fluxb',Fluxc',theta');
            Fluxd = Fluxdq(1,:);
            Fluxq = Fluxdq(2,:);
        end
    end
    % '------------------------------------------------------------------------
    %% Results Graphs of d-q Circuit Currents
    % iq = 2/3*(+(ia-0.5*ib-0.5*ic).*sin(theta)-sqrt(3)/2*(ib-ic).*cos(theta));
    % id = -2/3*((ia-0.5*ib-0.5*ic).*cos(theta)+sqrt(3)/2*(ib-ic).*sin(theta));
    
    if(flagCharger)
         if(isnan(n3ph))
            for i=1:NTL
                idq(:,i) = abcde2dqd3q30([ia(i);ib(i);ic(i);id(i);ie(i)],theta);
            end
            id1 = idq(1,:);
            iq1 = idq(2,:);
            id3 = idq(3,:);
            iq3 = idq(4,:);
            i0 = idq(5,:);
         else
            for i=1:NTL
                idq(:,i) = abc2dq0_g(n3ph, ia(i,:),ib(i,:),ic(i,:),theta,win_delta);
            end
            for i=1:n3ph
                id(:,i) = idq(3*i-2,:);
                iq(:,i) = idq(3*i-1,:);
                i0(:,i) = idq(3*i,:);
            end
        end     
    else
        if(isnan(n3ph))
            for i=1:NTL
                idq(:,i) = abcde2dqd3q30([ia(i);ib(i);ic(i);id(i);ie(i)],theta);
            end
            id1 = idq(1,:);
            iq1 = idq(2,:);
            id3 = idq(3,:);
            iq3 = idq(4,:);
        else
            [idq] = abc2dq(ia',ib',ic',theta');
            id = idq(1,:);
            iq = idq(2,:);
        end
    end
    
    if(flagCharger)
        %% output
        % '------------------------------------------------------------------------
        %electrical angle in deg
        th = tempo*2*pi*50*(180/pi); 
    else
        %% output
        % '------------------------------------------------------------------------
        %electrical angle in deg
        th = tempo.*per.EvalSpeed*(180/pi)*geo.p*2*pi/60; 
    end
    % flux, currents, resitance
    SOL.th = th';
    
    if(flagCharger)
        if(isnan(n3ph))
            SOL.fa = Fluxa';
            SOL.fb = Fluxb';
            SOL.fc = Fluxc';
            SOL.fd = Fluxd';
            SOL.fe = Fluxe';
            SOL.ia = ia';
            SOL.ib = ib';
            SOL.ic = ic';
            SOL.id = id';
            SOL.ie = ie';
            SOL.ua = ua';
            SOL.ub = ub';
            SOL.uc = uc';
            SOL.ud = ud';
            SOL.ue = ue';
            SOL.id1 = id1;
            SOL.iq1 = iq1;
            SOL.id3 = id3;
            SOL.iq3 = iq3;
            SOL.i0 = i0;
            SOL.fd1 = Fluxd1;
            SOL.fq1 = Fluxq1;
            SOL.fd3 = Fluxd3;
            SOL.fq3 = Fluxq3;
            SOL.f0 = Flux0;
        
            if(strcmp(study_name,'PWM'))
                
                T(1) = mean(electromagnetic_torque(1));
    
                endp = 1;
        
                for i=1:CalcDivision
                    inp = endp+1;
                    endp = inp -1 + pwm_pts(i);
                    irip(1,i) = abs(max(ia(inp:endp)) - min(ia(inp:endp)));
                    irip(2,i) = abs(max(ib(inp:endp)) - min(ib(inp:endp)));
                    irip(3,i) = abs(max(ic(inp:endp)) - min(ic(inp:endp)));
                    irip(4,i) = abs(max(id(inp:endp)) - min(id(inp:endp)));
                    irip(5,i) = abs(max(ie(inp:endp)) - min(ie(inp:endp)));
                    T(i+1) = mean(electromagnetic_torque(inp:endp));
                end
        
                SOL.Tpwm = T;
                SOL.irip = irip;
                SOL.tpwm = cell2mat(parray_data(:,1));
                % SOL.tsin = linspace(0,1/(50*CalcDivision),CalcDivision+1);
                SOL.tsin = linspace(0,1/(50),CalcDivision+1);
            end
        else
            SOL.fa = Fluxa;
            SOL.fb = Fluxb;
            SOL.fc = Fluxc;

            SOL.ia = ia;
            SOL.ib = ib;
            SOL.ic = ic;

            SOL.ua = ua;
            SOL.ub = ub;
            SOL.uc = uc;

            SOL.id = id;
            SOL.iq = iq;
            SOL.i0 = i0;

            SOL.fd = Fluxd;
            SOL.fq = Fluxq;
            SOL.f0 = Flux0;

            if(strcmp(study_name,'PWM'))
                
                T(1) = mean(electromagnetic_torque(1));
    
                endp = 1;
        
                for i=1:CalcDivision
                    inp = endp+1;
                    endp = inp -1 + pwm_pts(i);

                    for j=1:1:n3ph
                        irip(3*j-2,i) = abs(max(ia(inp:endp,j)) - min(ia(inp:endp,j)));
                        irip(3*j-1,i) = abs(max(ib(inp:endp,j)) - min(ib(inp:endp,j)));
                        irip(3*j,i) = abs(max(ic(inp:endp,j)) - min(ic(inp:endp,j)));
                    end

                    T(i+1) = mean(electromagnetic_torque(inp:endp));
                end
        
                SOL.Tpwm = T;
                SOL.irip = irip;
                SOL.tpwm = cell2mat(parray_data(:,1));
                SOL.tsin = linspace(0,1/50,CalcDivision+1);
            end
        end
    else
        if(isnan(n3ph))
            SOL.fa = Fluxa';
            SOL.fb = Fluxb';
            SOL.fc = Fluxc';
            SOL.fd = Fluxd';
            SOL.fe = Fluxe';
            SOL.ia = ia';
            SOL.ib = ib';
            SOL.ic = ic';
            SOL.id = id';
            SOL.ie = ie';
            SOL.ua = ua';
            SOL.ub = ub';
            SOL.uc = uc';
            SOL.ud = ud';
            SOL.ue = ue';
            SOL.id1 = id1;
            SOL.iq1 = iq1;
            SOL.id3 = id3;
            SOL.iq3 = iq3;
            SOL.fd1 = Fluxd1;
            SOL.fq1 = Fluxq1;
            SOL.fd3 = Fluxd3;
            SOL.fq3 = Fluxq3;
        
            if(strcmp(study_name,'PWM'))
                
                T(1) = mean(electromagnetic_torque(1));
    
                endp = 1;
        
                for i=1:CalcDivision
                    inp = endp+1;
                    endp = inp -1 + pwm_pts(i);
                    irip(1,i) = abs(max(ia(inp:endp)) - min(ia(inp:endp)));
                    irip(2,i) = abs(max(ib(inp:endp)) - min(ib(inp:endp)));
                    irip(3,i) = abs(max(ic(inp:endp)) - min(ic(inp:endp)));
                    irip(4,i) = abs(max(id(inp:endp)) - min(id(inp:endp)));
                    irip(5,i) = abs(max(ie(inp:endp)) - min(ie(inp:endp)));
                    T(i+1) = mean(electromagnetic_torque(inp:endp));
                end
        
                SOL.Tpwm = T;
                SOL.irip = irip;
                SOL.tpwm = cell2mat(parray_data(:,1));
                SOL.tsin = linspace(0,1/(50*CalcDivision),CalcDivision+1);
            end
        else
            SOL.fa = Fluxa';
            SOL.fb = Fluxb';
            SOL.fc = Fluxc';   
            SOL.ia = ia';
            SOL.ib = ib';
            SOL.ic = ic';
            SOL.ua = ua';
            SOL.ub = ub';
            SOL.uc = uc';
            SOL.id = id;
            SOL.iq = iq;
            SOL.fd = Fluxd;
            SOL.fq = Fluxq;
        end
    end

    % torque
    SOL.T = electromagnetic_torque';

    % losses
    SOL.Pfes = IronLoss_stator;
    SOL.Pfer = IronLoss_rotor;
    SOL.Ppm  = Loss_PM;
    SOL.Pfes_h = HystersisLoss_stator; 
    SOL.Pfer_h = HystersisLoss_rotor; 
    SOL.Pfes_c = EddyLoss_stator;
    SOL.Pfer_c = EddyLoss_rotor;

    SOL_old = SOL;
end

if(flagCharger)
    [parms, perf] = param_charger_JMAG(geo, per, theta, win_delta, th0, rotorPos, SOL, pathname, filename, JDesigner);
    SOL.parms{curr_pos} = parms;
    SOL.perf{curr_pos} = perf;
end

if strcmp(geo.RotType,'Hybrid')
    dispHybridResults(study_name, study, JDesigner);
end

% % Save and Exits JMAG Designer
JDesigner.Save();
JDesigner.Quit();
