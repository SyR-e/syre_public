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


%'----------------------------Open JMAG Designer---------------------------
 JMAGversion='222'; 
% '# Create an "app" application object to launch JMAG-Designer
JDesigner = actxserver(strcat('designer.Application.',JMAGversion));
JDesigner.Show(); % Show the JMAG interface
% '------------------------------------------------------------------------
% 'Motor Winding diagram definitions :
Ncoilseries = geo.win.Ns;%number of series turns
coilTurn = Ncoilseries/((geo.Qs/3)*(geo.p*2));%stator coil number set in external circuit
coilRes = per.Rs;% 'Stator coil resistance in ohms
coilLind = geo.win.Lend;% 'Stator coil LeakageInductance in Henry
phase_order = 1;% 'Phase order: 1--> UWV and 0 --> UVW
% '------------------------------------------------------------------------
% 'Motion Setting definitions :
th0 = (-geo.th0)/(geo.p);% 'Initial position of rotor in mechanical degrees
% '------------------------------------------------------------------------
% 'Study property settings: Step Control and Resolution
% 'Resolution or number of Step division per cycle of calculations
CalcDivision = per.nsim_singt*geo.p*2; 
%'#Number of electric periods to be calculed
% ncylce=per.delta_sim_singt*geo.p*2/360;
ncylce = 1;
% '------------------------------------------------------------------------
MagnetCG=geo.BLKLABELS.rotore.xy(geo.BLKLABELS.rotore.xy(:,3)==6,1:2);%center of gravity(G) of Magnet
Magnetxy=geo.rotor(geo.rotor(:,end-1)==6,:);
% '________________________________________________________________________
% '# Open jproj (Motor Drawing File)
JDesigner.NewProject('Untitled');
JDesigner.Load (fullfile(pathname, strcat(strrep(filename,'.mat','.jmag'),'.jproj')));
study = JDesigner.GetCurrentStudy();
% '------------------------------------------------------------------------
% 'Create the Design Parameters:
Design_parameters = study.GetDesignTable();

% 'initial position(motion) of the rotor in degree
createDesignParameters_JMAG(Design_parameters,'InitPosition',0,th0,'Initial position(motion) of the rotor (degrees)')
% 'rotation speed in rpm
createDesignParameters_JMAG(Design_parameters,'rspeed',0,per.EvalSpeed,'rotation speed (rpm)')
% 'Stator coil turn number
createDesignParameters_JMAG(Design_parameters,'coilTurn',0,coilTurn,'Coil turn number')
% 'Stator coil resistance in ohms
createDesignParameters_JMAG(Design_parameters,'coilRes',0,coilRes,'Coil resistance (ohms)')
% 'Stator coil LeakageInductance in Henry
createDesignParameters_JMAG(Design_parameters,'coilLind',0,coilLind,'coil LeakageInductance (Henry)')
% % 'Aumont of drive in volt for voltage type and in ampere for current type for the stator coils
createDesignParameters_JMAG(Design_parameters,'Aumont_Source',0,per.i0*per.overload,'Amplitude of motor drive (A or V)')
% 'Phase angle of drive for voltage type or current type for the stator coils
createDesignParameters_JMAG(Design_parameters,'Phase_Source',0,per.gamma+90,'Phase angle of motor drive (degrees)')
% 'Frequency of drive for voltage type or current type for the stator coils
createDesignParameters_JMAG(Design_parameters,'Freq_Source',0,per.EvalSpeed*geo.p/60,'Frequncy of motor drive (Hz)')
% 'Resolution or number of Step division per cycle of calculations
createDesignParameters_JMAG(Design_parameters,'CalcDivision',0,CalcDivision,'Resolution or Step division per cycle of calculations')
%'#Number of electric periods to be calculed
createDesignParameters_JMAG(Design_parameters,'ncylce',0,ncylce,'Number of electric periods to be calculed')

%% 'Study property settings
% '------------------------------------------------------------------------
%'Calculation time for one electric period / poles:
CalcTime='ncylce/(rspeed*PolePair/60)';   % in second
% ''#Number of calculation steps for Step Control
CalcStep='CalcDivision+1';% in order to have ncycle of electric period

% '# 'Step Control settings
StepControl = study.GetStep();
StepControl.SetValue('Step', strcat(CalcStep));
StepControl.SetValue('StepType', 1); % 0 = 'Simple', 1 = 'Uniform', 2 = 'Time'
StepControl.SetValue('StepDivision ', 'CalcDivision');
StepControl.SetValue('EndPoint', strcat(CalcTime)); %for 'StepType' = 1 

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
%#create the FEM Coil components
%#FEM Coils U, V, W
Gu(1,1)=X0; Gu(1,2)=Y0;
Gv(1,1)=X0; Gv(1,2)=Y0-duvw;
Gw(1,1)=X0; Gw(1,2)=Y0-2*duvw;
phaseLabel = {'U', 'V', 'W'}; % or any other label
positionU = [Gu(1,1), Gu(1,2)];
positionV = [Gv(1,1), Gv(1,2)];
positionW = [Gw(1,1), Gw(1,2)];
FEMCoil_Creation_JMAG(circuit, phaseLabel{1}, positionU);
FEMCoil_Creation_JMAG(circuit, phaseLabel{2}, positionV);
FEMCoil_Creation_JMAG(circuit, phaseLabel{3}, positionW);


%# Add Electric source components into the circuit
dx2=Gv(1,1)-4*dx; dy2=Gv(1,2);
circuit.CreateComponent('3PhaseCurrentSource', 'CS1');
circuit.CreateInstance('CS1', dx2, dy2);
circuit.GetComponent('CS1').SetName(strcat('3Phase-','Current'));
circuit.GetComponent(strcat('3Phase-','Current')).SetValue('XType ', 'Time');
circuit.GetComponent(strcat('3Phase-','Current')).SetValue('CommutatingSequence', phase_order);
circuit.GetComponent(strcat('3Phase-','Current')).SetValue('Amplitude', 'Aumont_Source');
circuit.GetComponent(strcat('3Phase-','Current')).SetValue('Frequency', 'Freq_Source');
circuit.GetComponent(strcat('3Phase-','Current')).SetValue('PhaseU', 'Phase_Source');

%Add Voltage probe for U-Phase
AddVoltageProbe_JMAG (circuit, 1, positionU+[-dx, 2], phaseLabel{1})
%Add Voltage probe for V-Phase
AddVoltageProbe_JMAG (circuit, 2, positionV+[-dx, 2], phaseLabel{2})
%Add Voltage probe for W-Phase
AddVoltageProbe_JMAG (circuit, 3, positionW+[-dx, 2], phaseLabel{3})
%Add a Ground component
circuit.CreateComponent('Ground', 'Ground');
circuit.CreateInstance('Ground', Gv(1,1)+3*dx, Gv(1,2)-dy);

%#create star connection of FEM Coils
circuit.CreateWire(Gu(1,1)+2, Gu(1,2), Gv(1,1)+2, Gv(1,2));
circuit.CreateWire(Gv(1,1)+2, Gv(1,2), Gw(1,1)+2, Gw(1,2));

%Connect Electric source to the Coils.
circuit.CreateWire(dx2+2, dy2+2, Gu(1,1)-2, Gu(1,2));
circuit.CreateWire(dx2+2, dy2, Gv(1,1)-2, Gv(1,2));
circuit.CreateWire(dx2+2, dy2-2, Gw(1,1)-2, Gw(1,2));

%Connect Ground
circuit.CreateWire(Gv(1,1)+2, Gv(1,2), Gv(1,1)+3*dx, Gv(1,2)-dy+2);




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
Motion_Condition.AddSelected(sel);
JDesigner.SetCurrentStudy('Load');

% '#'_______________________________________________________________
% '# Apply Torque Condition
study.CreateCondition('Torque', 'Electromagnetic_Torque');
Torque_Condition=study.GetCondition('Electromagnetic_Torque');
Torque_Condition.SetValue('TargetType', 0);
Torque_Condition.ClearParts();
Torque_Condition.SetValue('TargetType', 1);
Torque_Condition.SetLinkWithType('LinkedMotion', 'Speed');

JDesigner.SetCurrentStudy('Load');
% '#'_______________________________________________________________
% '# Apply Iron loss condition on Stator Core
Iron_losscondition_JMAG (study,'Stator Core','Stator_Core_Loss',1,2,1)
% '#'_______________________________________________________________
% '# Apply Iron loss condition on Rotor Core
Iron_losscondition_JMAG (study,'Rotor Core','Rotor_Core_Loss',1,2,1)
% '#'_______________________________________________________________
% '##FEM Coil conditions for Winding Connections
% '_______________________________________________________________
% coil direction Definition: Upward or Downward
% coildirectionp = 0; % upward coil direction, coildirectionn = 1; % downward coil direction

% 'Create FEM Coil conditions for U-phase
Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{1},+1,-1)
% '_______________________________________________________________
% 'Create FEM Coil conditions for W-phase
Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{3},+3,-3)
% '_______________________________________________________________
% 'Create FEM Coil conditions for V-phase
Coil_DirectionDefinition_JMAG (geo,study,phaseLabel{2},+2,-2)


JDesigner.SetCurrentStudy('Load');

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
sel.SelectPart(strcat('Coil','_U'));
sel.SelectPart(strcat('Coil','_V'));
sel.SelectPart(strcat('Coil','_W'));
Coil_mesh.AddSelected(sel);

JDesigner.SetCurrentStudy('load');

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
Uphase_name = 'U-Phase Coil';Vphase_name = 'V-Phase Coil';Wphase_name = 'W-Phase Coil';
% Results Graphs of Coil Flux-Linkage
CoilFluxLinkage_result = Data_Manager.GetDataSet('Coil Flux-Linkage');
Data_Manager.CreateAllCasesGraphModel(CoilFluxLinkage_result,'Coil Flux-Linkage');
% '------------------------------------------------------------------------
%% Results Graphs of CoilFluxLinkages
NTL = CoilFluxLinkage_result.GetRows();

Fluxa = zeros(NTL,1);
Fluxb = zeros(NTL,1);
Fluxc = zeros(NTL,1);

for row =1:1:NTL
    if strcmp(CoilFluxLinkage_result.GetColumnName(1), Uphase_name) == 1
        Fluxa(row,1) = CoilFluxLinkage_result.GetValue(row-1, 1);
    end
    if strcmp(CoilFluxLinkage_result.GetColumnName(2), Vphase_name) == 1
        Fluxb(row,1) = CoilFluxLinkage_result.GetValue(row-1, 2);
    end
    if strcmp(CoilFluxLinkage_result.GetColumnName(3), Wphase_name) == 1
        Fluxc(row,1) = CoilFluxLinkage_result.GetValue(row-1, 3);
    end
end
% '------------------------------------------------------------------------
%% Results Graphs of Circuit Current
ia = zeros(NTL,1);
ib = zeros(NTL,1);
ic = zeros(NTL,1);

CircuitCurrent_result = Data_Manager.GetDataSet('Circuit Current');
Data_Manager.CreateAllCasesGraphModel(CircuitCurrent_result,'Circuit Current Graph');
for row =1:1:NTL
    if strcmp(CircuitCurrent_result.GetColumnName(2), Uphase_name) == 1
        ia(row,1) = CircuitCurrent_result.GetValue(row-1, 2);
    end
    if strcmp(CircuitCurrent_result.GetColumnName(3), Vphase_name) == 1
        ib(row,1) = CircuitCurrent_result.GetValue(row-1, 3);
    end
    if strcmp(CircuitCurrent_result.GetColumnName(4), Wphase_name) == 1
        ic(row,1) = CircuitCurrent_result.GetValue(row-1, 4);
    end
end
% '------------------------------------------------------------------------
%% Results Graphs of Torque
electromagnetic_torque = zeros(NTL,1);
torque_result = Data_Manager.GetDataSet('Torque');
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
IronLoss_result = Data_Manager.GetDataSet('Iron Loss (Iron loss)');
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
HystersisLoss_result = Data_Manager.GetDataSet('Hysteresis Loss (Iron loss)');
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
EddyLoss_result = Data_Manager.GetDataSet('Joule Loss (Iron loss)');
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

PMLoss_result = Data_Manager.GetDataSet('Joule Loss');
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
% '------------------------------------------------------------------------
%% Results Graphs of rotor electrical position in deg
pos_Mov = tempo.*per.EvalSpeed*2*pi/60; % rotor position (rad mec)
theta = ( pos_Mov) .* geo.p; 

%% Results Graphs of d-q Flux Linkages
% Fluxq = 2/3*(+(Fluxa-0.5*Fluxb-0.5*Fluxc).*sin(theta)-sqrt(3)/2*(Fluxb-Fluxc).*cos(theta));
% Fluxd = -2/3*((Fluxa-0.5*Fluxb-0.5*Fluxc).*cos(theta)+sqrt(3)/2*(Fluxb-Fluxc).*sin(theta));

[Fluxdq] = abc2dq(Fluxa',Fluxb',Fluxc',theta');
Fluxd = Fluxdq(1,:);
Fluxq = Fluxdq(2,:);
% '------------------------------------------------------------------------
%% Results Graphs of d-q Circuit Currents
% iq = 2/3*(+(ia-0.5*ib-0.5*ic).*sin(theta)-sqrt(3)/2*(ib-ic).*cos(theta));
% id = -2/3*((ia-0.5*ib-0.5*ic).*cos(theta)+sqrt(3)/2*(ib-ic).*sin(theta));

[idq] = abc2dq(ia',ib',ic',theta');
id = idq(1,:);
iq = idq(2,:);


%% output
% '------------------------------------------------------------------------
%electrical angle in deg
th = tempo.*per.EvalSpeed*(180/pi)*geo.p*2*pi/60; 
% flux, currents, resitance
SOL.th = th';
SOL.fa = Fluxa';
SOL.fb = Fluxb';
SOL.fc = Fluxc';   
SOL.ia = ia';
SOL.ib = ib';
SOL.ic = ic';
SOL.id = id;
SOL.iq = iq;
SOL.fd = Fluxd;
SOL.fq = Fluxq;

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
% % Save and Exits JMAG Designer
JDesigner.Save();
JDesigner.Quit();
end