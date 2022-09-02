%% OpenSave.m

%%%%%%%%%%%%%%%%%%%%%% Load existing motor %%%%%%%%%%%%%%%%%%%%%%

[filemot, pathname, fltidx]=uigetfile('motorExamples\mot_01.fem', 'Pick a motor');
load(strrep(filemot,'.fem','.mat'));
[dataSet,geo,per] = back_compatibility(dataSet,geo,per,1);

%%%%%%%%%%%%%%%%%%%%%% Edit Main Data %%%%%%%%%%%%%%%%%%%%%%

switch dataSet.TypeOfRotor
    case 'Circular'
        Rotor = 1;
    case 'Seg'
        Rotor = 2;
    case 'ISeg'
        Rotor = 3;
    case 'Fluid'
        Rotor = 4;
    case 'SPM'
        Rotor = 5;
end

motor_type=listdlg('ListString', {'Circular', 'Seg', 'ISeg','Fluid','SPM'}, ...
    'ListSize',[250 150], ...
    'PromptString', 'Select motor type:', ...
    'InitialValue', Rotor, ...
    'SelectionMode', 'Single');

switch motor_type
    case 1
        dataSet.TypeOfRotor = 'Circular';
    case 2
        dataSet.TypeOfRotor = 'Seg';
    case 3
        dataSet.TypeOfRotor = 'ISeg';
    case 4
        dataSet.TypeOfRotor = 'Fluid';
    case 5
        dataSet.TypeOfRotor = 'SPM';
end

if motor_type == 5     %SPM
    if isoctave()
        parameters = inputdlg({'Number of pole pairs';'Number of slot per pole per phase';'Airgap thichness [mm]';
            'Stator outer radius [mm]';'Rotor outer radius [mm]';'Shaft radius [mm]';
            'Stack length [mm]'},'Edit main parameters',1, {dataSet.NumOfPolePairs;
            dataSet.NumOfSlots; dataSet.AirGapThickness; dataSet.StatorOuterRadius;
            dataSet.AirGapRadius; dataSet.ShaftRadius; dataSet.StackLength});
    else
        parameters = inputdlg({'Number of pole pairs';'Number of slot per pole per phase';'Airgap thichness [mm]';
            'Stator outer radius [mm]';'Rotor outer radius [mm]';'Shaft radius [mm]';
            'Stack length [mm]'},'Edit main parameters',[1 50; 1 50; 1 50; 1 50; 1 50; 1 50; 1 50],{num2str(dataSet.NumOfPolePairs);
            num2str(dataSet.NumOfSlots); num2str(dataSet.AirGapThickness); num2str(dataSet.StatorOuterRadius);
            num2str(dataSet.AirGapRadius); num2str(dataSet.ShaftRadius); num2str(dataSet.StackLength)});
    end
    dataSet.NumOfPolePairs    = eval(cell2mat(parameters(1)));        % number of pole pairs
    dataSet.NumOfSlots        = eval(cell2mat(parameters(2)));        % number of slot per pole per phase
    dataSet.AirGapThickness   = eval(cell2mat(parameters(3)));        % airgap thichness [mm]
    dataSet.StatorOuterRadius = eval(cell2mat(parameters(4)));        % stator outer radius [mm]
    dataSet.AirGapRadius      = eval(cell2mat(parameters(5)));        % rotor outer radius [mm]
    dataSet.ShaftRadius       = eval(cell2mat(parameters(6)));        % shaft radius [mm]
    dataSet.StackLength       = eval(cell2mat(parameters(7)));        % stack length [mm]
    
else
    if isoctave()
        parameters = inputdlg({'Number of pole pairs';'Number of slot per pole per phase';'Airgap thichness [mm]';
            'Stator outer radius [mm]';'Rotor outer radius [mm]';'Shaft radius [mm]';
            'Stack length [mm]';'Number of layers'},'Edit main parameters',1, {dataSet.NumOfPolePairs;
            dataSet.NumOfSlots; dataSet.AirGapThickness; dataSet.StatorOuterRadius;
            dataSet.AirGapRadius; dataSet.ShaftRadius; dataSet.StackLength; dataSet.NumOfLayers});
    else
        parameters = inputdlg({'Number of pole pairs';'Number of slot per pole per phase';'Airgap thichness [mm]';
            'Stator outer radius [mm]';'Rotor outer radius [mm]';'Shaft radius [mm]';
            'Stack length [mm]';'Number of layers'},'Edit main parameters',[1 50; 1 50; 1 50; 1 50; 1 50; 1 50; 1 50; 1 50],{num2str(dataSet.NumOfPolePairs);
            num2str(dataSet.NumOfSlots); num2str(dataSet.AirGapThickness); num2str(dataSet.StatorOuterRadius);
            num2str(dataSet.AirGapRadius); num2str(dataSet.ShaftRadius); num2str(dataSet.StackLength); num2str(dataSet.NumOfLayers)});
    end
    
    dataSet.NumOfPolePairs    = eval(cell2mat(parameters(1)));        % number of pole pairs
    dataSet.NumOfSlots        = eval(cell2mat(parameters(2)));        % number of slot per pole per phase
    dataSet.AirGapThickness   = eval(cell2mat(parameters(3)));        % airgap thichness [mm]
    dataSet.StatorOuterRadius = eval(cell2mat(parameters(4)));        % stator outer radius [mm]
    dataSet.AirGapRadius      = eval(cell2mat(parameters(5)));        % rotor outer radius [mm]
    dataSet.ShaftRadius       = eval(cell2mat(parameters(6)));        % shaft radius [mm]
    dataSet.StackLength       = eval(cell2mat(parameters(7)));        % stack length [mm]
    dataSet.NumOfLayers       = eval(cell2mat(parameters(8)));        % number of layers
    
end


%%%%%%%%%%%%%%%%%%%%%% Edit Boundary Data %%%%%%%%%%%%%%%%%%%%%%

if motor_type == 1 || motor_type == 2 || motor_type == 4     %Circ, Seg and Fluid
    if isoctave()
        parameters_bound = inputdlg({'Outer barrier angle (lower limit) [PU]';'Outer barrier angle (upper limit) [PU]';
            'Other barrier angle (lower limit) [PU]'; 'Other barrier angle (upper limit) [PU]';
            'Flux barrier width (lower limit) [PU]'; 'Flux barrier width (upper limit) [PU]';
            'Flux barrier translation (lower limit) [PU]'; 'Flux barrier translation (upper limit) [PU]'; 'Current angle phase (lower limit) [electr. deg.]'; 'Current angle phase (upper limit) [electr. deg.]'}, 'Edit boundary parameters',1, {dataSet.Alpha1Bou(1);
            dataSet.Alpha1Bou(2); dataSet.DeltaAlphaBou(1); dataSet.DeltaAlphaBou(2);
            dataSet.hcBou(1); dataSet.hcBou(2); dataSet.DfeBou(1); dataSet.DfeBou(2); dataSet.PhaseAngleCurrBou(1); dataSet.PhaseAngleCurrBou(2)});
    else
        parameters_bound = inputdlg({'Outer barrier angle (lower limit) [PU]';'Outer barrier angle (upper limit) [PU]';
            'Other barrier angle (lower limit) [PU]'; 'Other barrier angle (upper limit) [PU]';
            'Flux barrier width (lower limit) [PU]'; 'Flux barrier width (upper limit) [PU]';
            'Flux barrier translation (lower limit) [PU]'; 'Flux barrier translation (upper limit) [PU]'; 'Current angle phase (lower limit) [electr. deg.]'; 'Current angle phase (upper limit) [electr. deg.]'}, 'Edit boundary parameters',[1 60; 1 60; 1 60; 1 60; 1 60; 1 60; 1 60; 1 60; 1 60; 1 60], {num2str(dataSet.Alpha1Bou(1));
            num2str(dataSet.Alpha1Bou(2)); num2str(dataSet.DeltaAlphaBou(1)); num2str(dataSet.DeltaAlphaBou(2));
            num2str(dataSet.hcBou(1)); num2str(dataSet.hcBou(2)); num2str(dataSet.DfeBou(1)); num2str(dataSet.DfeBou(2)); num2str(dataSet.PhaseAngleCurrBou(1)); num2str(dataSet.PhaseAngleCurrBou(2))});
    end
    dataSet.Alpha1Bou         = [eval(cell2mat(parameters_bound(1))) eval(cell2mat(parameters_bound(2)))];    % first (outer) barrier angle
    dataSet.DeltaAlphaBou     = [eval(cell2mat(parameters_bound(3))) eval(cell2mat(parameters_bound(4)))];    % other barrier angle
    dataSet.hcBou             = [eval(cell2mat(parameters_bound(5))) eval(cell2mat(parameters_bound(6)))];    % flux barrier width
    dataSet.DfeBou            = [eval(cell2mat(parameters_bound(7))) eval(cell2mat(parameters_bound(8)))];    % flux barrier translation
    dataSet.PhaseAngleCurrBou = [eval(cell2mat(parameters_bound(9))) eval(cell2mat(parameters_bound(10)))];   % current phase angle
    dataSet.Dalpha1BouCheck = 1;
    dataSet.DalphaBouCheck  = 1;
    dataSet.hcBouCheck      = 1;
    dataSet.DxBouCheck      = 1;
    dataSet.GammaBouCheck    =1;
    
elseif motor_type == 3     %ISeg
    if isoctave()
        parameters_bound = inputdlg({'Outer barrier angle (lower limit) [PU]';'Outer barrier angle (upper limit) [PU]';
            'Other barrier angle (lower limit) [PU]'; 'Other barrier angle (upper limit) [PU]';
            'Flux barrier width (lower limit) [PU]'; 'Flux barrier width (upper limit) [PU]'; 'Current angle phase (lower limit) [electr. deg.]'; 'Current angle phase (upper limit) [electr. deg.]'}, 'Edit boundary parameters',1, {dataSet.Alpha1Bou(1);
            dataSet.Alpha1Bou(2); dataSet.DeltaAlphaBou(1); dataSet.DeltaAlphaBou(2);
            dataSet.hcBou(1); dataSet.hcBou(2); dataSet.PhaseAngleCurrBou(1); dataSet.PhaseAngleCurrBou(2)});
        
    else
        parameters_bound = inputdlg({'Outer barrier angle (lower limit) [PU]';'Outer barrier angle (upper limit) [PU]';
            'Other barrier angle (lower limit) [PU]'; 'Other barrier angle (upper limit) [PU]';
            'Flux barrier width (lower limit) [PU]'; 'Flux barrier width (upper limit) [PU]'; 'Current angle phase (lower limit) [electr. deg.]'; 'Current angle phase (upper limit) [electr. deg.]'}, 'Edit boundary parameters',[1 60; 1 60; 1 60; 1 60; 1 60; 1 60; 1 60; 1 60], {num2str(dataSet.Alpha1Bou(1));
            num2str(dataSet.Alpha1Bou(2)); num2str(dataSet.DeltaAlphaBou(1)); num2str(dataSet.DeltaAlphaBou(2));
            num2str(dataSet.hcBou(1)); num2str(dataSet.hcBou(2)); num2str(dataSet.PhaseAngleCurrBou(1)); num2str(dataSet.PhaseAngleCurrBou(2))});
    end
    dataSet.Alpha1Bou         = [eval(cell2mat(parameters_bound(1))) eval(cell2mat(parameters_bound(2)))];    % first (outer) barrier angle
    dataSet.DeltaAlphaBou     = [eval(cell2mat(parameters_bound(3))) eval(cell2mat(parameters_bound(4)))];    % other barrier angle
    dataSet.hcBou             = [eval(cell2mat(parameters_bound(5))) eval(cell2mat(parameters_bound(6)))];    % flux barrier width
    dataSet.PhaseAngleCurrBou = [eval(cell2mat(parameters_bound(7))) eval(cell2mat(parameters_bound(8)))];    % current phase angle
    dataSet.Dalpha1BouCheck = 1;
    dataSet.DalphaBouCheck  = 1;
    dataSet.hcBouCheck      = 1;
    dataSet.GammaBouCheck    =1;
    
elseif motor_type == 5     %SPM
    if isoctave()
        parameters_bound = inputdlg({'Flux barrier width (lower limit) [PU]'; 'Flux barrier width (upper limit) [PU]'; 'Current angle phase (lower limit) [electr. deg.]'; 'Current angle phase (upper limit) [electr. deg.]'}, 'Edit boundary parameters',1, {dataSet.hcBou(1); dataSet.hcBou(2); dataSet.PhaseAngleCurrBou(1); dataSet.PhaseAngleCurrBou(2)});
    else
        parameters_bound = inputdlg({'Flux barrier width (lower limit) [PU]'; 'Flux barrier width (upper limit) [PU]'; 'Current angle phase (lower limit) [electr. deg.]'; 'Current angle phase (upper limit) [electr. deg.]'}, 'Edit boundary parameters',[1 60; 1 60; 1 60; 1 60], {num2str(dataSet.hcBou(1)); num2str(dataSet.hcBou(2)); num2str(dataSet.PhaseAngleCurrBou(1)); num2str(dataSet.PhaseAngleCurrBou(2))});
    end
    dataSet.hcBou             = [eval(cell2mat(parameters_bound(1))) eval(cell2mat(parameters_bound(2)))];    % flux barrier width
    dataSet.PhaseAngleCurrBou = [eval(cell2mat(parameters_bound(3))) eval(cell2mat(parameters_bound(4)))];    % current phase angle
    
    dataSet.hcBouCheck      = 1;
    dataSet.GammaBouCheck    =1;
end


%%%%%%%%%%%%%%%%%%%%%% Edit Performance Data %%%%%%%%%%%%%%%%%%%%%%

if isoctave()
    parameters_per = inputdlg({'Minimum expected torque'; 'Max admitted torque ripple'}, 'Edit performance parameters',1, {dataSet.MinExpTorque; dataSet.MaxRippleTorque});
else
    parameters_per = inputdlg({'Minimum expected torque'; 'Max admitted torque ripple'}, 'Edit performance parameters',[1 60; 1 60], {num2str(dataSet.MinExpTorque); num2str(dataSet.MaxRippleTorque)});
end

dataSet.MinExpTorque      = eval(cell2mat(parameters_per(1)));      % minimum expected torque
dataSet.MaxRippleTorque   = eval(cell2mat(parameters_per(2)));      % max admitted torque ripple
dataSet.TorqueOptCheck = 1;
dataSet.TorRipOptCheck = 1;


%%%%%%%%%%%%%%%%%%%%%% Edit Optimization parameters%%%%%%%%%%%%%%%%%%%%%%

if isoctave()
    parameters_opt = inputdlg({'Max number of generation'; 'Population size'}, 'Edit optimization parameters',1, {dataSet.MaxGen; dataSet.XPop});
    
else
    parameters_opt = inputdlg({'Max number of generation'; 'Population size'}, 'Edit optimization parameters',[1 60; 1 60], {num2str(dataSet.MaxGen); num2str(dataSet.XPop)});
end

dataSet.MaxGen = eval(cell2mat(parameters_opt(1)));      % Max number of generations
dataSet.XPop   = eval(cell2mat(parameters_opt(2)));      % Population size


%%%%%%%%%%%%%%%%%%%%%% Other parameters%%%%%%%%%%%%%%%%%%%%%%

dataSet.ALPHApu = ones(1,dataSet.NumOfLayers)*round(1/(dataSet.NumOfLayers+0.5)*100)/100;
dataSet.HCpu = ones(1,dataSet.NumOfLayers)*0.5;
dataSet.HCpu = round(dataSet.HCpu .*100) ./100;
dataSet.DepthOfBarrier = zeros(1,dataSet.NumOfLayers);
if  strcmp(dataSet.TypeOfRotor,'Seg')||strcmp(dataSet.TypeOfRotor,'ISeg')  %mod walter
    dataSet.Areavert0=zeros(1,4);
    dataSet.Areaob0=zeros(1,4);
    dataSet.Areatot=zeros(1,4);
    dataSet.dob=ones(1,4);
    dataSet.dvert=ones(1,4);
end

clear motor_type parameters parameters_bound parameters_per parameters_opt filemot fltidx pathname


%%%%%%%%%%%%%%%%%%%%%% Save dialog %%%%%%%%%%%%%%%%%%%%%%

save_yn = questdlg('Save new motor?', '', 'Yes');
if (strcmp(save_yn, 'Yes'))    % Save and plot
    dataSet.currentfilename = inputdlg({'Insert motor name'}, ' ',1, {'newmachine'});
    
    [bounds,objs,geo,per,mat] = data0(dataSet);
    dataSet.RQ = buildDefaultRQ(bounds);
    [geo,gamma,mat] = interpretRQ(geo.RQ,geo,mat);
    
%     openfemm(1)
    eval_type='MO_OA';
    [geo,mat] = draw_motor_in_FEMM(geo,mat,[cd '\'],[char(dataSet.currentfilename) '.fem']);
%     mi_zoomnatural
%     mi_saveas(strcat(char(dataSet.currentfilename),'.fem'));
    mi_close, closefemm
    
    clear filemot fltidx pathname psCalc Q QsCalc save_yn t2 Rotor ans
    
    if isoctave()
        save ('-mat7-binary', strcat(char(dataSet.currentfilename),'.mat'));
    else
        save (strcat(char(dataSet.currentfilename),'.mat'));
    end
    
else     % Plot only

    [bounds,objs,geo,per,mat] = data0(dataSet);
    dataSet.RQ = buildDefaultRQ(bounds);
    [geo,gamma,mat] = interpretRQ(geo.RQ,geo,mat);
    
end

Mat = [geo.stator;geo.rotor];
figure(); h = gca;
GUI_Plot_Machine(h,Mat);

clear all
