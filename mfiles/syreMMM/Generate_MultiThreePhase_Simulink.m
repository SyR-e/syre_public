function Generate_MultiThreePhase_Simulink(motorModel,n_set)
    
    ctrlFolder_name = [motorModel.data.motorName '_ctrl_INST'];
    syrePath = fileparts(which('GUI_Syre.mlapp'));

    ctrlFolder_path = [motorModel.data.pathname ctrlFolder_name];

    copyfile([syrePath '\syreDrive\MultiThreePhase'], ctrlFolder_path);

    motorModel.data.T_VSD = ComputeVSD(n_set);
    
    save([ctrlFolder_path '\motorModel.mat'],'motorModel');

    MMM_print_MotorDataH(motorModel);

    PrintClarke(motorModel,n_set);
    ComputeDecouplingMatrix(motorModel,n_set);

    %% -------------------Generate Simulink model------------------------%%
   
    
    motorName = motorModel.data.motorName;
    Simulink_path = [ctrlFolder_path '\' motorModel.data.motorName '_ctrl_INST.slx'];
    Simulink_Name = [motorModel.data.motorName '_ctrl_INST'];
    new_model(Simulink_Name);
    
    if strcmp(motorModel.data.axisType,'SR') && strcmp(motorModel.data.motorType,'SR')
        Quad_Maps = 0; %SyR Convention - 1st quadrant maps
    elseif strcmp(motorModel.data.axisType,'SR') && strcmp(motorModel.data.motorType,'PM')
        Quad_Maps = 1; %PM-SyR - 1st and 4th quadrant maps
    elseif strcmp(motorModel.data.axisType,'PM') && strcmp(motorModel.data.motorType,'PM')
        Quad_Maps = 2; %IPM - 1st and 2st quadrant maps
    end


    path = getComponentsPaths();
    SimulinkRefBlocks = 'RefBlocks';
    open([syrePath '\syreDrive\MultiThreePhase\RefBlocks.slx']);
    
    
    print_control_script(ctrlFolder_path,n_set,Quad_Maps);
    add_control(SimulinkRefBlocks,Simulink_Name,path,n_set);
    add_inverters(Simulink_Name,SimulinkRefBlocks,n_set);
    create_circuital_model(Simulink_Name,path,SimulinkRefBlocks,n_set);

    close_system([syrePath '\syreDrive\MultiThreePhase\RefBlocks.slx']);
    %--------------------------Set Solver Configuration-------------------%
    cs = getActiveConfigSet(Simulink_Name);
    cs.set_param('StartTime', '0.0');   % Start time
    cs.set_param('StopTime', '1.0');   % Stop time
    cs.set_param('SolverType', 'Variable-step');   % Type
    cs.set_param('RelTol', '1e-3');   % Relative Tolerance
    cs.set_param('AbsTol', '1e-3');   % Absolute tolerance
    cs.set_param('SolverName', 'ode23tb');   
    cs.set_param('MaxStep', 'Tstep');
    
    Init = sprintf(" run init_sim.m; mex Motor_ctrl.c User_functions\\Src\\*.c");
    set_param(Simulink_Name,'InitFcn',Init);

    save_system(Simulink_path);
    close_system(Simulink_Name);
    delete([ctrlFolder_path  '\RefBlocks.slx']);

    clc
    disp('Simulink model created!');
    disp(['pathname:'])
    disp(['  ' ctrlFolder_path '\'])
    disp(['filename:'])
    disp(['  ' motorModel.data.motorName '_ctrl_INST.slx'])
    
end