function create_circuital_model(filepath,path,SimulinkRefBlocks,n_set)

%% ------------------------------------------------------------------------%
%                      ADD MOTOR MODEL SUBSYSTEM
%--------------------------------------------------------------------------%

% Add Motor Model Subsystem whose position is shifted respect to Inverter 1
% Port 1 of Motor Model Subsystem and Port 1 of Inverter 1 are taken for
% coordinates computations 

deltaX_motorModel = 450;
deltaY_motorModel = 50;

PlaceNewBlock(filepath,'Inverter Model 1',path.Subsystem,'Motor Model',deltaX_motorModel,deltaY_motorModel,1,1);
motor_model_path = [filepath '/Motor Model'];

% The motor model subsystem is enlarged respect to the number of sets 
width_motorModel = 50;
deltaY_motorModel= 150; 
dim(motor_model_path,0,0,width_motorModel,deltaY_motorModel*n_set);

% The Motor Model subsystem content is deleted 
Simulink.SubSystem.deleteContents([filepath '/Motor Model']);


%% ------------------------------------------------------------------------%
%             ADD CIRCUITAL AND MECHANICAL MODEL SUBSYSTEMS
%--------------------------------------------------------------------------%


% The Circuital Model Subsytem is added and enlarged inside the Motor Model
% Subsystem

circuit_model_path = [motor_model_path '/Circuital Model'];
add_block(path.Subsystem,circuit_model_path);
deltaY_circuitalModel = 150;
width_circuitalModel = 75;
dim(circuit_model_path,0,0,width_circuitalModel,deltaY_circuitalModel*n_set);

% Add the Mechanical equation subsystem respect to the Circuital Model
% Subsystem. The vertical displacemenet is computed between port 1 of both
% systems

deltaY_mechanicalEquations = deltaY_circuitalModel*n_set*1.2;
PlaceNewBlock(motor_model_path,'Circuital Model',[SimulinkRefBlocks '/Mechanical equations'],'Mechanical equations',width_circuitalModel*0.5,deltaY_mechanicalEquations,1,1);

% Add input and outports of the Motor Model connected to the mechanical
% Equations Subsytem
PlaceConnectNewBlock(motor_model_path,'Mechanical equations',path.Inport,'Tload',-50,0,1,1,'connect');
PlaceConnectNewBlock(motor_model_path,'Mechanical equations',path.Inport,'nload',-50,0,3,1,'connect');
PlaceConnectNewBlock(motor_model_path,'Mechanical equations',path.Outport,'n_m',50,0,4,1,'connect');
PlaceConnectNewBlock(motor_model_path,'Mechanical equations',path.Outport,'thetar',50,0,6,1,'connect');

% The circuital model content is deleted
Simulink.SubSystem.deleteContents(circuit_model_path);
ports_CCS = [];
% For cycle is exploited in order to copy a single phase n times
% For cycle is exploited to repeat the 3-phase circuits n_set times
% Blocks are added and placed with the function addpos and strcat

%% ------------------------------------------------------------------------%
%             FILL THE CIRCUITAL MODEL SUBSYSTEM
%--------------------------------------------------------------------------%


%% Add electrical inports of the Circuital Model Subsystem
phase = [];
for i = 1:n_set
    a = ['a' num2str(i)];
    b = ['b' num2str(i)];
    c = ['c' num2str(i)];
    phase = [phase; a; b; c];  % vector of electric ports names
    if i == 1
        %reference = 'a1';
        % Add a1, b1 refered to a1 and c1 refered to b1
        add_block(path.PMCport,[circuit_model_path '/' a]);
        PlaceNewBlock(circuit_model_path,a,path.PMCport,b,0,350,1,1);
        PlaceNewBlock(circuit_model_path,b,path.PMCport,c,0,350,1,1);
    else
        % Add ak refered to ck-1 
        reference = ['c' num2str(i-1)];
        PlaceNewBlock(circuit_model_path,reference,path.PMCport,a,0,350,1,1);
        PlaceNewBlock(circuit_model_path,a,path.PMCport,b,0,350,1,1);
        PlaceNewBlock(circuit_model_path,b,path.PMCport,c,0,350,1,1);
    end   
end

%% Add Electric Circuit inside the Circuital Model Subsystem 

for i = 1:(3*n_set)
    istring = num2str(i);
    PlaceNewBlock(circuit_model_path,phase(i,:),path.Resistor,['Ra',istring],200,0,1,1);
    PlaceNewBlock(circuit_model_path,phase(i,:),path.ControlledCurrentSource,['Controlled Current Source',istring],400,0,1,3,'left');
    PlaceNewBlock(circuit_model_path,['Controlled Current Source',istring],path.From,['From_',istring],-50,0,2,1);
    PlaceNewBlock(circuit_model_path,['Controlled Current Source',istring],path.SimulinkPSConverter,['Simulink-PS Converter',istring],-40,0,2,1);
    PlaceNewBlock(circuit_model_path,['Controlled Current Source',istring],path.Resistor,['Rea',istring],0,50,3,1);
    PlaceNewBlock(circuit_model_path,['Controlled Current Source',istring],path.VoltageSensor,['Ea sensor',istring],0,100,3,1,'right');
    PlaceNewBlock(circuit_model_path,['Controlled Current Source',istring],path.VoltageSensor,['Va sensor',istring],0,-100,3,1,'right');
    PlaceNewBlock(circuit_model_path,['Va sensor',istring],path.PSSimulinkConverter,['PS-Simulink Converter1',istring],50,0,2,1);
    PlaceNewBlock(circuit_model_path,['Va sensor',istring],path.Goto,['Goto1',istring],80,0,2,1);
    PlaceNewBlock(circuit_model_path,['Ea sensor',istring],path.PSSimulinkConverter,['PS-Simulink Converter2',istring],50,0,2,1);
    PlaceNewBlock(circuit_model_path,['Ea sensor',istring],path.Goto,['Goto2',istring],80,0,2,1);
    PlaceNewBlock(circuit_model_path,'Goto11',path.From,['From2',istring],200,i*50,1,1);
    if i <= 3
        PlaceNewBlock(circuit_model_path,['From2',istring],path.From,['FromV',istring],700,800,1,1);
        set_param([circuit_model_path,'/FromV',istring],'GotoTag',['V',phase(i,:)]);
    end

% The parameters of the blocks of each phase are set
    set_param([circuit_model_path,'/Ra',istring],'R','Rs');
    set_param([circuit_model_path,'/Rea',istring],'R','1e9');
    set_param([circuit_model_path,'/From_',istring],'GotoTag',['I',phase(i,:)]);
    set_param([circuit_model_path,'/Goto1',istring],'GotoTag',['V',phase(i,:)]);
    set_param([circuit_model_path,'/Goto2',istring],'GotoTag',['E',phase(i,:)]);
    set_param([circuit_model_path,'/From2',istring],'GotoTag',['E',phase(i,:)]);
%     set_param([circuit_model_path,'/From1',istring],'GotoTag',['V',phase(i,:)]);

% The data of the ports of the simscape blocks are collected to connect them
    ports_PMC = get_param([circuit_model_path,'/',phase(i,:)],'PortHandles');
    ports_Ra = get_param([circuit_model_path,'/Ra',istring],'PortHandles');
    ports_Rea = get_param([circuit_model_path,'/Rea',istring],'PortHandles');
    ports_Vasensor = get_param([circuit_model_path,'/Va sensor',istring],'PortHandles');
    ports_Easensor = get_param([circuit_model_path,'/Ea sensor',istring],'PortHandles');
    ports_ControlledCurrentSource = get_param([circuit_model_path,'/Controlled Current Source',istring],'PortHandles');
    ports_CCS = [ports_CCS; ports_ControlledCurrentSource];
    ports_SimulinkPSConverter = get_param([circuit_model_path,'/Simulink-PS Converter',istring],'PortHandles');
    ports_PSSimulinkConverter1 = get_param([circuit_model_path,'/PS-Simulink Converter1',istring],'PortHandles');
    ports_PSSimulinkConverter2 = get_param([circuit_model_path,'/PS-Simulink Converter2',istring],'PortHandles');       
    
    % Connect the electrical ports between them 
    add_line(circuit_model_path,ports_PMC.RConn,ports_Vasensor.LConn,'autorouting','smart');
    add_line(circuit_model_path,ports_PMC.RConn,ports_Ra.LConn,'autorouting','smart');
    add_line(circuit_model_path,ports_Ra.RConn,ports_ControlledCurrentSource.RConn(2),'autorouting','smart');
    add_line(circuit_model_path,ports_Ra.RConn,ports_Rea.LConn,'autorouting','smart');
    add_line(circuit_model_path,ports_Ra.RConn,ports_Easensor.LConn,'autorouting','smart');
    add_line(circuit_model_path,ports_Vasensor.RConn(1),ports_PSSimulinkConverter1.LConn,'autorouting','smart');
    add_line(circuit_model_path,ports_Easensor.RConn(1),ports_PSSimulinkConverter2.LConn,'autorouting','smart');
    add_line(circuit_model_path,ports_SimulinkPSConverter.RConn,ports_ControlledCurrentSource.RConn(1),'autorouting','smart');
    add_line(circuit_model_path,ports_Vasensor.RConn(2),ports_ControlledCurrentSource.LConn,'autorouting','smart');
    add_line(circuit_model_path,ports_Rea.RConn,ports_ControlledCurrentSource.LConn,'autorouting','smart');
    add_line(circuit_model_path,ports_Easensor.RConn(2),ports_ControlledCurrentSource.LConn,'autorouting','smart');

    % simulink blocks are connected differentyly and easily as always
    add_line(circuit_model_path,['From_',istring,'/1'],['Simulink-PS Converter',istring,'/1']);
    add_line(circuit_model_path,['PS-Simulink Converter1',istring,'/1'],['Goto1',istring,'/1']);
    add_line(circuit_model_path,['PS-Simulink Converter2',istring,'/1'],['Goto2',istring,'/1']);
    if (i/3) == ceil(i/3)
        add_line(circuit_model_path,ports_CCS(i-2).LConn,ports_CCS(i-1).LConn,'autorouting','smart');
        add_line(circuit_model_path,ports_CCS(i).LConn,ports_CCS(i-1).LConn,'autorouting','smart');
    end
end

%% Add from and goto references

% Create Voltage mux for the first set
PlaceConnectNewBlock(circuit_model_path,'From21',path.Mux,'MuxE',50,0,1,1,'connect');
PlaceConnectNewBlock(circuit_model_path,'FromV1',path.Mux,'MuxV',50,40,1,1,'connect');
set_param([circuit_model_path,'/MuxE'],'Inputs',num2str(3*n_set));
dim([circuit_model_path,'/MuxE'],0,0,0,n_set*130);
set_param([circuit_model_path,'/MuxV'],'Inputs','3');

% Create Voltage mux for the other sets
for i = 2:(n_set*3)
    add_line(circuit_model_path,['From2' num2str(i) '/1'],['MuxE/' num2str(i)],'autorouting','smart');
    if i <= 3
        add_line(circuit_model_path,['FromV' num2str(i) '/1'],['MuxV/' num2str(i)],'autorouting','smart');
    end
end

% Add Tm tag
PlaceNewBlock(circuit_model_path,'FromV2',path.From,'FromTm',0,-250-50*n_set,1,1);
set_param([circuit_model_path,'/FromTm'],'GotoTag','Tm');
PlaceConnectNewBlock(circuit_model_path,'FromTm',path.Outport,'Tm',100,0,1,1,'connect');

% Add Iabc tag respect to Tm tag
PlaceNewBlock(circuit_model_path,'FromTm',path.From,'FromIabc1',0,50,1,1);
set_param([circuit_model_path,'/FromIabc1'],'GotoTag','Iabc1');
PlaceConnectNewBlock(circuit_model_path,'FromIabc1',path.Outport,'Iabc1',100,0,1,1,'connect');

% Add Fdq tag respect to Tm tag
PlaceNewBlock(circuit_model_path,'FromTm',path.From,'From_Fdq',0,100,1,1);
set_param([circuit_model_path,'/From_Fdq'],'GotoTag','Fdq');
PlaceConnectNewBlock(circuit_model_path,'From_Fdq',path.Outport,'Fdq',100,0,1,1,'connect');

% Add Idq_m tag respect to Tm tag
PlaceNewBlock(circuit_model_path,'FromTm',path.From,'From_Idq_m',0,150,1,1);
set_param([circuit_model_path,'/From_Idq_m'],'GotoTag','Idq_m');
PlaceConnectNewBlock(circuit_model_path,'From_Idq_m',path.Outport,'Idq_m',100,0,1,1,'connect');

PlaceConnectNewBlock(circuit_model_path,'MuxV',path.Outport,'Vabc1',35,0,4,1,'connect');
 
if(n_set==1)
    circuital_model_outputs = ["Tm"; "Iabc"; "Fdq"; "Idq_m"; "Vabc";];
else
    circuital_model_outputs = ["Tm"; "Iabc1"; "Fdq"; "Idq_m"; "Vabc1";];
end

if n_set > 1
        % Add Ixy tag respect to Tm tag
        %PlaceNewBlock(circuit_model_path,'FromTm',path.From,'From_Ixy_m1',0,250,1,1);
        %set_param([circuit_model_path,'/From_Ixy_m1'],'GotoTag','Ixy_m1');
        %PlaceConnectNewBlock(circuit_model_path,'From_Ixy_m1',path.Outport,'Ixy_m1',100,0,1,1,'connect');    
    for i = 2:n_set
    %     if i == 2
    %         k = -50;
    %     else
    %         k = 0;
    %     end

        % Add Iabc tags respect to Tm tag and number n_phase
        PlaceNewBlock(circuit_model_path,'FromTm',path.From,['FromIabc' num2str(i)],0,100*i,1,1);   %k +50*(i-2)
        set_param([circuit_model_path,['/FromIabc' num2str(i)]],'GotoTag',['Iabc' num2str(i)]);
        PlaceConnectNewBlock(circuit_model_path,['FromIabc' num2str(i)],path.Outport,['Iabc' num2str(i)],100,0,1,1,'connect');
        circuital_model_outputs = [circuital_model_outputs; strcat('Iabc',num2str(i))];
    end
end

% Add sinCos tag
PlaceNewBlock(circuit_model_path,['From2' num2str(3*n_set)],path.From,'From_SinCos',0,50,1,1);
set_param([circuit_model_path,'/From_SinCos'],'GotoTag','SinCos');

%%  FLUX INTEGRATION subsystem

FIpath = [circuit_model_path '/Flux Integration'];
PlaceNewBlock(circuit_model_path,'MuxE',path.Subsystem,'Flux Integration',50,10,3*n_set + 1,1);
FIoutputs = [];
FIsuboutputs = [];
IFMoutputs = ["Idq_m"; "Tm"];
CCoutputs = ["Iabc1"; "Ia1"; "Ib1"; "Ic1"];
if n_set == 1
    FIoutputs = ["Fdq"; "F0"];
    FIsuboutputs = ["eab"; "e0"];
else
    FIoutputs = [FIoutputs; "Fdq"];
    FIsuboutputs = [FIsuboutputs; "eab"];
    for i = 1:(n_set-1)
        Fxy = strcat("Fxy",num2str(i));
        FIoutputs = [FIoutputs; Fxy];
        exy = strcat("exy",num2str(i));
        FIsuboutputs = [FIsuboutputs; exy];
        Ixy = strcat("Ixy_m",num2str(i));
        IFMoutputs = [IFMoutputs; Ixy];
        CCoutputs = [CCoutputs; strcat("Iabc",num2str(i+1)); strcat("Ia",num2str(i+1)); strcat("Ib",num2str(i+1)); strcat("Ic",num2str(i+1))];
    end
    FIoutputs = [FIoutputs; "F0"];
    FIsuboutputs = [FIsuboutputs; "e0"];
end

subadjust(circuit_model_path,FIpath,["e_abc" "SinCos"],FIoutputs,["no" "yes"],50,-25*n_set);
add_line(circuit_model_path,'MuxE/1','Flux Integration/1','autorouting','smart');
add_line(circuit_model_path,'From_SinCos/1','Flux Integration/2','autorouting','smart');

for i = 1:length(FIoutputs)
    pos(FIpath,'e_abc',FIoutputs(i),550,100*(i-1),1,1);
    if n_set > 1
        if (1 < i) && (i < length(FIoutputs))
            PlaceConnectNewBlock(FIpath,FIoutputs(i),path.Integrator,['Integrator' num2str(i)],-250,0,1,2,'connect');
        end
    end
end
pos(FIpath,'Fdq','SinCos',-130,15,1,1);
PlaceConnectNewBlock(FIpath,'Fdq',[SimulinkRefBlocks '/ab->dq'],'ab->dq',-50,0,1,3,'connect');
add_line(FIpath,'SinCos/1','ab->dq/2','autorouting','smart');
PlaceConnectNewBlock(FIpath,'ab->dq',path.Mux,'Mux',-80,0,1,3,'connect');
PlaceConnectNewBlock(FIpath,'Mux',path.Integrator,'IntegratorA',-50,-20,1,2,'connect');
PlaceConnectNewBlock(FIpath,'Mux',path.Integrator,'IntegratorB',-50,20,2,2,'connect');
PlaceNewBlock(FIpath,'Mux',path.Demux,'Demux',-170,0,3,1);
add_line(FIpath,'Demux/1','IntegratorA/1','autorouting','smart');
add_line(FIpath,'Demux/2','IntegratorB/1','autorouting','smart');
set_param([FIpath '/IntegratorA'],'InitialCondition','InitIntg_d');
set_param([FIpath '/IntegratorB'],'InitialCondition','InitIntg_q');

% abc->abxy0 subsystem of the FI subsystem
FIsubpath = [FIpath '/abc->abxy0'];
PlaceNewBlock(FIpath,'e_abc',path.Subsystem,'abc->abxy0',50,0,1,1);
subadjust(FIpath,FIsubpath,"e_abc",FIsuboutputs,["no" "no"]);
add_line(FIpath,'e_abc/1','abc->abxy0/1','autorouting','smart');
add_line(FIpath,'abc->abxy0/1','Demux/1','autorouting','smart');
add_line(FIpath,['abc->abxy0/' num2str(length(FIsuboutputs))],'F0/1','autorouting','smart');
for i = 2:(length(FIsuboutputs) - 1)
    add_line(FIpath,['abc->abxy0/' num2str(i)],['Integrator' num2str(i) '/1'],'autorouting','smart');
end
Mux_ports = [];
for i = 1:length(FIsuboutputs)
    pos(FIsubpath,'e_abc',FIsuboutputs(i),300,80*(i-1) - 40*(length(FIsuboutputs)-1),1,1);
    PlaceConnectNewBlock(FIsubpath,FIsuboutputs(i),path.Mux,['Mux' num2str(i)],-50,0,1,3,'connect');
    if i < length(FIsuboutputs)
        port1 = strcat("Mux",num2str(i),"/1");
        port2 = strcat("Mux",num2str(i),"/2");
        Mux_ports = [Mux_ports; port1; port2];
    elseif i == length(FIsuboutputs)
        set_param([FIsubpath '/Mux' num2str(i)],'Inputs',num2str(n_set));
        for k = 1:n_set
            portk = strcat("Mux",num2str(i),'/',num2str(k));
            Mux_ports = [Mux_ports; portk];
        end
    end
end
PlaceConnectNewBlock(FIsubpath,'e_abc',path.Gain,'Gain',50,0,1,1,'connect');
set_param([FIsubpath '/Gain'],'Gain','T','Multiplication','Matrix(K*u) (u vector)');
PlaceConnectNewBlock(FIsubpath,'Gain',path.Demux,'Demux',50,0,2,1,'connect');
set_param([FIsubpath '/Demux'],'Outputs',num2str(2*(length(FIsuboutputs) - 1) + n_set));
dim([FIsubpath '/Gain'],-20,-20,20,20);
dim([FIsubpath '/Demux'],0,-80,0,80);
for i = 1:(2*(length(FIsuboutputs) - 1) + n_set)
    add_line(FIsubpath,['Demux/' num2str(i)],Mux_ports(i),'autorouting','smart');
end

%% Rotation transformation subsystem
subpath = [circuit_model_path '/Subsystem'];
PlaceNewBlock(circuit_model_path,'From_SinCos',path.Inport,'theta_r',0,100,1,1);
PlaceConnectNewBlock(circuit_model_path,'theta_r',[SimulinkRefBlocks '/Subsystem'],'Subsystem',50,0,1,1,'connect');
PlaceConnectNewBlock(circuit_model_path,'Subsystem',path.Goto,'GotoSinCos',50,0,2,1,'connect');
set_param([circuit_model_path '/GotoSinCos'],'GotoTag','SinCos');



%% -----------------------------------------------------------------------%
%                   INVERSE FLUX MAPS SUBSYSTEM
% ------------------------------------------------------------------------%


IFMpath = [circuit_model_path '/Inverse Flux Maps'];
PlaceNewBlock(circuit_model_path,'Subsystem',[SimulinkRefBlocks '/Inverse Flux Maps'],'Inverse Flux Maps',0,80,1,1);
add_line(circuit_model_path,'theta_r/1','Inverse Flux Maps/1','autorouting','smart');
for i = 2:(length(FIoutputs) - 1)
    PlaceNewBlock(IFMpath,'Fdq',path.Inport,FIoutputs(i),0,250 + 50*(i-1),1,1);
    PlaceConnectNewBlock(IFMpath,FIoutputs(i),path.Gain,['Gain' num2str(i-1)],200,0,1,1,'connect');
    PlaceConnectNewBlock(IFMpath,['Gain' num2str(i-1)],path.Outport,IFMoutputs(i+1),200,0,2,1,'connect');
    set_param([IFMpath '/Gain' num2str(i-1)],'Gain','1/L_sigma');
    dim([IFMpath '/Gain' num2str(i-1)],-20,-10,20,10);
end
for i = 1:(length(FIoutputs) - 1)
    PlaceConnectNewBlock(circuit_model_path,'Inverse Flux Maps',path.From,strcat('From',FIoutputs(i)),-50,35*(i-1),i+1,1,'connect');
    set_param(strcat(circuit_model_path,'/From',FIoutputs(i)),'GotoTag',FIoutputs(i));
end
for i = 1:length(IFMoutputs)
    PlaceConnectNewBlock(circuit_model_path,'Inverse Flux Maps',path.Goto,strcat("Goto",IFMoutputs(i)),50,50*(i-1),length(FIoutputs)+i,1,'connect');
    set_param(strcat(circuit_model_path,'/Goto',IFMoutputs(i)),'GotoTag',IFMoutputs(i));
end



%% -----------------------------------------------------------------------%
%                   IRON LOSSES MAP SUBSYSTEM
% ------------------------------------------------------------------------%


ILMpath = [circuit_model_path '/Iron Losses Maps'];
PlaceNewBlock(circuit_model_path,'Inverse Flux Maps',[SimulinkRefBlocks '/Iron Losses Maps'],'Iron Losses Maps',0,100*n_set,1,1);
PlaceConnectNewBlock(circuit_model_path,'Iron Losses Maps',path.From,'FromIdq_m',-50,0,1,1,'connect');
set_param([circuit_model_path '/FromIdq_m'],'GotoTag','Idq_m');
PlaceConnectNewBlock(circuit_model_path,'Iron Losses Maps',path.Inport,'n_m',-50,20,2,1,'connect');
PlaceConnectNewBlock(circuit_model_path,'Iron Losses Maps',path.From,'From_Fdq_m',-50,35,3,1,'connect');
set_param([circuit_model_path '/From_Fdq_m'],'GotoTag','Fdq');
PlaceConnectNewBlock(circuit_model_path,'Iron Losses Maps',path.From,'From_Tm',-50,60,4,1,'connect');
set_param([circuit_model_path '/From_Tm'],'GotoTag','Tm');
PlaceConnectNewBlock(circuit_model_path,'Iron Losses Maps',path.Goto,'Goto_Idq',50,0,5,1,'connect');
set_param([circuit_model_path '/Goto_Idq'],'GotoTag','Idq');


%% -----------------------------------------------------------------------%
%                   CURRENTS COMPUTATION SUBSYTEM 
% ------------------------------------------------------------------------%


CCpath = [circuit_model_path '/Currents Calculation'];
PlaceNewBlock(circuit_model_path,'Iron Losses Maps',path.Subsystem,'Currents Calculation',0,150,4,1);
dim(CCpath,0,0,0,150*n_set);
CCinputs = ["Idq"; "SinCos"; IFMoutputs(3:end); "F0"];
subadjust(circuit_model_path,CCpath,CCinputs,CCoutputs,["yes" "yes"],0);
pos(CCpath,'Idq','SinCos',0,35,1,1);
if n_set > 1
    for i = 3:(length(CCinputs) - 1)
        pos(CCpath,'SinCos',CCinputs(i),0,80+40*(i-3),1,1);
    end
    pos(CCpath,CCinputs(i),'F0',0,80,1,1);
else
    pos(CCpath,'SinCos','F0',0,120,1,1);
end

pos(CCpath,'Idq','Iabc1',700,-100*(n_set-1),1,1); %%%
PlaceConnectNewBlock(CCpath,'Iabc1',path.Mux,'Mux1',-50,0,1,3,'connect');
set_param([CCpath '/Mux1'],'Inputs','3');
for i = 1:(length(CCoutputs) - 1)
    delta = 50;
    if (i/4 == ceil(i/4)) 
        delta = 100;
    end
    pos(CCpath,CCoutputs(i),CCoutputs(i+1),0,delta,1,1);
    if (i/4 == ceil(i/4)) 
        PlaceConnectNewBlock(CCpath,CCoutputs(i+1),path.Mux,['Mux' num2str(i+1)],-50,0,1,3,'connect');
        set_param([CCpath '/Mux' num2str(i+1)],'Inputs','3');
    end
end
PlaceConnectNewBlock(CCpath,'Idq',[SimulinkRefBlocks '/dq->ab1'],'dq->ab1',50,0,1,1,'connect');
add_line(CCpath,'SinCos/1','dq->ab1/2','autorouting','smart');
PlaceConnectNewBlock(CCpath,'dq->ab1',path.Mux,'Mux',100,100,3,1,'connect');
set_param([CCpath '/Mux'],'Inputs',num2str(n_set+1));
dim([CCpath '/Mux'],0,-15*(n_set-1),0,15*(n_set-1));
PlaceConnectNewBlock(CCpath,'F0',path.Gain,'Gain1',50,0,1,1,'connect');
set_param([CCpath '/Gain1'],'Gain','0.0');
PlaceConnectNewBlock(CCpath,'Gain1',path.Gain,'Gain2',50,0,2,1,'connect');
set_param([CCpath '/Gain2'],'Gain','1/3');
add_line(CCpath,'Gain2/1',['Mux/' num2str(n_set+1)],'autorouting','smart')
if n_set > 1
    for i = 3:(n_set+1)
        add_line(CCpath,strcat(CCinputs(i),'/1'),['Mux/' num2str(i-1)],'autorouting','smart')
    end
end
PlaceConnectNewBlock(CCpath,'Mux',path.Gain,'GainT',50,0,(n_set+1)+1,1,'connect');
set_param([CCpath '/GainT'],'Gain','T_inv','Multiplication','Matrix(K*u) (u vector)');
dim([CCpath '/GainT'],-20,-10,20,10);
PlaceConnectNewBlock(CCpath,'GainT',path.Memory,'Memory',50,0,2,1,'connect');
PlaceConnectNewBlock(CCpath,'Memory',path.Demux,'Demux',50,0,2,1,'connect');
set_param([CCpath '/Demux'],'Outputs',num2str(3*n_set));
dim([CCpath '/Demux'],0,-30*(n_set-1),0,30*(n_set-1));
for i = 1:n_set
    add_line(CCpath,['Demux/' num2str(3*i-2)],['Ia' num2str(i) '/1'],'autorouting','smart');
    add_line(CCpath,['Demux/' num2str(3*i-1)],['Ib' num2str(i) '/1'],'autorouting','smart');
    add_line(CCpath,['Demux/' num2str(3*i)],['Ic' num2str(i) '/1'],'autorouting','smart');

    add_line(CCpath,['Demux/' num2str(3*i-2)],['Mux' num2str(3*i+(i-3)) '/1'],'autorouting','smart');
    add_line(CCpath,['Demux/' num2str(3*i-1)],['Mux' num2str(3*i+(i-3)) '/2'],'autorouting','smart');
    add_line(CCpath,['Demux/' num2str(3*i)],['Mux' num2str(3*i+(i-3)) '/3'],'autorouting','smart');
end

%% MOTOR MODEL connections

add_line(motor_model_path,'Circuital Model/1','Mechanical equations/2','autorouting','smart');
add_line(motor_model_path,'Mechanical equations/1','Circuital Model/2','autorouting','smart');
add_line(motor_model_path,'Mechanical equations/3','Circuital Model/1','autorouting','smart');

for i = 1:(3*n_set)
    if(n_set==1)
    PlaceNewBlock(motor_model_path,'Circuital Model',path.PMCport,phase(i,:),-200,0,1+length(circuital_model_outputs)+i,1);
    else
        PlaceNewBlock(motor_model_path,'Circuital Model',path.PMCport,phase(i,:),-200,0,2+length(circuital_model_outputs)+i,1);
    end

    ports_PMC = get_param(strcat(motor_model_path,'/',phase(i,:)),'PortHandles');
    ports_CircuitalModel = get_param(circuit_model_path,'PortHandles');
    add_line(motor_model_path,ports_PMC.RConn,ports_CircuitalModel.LConn(i),'autorouting','smart');
%     add_line(motor_model_path,strcat(phase(i,:),'/1'),[circuit_model_path '/' num2str(i+2)],'autorouting','smart');

end
for i = 1:length(circuital_model_outputs)
    PlaceNewBlock(motor_model_path,'Circuital Model',path.Outport,circuital_model_outputs(i),100,0,2+i,1);
    add_line(motor_model_path,['Circuital Model/' num2str(i)],strcat(circuital_model_outputs(i),"/1"),'autorouting','smart');
end

%% Connections between motor and inverter

port_motor_model = get_param([filepath '/Motor Model'],'PortConnectivity');

for i=1:length(port_motor_model)
    if(strcmp(port_motor_model(i).Type,'LConn1'))
        index = i;  % indice da cui partono le connessioni a sinitra del motor model
    end
end
    
for i=1:n_set
    add_block('nesl_utility/Connection Label',[filepath '/aa' num2str(i)]);
    set_param([filepath '/aa' num2str(i)],'Orientation','left');
    set_param([filepath '/aa' num2str(i)],'label',['a' num2str(i)]);
    set_param([filepath '/aa' num2str(i)],'ShowName','off');
    tmp = get_param([filepath '/aa' num2str(i)],'position');
    width = (tmp(3)-tmp(1));
    height = tmp(4)-tmp(2);
    tmp = [port_motor_model(index).Position(1)-80 port_motor_model(index).Position(2)-height*0.5 port_motor_model(index).Position(1)+width-80 port_motor_model(index).Position(2)+height*0.5];
    set_param([filepath '/aa' num2str(i)],'Position',tmp);
    set_param([filepath '/aa' num2str(i)],'ShowName','off');
    add_line(filepath,[tmp(3) tmp(2)+height*0.5; port_motor_model(index).Position]);
    
    
    add_block('nesl_utility/Connection Label',[filepath '/bb' num2str(i)]);
    set_param([filepath '/bb' num2str(i)],'Orientation','left');
    set_param([filepath '/bb' num2str(i)],'label',['b' num2str(i)]);
    set_param([filepath '/bb' num2str(i)],'ShowName','off');
    tmp = get_param([filepath '/bb' num2str(i)],'position');
    width = (tmp(3)-tmp(1));
    height = tmp(4)-tmp(2);
    tmp = [port_motor_model(index+1).Position(1)-80 port_motor_model(index+1).Position(2)-height*0.5 port_motor_model(index+1).Position(1)+width-80 port_motor_model(index+1).Position(2)+height*0.5];
    set_param([filepath '/bb' num2str(i)],'Position',tmp);
    set_param([filepath '/bb' num2str(i)],'ShowName','off');
    add_line(filepath,[tmp(3) tmp(2)+height*0.5; port_motor_model(index+1).Position]);
    

    add_block('nesl_utility/Connection Label',[filepath '/cc' num2str(i)]);
    set_param([filepath '/cc' num2str(i)],'Orientation','left');
    set_param([filepath '/cc' num2str(i)],'label',['c' num2str(i)]);
    set_param([filepath '/cc' num2str(i)],'ShowName','off');
    tmp = get_param([filepath '/cc' num2str(i)],'position');
    width = (tmp(3)-tmp(1));
    height = tmp(4)-tmp(2);
    tmp = [port_motor_model(index+2).Position(1)-80 port_motor_model(index+2).Position(2)-height*0.5 port_motor_model(index+2).Position(1)+width-80 port_motor_model(index+2).Position(2)+height*0.5];
    set_param([filepath '/cc' num2str(i)],'Position',tmp);
    set_param([filepath '/cc' num2str(i)],'ShowName','off');
    add_line(filepath,[tmp(3) tmp(2)+height*0.5; port_motor_model(index+2).Position]);

    index = index+3;
    
end


%% -----------------------------------------------------------------------%
%                               ADD BUS CREATOR 
%%------------------------------------------------------------------------%

PlaceNewBlock(filepath,'Motor Model',path.BusCreator,'BusCreator',300,0,1,1);

if(n_set==1)
    set_param([filepath '/BusCreator'],'inputs',num2str(7));
else
    set_param([filepath '/BusCreator'],'inputs',num2str(7+n_set-1));
end    

dim([filepath '/BusCreator'],0,0,0,50*11);

Bus_inputs = get_param([filepath '/BusCreator'],'PortHandles');
Motor_outputs = get_param([filepath '/Motor Model'],'PortHandles');


% Add connection and signal name
tmp = add_line(filepath,Motor_outputs.Outport(1),Bus_inputs.Inport(1),'autorouting','smart');
set_param(tmp, 'Name', 'n_m');
tmp = add_line(filepath,Motor_outputs.Outport(2),Bus_inputs.Inport(2),'autorouting','smart');
set_param(tmp, 'Name', 'theta_r');
tmp = add_line(filepath,Motor_outputs.Outport(3),Bus_inputs.Inport(3),'autorouting','smart');
set_param(tmp, 'Name', 'T_m');
tmp = add_line(filepath,Motor_outputs.Outport(4),Bus_inputs.Inport(4),'autorouting','smart');
set_param(tmp, 'Name', 'Iabc1');
tmp = add_line(filepath,Motor_outputs.Outport(5),Bus_inputs.Inport(5),'autorouting','smart');
set_param(tmp, 'Name', 'lambda_dq');
tmp = add_line(filepath,Motor_outputs.Outport(6),Bus_inputs.Inport(6),'autorouting','smart');
set_param(tmp, 'Name', 'Idq_m');
tmp = add_line(filepath,Motor_outputs.Outport(7),Bus_inputs.Inport(7),'autorouting','smart');
set_param(tmp, 'Name', 'Vabc1');


for i=2:n_set
    tmp=add_line(filepath,Motor_outputs.Outport(7+i-1),Bus_inputs.Inport(7+i-1),'autorouting','smart');
    set_param(tmp, 'Name', ['Iabc' num2str(i)]);

end    

PlaceNewBlock(filepath,'BusCreator',path.ToWorkspace,'MotorOutputs',100,n_set*50,1,1);
port_toworkspace = get_param([filepath '/MotorOutputs'],'PortHandles');
PlaceNewBlock(filepath,'BusCreator',path.Goto,'Out_M',100,n_set*70,1,1);
port_out_ctrl = get_param([filepath '/Out_M'],'PortHandles');
set_param([filepath '/Out_M'],'GotoTag','Out_M');
set_param([filepath '/Out_M'],'TagVisibility','global')


add_line(filepath,Bus_inputs.Outport,port_toworkspace.Inport,'autorouting','smart');
add_line(filepath,Bus_inputs.Outport,port_out_ctrl.Inport,'autorouting','smart');

%% Add current goto references

% corrente iab1
%PlaceConnectNewBlock(filepath,'Motor Model ',path.Goto,'iabc_1',300,-100,1,1,'connect');
PlaceNewBlock(filepath,'BusCreator',path.Goto,'iabc_1',0,-60,1,1);
set_param([filepath '/iabc_1' ],'GotoTag','iabc1');
port_goto_iabc_1 = get_param([filepath '/iabc_1'],'PortHandles');
PlaceNewBlock(filepath,'BusCreator',path.Goto,'theta_r',0,-120,1,1);
set_param([filepath '/theta_r' ],'GotoTag','theta_mec_meas');
port_theta_r = get_param([filepath '/theta_r'],'PortHandles');

add_line(filepath,Motor_outputs.Outport(4),port_goto_iabc_1.Inport,'autorouting','smart');
add_line(filepath,Motor_outputs.Outport(2),port_theta_r.Inport,'autorouting','smart');

% altre correnti
% 
for i=2:n_set
    %addpos(filepath,'BusCreator',path.Goto,['iabc_' num2str(i)],0,50*(i-1),11,1);
    y = 80+50*i;
    z = 8+i;
    PlaceConnectNewBlock(filepath,'Motor Model',path.Goto,['iabc_' num2str(i)],150,y,z,1,'connect');
    set_param([filepath '/iabc_' num2str(i)],'GotoTag',['iabc' num2str(i)]);
end

%% Add load torque and speed

add_block([SimulinkRefBlocks '/load torque and speed'],[filepath '/load torque speed']);
dim([filepath '/load torque speed'],900,0,900,0);
load_ports = get_param([filepath '/load torque speed'],'PortHandles');
port_motor_model = get_param([filepath '/Motor Model'],'PortHandles');
add_line(filepath,load_ports.Outport(1),port_motor_model.Inport(1),'autorouting','smart');
add_line(filepath,load_ports.Outport(2),port_motor_model.Inport(2),'autorouting','smart');


%% -----------------------------------------------------------------------%
%                           ADD SCOPES SUBSYTEM
% ------------------------------------------------------------------------%

add_block([SimulinkRefBlocks '/Scopes'],[filepath '/Scopes']);
dim([filepath '/Scopes'],0,500,0,500);


end