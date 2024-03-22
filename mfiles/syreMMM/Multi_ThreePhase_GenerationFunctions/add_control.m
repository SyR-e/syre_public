function add_control(SimulinkRefBlock,filepath,path,n_set)

Simulink.SubSystem.copyContentsToBlockDiagram([SimulinkRefBlock '/Digital Control'],filepath);

% Change subsytem size
pos = get_param([filepath '/Digital Control'],'position');
delta_y = 50;   
new_pos = [pos(1) pos(2) pos(3) pos(4)+delta_y*n_set];
set_param([filepath '/Digital Control'],'position',new_pos);


% Resize mux 
demux_outputs = str2num(get_param([filepath '/Digital Control/Demux2'],'Outputs'));
set_param([filepath '/Digital Control/Demux2'],'Outputs',num2str(demux_outputs+(n_set-1)*4));
demux_position = get_param([filepath '/Digital Control/Demux2'],'position');
demux_position(4) = demux_position(4)+(n_set-1)*4*40;
set_param([filepath '/Digital Control/Demux2'],'position',demux_position);

mux_inputs = str2num(get_param([filepath '/Digital Control/Mux'],'Inputs'));
set_param([filepath '/Digital Control/Mux'],'Inputs',num2str(mux_inputs+(n_set-1)));

mux_position = get_param([filepath '/Digital Control/Mux'],'position');
mux_position(4) = mux_position(4)+(n_set-1)*50;
set_param([filepath '/Digital Control/Mux'],'position',mux_position);

demux_ports = get_param([filepath '/Digital Control/Demux2'],'PortConnectivity');
mux_ports = get_param([filepath '/Digital Control/Mux'],'PortConnectivity');
demux_k = 36;
mux_k = 10;

if(n_set>1)
    for i=2:n_set
        add_block(path.Inport,[filepath '/Digital Control/isabc_' num2str(i-1)]);
        tmp = get_param([filepath '/Digital Control/isabc_' num2str(i-1)],'Position');
        width = tmp(3)-tmp(1);
        height = tmp(4)-tmp(2);
        tmp = [mux_ports(mux_k).Position(1)-100 mux_ports(mux_k).Position(2)-height*0.5 mux_ports(mux_k).Position(1)+width-100 mux_ports(mux_k).Position(2)+height*0.5];
        set_param([filepath '/Digital Control/isabc_' num2str(i-1)],'Position',tmp);
        add_line([filepath '/Digital Control'],[tmp(3) tmp(2)+height*0.5; mux_ports(mux_k).Position]);
        mux_k = mux_k+1; 
    end
end

if(n_set>1)
    for i=2:n_set
    add_block(path.Mux,[filepath '/Digital Control/mux' num2str(i)])
    set_param([filepath '/Digital Control/mux' num2str(i)],'Inputs','3');
    tmp = get_param([filepath '/Digital Control/mux' num2str(i)],'position');
    width = tmp(3)-tmp(1);
    height = tmp(4)-tmp(2);
    height = height+3*20; 
    tmp = [demux_ports(demux_k+1).Position(1)+100 demux_ports(demux_k+1).Position(2)-height*0.5 demux_ports(demux_k+1).Position(1)+width+100 demux_ports(demux_k+1).Position(2)+height*0.5];
    set_param([filepath '/Digital Control/mux' num2str(i)],'position',tmp);
    tmp_mux_port = get_param([filepath '/Digital Control/mux' num2str(i)],'PortConnectivity');
    add_line([filepath '/Digital Control'],[demux_ports(demux_k).Position; tmp_mux_port(1).Position]);
    add_line([filepath '/Digital Control'],[demux_ports(demux_k+1).Position; tmp_mux_port(2).Position]);
    add_line([filepath '/Digital Control'],[demux_ports(demux_k+2).Position; tmp_mux_port(3).Position]);
    
    add_block(path.Outport,[filepath '/Digital Control/duty_abc_' num2str(i)]);
    tmp = get_param([filepath '/Digital Control/duty_abc_' num2str(i)],'Position');
    width = tmp(3)-tmp(1);
    height = tmp(4)-tmp(2);
    tmp = [tmp_mux_port(4).Position(1)+100 tmp_mux_port(4).Position(2)-height*0.5 tmp_mux_port(4).Position(1)+width+100 tmp_mux_port(4).Position(2)+height*0.5];
    set_param([filepath '/Digital Control/duty_abc_' num2str(i)],'Position',tmp);
    add_line([filepath '/Digital Control'],[tmp_mux_port(4).Position; tmp(1) tmp(2)+height*0.5]);

    
    add_block(path.Outport,[filepath '/Digital Control/pwm_stop_' num2str(i)]);
    tmp = get_param([filepath '/Digital Control/pwm_stop_' num2str(i)],'Position');
    width = tmp(3)-tmp(1);
    height = tmp(4)-tmp(2);
    tmp = [demux_ports(demux_k+3).Position(1)+100 demux_ports(demux_k+3).Position(2)-height*0.5 demux_ports(demux_k+3).Position(1)+width+100 demux_ports(demux_k+3).Position(2)+height*0.5];
    set_param([filepath '/Digital Control/pwm_stop_' num2str(i)],'Position',tmp);
    add_line([filepath '/Digital Control'],[demux_ports(demux_k+3).Position; tmp(1) tmp(2)+height*0.5]);
    demux_k = demux_k+4; 
    
    end



 % Add out froms 

 ports_control = get_param([filepath '/Digital Control'],'PortConnectivity');

 k = 8;
 for i=1:n_set
    add_block('simulink/Signal Routing/From',[filepath '/iabc' num2str(i) ]);
    set_param([filepath '/iabc' num2str(i) ],'GotoTag',['iabc' num2str(i)]);
    set_param([filepath '/iabc' num2str(i) ],'GotoTag',['iabc' num2str(i)]);
    tmp = get_param([filepath '/iabc' num2str(i) ],'position');
    width = (tmp(3)-tmp(1))*1.5;
    height = tmp(4)-tmp(2);
    
    if i==1
        tmp = [ports_control(1).Position(1)-100 ports_control(1).Position(2)-height*0.5 ports_control(1).Position(1)+width-100 ports_control(1).Position(2)+height*0.5];
    else 
        tmp = [ports_control(k).Position(1)-100 ports_control(k).Position(2)-height*0.5 ports_control(k).Position(1)+width-100 ports_control(k).Position(2)+height*0.5];
    end    
    
    set_param([filepath '/iabc' num2str(i) ],'Position',tmp);
    if i==1
        add_line(filepath,[tmp(3) tmp(2)+height*0.5; ports_control(1).Position]);
    else
        add_line(filepath,[tmp(3) tmp(2)+height*0.5; ports_control(k).Position]);
    end
    k=k+1;
 end

end


