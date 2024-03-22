
function add_inverters(filepath,SimulinkRefBlocks,n_set)

add_block([SimulinkRefBlocks '/Inverter Model'],[filepath '/Inverter Model 1']);
inv1_pos = get_param([filepath '/Inverter Model 1'],'position');
set_param([filepath '/Inverter Model 1'],'position',[inv1_pos(1)+350 inv1_pos(2) inv1_pos(3)+350 inv1_pos(4) ]);
if(n_set>1)
    for i=2:n_set
        refblock = [filepath '/Inverter Model ' num2str(i-1)];
        path_newblock = [filepath '/Inverter Model ' num2str(i)];        
        delta_y = 200;
        refblock_position = get_param(refblock,'position');
        newblock_coordinates = [refblock_position(1);refblock_position(2)+delta_y; refblock_position(3); refblock_position(4)+delta_y];
        add_block([SimulinkRefBlocks '/Inverter Model'],path_newblock);
        set_param(path_newblock,'position',newblock_coordinates);
    end
end
ports_DigitalControl = get_param([filepath '/Digital Control'],'PortHandles');

k =1;
for i=1:n_set
    tmp = get_param([filepath '/Inverter Model ' num2str(i)],'PortHandles');
    eval([ 'ports_Inverter' num2str(i) '= tmp']);
    add_line(filepath,ports_DigitalControl.Outport(k),eval([ 'ports_Inverter' num2str(i) '.Inport(1)']),'autorouting','smart'); %duty_abc
    add_line(filepath,ports_DigitalControl.Outport(k+1),eval([ 'ports_Inverter' num2str(i) '.Inport(2)']),'autorouting','smart'); %pmw_stop
    if(i==1)
        k=k+3;
    else
        k=k+2;
    end
end

for i=1:n_set
    ports_inv = get_param([filepath '/Inverter Model ' num2str(i)],'PortConnectivity');
    if(i==1) % Add vdc measurement only for inverter 1        
        add_block('simulink/Signal Routing/Goto',[filepath '/Vdc' num2str(i) ]);
        set_param([filepath '/Vdc' num2str(i) ],'GotoTag',['vdc' num2str(i)]);
        set_param([filepath '/Vdc' num2str(i) ],'ShowName','off');
        tmp = get_param([filepath '/Vdc' num2str(i) ],'position');
        width = (tmp(3)-tmp(1))*1.5;
        height = tmp(4)-tmp(2);
        tmp = [ports_inv(4).Position(1)+80 ports_inv(4).Position(2)-height*0.5 ports_inv(4).Position(1)+width+80 ports_inv(4).Position(2)+height*0.5];
        set_param([filepath '/Vdc' num2str(i) ],'Position',tmp);
        add_line(filepath,[ports_inv(4).Position;tmp(1) tmp(2)+height*0.5]);
    else
        delete_block([filepath '/Inverter Model ' num2str(i) '/Vdc']);
    end

    % Add From iabc
    add_block('simulink/Signal Routing/From',[filepath '/From_iabc' num2str(i) ]);
    set_param([filepath '/From_iabc' num2str(i) ],'GotoTag',['iabc' num2str(i)]);
    tmp = get_param([filepath '/From_iabc' num2str(i) ],'position');
    width = (tmp(3)-tmp(1))*1.5;
    height = tmp(4)-tmp(2);
    tmp = [ports_inv(3).Position(1)-80 ports_inv(3).Position(2)-height*0.5 ports_inv(3).Position(1)+width-80 ports_inv(3).Position(2)+height*0.5];
    set_param([filepath '/From_iabc' num2str(i) ],'Position',tmp);
    set_param([filepath '/From_iabc' num2str(i) ],'ShowName','off');
    add_line(filepath,[tmp(3) tmp(2)+height*0.5; ports_inv(3).Position]);
end

for i=1:n_set
    % Add electrical connection labels
    ports_inv = get_param([filepath '/Inverter Model ' num2str(i)],'PortConnectivity');
    
    if(i==1)
        x = 5;
        y = 6;
        z = 7;
    else
        x = 4;
        y = 5;
        z = 6;
    end

    add_block('nesl_utility/Connection Label',[filepath '/a' num2str(i)]);
    tmp = get_param([filepath '/a' num2str(i)],'position');
    set_param([filepath '/a' num2str(i)],'label',['a' num2str(i)]);
    set_param([filepath '/a' num2str(i)],'ShowName','off');
    width =  tmp(3)-tmp(1);
    height = tmp(4)-tmp(2);
    tmp = [ports_inv(x).Position(1)+80 ports_inv(x).Position(2)-height*0.5 ports_inv(x).Position(1)+width+80 ports_inv(x).Position(2)+height*0.5];
    set_param([filepath '/a' num2str(i) ],'Position',tmp);
    add_line(filepath,[ports_inv(x).Position;tmp(1) tmp(2)+height*0.5]);

    add_block('nesl_utility/Connection Label',[filepath '/b' num2str(i)]);
    tmp = get_param([filepath '/b' num2str(i)],'position');
    set_param([filepath '/b' num2str(i)],'label',['b' num2str(i)]);
    set_param([filepath '/b' num2str(i)],'ShowName','off');
    width = (tmp(3)-tmp(1));
    height = tmp(4)-tmp(2);
    tmp = [ports_inv(y).Position(1)+80 ports_inv(y).Position(2)-height*0.5 ports_inv(y).Position(1)+width+80 ports_inv(y).Position(2)+height*0.5];
    set_param([filepath '/b' num2str(i) ],'Position',tmp);
    add_line(filepath,[ports_inv(y).Position;tmp(1) tmp(2)+height*0.5]);
    
    add_block('nesl_utility/Connection Label',[filepath '/c' num2str(i)]);
    tmp = get_param([filepath '/c' num2str(i)],'position');
    set_param([filepath '/c' num2str(i)],'label',['c' num2str(i)]);
    set_param([filepath '/c' num2str(i)],'ShowName','off');
    width = (tmp(3)-tmp(1));
    height = tmp(4)-tmp(2);
    tmp = [ports_inv(z).Position(1)+80 ports_inv(z).Position(2)-height*0.5 ports_inv(z).Position(1)+width+80 ports_inv(z).Position(2)+height*0.5];
    set_param([filepath '/c' num2str(i) ],'Position',tmp);
    add_line(filepath,[ports_inv(z).Position;tmp(1) tmp(2)+height*0.5]);
    
end


