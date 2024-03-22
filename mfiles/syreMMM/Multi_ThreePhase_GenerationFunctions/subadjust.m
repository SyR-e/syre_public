function subadjust(filepath,subpath,inputs,outputs,activation,delta,y)
% activation must be a vector containing combinations of yes or no.
% Ex: ["no" "yes"] in order to only connect the output ports of the subsystem to Goto blocks

Simulink.SubSystem.deleteContents(subpath);
inpquantity = length(inputs);
outquantity = length(outputs);
path = getComponentsPaths();

subsystemstring = strrep(subpath,[filepath '/'],'');
if nargin < 6
    delta = 50;
    y = 0;
elseif nargin < 7
    y = 0;
end

for i = 1:inpquantity
    add_block(path.Inport,strcat(subpath,'/',inputs(i)));
end

for i = 1:outquantity
    add_block(path.Outport,strcat(subpath,'/',outputs(i)));
%     addposlink(filepath,subsystemstring,path.Goto,strcat('Goto',outputs(i)),50,delta*(i-1) + y,inpquantity + i,1,'connect')
%     set_param(strcat(filepath,'/Goto',outputs(i)),'GotoTag',outputs(i));
end

% for i = 1:outquantity
%     addposlink(filepath,subsystemstring,path.Goto,strcat('Goto',outputs(i)),50,delta*(i-1) + y,inpquantity + i,1,'connect')
%     set_param(strcat(filepath,'/Goto',outputs(i)),'GotoTag',outputs(i));
% end

if nargin > 4
    logic_vector = (activation == ["yes" "yes"]);

    if logic_vector(1) == 1
        for i = 1:inpquantity
            PlaceConnectNewBlock(filepath,subsystemstring,path.From,strcat('From',inputs(i)),-50,delta*(i-1) + y,i,1,'connect')
            set_param(strcat(filepath,'/From',inputs(i)),'GotoTag',inputs(i));
        end
    end

    if logic_vector(2) == 1
        for i = 1:outquantity
            PlaceConnectNewBlock(filepath,subsystemstring,path.Goto,strcat('Goto',outputs(i)),50,delta*(i-1) + y,inpquantity + i,1,'connect')
            set_param(strcat(filepath,'/Goto',outputs(i)),'GotoTag',outputs(i));
        end
    end
end

end

