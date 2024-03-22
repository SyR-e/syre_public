function PlaceConnectNewBlock(filepath,refblock,path_newblock,newblock,x,y,port1,port2,flip,link) 
% refblock/newblock must be the string referred to the path of the block in the simulink file
% path_newblock must be the string which allows to take the block from the libraries of simulink
% x and y are the displacments of the new block with respect to the reference one
% port1 and port2 are the ports of the reference and the new block respectively, that set the frame of reference of x and y (on port1) in order to establish were the port2 will be placed    
% flip must be a string that, if necessary, allows to rotate the new block up, down, left or right


    ref_block_path = strcat(filepath,'/',refblock);
    new_block_path = strcat(filepath,'/',newblock);

    % Get reference block type
    refport = get_param(ref_block_path,'PortConnectivity');
    link_string = 'connect';
    add_block(path_newblock,new_block_path);

    a = strcmp(flip,link_string);
    if nargin >= 9
        if a == 0
            set_param(new_block_path,'Orientation',flip);        
        end
    end

    % Get new block dimensions 
    newdim = get_param(new_block_path,'position');
    wn = newdim(3)-newdim(1);   % width of the new block
    hn = newdim(4)-newdim(2);   % height of the new block
    centrenew = [(newdim(1)+newdim(3))/2; (newdim(2)+newdim(4))/2]; % coordinates of the new block centre
        
% Get ports of new block    
    newport = get_param(new_block_path,'PortConnectivity');    
    deltax = wn/2 - (centrenew(1) - newport(port2).Position(1));
    deltay = hn/2 - (centrenew(2) - newport(port2).Position(2));

% Set coordinates of new block, computing the displacement x and y taking the position of port1 of reference block     
    coordinates = [refport(port1).Position(1) + x - deltax; refport(port1).Position(2) + y - deltay; refport(port1).Position(1) + x + (wn - deltax); refport(port1).Position(2) + y + (hn - deltay)];
    set_param(new_block_path,'position',coordinates);

%% Connecting simulink blocks

    refPortHandles = get_param(ref_block_path,'PortHandles');
    refnumberofportsinput = length(refPortHandles.Inport);
    refoutputport = num2str(port1-refnumberofportsinput);
    refinputport = num2str(port1);

    newPortHandles = get_param(new_block_path,'PortHandles');
    newnumberofportsinput = length(newPortHandles.Inport);     
    newinputport = num2str(port2);
    newoutputport = num2str(port2-newnumberofportsinput);

    if x >= 0
        add_line(filepath,strcat(refblock,'/',refoutputport),strcat(newblock,'/',newinputport),'autorouting','smart');
    else
        add_line(filepath,strcat(newblock,'/',newoutputport),strcat(refblock,'/',refinputport),'autorouting','smart');
    end
end
    