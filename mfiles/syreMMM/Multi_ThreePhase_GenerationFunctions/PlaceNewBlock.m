function PlaceNewBlock(filepath,refblock,path_newblock,newblock,x,y,port1,port2,flip) 

% refblock/newblock must be the string referred to the path of the block in the simulink file
% path_newblock must be the string which allows to take the block from the libraries of simulink
% x and y are the displacments of the new block with respect to the reference one
% port1 and port2 are the ports of the reference and the new block respectively, that set the frame of reference of x and y (on port1) in order to establish were the port2 will be placed    
% flip must be a string that, if necessary, allows to rotate the new block up, down, left or right

% Define paths 
ref_block_path = strcat(filepath,'/',refblock);
new_block_path = strcat(filepath,'/',newblock);


% Get ports of reference block
    refport = get_param(ref_block_path,'PortConnectivity');
    add_block(path_newblock,new_block_path);

% Flip block if necessary
    if nargin == 9
        set_param(new_block_path,'Orientation',flip);
    end

% Get dimensions of new block     
    newdim = get_param(new_block_path,'position');
    wn = newdim(3)-newdim(1); % width of the new block
    hn = newdim(4)-newdim(2); % height of the new block
    centrenew = [(newdim(1)+newdim(3))/2; (newdim(2)+newdim(4))/2]; % coordinates of the new block centre
    
% Get ports of new block    
    newport = get_param(new_block_path,'PortConnectivity');    
    deltax = wn/2 - (centrenew(1) - newport(port2).Position(1));
    deltay = hn/2 - (centrenew(2) - newport(port2).Position(2));
    
% Set coordinates of new block, computing the displacement x and y taking the position of port1 of reference block     
    coordinates = [refport(port1).Position(1) + x - deltax; refport(port1).Position(2) + y - deltay; refport(port1).Position(1) + x + (wn - deltax); refport(port1).Position(2) + y + (hn - deltay)];
    set_param(new_block_path,'position',coordinates);
  
end