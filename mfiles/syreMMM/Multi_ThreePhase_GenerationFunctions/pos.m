function pos(filepath,refblock,newblock,x,y,port1,port2,flip) 
% refblock/newblock must be the string referred to the path of the block in the simulink file
% x and y are the displacments of the new block with respect to the reference one
% port1 and port2 are the ports of the reference and the new block respectively, that set the frame of reference of x and y (on port1) in order to establish were the port2 will be placed    
% flip must be a string that, if necessary, allows to rotate the new block up, down, left or right

    refport = get_param(strcat(filepath,'/',refblock),'PortConnectivity');
    refmatrix = [];

    if nargin == 8
        set_param(strcat(filepath,'/',newblock),'Orientation',flip);
    end

    newdim = get_param(strcat(filepath,'/',newblock),'position');
    wn = newdim(3)-newdim(1); % width of the new block
    hn = newdim(4)-newdim(2); % height of the new block
    centrenew = [(newdim(1)+newdim(3))/2; (newdim(2)+newdim(4))/2]; % coordinates of the new block centre
    newport = get_param(strcat(filepath,'/',newblock),'PortConnectivity');
    newmatrix = [];

    for i = 1:numel(refport)       
        % position of port connector
        xpref = refport(i).Position(1);
        ypref = refport(i).Position(2);            
        % port block datas      
        refmatrix = [refmatrix; xpref ypref];
    end

    for i = 1:numel(newport)       
        % position of port connector
        xpnew = newport(i).Position(1);
        ypnew = newport(i).Position(2); 
        % port block datas
        newmatrix = [newmatrix; xpnew ypnew]; 
    end
    
    deltax = wn/2 - (centrenew(1) - newmatrix(port2,1));
    deltay = hn/2 - (centrenew(2) - newmatrix(port2,2));    
    coordinates = [refmatrix(port1,1) + x - deltax; refmatrix(port1,2) + y - deltay; refmatrix(port1,1) + x + (wn - deltax); refmatrix(port1,2) + y + (hn - deltay)];
    set_param(strcat(filepath,'/',newblock),'position',coordinates);
  
end