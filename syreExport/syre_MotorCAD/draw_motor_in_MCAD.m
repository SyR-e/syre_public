function draw_motor_in_MCAD(filename, pathname)

load([pathname filename])
mcad=actxserver('MotorCAD.AppAutomation');
 
%% Materials
materials = geo.BLKLABELS.materials;

%%% note: the material names in syre and motorcad must be the same

%Iron
tmp = materials(4);
%tmp = strcat(tmp,'_ImportedFromSyR-e');
tmp = convertCharsToStrings(tmp);
invoke(mcad,'SetComponentMaterial','Stator Lam (Back Iron)',tmp);
invoke(mcad,'SetComponentMaterial','Stator Lam (Tooth)',tmp);      
invoke(mcad,'SetComponentMaterial','Rotor Lam (Back Iron)',tmp);    
invoke(mcad,'SetComponentMaterial','Shaft [Active]',tmp); 

%Magnet
tmp = materials(6);
tmp = convertCharsToStrings(tmp);
invoke(mcad,'SetComponentMaterial','Magnet',tmp);

%% Stator and Rotor
invoke(mcad,'SetVariable','StatorRotation',0);
invoke(mcad,'SetVariable','RotorRotation',(90/geo.p));

draw_stator_in_MCAD(mcad,geo,0,0,dataSet)
draw_rotor_in_MCAD(mcad,geo,0,0)

% %Lamination stacking factor
% invoke(mcad,'SetVariable','Stacking_Factor_[Stator]',1);
% invoke(mcad,'SetVariable','Stacking_Factor_[Rotor]',1);

file_mot=strrep(filename,'.mat','.mot');
invoke(mcad,'SaveToFile',[pathname file_mot]);

%% Winding
windingSyreToMCAD(mcad,pathname,filename,file_mot)  

%% dxf
syreToDxfMCAD(pathname,filename)              

%.dxf MCAD settings
invoke(mcad,'LoadDXFFile',[pathname filename(1:end-4),'.dxf']); 
invoke(mcad,'SetVariable','UseDXFImportForFEA_Magnetic', true);
invoke(mcad,'SetVariable','UseDXFImportForFEA_Mechanical',true);
invoke(mcad,'SetVariable','DXFImportType',1);

%% Output 
%Save MCAD model
invoke(mcad,'SaveToFile',[pathname file_mot]);

%Save workspace
save([pathname,filename],'dataSet','geo','per','mat');

disp('Motor-CAD file saved in:')
disp([pathname file_mot])
disp(' ')
disp('Syr-e file saved in:')
disp([pathname filename])
disp(' ')
end