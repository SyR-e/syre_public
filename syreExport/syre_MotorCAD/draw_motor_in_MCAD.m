function draw_motor_in_MCAD(filename, pathname)

load([pathname filename])
mcad=actxserver('MotorCAD.AppAutomation');
 
% %%%%%%%%%% Materials
% % tmp=geo.BLKLABELS.materials(4);
% % tmp = convertCharsToStrings(tmp);
% % invoke(mcad,'SetComponentMaterial','Stator Lam (Back Iron)',tmp);
% % invoke(mcad,'SetComponentMaterial','Stator Lam (Tooth)',tmp);       %stator
% % invoke(mcad,'SetComponentMaterial','Rotor Lam (Back Iron)',tmp);    %rotor
% % invoke(mcad,'SetComponentMaterial','Shaft [Active]',tmp);           %shaft

%Stator and Rotor angle
invoke(mcad,'SetVariable','StatorRotation',0);
invoke(mcad,'SetVariable','RotorRotation',(90/geo.p));
%invoke(mcad,'SetVariable','BPM_Rotor','13');

% stator parameters 
draw_stator_in_MCAD(mcad,geo,per,mat,dataSet)
% rotor parameters
draw_rotor_in_MCAD(mcad,geo,per,mat)

%%%%%%Lamination stacking factor
invoke(mcad,'SetVariable','Stacking_Factor_[Stator]',1);
invoke(mcad,'SetVariable','Stacking_Factor_[Rotor]',1);

file_mot=strrep(filename,'.mat','.mot');
invoke(mcad,'SaveToFile',[pathname file_mot]);

windingSyreToMCAD(mcad,pathname,filename,file_mot)          %export winding to MotorCAD
syreToDxfMCAD(pathname,filename)             %create a proper .dxf for MotorCAD with 1 slot and 1 pole

%.dxf MCAD settings
invoke(mcad,'LoadDXFFile',[pathname filename(1:end-4),'.dxf']);        %load the .dxf to MotorCAD
invoke(mcad,'SetVariable','UseDXFImportForFEA_Magnetic', true);
invoke(mcad,'SetVariable','UseDXFImportForFEA_Mechanical',true);
invoke(mcad,'SetVariable','DXFImportType',1);

%Save MCAD model
invoke(mcad,'SaveToFile',[pathname file_mot]);
% invoke(mcad,'Quit');

%Save workspace
save([pathname,filename],'dataSet','geo','per','mat');

disp('Motor-CAD file saved in:')
disp([pathname file_mot])
disp(' ')
disp('Syr-e file saved in:')
disp([pathname filename])
disp(' ')

end