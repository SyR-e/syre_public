% Copyright 2024
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [dataSet,model] = DrawAndSaveMachineCOMSOL(dataSet,filename,pathname) %FILE NAME E PATH NAME da specificare quando inserisco procedura nel piano xb
% function dataSet = DrawPushMachine(handles,filename,pathname)

% if ~isfield(app,'dataSet')  % not launched from GUI_Syre
%     flagGUI=0;
%     dataSet=app;
% else                            % launched from GUI_Syre
%     flagGUI=1;
%     dataSet=app.dataSet;
% end
clc
disp('Comsol Motor Model creation...')

nameIn  = dataSet.currentfilename;
nameIn  = strrep(nameIn,'.mat','.mph');
pathIn  = dataSet.currentpathname;
fileIn  = [pathIn nameIn];

if dataSet.custom
    button = questdlg('Save custom machine?','SELECT','Yes','Cancel','Yes');
    if strcmp(dataSet.TypeOfRotor,'EESM')
        warning('Mechanical analysis not working with custom geometry')
    end
end

if dataSet.custom==0 || (isequal(button,'Yes') && (dataSet.custom)) 
    if nargin==1
        %[filename,pathname] = uiputfile([pathIn nameIn],'Input machine name and location');
        % if ~filename
        %     error('No file name selected');
        % end
        filename = nameIn;
        pathname = pathIn;
    else
        filename = strrep(filename,'.mat','.mph');
    end
    if ~isfolder(strcat(pathIn, strrep(nameIn,'.mph','_Comsol')))
        mkdir(pathname, strrep(nameIn,'.mph','_Comsol'));
    end

    [~, ~, geo, per, mat] = data0(dataSet);
    RQ = dataSet.RQ;
    
    fileans   = strrep(filename,'.mph','.ans');
    
    geo.custom = dataSet.custom;

    eval_type = 'singt';
end

fem = dimMesh(geo,eval_type);
[~,BLKLABLESrotTemp,geo,~] = ROTmatr(geo,fem,mat);
[geo,~,BLKLABLESstatTemp] = STATmatr(geo,fem);
geo.BLKLABELS.rotore  = BLKLABLESrotTemp;
geo.BLKLABELS.statore = BLKLABLESstatTemp;

if dataSet.custom
    if isequal(button,'Yes')
        [model,geo] = draw_motor_in_COMSOL(geo,mat,pathIn,nameIn);
                
%         if ~strcmp(fileIn,[pathname filename])
            fileTmp = [cd '\tmp\' filename];  
            copyfile(fileIn , fileTmp);
            copyfile(fileTmp,[pathname filename])
            delete (fileTmp)
%         end
                
        if isfile([pathname fileans])
            delete ([pathname fileans])
        end
         
        [model,geo] = draw_motor_in_COMSOL(geo,mat,pathIn,nameIn);

    else
        disp('Custom machine not saved')
    end
    
else
    [model,geo] = draw_motor_in_COMSOL(geo,mat,pathIn,nameIn);
end

if dataSet.custom==0 || (isequal(button,'Yes') && (dataSet.custom)) 
    geo.RQ = RQ;
    
    filename = strrep(filename,'mph','mat');
    dataSet.currentpathname = [pathname '\'];
    dataSet.currentfilename = filename;
    dataSet.slidingGap      = 1; % R347
    if dataSet.Qs==6*dataSet.Num3PhaseCircuit*dataSet.NumOfSlots*dataSet.NumOfPolePairs
        dataSet.slidingGap = 0;
    end
    
    % refresh GUI display data
    % if flagGUI
    %     set(app.currentMotFileName,'Value',filename);  % update display
    % %     load([pathname filename]);
    
    
    % end
    dataSet.RQ = round(dataSet.RQ,4);
    dataSet.currentpathname = pathname;
    dataSet.currentfilename = filename;
    
    geo = orderfields(geo);
    per = orderfields(per);
    dataSet = orderfields(dataSet);
    mat = orderfields(mat);
    save([pathname filename],'geo','per','dataSet','mat');
end
% cd(currentDir);

dataSet.infoComsol = geo.infoComsol;

disp('Done!')
