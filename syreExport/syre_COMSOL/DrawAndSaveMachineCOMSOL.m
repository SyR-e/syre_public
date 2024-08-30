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

function dataSet = DrawAndSaveMachineCOMSOL(dataSet,filename,pathname)
% function dataSet = DrawPushMachine(handles,filename,pathname)

% if ~isfield(app,'dataSet')  % not launched from GUI_Syre
%     flagGUI=0;
%     dataSet=app;
% else                            % launched from GUI_Syre
%     flagGUI=1;
%     dataSet=app.dataSet;
% end

nameIn  = dataSet.currentfilename;
nameIn  = strrep(nameIn,'.mat','.mph');
pathIn  = dataSet.currentpathname;
fileIn  = [pathIn nameIn];

if dataSet.custom
    button = questdlg('Save custom machine?','SELECT','Yes','Cancel','Yes');

end

if dataSet.custom==0 || (isequal(button,'Yes') && (dataSet.custom)) 
    if nargin==1
        [filename,pathname] = uiputfile([dataSet.currentpathname 'newmachine.mph'],'Input machine name and location');
        if ~filename
            error('No file name selected');
        end
    else
        filename = strrep(filename,'.mat','.mph');
    end
    
    [~, ~, geo, per, mat] = data0(dataSet);
    RQ = dataSet.RQ;
    
    fileans   = strrep(filename,'.mph','.ans');
    
    geo.custom = dataSet.custom;

    eval_type = 'singt';
end

if dataSet.custom
    if isequal(button,'Yes')
        [geo,mat] = draw_motor_in_COMSOL(geo,mat, pathname, filename);
                
%         if ~strcmp(fileIn,[pathname filename])
            fileTmp = [cd '\tmp\' filename];  
            copyfile(fileIn , fileTmp);
            copyfile(fileTmp,[pathname filename])
            delete (fileTmp)
%         end
                
        if isfile([pathname fileans])
            delete ([pathname fileans])
        end
         
        [geo,mat] = draw_motor_in_COMSOL(geo,mat, pathname, filename);

    else
        disp('Custom machine not saved')
    end
    
else
    [geo,mat] = draw_motor_in_COMSOL(geo,mat, pathname, filename);
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
