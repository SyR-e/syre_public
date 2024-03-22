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

function dataSet = DrawAndSaveMachine_JMAG(dataSet,filename,pathname)
% dataSet.currentpathname='D:\syre_maedeh\motorExamples\';
% dataSet.currentfilename='THOR.mat';%syreDefaultMotor,THOR,SPM1

% ==== TAKE THE DIR AND THE STRUCT DATASET ===============================
pathname = dataSet.currentpathname;
filename = dataSet.currentfilename;
% warning('Warning!!! Debug mode active!!!!!!!')
% pathname = [cd '\motorExamples\'];
load([pathname filename]);
% '#model name
model_name= strrep(dataSet.currentfilename,'.mat','.jmag');

%% % ########################Open JMAG Designer############################
 JMAGversion='222'; 
% '# Create an "app" application object to launch JMAG-Designer
JDesigner = actxserver(strcat('designer.Application.',JMAGversion));
JDesigner.Show(); % Show the JMAG interface
%'Get the installation folder of JMAG-Designer
JDesigner.GetAppDir();
% '# Create a new project.
JDesigner.NewProject('Untitled');
% '________________________________________________________________________
% '# Set name of project and save it.
JDesigner.NewProject(strrep(dataSet.currentfilename,'.mat','.jmag'));
JDesigner.SaveAs(strcat(pathname,'\',strrep(dataSet.currentfilename,'.mat','.jmag'),'.jproj'));
%% % ########################Create Motor Geometry#########################
draw_motor_in_JMAG(geo,mat,model_name,JDesigner)

dataSet.currentpathname = [pathname];
dataSet.currentfilename = filename;

disp('JMAG model saved!')