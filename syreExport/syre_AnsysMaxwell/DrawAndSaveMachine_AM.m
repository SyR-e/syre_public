% Copyright 2021
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


function dataSet = DrawAndSaveMachine_AM(dataSet)
%path motor
% [filename, filepath] = uigetfile('\*.mat');
% addpath(filepath);
% load([filepath filename])

filepath = dataSet.currentpathname;
filename = dataSet.currentfilename;
load([filepath filename]);

syreToDxfansys(geo.stator,geo.rotor,filepath,filename) %create dxf of selected motor

filepath      = strrep(filepath,'\','/');
currentFolder = pwd;
currentFolder = strrep(currentFolder,'\','/');
%iron python inside ansys path
%ipypath='"C:\Program Files\AnsysEM\AnsysEM20.1\Win64\common\IronPython\ipy64.exe" "';
%defipypath = 'C:\Program Files\AnsysEM\AnsysEM20.1\Win64\common\IronPython\ipy64.exe';
ipypath    = 'C:\Program Files\AnsysEM\v222\Win64\common\IronPython\';
ipy64exe   = 'ipy64.exe';
% [ipy64exe, ipypath] = uigetfile(defipypath,'Select ipy64.exe (Ansys) Directory');
ipypath=strcat('"',ipypath,ipy64exe,'" "');

%% Savings Variabiles 
n_pontR           = nnz(geo.pontR);                        %number of radial gap
R                 = geo.R;                                 %stator radius
r                 = geo.r;                                 %rotor radius
p                 = geo.p;                                 %pole pair
g                 = geo.g;                                 %air gap thickness
nlay              = geo.nlay;                              %layer of rotor barriers
q                 = geo.q;                                 %slot per cave per phase
l                 = geo.l;                                 %stack length
radial_ribs_split = nnz(geo.radial_ribs_split.*geo.pontR); %half added barrier area using split

Rotor_MatName    = mat.Rotor.MatName;    %rotor material
Stator_MatName   = mat.Stator.MatName;   %stator material
Shaft_MatName    = mat.Shaft.MatName;    %shaft material
SlotCond_MatName = mat.SlotCond.MatName; %slot conductor material
LayerMag_MatName = mat.LayerMag.MatName; %magnet material
LayerAir_MatName = mat.LayerAir.MatName; %air material

LayerMag_Hc      = mat.LayerMag.Hc;      %magnet Hc
LayerMag_mu      = mat.LayerMag.mu;      %magnet permeability
LayerMag_sigmaPM = mat.LayerMag.sigmaPM; %magnet conductivity
LayerMag_kgm3    = mat.LayerMag.kgm3;    %magnet density


SlotCond_sigma = mat.SlotCond.sigma; %cond conductivity
SlotCond_alpha = mat.SlotCond.alpha; %cond thermal coefficient
SlotCond_kgm3  = mat.SlotCond.kgm3;  %cond density

Shaft_alpha = mat.Shaft.alpha; %alfa loss coefficient shaft
Shaft_beta  = mat.Shaft.beta;  %beta loss coefficient shaft
Shaft_kh    = mat.Shaft.kh;    %hysteresis loss constant
Shaft_ke    = mat.Shaft.ke;    %eddy current loss constant
Shaft_BH    = mat.Shaft.BH;    %BH curve shaft
Shaft_kgm3  = mat.Shaft.kgm3;  %Shaft density;

Rotor_alpha = mat.Rotor.alpha; %alfa loss coefficient rotor
Rotor_beta  = mat.Rotor.beta;  %beta loss coefficient rotor
Rotor_kh    = mat.Rotor.kh;    %hysteresis loss constant
Rotor_ke    = mat.Rotor.ke;    %eddy current loss constant
Rotor_BH    = mat.Rotor.BH;    %BH curve Rotor
Rotor_kgm3  = mat.Rotor.kgm3;  %Rotor density;

Stator_alpha = mat.Stator.alpha; %alfa loss coefficient rotor
Stator_beta  = mat.Stator.beta;  %beta loss coefficient rotor
Stator_kh    = mat.Stator.kh;    %hysteresis loss constant
Stator_ke    = mat.Stator.ke;    %eddy current loss constant
Stator_BH    = mat.Stator.BH;    %BH curve Stator
Stator_kgm3  = mat.Stator.kgm3;  %Stator density;

BLKLABELS_rotore_xy = geo.BLKLABELS.rotore.xy; %coordinate materiali rotore


%% Export Variables 

fclose('all');
pyfile = fopen([currentFolder '\syreExport\syre_AnsysMaxwell\setup_draw_motor_in_ansys.txt'], 'wt' );

% coordinate system permanent magnets
[k,]     = find(BLKLABELS_rotore_xy(:,3)==6);
geo.n_PM = length(k);   %number Permanent Magnet
n_PM     = geo.n_PM;
PM_CS    = '"PM_CS":['; %init. coordinate system PMs


for ii = [1,2,6,7]
    PM_CS = strcat(PM_CS,'[');
    for jj = 1:length(k)
        PM_CS = strcat(PM_CS,num2str(BLKLABELS_rotore_xy(jj,ii)),',');    
    end
    PM_CS = strcat(PM_CS,'],');
end
PM_CS = strcat(PM_CS,'],');
 
%radial_split_ribs
export={
    'import sys,os'
    'import pickle'
    
    'geodata={'
    '"R":%f, #stator radius'
    '"r":%f, #rotor radius'
    '"p" :%d,   #pole pair'
    '"g" :%f,   #air gap thickness'
    '"nlay" :%d,   #layer of rotor barrier'
    '%s #Coordinate system magnet: 0-x 1-y 2-xvect 3-yvect'
    '"q":%d,  #slot per cave per phase'
    '"radial_ribs_split":%d, #area added by split barrier'
    '"n_PM":%d, #number of permanent magnet (number of magnetic segments)'
    '"l":%f, #stack length'
    '"filepath":"%s", #motor"s folder'
    '"filename":"%s", #motor"s file mat name'
    '}'
    '#struct material names'
    'material={'
    '"rotor":"%s",'
    '"stator":"%s",'
    '"shaft":"%s",'
    '"slotcond":"%s",'
    '"magnet":"%s",'
    '}'    
    
    "with open('%s/syreExport/syre_AnsysMaxwell/temp/temp.pkl', 'wb') as export:"
    '    pickle.dump([geodata,material],export,protocol=2)'
};

str = sprintf('%s\n',export{:});
fprintf(pyfile,str,R,r,p,g,nlay,PM_CS,q,radial_ribs_split,...
    n_PM,l,filepath,filename,Rotor_MatName,Stator_MatName,...
    Shaft_MatName,SlotCond_MatName,LayerMag_MatName,currentFolder);
clear export


%% Programm Start

start={
    'sys.path.append(r"C:/Program Files/AnsysEM/v222/Win64")'
    'sys.path.append(r"C:/Program Files/AnsysEM/v222/Win64/PythonFiles/DesktopPlugin")'
    'import ScriptEnv'
    ''
    'ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")'
    'oDesktop.RestoreWindow()'
    'oProject = oDesktop.NewProject()'
    'oProject.InsertDesign("Maxwell 2D", "Maxwell2DDesign1", "Transient", "")'
    'oDesign = oProject.SetActiveDesign("Maxwell2DDesign1")'
    'oEditor = oDesign.SetActiveEditor("3D Modeler")'
    'oProject.SaveAs("%s%s.aedt", True)'
    'oDesktop.ClearMessages("", "",3)'};
str=sprintf('%s\n',start{:});
fprintf(pyfile,str,filepath,erase(filename,".mat"));
clear start

%%%%%%%%%%%%%%%%Air%%%%%%%%%%%%%%

air={
    '#########Add materials to library'
    '#########Air##########'
    'oDefinitionManager = oProject.GetDefinitionManager()'
    'oDefinitionManager.AddMaterial('
    '	['
    '		"NAME:%s ",'
    '		"CoordinateSystemType:=", "Cartesian",'
    '		"BulkOrSurfaceType:="	, 1,'
    '		['
    '			"NAME:PhysicsTypes",'
    '			"set:="			, ["Electromagnetic"]'
    '		]'
    '	])'};
str=sprintf('%s\n',air{:});
fprintf(pyfile,str,LayerAir_MatName);
clear air

%%%%%%%%%%%%%Magnet%%%%%%%%%%%%%%%%%%%%

if LayerMag_Hc == 0
mag={''};
else
    if isfield(mat.LayerMag,'temp')   
        [LayerMag_Tcoeff,LayerMag_Href] = MagHc_temp_fitting(mat);
    else
        LayerMag_Tcoeff=0;
        LayerMag_Href=LayerMag_Hc;
    end
mag={
    '############Magnet############'
    'oDefinitionManager.AddMaterial('
    '	['
	'	"NAME:%s",'
	'	"CoordinateSystemType:=", "Cartesian",'
	'	"BulkOrSurfaceType:="	, 1,'
	'	['
	'		"NAME:PhysicsTypes",'
	'		"set:="			, ["Electromagnetic"]'
	'	],'
	'	['
	'		"NAME:ModifierData",'
	'		['
	'			"NAME:ThermalModifierData",'
	'			"modifier_data:="	, "thermal_modifier_data",'
	'			['
	'				"NAME:all_thermal_modifiers",'
	'				['
	'					"NAME:one_thermal_modifier",'
	'					"Property::="		, "magnetic_coercivity",'
	'					"Index::="		, 0,'
	'					"prop_modifier:="	, "thermal_modifier",'
	'					"use_free_form:="	, False,'
	'					"Tref:="		, "20cel",'
	'					"C1:="			, "%f",'
	'					"C2:="			, "0",'
	'					"TL:="			, "-273.15cel",'
	'					"TU:="			, "1000cel",'
	'					"auto_calculation:="	, True'
	'				]'
	'			]'
	'		]'
	'	],'
	'	"permeability:="	, "%f",'
	'	"conductivity:="	, "%f",'
	'	['
	'		"NAME:magnetic_coercivity",'
	'		"property_type:="	, "VectorProperty",'
	'		"Magnitude:="		, "-%fA_per_meter",'
	'		"DirComp1:="		, "1",'
	'		"DirComp2:="		, "0",'
	'		"DirComp3:="		, "0"'
	'	],'
	'	"mass_density:="	, "%f"'
	'])'
    };


str = sprintf('%s\n',mag{:});

fprintf(pyfile,str,LayerMag_MatName,LayerMag_Tcoeff,LayerMag_mu,LayerMag_sigmaPM,LayerMag_Href,...
        LayerMag_kgm3);
    
clear mag
end
%%%%%%%%%%%%%%%%%%%% Slot Conductor %%%%%%%%%%%%%%%%%%%%%%%%%

slotcond={
    '########## Slot conductor ##########'
    'oDefinitionManager.AddMaterial('
	'['
	'	"NAME:%s",'
	'	"CoordinateSystemType:=", "Cartesian",'
	'	"BulkOrSurfaceType:="	, 1,'
	'	['
	'		"NAME:PhysicsTypes",'
	'		"set:="			, ["Electromagnetic","Thermal"]'
	'	],'
	'	"conductivity:="	, "%f",'
	'	"thermal_conductivity:=", "%f",'
	'	"mass_density:="	, "%f"'
	'])'
    };
str = sprintf('%s\n',slotcond{:});
fprintf(pyfile,str,SlotCond_MatName,SlotCond_sigma,SlotCond_alpha*1e5,...
    SlotCond_kgm3);
clear slotcond

%%%%%%%%%%%%%%%%%%% Shaft %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Shaft_MatName,'ShaftAir')
    shaft={
    '#########Add materials to library'
    '#########ShaftAir##########'
    'oDefinitionManager = oProject.GetDefinitionManager()'
    'oDefinitionManager.AddMaterial('
    '	['
    '		"NAME:%s ",'
    '		"CoordinateSystemType:=", "Cartesian",'
    '		"BulkOrSurfaceType:="	, 1,'
    '		['
    '			"NAME:PhysicsTypes",'
    '			"set:="			, ["Electromagnetic"]'
    '		]'
    '	])'};
str = sprintf('%s\n',shaft{:});
fprintf(pyfile,str,Shaft_MatName);
clear shaft

else    
[Kh,Kadd,Ke]=losscoefficients(Shaft_alpha,Shaft_beta,Shaft_kh,Shaft_ke,Shaft_kgm3);
shaft={
    '############ Shaft plates ###########'
    'oDefinitionManager.AddMaterial('
	'['
	'	"NAME:%s",'
	'	"CoordinateSystemType:=", "Cartesian",'
	'	"BulkOrSurfaceType:="	, 1,'
	'	['
	'		"NAME:PhysicsTypes",'
	'		"set:="			, ["Electromagnetic"]'
	'	],'
	'	['
	'		"NAME:permeability",'
	'		"property_type:="	, "nonlinear",'
	'		"BTypeForSingleCurve:="	, "normal",'
	'		"HUnit:="		, "A_per_meter",'
	'		"BUnit:="		, "tesla",'
	'		"IsTemperatureDependent:=", False,'
	'		['
	'			"NAME:BHCoordinates",'
	'			['
	'				"NAME:DimUnits", '
	'				"", '
	'				""'
	'			],'
    };
str = sprintf('%s\n',shaft{:});
fprintf(pyfile,str,Shaft_MatName);
clear shaft

for ii=1:length(Shaft_BH(:,1))
shaft={
	'			["NAME:Coordinate",["NAME:CoordPoint",%f,%f]],'			
    };
str = sprintf('%s\n',shaft{:});
fprintf(pyfile,str,Shaft_BH(ii,2),Shaft_BH(ii,1));
clear shaft
end
    clear shaft
    shaft={
	'		],'
	'		['
	'			"NAME:Temperatures"'
	'		]'
	'	],'
	'	['
	'		"NAME:magnetic_coercivity",'
	'		"property_type:="	, "VectorProperty",'
	'		"Magnitude:="		, "0A_per_meter",'
	'		"DirComp1:="		, "1",'
	'		"DirComp2:="		, "0",'
	'		"DirComp3:="		, "0"'
	'	],'
	'	['
	'		"NAME:core_loss_type",'
	'		"property_type:="	, "ChoiceProperty",'
	'		"Choice:="		, "Electrical Steel"'
	'	],'
	'	"core_loss_kh:="	, "%f",'
	'	"core_loss_kc:="	, "%f",'
	'	"core_loss_ke:="	, "%f",'
	'	"core_loss_kdc:="	, "0",'
	'	"mass_density:="	, "%f",'
	'	"core_loss_equiv_cut_depth:=", "0.001meter"'
	'])'
    
};
str=sprintf('%s\n',shaft{:});
fprintf(pyfile,str,Kh,Ke,Kadd,Shaft_kgm3);
clear shaft
end

%%%%%%%%%%%%%%%%%%% Rotor's plates %%%%%%%%%%%%%%%%%%%
if strcmp(Shaft_MatName,Rotor_MatName)==0
[Kh,Kadd,Ke] = losscoefficients(Rotor_alpha,Rotor_beta,Rotor_kh,Rotor_ke,Rotor_kgm3);
rotor={
    '############ Rotor"s plates ###########'
    'oDefinitionManager.AddMaterial('
	'['
	'	"NAME:%s",'
	'	"CoordinateSystemType:=", "Cartesian",'
	'	"BulkOrSurfaceType:="	, 1,'
	'	['
	'		"NAME:PhysicsTypes",'
	'		"set:="			, ["Electromagnetic"]'
	'	],'
	'	['
	'		"NAME:permeability",'
	'		"property_type:="	, "nonlinear",'
	'		"BTypeForSingleCurve:="	, "normal",'
	'		"HUnit:="		, "A_per_meter",'
	'		"BUnit:="		, "tesla",'
	'		"IsTemperatureDependent:=", False,'
	'		['
	'			"NAME:BHCoordinates",'
	'			['
	'				"NAME:DimUnits", '
	'				"", '
	'				""'
	'			],'
    };
str = sprintf('%s\n',rotor{:});
fprintf(pyfile,str,Rotor_MatName);
clear rotor
for ii=1:length(Rotor_BH(:,1))
rotor={
	'			["NAME:Coordinate",["NAME:CoordPoint",%f,%f]],'		
    };
str=sprintf('%s\n',rotor{:});
fprintf(pyfile,str,Rotor_BH(ii,2),Rotor_BH(ii,1));
clear rotor
end
    clear rotor
    rotor={
	'		],'
	'		['
	'			"NAME:Temperatures"'
	'		]'
	'	],'
	'	['
	'		"NAME:magnetic_coercivity",'
	'		"property_type:="	, "VectorProperty",'
	'		"Magnitude:="		, "0A_per_meter",'
	'		"DirComp1:="		, "1",'
	'		"DirComp2:="		, "0",'
	'		"DirComp3:="		, "0"'
	'	],'
	'	['
	'		"NAME:core_loss_type",'
	'		"property_type:="	, "ChoiceProperty",'
	'		"Choice:="		, "Electrical Steel"'
	'	],'
	'	"core_loss_kh:="	, "%f",'
	'	"core_loss_kc:="	, "%f",'
	'	"core_loss_ke:="	, "%f",'
	'	"core_loss_kdc:="	, "0",'
	'	"mass_density:="	, "%f",'
	'	"core_loss_equiv_cut_depth:=", "0.001meter"'
	'])'
    
};
str = sprintf('%s\n',rotor{:});
fprintf(pyfile,str,Kh,Ke,Kadd,Rotor_kgm3);
clear rotor
end
%%%%%%%%%%%%%%%%%%% Stator's plates %%%%%%%%%%%%%%
if strcmp(Stator_MatName,Rotor_MatName)==0 &&...
        strcmp(Stator_MatName,Shaft_MatName)==0
[Kh,Kadd,Ke] = losscoefficients(Stator_alpha,Stator_beta,Stator_kh,Stator_ke,Stator_kgm3);
stator={
    '############ Stator"s plates ###########'
    'oDefinitionManager.AddMaterial('
	'['
	'	"NAME:%s",'
	'	"CoordinateSystemType:=", "Cartesian",'
	'	"BulkOrSurfaceType:="	, 1,'
	'	['
	'		"NAME:PhysicsTypes",'
	'		"set:="			, ["Electromagnetic"]'
	'	],'
	'	['
	'		"NAME:permeability",'
	'		"property_type:="	, "nonlinear",'
	'		"BTypeForSingleCurve:="	, "normal",'
	'		"HUnit:="		, "A_per_meter",'
	'		"BUnit:="		, "tesla",'
	'		"IsTemperatureDependent:=", False,'
	'		['
	'			"NAME:BHCoordinates",'
	'			['
	'				"NAME:DimUnits", '
	'				"", '
	'				""'
	'			],'
    };

str = sprintf('%s\n',stator{:});
fprintf(pyfile,str,Stator_MatName);
clear stator
for ii=1:length(Stator_BH(:,1))
stator = {
	'			["NAME:Coordinate",["NAME:CoordPoint",%f,%f]],'		
    };
str = sprintf('%s\n',stator{:});
fprintf(pyfile,str,Stator_BH(ii,2),Stator_BH(ii,1));
clear stator
end
    clear stator
    stator={
	'		],'
	'		['
	'			"NAME:Temperatures"'
	'		]'
	'	],'
	'	['
	'		"NAME:magnetic_coercivity",'
	'		"property_type:="	, "VectorProperty",'
	'		"Magnitude:="		, "0A_per_meter",'
	'		"DirComp1:="		, "1",'
	'		"DirComp2:="		, "0",'
	'		"DirComp3:="		, "0"'
	'	],'
	'	['
	'		"NAME:core_loss_type",'
	'		"property_type:="	, "ChoiceProperty",'
	'		"Choice:="		, "Electrical Steel"'
	'	],'
	'	"core_loss_kh:="	, "%f",'
	'	"core_loss_kc:="	, "%f",'
	'	"core_loss_ke:="	, "%f",'
	'	"core_loss_kdc:="	, "0",'
	'	"mass_density:="	, "%f",'
	'	"core_loss_equiv_cut_depth:=", "0.001meter"'
	'])'
    
};

str = sprintf('%s\n',stator{:});
fprintf(pyfile,str,Kh,Ke,Kadd,Stator_kgm3);
clear stator
end


runpy = {'sys.exit()'};
str = sprintf('\n%s\n',runpy{:});
fprintf(pyfile,str,[currentFolder '\syreExport\syre_AnsysMaxwell']);
clear runpy

fclose(pyfile);
clear pyfile

%from txt to py
movefile([currentFolder '\syreExport\syre_AnsysMaxwell\setup_draw_motor_in_ansys.txt'],[currentFolder '\syreExport\syre_AnsysMaxwell\setup_draw_motor_in_ansys.py'],'f');

%% Create closing project script
fclose("all");

pyfile = fopen([currentFolder '\syreExport\syre_AnsysMaxwell\close_projectAM.txt'], 'wt');

closer = {
    'import sys'
    'sys.path.append(r"C:/Program Files/AnsysEM/v222/Win64")'
    'sys.path.append(r"C:/Program Files/AnsysEM/v222/Win64/PythonFiles/DesktopPlugin")'
    'import ScriptEnv'
    'ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")'
    'oDesktop.RestoreWindow()'
    'oDesktop.CloseProject("%s")'
    };

str = sprintf('%s\n',closer{:});
fprintf(pyfile,str,erase(filename,".mat"));
clear closer

fclose(pyfile);

movefile([currentFolder '\syreExport\syre_AnsysMaxwell\close_projectAM.txt'],[currentFolder '\syreExport\syre_AnsysMaxwell\close_projectAM.py'],'f');

%Ansys start e py run
command = strcat(ipypath,strrep([currentFolder '\syreExport\syre_AnsysMaxwell'],'/','\'),'\setup_draw_motor_in_ansys.py"');
[status,cmdout] = dos(command);
if status ~= 0
    fprintf('An error occurred - program aborted - \n %s',cmdout)
end

command = strcat(ipypath,strrep([currentFolder '\syreExport\syre_AnsysMaxwell'],'/','\'),'\draw_motor_in_ansys.py"');

[status,cmdout] = dos(command);

if status ~= 0
    fprintf('An error occurred - program aborted - \n %s',cmdout)
end
