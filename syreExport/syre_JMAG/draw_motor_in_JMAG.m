
%  Copyright 2024
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function draw_motor_in_JMAG(geo,mat,model_name,JDesigner)

% '######################## Draw  CAD Geometry ############################
% AreaCodes %load the area codes for stator geometry creation
indCoilHalf  = 0;
indCoil = 1;
indWedge  = geo.Qs+1;
% % CAD Parameters
sid=2*(geo.r+geo.g);
sod=2*geo.R;
nslot=geo.Qs*geo.p*2;
rsh=geo.Ar;
rrout=geo.r;
alpha_s=2*pi/(2*nslot);% angle betwee two slots
alpha_r=2*pi/(2*geo.p);%pole pitch 
Wedge=geo.stator(geo.stator(:,end)==indWedge,:);
Coil=geo.stator(geo.stator(:,end)==indCoil,:);
CoilHalf=geo.stator(geo.stator(:,end)==indCoilHalf,:);
% '------------------------------------------------------------------------
% 'Mesh Size Control (mm):
% ' Motor Part mesh size
stator_mesh_size=1.0;   % 'Stator mesh size
coil_mesh_size=3.0;     % 'Coil mesh size
rotor_mesh_size=1.0;    % 'Rotor mesh size
magnet_mesh_size=1.50;  % 'Magnet mesh size
% 'Airgap mesh size control using Slide Mesh Division
RDivision = 3;    % 'Radial Division in the airgap 
mech_angle=2*360/(geo.p*2);
CDivision = mech_angle*gcd(nslot,geo.p*2);  % 'Circumferential Division in the airgap (to consider the lowest subharmonic)
% '------------------------------------------------------------------------
% 'Open Geometry Editor interface
JDesigner.LaunchGeometryEditor();
geomApp = JDesigner.CreateGeometryEditor('true');
% '------------------------------------------------------------------------
% 'Create the CAD Parameters:
CADEquationTable = geomApp.GetDocument().GetDesignTable();
% '#stator slots per pole
CreateCADParameters_JMAG(CADEquationTable,'Qs', 0, strcat(num2str(geo.Qs)), 'stator slot Number per pole');
% '#rotor poles number: Poles
CreateCADParameters_JMAG(CADEquationTable,'P', 0, strcat(num2str(geo.p*2)), 'rotor pole number');
% '#Angle between 2 stator slots: alpha_s
CreateCADParameters_JMAG(CADEquationTable,'alpha_s', 1, '360.0/(P*Qs)', 'Mechanical Angle between 2 stator slots (Degrees)');
% '#_______________________________________________________________________
% '#Create the sketch of the Stator
%part indexes for stator,coil,rotor and magnet are 0,1,2,3
ref1 = geomApp.GetDocument().GetAssembly().GetItem(0);
ref2 = geomApp.GetDocument().CreateReferenceFromItem(ref1);
sketch_stator = geomApp.GetDocument().GetAssembly().CreateSketch(ref2);
sketch_stator.OpenSketch();
sketch_stator.SetProperty('Name', 'Stator');
sketch_stator.SetProperty('Color', 'fuchsia');
sketch_stator.SetProperty('UseParentColor', 0);
% '------------------------------------------------------------------------
% 'Define dimensions of the Stator 
% 'S0 = R0 origin of cordinate system (0,0)
S0x=0.0;
S0y=0.0;
% 'S1
S1x=sid/2.0;
S1y=0.0;
% 'S2
S2x=sod/2.0;
S2y=0.0;
% 'S3
S3x=(sod/2.0)*cos(alpha_s);
S3y=(sod/2.0)*sin(alpha_s);
% 'S4
S4x=Coil(6,1);
S4y=Coil(6,2);
% 'S5
S5x=Coil(7,3);
S5y=Coil(7,4);
% 'S6
S6x=Coil(7,5);
S6y=Coil(7,6);
% 'S7
S7x=Coil(7,1);
S7y=Coil(7,2);
% 'S8
S8x=Coil(8,3);
S8y=Coil(8,4);
% 'S9
S9x=Coil(9,3);
S9y=Coil(9,4);
% 'S10
S10x=Wedge(2,3);
S10y=Wedge(2,4);
 % '------------------------------------------------------------------------
% '#Create a vertex
% S0
S0=sketch_stator.CreateVertex(S0x, S0y);
S0.SetProperty('Name', 'S0');
ref_S0 = geomApp.GetDocument().CreateReferenceFromItem(S0);
geomApp.GetDocument().GetSelection().Clear();
geomApp.View().Xy();

% '#Create the Lines and Arc
% '#Line S1S2
Line_S1S2=sketch_stator.CreateLine(S1x, S1y, S2x, S2y);
Line_S1S2.SetProperty('Name', 'S1S2');
% '#Arc S2S3
Arc_S2S3=sketch_stator.CreateArc(S0x, S0y, S2x, S2y, S3x, S3y);
Arc_S2S3.SetProperty('Name', 'S2S3');
% '#Line S3S4
Line_S3S4=sketch_stator.CreateLine(S3x, S3y, S4x, S4y);
Line_S3S4.SetProperty('Name', 'S3S4');
ref_S3S4 = geomApp.GetDocument().CreateReferenceFromItem(Line_S3S4);
% '#Line S4S5
Line_S4S5=sketch_stator.CreateLine(S4x, S4y, S5x, S5y);
Line_S4S5.SetProperty('Name', 'S4S5');
% '#Arc S5S6
Arc_S5S6=sketch_stator.CreateArc(S7x, S7y, S6x, S6y, S5x, S5y);
Arc_S5S6.SetProperty('Name', 'S5S6');
% '#Line S6S8
Line_S6S8=sketch_stator.CreateLine(S6x, S6y, S8x, S8y);
Line_S6S8.SetProperty('Name', 'S6S8');
% '#Line S8S9
Line_S8S9=sketch_stator.CreateLine(S8x, S8y, S9x, S9y);
Line_S8S9.SetProperty('Name', 'S8S9');
% '#Line S9S10
Line_S9S10=sketch_stator.CreateLine(S9x, S9y, S10x, S10y);
Line_S9S10.SetProperty('Name', 'S9S10');
% '#Arc S10S1
Arc_S10S1=sketch_stator.CreateArc(S0x, S0y, S1x, S1y, S10x, S10y);
Arc_S10S1.SetProperty('Name', 'S10S1');

geomApp.GetDocument().GetSelection().Clear();

geomApp.View().Xy();

% '############### Regions of the Stator ##################################
% 
% '#Stator Core Region (S1-S2-S3-S4-S5-S6-S7-S8-S1)
geomApp.GetDocument().GetSelection().Add(Line_S1S2);
geomApp.GetDocument().GetSelection().Add(Arc_S2S3);
geomApp.GetDocument().GetSelection().Add(Line_S3S4);
geomApp.GetDocument().GetSelection().Add(Line_S4S5);
geomApp.GetDocument().GetSelection().Add(Arc_S5S6);
geomApp.GetDocument().GetSelection().Add(Line_S6S8);
geomApp.GetDocument().GetSelection().Add(Line_S8S9);
geomApp.GetDocument().GetSelection().Add(Line_S9S10);
geomApp.GetDocument().GetSelection().Add(Arc_S10S1);

sketch_stator.CreateRegions();

geomApp.GetDocument().GetSelection().Clear();

sketch_RegionItem = geomApp.GetDocument().GetAssembly().GetItem('Stator');
%Obtains a Region by selecting a name
stator_region = sketch_RegionItem.GetItem('Region');
stator_region.SetProperty('Name', 'Stator');
stator_region.SetProperty('Color', 'fuchsia');
stator_region.SetProperty('UseParentColor', 0);
%‘Creates an object for reference using a Region object
Ref_statorRegion = geomApp.GetDocument().CreateReferenceFromItem(stator_region);
%Adds an object ID for reference to an array.
refarray = Ref_statorRegion.GetIdentifier();
%refarray = 'faceregion(TRegionItem105)';
 
 
% '#Create Partial of full model using the symmetry, mirror... functions
% '# Partial or full stator model 
 
% '#Region Mirror Copy of stator core region
RegionMirror=sketch_stator.CreateRegionMirrorCopy();
RegionMirror.SetProperty('Name', 'RegionMirror_Stator');
RegionMirror.SetProperty('Merge', 1);
RegionMirror.SetProperty('SymmetryType', 0);
RegionMirror.SetPropertyByReference('Symmetry', ref_S3S4);
RegionMirror.SetProperty('Region', refarray);

geomApp.GetDocument().GetSelection().Clear();

% '#Region Radial Pattern of stator core region
RegionRadialPatt=sketch_stator.CreateRegionRadialPattern();
RegionRadialPatt.SetProperty('Name', 'RadialPattern_stator');
RegionRadialPatt.SetProperty('Merge', 1);
RegionRadialPatt.SetProperty('Region', refarray);
RegionRadialPatt.SetPropertyByReference('Center', ref_S0);
RegionRadialPatt.SetProperty('Angle', 'alpha_s');
RegionRadialPatt.SetProperty('Instance', 'Qs');
RegionRadialPatt.SetProperty('MergeAxis', 1);

geomApp.GetDocument().GetSelection().Clear();

sketch_stator.CloseSketch()
% '#_______________________________________________________________________
% '#Create the sketch of the Coil
ref1 = geomApp.GetDocument().GetAssembly().GetItem(1);
ref2 = geomApp.GetDocument().CreateReferenceFromItem(ref1);
sketch_coil = geomApp.GetDocument().GetAssembly().CreateSketch(ref2);
sketch_coil.OpenSketch();
sketch_coil.SetProperty('Name', 'Coil');
sketch_coil.SetProperty('Color', 'green');
sketch_coil.SetProperty('UseParentColor', 0);
% Lower Coils
% '------------------------------------------------------------------------
% 'Define dimensions of the coil
% 'C0 = S0 origin of cordinate system (0,0)
C0x=0.0;
C0y=0.0; 
% 'C1
C1x=CoilHalf(1,3);
C1y=CoilHalf(1,4);
% 'C2
C2x=(CoilHalf(1,1)+CoilHalf(1,3))/2;
C2y=(CoilHalf(1,2)+CoilHalf(1,4))/2;
% 'C3
C3x=S4x;
C3y=S4y;
% 'C4
C4x=S5x;
C4y=S5y;
% 'C5
C5x=S6x;
C5y=S6y;
% 'C6
C6x=S7x;
C6y=S7y;
% 'C7
C7x=Coil(10,3);
C7y=Coil(10,4);
% 'C8
C8x=S9x;
C8y=S9y;
% 'C9
C9x=S8x;
C9y=S8y;
% '#Create a vertex
% C0
% C0
C0=sketch_coil.CreateVertex(C0x, C0y);
C0.SetProperty('Name', 'C0');
ref_C0 = geomApp.GetDocument().CreateReferenceFromItem(C0);
geomApp.GetDocument().GetSelection().Clear();
geomApp.View().Xy()

% '#Create the Lines and Arc
%for lower layer
% '#Line C1C2
Line_C1C2=sketch_coil.CreateLine(C1x, C1y, C2x, C2y);
Line_C1C2.SetProperty('Name', 'C1C2');
% '#Line C2C3
Line_C2C3=sketch_coil.CreateLine(C2x, C2y, C3x, C3y);
Line_C2C3.SetProperty('Name', 'C2C3');
ref_C2C3 = geomApp.GetDocument().CreateReferenceFromItem(Line_C2C3);
% '#Line C3C4
Line_C3C4=sketch_coil.CreateLine(C3x, C3y, C4x, C4y);
Line_C3C4.SetProperty('Name', 'C3C4');
% '#Arc C4C5
Arc_C4C5=sketch_coil.CreateArc(C6x, C6y, C5x, C5y, C4x, C4y);
Arc_C4C5.SetProperty('Name', 'C4C5');
% '#Line C5C1
Line_C5C1=sketch_coil.CreateLine(C5x, C5y, C1x, C1y);
Line_C5C1.SetProperty('Name', 'C5C1');
%for upper layer
% '#Line C1C2p
Line_C1C2p=sketch_coil.CreateLine(C1x, C1y, C2x, C2y);
Line_C1C2p.SetProperty('Name', 'C1C2p');
% '#Line C2C7
Line_C2C7=sketch_coil.CreateLine(C2x, C2y, C7x, C7y);
Line_C2C7.SetProperty('Name', 'C2C7');
ref_C2C7 = geomApp.GetDocument().CreateReferenceFromItem(Line_C2C7);
% '#Line C7C8
Line_C7C8=sketch_coil.CreateLine(C7x, C7y, C8x, C8y);
Line_C7C8.SetProperty('Name', 'C7C8');
% '#Line C8C9
Line_C8C9=sketch_coil.CreateLine(C8x, C8y, C9x, C9y);
Line_C8C9.SetProperty('Name', 'C8C9');
% '#Line C1C9
Line_C1C9=sketch_coil.CreateLine(C1x, C1y, C9x, C9y);
Line_C1C9.SetProperty('Name', 'C1C9');
% '########################################################################

 
% '########################################################################
% '####################################### Regions of the Coil ############
% '#lower Coil Region (C1-C2-C3-C4-C5-C1)
geomApp.GetDocument().GetSelection().Add(Line_C1C2);
geomApp.GetDocument().GetSelection().Add(Line_C2C3);
geomApp.GetDocument().GetSelection().Add(Line_C3C4);
geomApp.GetDocument().GetSelection().Add(Arc_C4C5);
geomApp.GetDocument().GetSelection().Add(Line_C5C1);
sketch_coil.CreateRegions();
% '#upper Coil Region (C1-C2-C7-C8-C9-C1)
geomApp.GetDocument().GetSelection().Add(Line_C1C2p);
geomApp.GetDocument().GetSelection().Add(Line_C2C7);
geomApp.GetDocument().GetSelection().Add(Line_C7C8);
geomApp.GetDocument().GetSelection().Add(Line_C8C9);
geomApp.GetDocument().GetSelection().Add(Line_C1C9);
sketch_coil.CreateRegions();

geomApp.GetDocument().GetSelection().Clear();

sketch_RegionItem = geomApp.GetDocument().GetAssembly().GetItem('Coil');
%Obtains a Region by selecting a name
coil_region = sketch_RegionItem.GetItem('Region');
coil_region.SetProperty('Name', 'LCoil');
coil_region.SetProperty('Color', 'green');
coil_region.SetProperty('UseParentColor', 0);
coil_region1 = sketch_RegionItem.GetItem('Region.2');
coil_region1.SetProperty('Name', 'UCoil');
%‘Creates an object for reference using a Region object
Ref_coilRegion = geomApp.GetDocument().CreateReferenceFromItem(coil_region);
Ref_coilRegion1 = geomApp.GetDocument().CreateReferenceFromItem(coil_region1);
%Adds an object ID for reference to an array.
refarray = Ref_coilRegion.GetIdentifier();
refarray1 = Ref_coilRegion1.GetIdentifier();

% '########################################################################
% '#Create Partial of full model using the symmetry, mirror... functions
% '#Region Mirror Copy of coil region
RegionMirror=sketch_coil.CreateRegionMirrorCopy();
RegionMirror.SetProperty('Name', 'RegionMirror_coil');
RegionMirror.SetProperty('Merge', 1);
RegionMirror.SetProperty('SymmetryType', 0);
RegionMirror.SetPropertyByReference('Symmetry', ref_C2C3);
RegionMirror.SetProperty('Region', refarray);

RegionMirror=sketch_coil.CreateRegionMirrorCopy();
RegionMirror.SetProperty('Name', 'RegionMirror_coil1');
RegionMirror.SetProperty('Merge', 1);
RegionMirror.SetProperty('SymmetryType', 0);
RegionMirror.SetPropertyByReference('Symmetry', ref_C2C7);
RegionMirror.SetProperty('Region', refarray1);

geomApp.GetDocument().GetSelection().Clear();
% '#Region Radial Pattern of coil region
RegionRadialPatt=sketch_coil.CreateRegionRadialPattern();
RegionRadialPatt.SetProperty('Name', 'RadialPattern_coil');
RegionRadialPatt.SetProperty('Merge', 1);
RegionRadialPatt.SetProperty('Region', refarray);
RegionRadialPatt.SetPropertyByReference('Center', ref_C0);
RegionRadialPatt.SetProperty('Angle', 'alpha_s');
RegionRadialPatt.SetProperty('Instance', 'Qs');
RegionRadialPatt.SetProperty('MergeAxis', 1);

RegionRadialPatt=sketch_coil.CreateRegionRadialPattern();
RegionRadialPatt.SetProperty('Name', 'RadialPattern_coil1');
RegionRadialPatt.SetProperty('Merge', 1);
RegionRadialPatt.SetProperty('Region', refarray1);
RegionRadialPatt.SetPropertyByReference('Center', ref_C0);
RegionRadialPatt.SetProperty('Angle', 'alpha_s');
RegionRadialPatt.SetProperty('Instance', 'Qs');
RegionRadialPatt.SetProperty('MergeAxis', 1);

sketch_coil.CloseSketch()

geomApp.GetDocument().GetSelection().Clear();
geomApp.View().Xy()


% '#_______________________________________________________________________
% '#Create the sketch of the Rotor
ref1 = geomApp.GetDocument().GetAssembly().GetItem(2);
ref2 = geomApp.GetDocument().CreateReferenceFromItem(ref1);
sketch_rotor = geomApp.GetDocument().GetAssembly().CreateSketch(ref2);
sketch_rotor.OpenSketch();
sketch_rotor.SetProperty('Name', 'Rotor');
sketch_rotor.SetProperty('Color', 'fuchsia');
sketch_rotor.SetProperty('UseParentColor', 0);
% '------------------------------------------------------------------------
% 'Define dimensions of the Rotor                    
% 'R0 = origin of cordinate system (0,0)
R0x=0.0;
R0y=0.0;
% 'R1
R1x=rsh;
R1y=0.0;
% 'R2
R2x=rrout;
R2y=0.0;
% 'R3
R3x=rrout*cos(alpha_r);
R3y=rrout*sin(alpha_r);
% 'R4
R4x=rsh*cos(alpha_r);
R4y=rsh*sin(alpha_r);
% '#Create the Lines and Arc
% '#Line R1R2
Line_R1R2=sketch_rotor.CreateLine(R1x, R1y, R2x, R2y);
Line_R1R2.SetProperty('Name', 'R1R2');
% '#Arc R2R3
Arc_R2R3=sketch_rotor.CreateArc(R0x, R0y, R2x, R2y, R3x, R3y);
Arc_R2R3.SetProperty('Name', 'R2R3');
% '#Line R3R4
Line_R3R4=sketch_rotor.CreateLine(R3x, R3y, R4x, R4y);
Line_R3R4.SetProperty('Name', 'R3R4');
% '#Arc R4R1
Arc_R4R1=sketch_rotor.CreateArc(R0x, R0y, R1x, R1y, R4x, R4y);
Arc_R4R1.SetProperty('Name', 'R4R1');


% case codMatBar   % AirRotorBars
% '#Create lines and arcs       
materialCodes % load the material codes
AirRotorxy=geo.rotor(geo.rotor(:,end-1)==codMatAirRot,:);
AirRotorArea=unique(AirRotorxy(:,9)); %the area numbers for air rotor regions
occurance=histcounts(AirRotorxy(:,9)',[AirRotorArea', inf]);% Count occurrences of each unique value
lines = cell(1, length(AirRotorxy(:,1)));

for kk=1:length(AirRotorxy(:,1))

            switch AirRotorxy(kk,7)
                case 0
            if ~(AirRotorxy(kk,1) == AirRotorxy(kk,3) &&  AirRotorxy(kk,2) == AirRotorxy(kk,4))      
                    lines{kk} = sketch_rotor.CreateLine(AirRotorxy(kk,1), AirRotorxy(kk,2), AirRotorxy(kk,3), AirRotorxy(kk,4));
                    lines{kk}.SetProperty('Name', ['LR', num2str(kk), num2str(AirRotorxy(kk,9))]);
            end
                case 1
                    lines{kk} = sketch_rotor.CreateArc(AirRotorxy(kk,1), AirRotorxy(kk,2), AirRotorxy(kk,3), AirRotorxy(kk,4), AirRotorxy(kk,5), AirRotorxy(kk,6));
                    lines{kk}.SetProperty('Name', ['LR', num2str(kk) num2str(AirRotorxy(kk,9))]);
                case -1
                    lines{kk} = sketch_rotor.CreateArc(AirRotorxy(kk,1), AirRotorxy(kk,2), AirRotorxy(kk,5), AirRotorxy(kk,6), AirRotorxy(kk,3), AirRotorxy(kk,4));
                    lines{kk}.SetProperty('Name', ['LR', num2str(kk) num2str(AirRotorxy(kk,9))]);

            end
end

if length(AirRotorxy(:,1)) > 8  %condition for Motors except SPMs
% '#Create regions
geomApp.GetDocument().GetSelection().Add(Line_R1R2);
geomApp.GetDocument().GetSelection().Add(Arc_R2R3);
geomApp.GetDocument().GetSelection().Add(Line_R3R4);
geomApp.GetDocument().GetSelection().Add(Arc_R4R1);
for kkk=1:length(AirRotorxy(:,1))
        geomApp.GetDocument().GetSelection().Add(lines{kkk});
end
sketch_rotor.CreateRegions();
geomApp.GetDocument().GetSelection().Clear();

sketch_RegionItem = geomApp.GetDocument().GetAssembly().GetItem('Rotor');
%Obtains a Region by selecting a name
rotor_regiondel=sketch_RegionItem.GetItem('Region');
%‘Creates an object for reference using a Region object
geomApp.GetDocument().GetSelection().Add(rotor_regiondel);
geomApp.GetDocument().GetSelection().Delete()

for i = 2:length(AirRotorArea(:,1))
    regionName = ['Region.', num2str(i)];
    rotor_regiondel = sketch_RegionItem.GetItem(regionName);
    % Check if the region exists before attempting to delete
    if ~isempty(rotor_regiondel)
        geomApp.GetDocument().GetSelection().Add(rotor_regiondel);
        geomApp.GetDocument().GetSelection().Delete();
    end
end
%Obtains a Region by selecting a name
regionName = ['Region.', num2str(length(AirRotorArea(:,1))+1)];
rotor_region = sketch_RegionItem.GetItem(regionName);
rotor_region.SetProperty('Name', 'Rotor');
rotor_region.SetProperty('Color', 'fuchsia');
rotor_region.SetProperty('UseParentColor', 0);

sketch_rotor.CloseSketch()
geomApp.View().Xy()
%%%%%%%%%%%%%%%%%%%SPM Motors%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    % '#Create regions
    final=cumsum(occurance);init=cumsum([1 occurance(1,1:end-1)]);
    for kkk=1:length(AirRotorArea)
    for j=init(kkk):final(kkk)
        geomApp.GetDocument().GetSelection().Add(lines{j});
    end
    sketch_rotor.CreateRegions();
    end

geomApp.GetDocument().GetSelection().Add(Line_R1R2);
geomApp.GetDocument().GetSelection().Add(Arc_R2R3);
geomApp.GetDocument().GetSelection().Add(Line_R3R4);
geomApp.GetDocument().GetSelection().Add(Arc_R4R1);
sketch_rotor.CreateRegions();
%%%%%%%%%%%%%%%%%%%Subtraction extra regions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sketch_RegionItem = geomApp.GetDocument().GetAssembly().GetItem('Rotor');
%Obtains a Region by selecting a name
rotor_regionsub=sketch_RegionItem.GetItem('Region.3');

%‘Creates an object for reference using a Region object
geomApp.GetDocument().GetSelection().Add(rotor_regionsub);
RegionBoolean=sketch_RegionItem.CreateBoolean();

Ref_totalregion = geomApp.GetDocument().CreateReferenceFromItem(rotor_regionsub);
refarray = Ref_totalregion.GetIdentifier();
RegionBoolean.SetProperty("Region1", refarray);

subtract_region1 = sketch_RegionItem.GetItem('Region');
Ref_subtRegion1 = geomApp.GetDocument().CreateReferenceFromItem(subtract_region1);
refarray1 = Ref_subtRegion1.GetIdentifier();
RegionBoolean.SetProperty("Region2", refarray1)
RegionBoolean.SetProperty('DeleteRegion2', 1)


RegionBoolean1=sketch_RegionItem.CreateBoolean();
subtract_region2 = sketch_RegionItem.GetItem('Region.2');
Ref_subtRegion2 = geomApp.GetDocument().CreateReferenceFromItem(subtract_region2);
refarray2 = Ref_subtRegion2.GetIdentifier();
RegionBoolean1.SetProperty("Region1", refarray);
RegionBoolean1.SetProperty("Region2", refarray2);
RegionBoolean1.SetProperty('DeleteRegion2', 1)

end

sketch_rotor.CloseSketch()
geomApp.GetDocument().GetSelection().Clear();
geomApp.View().Xy()
% '#_______________________________________________________________________
Magnetxy=geo.rotor(geo.rotor(:,end-1)==codMatBar,:);
if any(unique(Magnetxy(:,9)) ~= 0)
% '#Create the sketch of the Magnet
ref1 = geomApp.GetDocument().GetAssembly().GetItem(3);
ref2 = geomApp.GetDocument().CreateReferenceFromItem(ref1);
sketch_magnet = geomApp.GetDocument().GetAssembly().CreateSketch(ref2);
sketch_magnet.OpenSketch();
sketch_magnet.SetProperty('Name', 'Magnet');
sketch_magnet.SetProperty('Color', 'red');
sketch_magnet.SetProperty('UseParentColor', 0);

% '------------------------------------------------------------------------
% 'Define dimensions of the Magnet
% '#Create lines and arcs       
% case codMatBar   % PM
Magnetxy=geo.rotor(geo.rotor(:,end-1)==codMatBar,:);
MagnetArea=unique(Magnetxy(:,9)); %the area numbers for magnet material
occurance=histcounts(Magnetxy(:,9)',[MagnetArea', inf]);% Count occurrences of each unique value

lines = cell(1, length(Magnetxy(:,1)));
for kk=1:length(Magnetxy(:,1))

            switch Magnetxy(kk,7)
                case 0
                    lines{kk} = sketch_magnet.CreateLine(Magnetxy(kk,1), Magnetxy(kk,2), Magnetxy(kk,3), Magnetxy(kk,4));
                    lines{kk}.SetProperty('Name', ['LM', num2str(kk), num2str(Magnetxy(kk,9))]);
                case 1
                    lines{kk} = sketch_magnet.CreateArc(Magnetxy(kk,1), Magnetxy(kk,2), Magnetxy(kk,3), Magnetxy(kk,4), Magnetxy(kk,5), Magnetxy(kk,6));
                    lines{kk}.SetProperty('Name', ['LM', num2str(kk) num2str(Magnetxy(kk,9))]);
                case -1
                    lines{kk} = sketch_magnet.CreateArc(Magnetxy(kk,1), Magnetxy(kk,2), Magnetxy(kk,5), Magnetxy(kk,6), Magnetxy(kk,3), Magnetxy(kk,4));
                    lines{kk}.SetProperty('Name', ['LM', num2str(kk) num2str(Magnetxy(kk,9))]);

            end
end

% '#Create regions
final=cumsum(occurance);init=cumsum([1 occurance(1,1:end-1)]);
for kkk=1:length(MagnetArea)
    for j=init(kkk):final(kkk)
        geomApp.GetDocument().GetSelection().Add(lines{j});
    end
    sketch_magnet.CreateRegions();

end

sketch_magnet.CloseSketch()
geomApp.View().Xy()
end
% % '# 'Import CAD Geometry into JMAG Designer
% '------------------------------------------------------------------------
JDesigner.ImportDataFromGeometryEditor();
% '------------------------------------------------------------------------
JDesigner.GetCurrentModel().SetName(model_name);
% Get the current model
model = JDesigner.GetCurrentModel();
JDesigner.Save();

% '#_______________________________________________________________________
% % '# 'Partname Assembly
% '########################### Get ID number of each part.
% Rotor
rowIndex = find(geo.BLKLABELS.rotore.xy(:, 3) == codMatFeRot, 1);
Protor = JDesigner.CreatePoint(geo.BLKLABELS.rotore.xy(rowIndex,1), geo.BLKLABELS.rotore.xy(rowIndex,2), 0);%shoulld be modified
rotor = model.GetPartByPosition(Protor);
rotor.SetName('Rotor');
%Stator
Pstator = JDesigner.CreatePoint(geo.BLKLABELS.statore.xy(4,1), geo.BLKLABELS.statore.xy(4,2), 0);
stator = model.GetPartByPosition(Pstator);
stator.SetName('Stator');

% '########################### Rotor Core group
model.GetGroupList().CreateGroup('Rotor Core');
model.GetGroupList().AddPartToGroup('Rotor Core', 'Rotor');
model.SetColor('Rotor Core', 'magenta');
% '########################### Magnet group
MagnetCG=geo.BLKLABELS.rotore.xy(geo.BLKLABELS.rotore.xy(:,3)==codMatBar,1:2);%center of gravity(G) of Magnet
% '########################### Stator Core group
model.GetGroupList().CreateGroup('Stator Core');
model.GetGroupList().AddPartToGroup('Stator Core', 'Stator');
model.SetColor('Stator Core', 'magenta');

% '########################### Coil groups
CoilCG=geo.BLKLABELS.statore.xy(geo.BLKLABELS.statore.xy(:,3)==codMatCu,1:2);%center of gravity(G) of Copper Coils
CoilCGX_reshape = reshape (CoilCG(:,1),2,geo.Qs);%for two-layer windings
CoilCGY_reshape = reshape (CoilCG(:,2),2,geo.Qs);%for two-layer windings

% '#"Coil ID Creation needed for further Post Processing
coil_ID=zeros(2*geo.Qs,1); %two-layer winding
for slot = 1:1:2*geo.Qs
    P = JDesigner.CreatePoint(CoilCGX_reshape(slot),CoilCGY_reshape(slot),0);
    coil = model.GetPartByPosition(P); 
    coil_ID(slot) = coil.ID(); 
end
% '########################### Create Groups for Coil U, W, and V
model.GetGroupList().CreateGroup(strcat('Coil', '_U'));
model.GetGroupList().CreateGroup(strcat('Coil', '_W'));
model.GetGroupList().CreateGroup(strcat('Coil', '_V'));
% '#"Coil U" groups
createCoilGroups_JMAG(model, JDesigner, CoilCGX_reshape, CoilCGY_reshape, geo, 'Coil','Up', +1, 'green');
createCoilGroups_JMAG(model, JDesigner, CoilCGX_reshape, CoilCGY_reshape, geo, 'Coil','Un', -1, 'green');
% '#"Coil W" groups
createCoilGroups_JMAG(model, JDesigner, CoilCGX_reshape, CoilCGY_reshape, geo, 'Coil','Wp', +3, 'maroon');
createCoilGroups_JMAG(model, JDesigner, CoilCGX_reshape, CoilCGY_reshape, geo, 'Coil','Wn', -3, 'maroon');
% '#"Coil V" groups
createCoilGroups_JMAG(model, JDesigner, CoilCGX_reshape, CoilCGY_reshape, geo, 'Coil','Vp', +2, 'yellow');
createCoilGroups_JMAG(model, JDesigner, CoilCGX_reshape, CoilCGY_reshape, geo, 'Coil','Vn', -2, 'yellow');

% '###########################SET for Rotor core
createSet_JMAG(model, 'Rotor Core', 'Rotor',MagnetCG);
% '###########################SET for Stator core
createSet_JMAG(model, 'Stator Core', 'Stator',MagnetCG);

% '#_______________________________________________________________________
%% % Close Geometry Editor
model.CloseCadLink(); 

% '------------------------------------------------------------------------
%% %'Create a new study into JMAG Designer
study_name = 'Load';% '#study name
study = model.CreateStudy('Transient2D', study_name);
JDesigner.SetCurrentStudy(study_name);
% '------------------------------------------------------------------------
% 'Study property settings: Full Circuit Conversion settings 
CircuitConnectType = 0;     % 0 = 'Series' : Use Series connection type
% '------------------------------------------------------------------------
% 'Study property settings: Linear Solver and Nonlinear Calculation
% 'Nonlinear Maximum iteration number
Nonlinear_MaxIteration=50;
% Nonlinear calculation convergence tolerance
Nonlinear_Tolerance=0.0001;
% Use High Speed Solver to short the nonlinear calculation time when using non-linear material for the transient analysis
Nonlinear_Speedup=1;

% 'Create the Design Parameters:
% '#'______________________________________________________________________
Design_parameters = study.GetDesignTable();
% 'rotor pole pair number
createDesignParameters_JMAG(Design_parameters,'PolePair',0,geo.p,'pole pair number')
% 'Periodic angle of the model
createDesignParameters_JMAG(Design_parameters,'model_periodic_angle',1,'(360.0/PolePair/2)','Periodic angle of the model (degrees)')
% '#'______________________________________________________________________
% 'Mesh Control parameters: Radial Division in the air gap
createDesignParameters_JMAG(Design_parameters,'RDivision',0,RDivision,'Radial mesh Division from airgap')

% 'Mesh Control parameters: Circumferential Division in the air gap
createDesignParameters_JMAG(Design_parameters,'CDivision',0,CDivision,'Circumferential mesh Division from airgap')

% 'Mesh Control parameters: Element Size of stator
createDesignParameters_JMAG(Design_parameters,'stator_mesh_size',0,stator_mesh_size,'Stator Mesh Element Size')

% 'Mesh Control parameters: Element Size of coil
createDesignParameters_JMAG(Design_parameters,'coil_mesh_size',0,coil_mesh_size,'Coil Mesh Element Size')

% 'Mesh Control parameters: Element Size of rotor
createDesignParameters_JMAG(Design_parameters,'rotor_mesh_size',0,rotor_mesh_size,'Rotor Mesh Element Size')

% 'Mesh Control parameters: Element Size of magnet
createDesignParameters_JMAG(Design_parameters,'magnet_mesh_size',0,magnet_mesh_size,'Magnet Mesh Element Size')
% '------------------------------------------------------------------------
% 'motor stack length (mm)
motor_stack_length=geo.l;% 'Motor stack length (mm)
createDesignParameters_JMAG(Design_parameters,'motor_stack_length',0,motor_stack_length,'motor stack length (mm)')
% '------------------------------------------------------------------------

% 'Study property settings
% '------------------------------------------------------------------------
StudyProperties = study.GetStudyProperties();
% '------------------------------------------------------------------------
% 'Full Model Conversion settings
FullModelConversionType = 1;  % 1 = 'On' : Convert in full CAD model
StudyProperties.SetValue('FullModelConversion', FullModelConversionType);
StudyProperties.SetValue('ConversionFactor', 1*1); % The conversion factor excluding periodic boundary

%Set the motor stack length
StudyProperties.SetValue('ModelThickness', 'motor_stack_length');
% '------------------------------------------------------------------------
% ' Full Circuit Conversion settings
CircuitConversionType = 1; % 1 = 'Periodic' : Convert in full circuit (Synchronize with periodic Boundary)
StudyProperties.SetValue('ConversionType ', CircuitConversionType);
StudyProperties.SetValue('CircuitConversionFactor ', 1*1); % Other than¨periodic boundary: The conversion factor excluding periodic boundary
StudyProperties.SetValue('CircuitConnect', CircuitConnectType); 
% '------------------------------------------------------------------------
% 'Nonlinear Calculation settings
StudyProperties.SetValue('NonlinearMaxIteration', Nonlinear_MaxIteration);
StudyProperties.SetValue('NonlinearTolerance', Nonlinear_Tolerance);
StudyProperties.SetValue('NonlinearSpeedup', Nonlinear_Speedup);
% '------------------------------------------------------------------------
% 'Iron Loss Activation settings
StudyProperties.SetValue('IronEddyCurrentDensityLamination', 1);
StudyProperties.SetValue('IronCurrentLossDensityLamination', 1);
StudyProperties.SetValue('IronHysteresisLossDensityLamination', 1);
StudyProperties.SetValue('EddyCurrentDensityLamination', 1);
StudyProperties.SetValue('CurrentLossDensityLamination', 1);
StudyProperties.SetValue('HysteresisLossDensityLamination', 1);

%% %% 'Define Material Characteristics
% '------------------------------------------------------------------------
BLKLABELS.materials=char(...
    mat.Stator.MatName,...
    mat.Rotor.MatName,...
    'N52M',...
    'Copper'...
    );  
% 'Magnet Operating Temperature
if any(unique(Magnetxy(:,9)) ~= 0)
PMTemp = mat.LayerMag.temp.temp(1,1);%20 C
end
% % '------------------------------------------------------------------------
% % 'Coil material
coil_mat=strtrim(BLKLABELS.materials(4,:));
% 'Electric resistivity of the coil
%Allow eddy current from coil : Yes(1) / No(0)
eddycurrent_c=1;    % 0= No eddy current 
                    % 1= Allow eddy current



 
% '------------------------------------------------------------------------
% 'stator & rotor material-customized material based on SYRe library
new_resistivity=1.673e-08;  %ohm.m   
stator_mat1 = BLKLABELS.materials(1,:);
stator_mat = strcat(stator_mat1,'Customized');
%'# New BH-curve data based on SyRe library
property = material_properties_iron(stator_mat1);
BH_curve = property.BH;
BH_curve(:, [1, 2]) = BH_curve(:, [2, 1]);
libPath = 'Custom Materials';   %Path in the [Custom Materials] folder from Material database  
JDesigner.GetMaterialLibrary().CreateCustomMaterial(stator_mat, libPath);
create_new_material=JDesigner.GetMaterialLibrary().GetUserMaterial(stator_mat);
create_new_material.SetValue('Manufacturer', 'Customer');
create_new_material.SetValue('Category', 'Semi-Magnetic');
%Magnetic properties
create_new_material.SetValue('MagneticPropertyType', 'MagneticSteel'); %0 or MagneticSteel / 1 or Magnet / 2 or MagnetizationMaterial 
create_new_material.SetValue('IsotropicType', 'Isotropic');% Material type : 0 or Isotropic / 1 or Anisotropic
create_new_material.SetValue('MagneticSteelPermeabilityType', 'BHCurve');   
% Magnetic steel permeability type: % 0 or LinearConstant
% 1 or LinearThermalDependency% 2 or BHCurve% 3 or BMCurve% 4 or BMurCurve
% 5 or StressDependency% 6 or MagneticSteelUserSubroutine% 7 or TemperatureDependencyBHCurve
create_new_material.GetTable('BhTable').SetName('BH_MyNewMaterial'); % BH curve name
create_new_material.GetTable('BhTable').SetTable([num2cell(BH_curve)]); % set BH curve data
%Electric properties
create_new_material.SetValue('ConductivityType', 1);
create_new_material.SetValue('Resistivity', new_resistivity);
%Loss properties
%IronLoss_type = 1 = IronLossFomula
create_new_material.SetValue('Loss_Type', 1);
create_new_material.SetValue('LossConstantKhX', property.kh);
create_new_material.SetValue('LossConstantKeX', property.ke);
create_new_material.SetValue('LossConstantAlphaX', property.alpha);
create_new_material.SetValue('LossConstantBetaX', property.beta);
create_new_material.SetValue('LossConstantGammaX', 2.0);
create_new_material.SetValue('LossConstantDeltaX', 2.0);
% '------------------------------------------------------------------------
% 'stator material
% stator_mat=BLKLABELS.materials(1,:);
% 'stator lanimation factor in %
stator_lamfactor=99.9;
% '#stator magnetic saturation factor in %
% stator_magSatfactor = 100;
% 'Electric conductivity of the stator
%Use Materials which contain Hysteresis Loops from stator: Yes(1) / No(0)
hysteresisloop_s=0; % 0= not use 
                    % 1=use
%Allow eddy current from stator : Yes(1) / No(0)
eddycurrent_s=0;    % 0= No eddy current 
                    % 1= Allow eddy current 
                    % 2= calculate distribution in steel plate 
% '------------------------------------------------------------------------
% 'Rotor Material
% rotor_mat=BLKLABELS.materials(2,:);
rotor_mat = stator_mat;
% 'Rotor lanimation factor in %
rotor_lamfactor=99.9;
% 'rotor magnetic saturation factor in %
% rotor_magSatfactor=100;
% 'Electric conductivity of the rotor
%Use Materials which contain Hysteresis Loops from rotor: Yes(1) / No(0)
hysteresisloop_r=0; % 0= not use 
                    % 1=use
%Allow eddy current from rotor : Yes(1) / No(0)
eddycurrent_r=0;    % 0= No eddy current 
                    % 1= Allow eddy current 
                    % 2= calculate distribution in steel plate
% '------------------------------------------------------------------------
% 'Magnet material
if any(unique(Magnetxy(:,9)) ~= 0)
magnet_mat=strtrim(BLKLABELS.materials(3,:));
% 'Magnetization Pattern for IPMs
magnet_patternI='Parallel';
% 'Magnetization Pattern for SPMs
magnet_patternS='RadialCircular';
% 'Electric resistivity of the magnet
% 'magnet_conductivity=xxxxxxx;
%Allow eddy current from magnet : Yes(1) / No(0)
eddycurrent_m=1;    % 0= No eddy current 
                    % 1= Allow eddy current


% 'Customized Magnet material
% new_material_name='FeNGen2';
% libPath = 'Custom Materials';   %Path in the [Custom Materials] folder from Material database  
% JDesigner.GetMaterialLibrary().CreateCustomMaterial(new_material_name, libPath);
% create_new_material=JDesigner.GetMaterialLibrary().GetUserMaterial(new_material_name);
% 
% %Magnetic properties
% create_new_material.SetValue('MagneticPropertyType', 1); %0 or MagneticSteel / 1 or Magnet / 2 or MagnetizationMaterial 
% create_new_material.SetValue('MagnetPermeabilityType', 3); 
% create_new_material.GetTable('DemagTemperatureTable').SetName('Demagnetization Table'); 
% curve = [25, 300000, 1.5, 1.1, 70, 180000];  %A/m
% create_new_material.GetTable('DemagTemperatureTable').SetTable([num2cell(curve)]); 
% magnet_mat=strtrim(new_material_name);                    
end
% '#'______________________________________________________________________
% '#'Set up Coil material
SetCoilMaterial_JMAG(study,coil_mat,eddycurrent_c,'_U')
SetCoilMaterial_JMAG(study,coil_mat,eddycurrent_c,'_V')
SetCoilMaterial_JMAG(study,coil_mat,eddycurrent_c,'_W')

% '#'_______________________________________________________________
% '#'Set up Stator core material
SetCoreMaterial_JMAG(study,'Stator Core',stator_mat,1,stator_lamfactor,hysteresisloop_s,eddycurrent_s)
% '#'_______________________________________________________________
% '#'Set up Rotor core material
SetCoreMaterial_JMAG(study,'Rotor Core',rotor_mat,1,rotor_lamfactor,hysteresisloop_r,eddycurrent_r)
% '#'_______________________________________________________________
% '#'Set up Magnet material
% % For grouped magnets
if any(unique(Magnetxy(:,9)) ~= 0)
MagnetCG=geo.BLKLABELS.rotore.xy(geo.BLKLABELS.rotore.xy(:,3)==codMatBar,1:2);%center of gravity(G) of Magnet
MagnetizationDirection=geo.BLKLABELS.rotore.xy(geo.BLKLABELS.rotore.xy(:,3)==codMatBar,6:7);
 %condition for Motors except SPMs
for ID = 1:1:length(MagnetCG)
        Pm = JDesigner.CreatePoint(MagnetCG(ID,1),MagnetCG(ID,2),0);
        part_m = model.GetPartByPosition(Pm);
        part_m.SetName(strcat('Magnet',num2str(ID)));
        study.SetMaterialByName(strcat('Magnet',num2str(ID)), magnet_mat);
        study.GetMaterial(strcat('Magnet',num2str(ID))).SetValue('EddyCurrentCalculation', eddycurrent_m);
        study.GetMaterial(strcat('Magnet',num2str(ID))).SetValue('UserConductivityType', 0);
        study.GetMaterial(strcat('Magnet',num2str(ID))).SetValue('TemperatureType', 1);
        study.GetMaterial(strcat('Magnet',num2str(ID))).SetValue('Temperature', PMTemp);
        study.GetMaterial(strcat('Magnet',num2str(ID))).SetOrientation(1); %' 1=Outward 0=Inward
            if length(AirRotorxy(:,1)) > 8  %condition for Motors except SPMs
        study.GetMaterial(strcat('Magnet',num2str(ID))).SetPattern(magnet_patternI);
        study.GetMaterial(strcat('Magnet',num2str(ID))).SetDirectionXYZ(MagnetizationDirection(ID,1), MagnetizationDirection(ID,2), 0);
            else
        study.GetMaterial(strcat('Magnet',num2str(ID))).SetPattern(magnet_patternS);
        study.GetMaterial(strcat('Magnet',num2str(ID))).SetValue('Poles', geo.p*2);
        study.GetMaterial(strcat('Magnet',num2str(ID))).SetValue('StartAngle', (alpha_r/2)*180/pi); %'variable name
        study.GetMaterial(strcat('Magnet',num2str(ID))).SetDirectionXYZ(1, 0, 0);
        study.GetMaterial(strcat('Magnet',num2str(ID))).SetOriginXYZ(0, 0, 0);
            end
end
createSet_JMAG(model, 'Rotor Magnet','Magnet',MagnetCG);
end
% '------------------------------------------------------------------------

%% %Apply Condictions
% '#'______________________________________________________________________
% '#'Apply the periodic boundary condition
JDesigner.SetCurrentStudy(study_name);
pbcrx=(R1x+R2x)/2.0; %'#rotor reference line
pbcry=(R1y+R2y)/2.0;
pbcsx=(S1x+S2x)/2.0; %'#stator reference line
pbcsy=(S1y+S2y)/2.0;
periodicboudary_name='Periodic';
study.CreateCondition('RotationPeriodicBoundary', 'Periodic Condition');
periodicity_type= 1; %'# 1=antiperiodic
study.GetCondition('Periodic Condition').SetName(periodicboudary_name);
periodic_boundary_condition = study.GetCondition(periodicboudary_name);%study.GetCondition(cond0)
periodic_boundary_condition.SetValue('BoundaryType', periodicity_type);
periodic_boundary_condition.SetValue('Angle', 'model_periodic_angle');
periodic_boundary_condition.SetValue('Axis2D', 0); %'#0=upward, 1=downward
periodic_boundary_condition.ClearParts();
sel = periodic_boundary_condition.GetSelection();
sel.SelectEdgeByPosition(pbcrx, pbcry, 0);
sel.SelectEdgeByPosition(pbcsx, pbcsy, 0);
periodic_boundary_condition.AddSelected(sel);

JDesigner.View().SelectEdge();
JDesigner.View().SelectPart();
JDesigner.SetCurrentStudy(study_name);

% '------------------------------------------------------------------------
 % % Save and Exits JMAG Designer
JDesigner.Save();
JDesigner.Quit();
