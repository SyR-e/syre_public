
% Copyright 2018
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

function [geo,mat] = draw_motor_in_MN(geo,mat,pathname,filename,h)

eval_type = 'singt';    % default
mesh_res = dimMesh(geo,eval_type);

% BLKLABELS.materials=char('AIR', 'AIR', 'Copper: 5.77e7 Siemens/meter','M600-50A','M600-50A','PM04: Br 0.4 mur 1.0','CR10: Cold rolled 1010 steel','Copper: 5.77e7 Siemens/meter');  %scelgo io i nomi dei materiali
BLKLABELS.materials=char('AIR', 'AIR', 'Copper: 5.77e7 Siemens/meter','M250-35A','M250-35A','F42SH_HREfree','CR10: Cold rolled 1010 steel','Copper: 5.77e7 Siemens/meter');  %scelgo io i nomi dei materiali

BLKLABELS.materials=char(...
    'AIR',...
    'AIR',...
    'Copper: 5.77e7 Siemens/meter',...
    mat.Stator.MatName,...
    mat.Rotor.MatName,...
    mat.LayerMag.MatName,...
    'CR10: Cold rolled 1010 steel',...
    'Copper: 5.77e7 Siemens/meter'...
    );



BLKLABELS.rotore = geo.BLKLABELS.rotore;
BLKLABELS.statore= geo.BLKLABELS.statore;

if ~isfile([pathname filename(1:end-4),'.dxf'])
    button='Yes';
else
    button = questdlg('DXF existing. Replace it?','SELECT','Yes','No','Yes');
end

if strcmp(button,'Yes')
    syreToDxf(geo.stator, geo.rotor,pathname, filename);
end

% inizio funzione vecchia
CGD='Call getDocument()';
CGDGV='Call getDocument().getView()';
IMCUS='infoMakeComponentUnionSurfaces';
IMCRV='infoMakeComponentRemoveVertices';
ITIS=' infoToggleInSelection';
IATS='infoAddToSelection';

Br = mat.LayerMag.Br(1);
Q = 6*geo.q*geo.p*geo.win.n3phase;
pc = 360/(Q)/2;         % mezzo passo cava
l = geo.l;

h = NewDocumentMagnet(h);

% Importa file dxf e crea geometria
mh = h.magnetHandler;
invoke(mh,'processCommand',['CALL getDocument().getView().importDXF("',[pathname filename(1:end-4),'.dxf'],'")']);

%% prova confronto materiali
% mh = h.magnetHandler;
% Command = ['CALL getDocument().getModelMaterialDatabase().getMaterials(materials)']
% invoke(mh, 'processCommand', Command);

% MN6 = actxserver('Magnet.application');
% Command = ['CALL getSystemMaterialDatabase().getMaterials(BMN-38,EH)']
% invoke(MN6, 'processCommand', Command);

% [traferro,BLKLABELStraf]=traferroMatr(geo.r,geo.g,geo.Qs,geo.ps,geo.p,fem.res_traf);
[traferro,BLKLABELStraf] = traferroMatr(geo,mesh_res);
[rig_traf, col_traf] = size(traferro);
% if strcmp(geo.RotType,'SPM')
%     % SF (30/08/2018): for the SPM machines the airgap cover also the air
%     % space on the PMs side, so the lines must be extended
%     tmp=traferro(5,:);
%     tmpAng=atan2(tmp(2),tmp(1));
%     tmpRad=geo.r-geo.lm;
%     tmp(1)=tmpRad*cos(tmpAng);
%     tmp(2)=tmpRad*sin(tmpAng);
%     traferro(5,:)=tmp;
%     tmp=traferro(6,:);
%     tmpAng=atan2(tmp(2),tmp(1));
%     tmpRad=geo.r-geo.lm;
%     tmp(1)=tmpRad*cos(tmpAng);
%     tmp(2)=tmpRad*sin(tmpAng);
%     traferro(6,:)=tmp;
% end

for ii=1:rig_traf
    if(traferro(ii,col_traf)==0)
        DrawLineMagnet(h,[traferro(ii,1) traferro(ii,2)],[traferro(ii,3) traferro(ii,4)]);
    else
        DrawArcMagnetperPunti(h,traferro(ii,1),traferro(ii,2),traferro(ii,3),traferro(ii,4),traferro(ii,5),traferro(ii,6));
    end
end

% Assegna materiali alle diverse parti della macchina
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

materialCodes;  % load the material codes for BLKLABELS decodification
% Stator
indexCu  = 1;
indexAir = 1;
indexFe  = 1;
xyStator = BLKLABELS.statore.xy;
namesStator = BLKLABELS.statore.names;

for kk=1:length(BLKLABELS.statore.xy(:,1))
    switch BLKLABELS.statore.xy(kk,3)
        case codMatCu       % slot conductor
            eleName = cell2mat(namesStator.Cu_slot{indexCu});
            indexCu = indexCu+1;
            h = MakeComponentMagnet(h,[xyStator(kk,1) xyStator(kk,2)],eleName,l,BLKLABELS.materials(xyStator(kk,3),:),'None',[0, 0, 1],xyStator(kk,4));
        case codMatAirSta   % slot air
            eleName = cell2mat(namesStator.air_slot{indexAir});
            indexAir=indexAir+1;
            h = MakeComponentMagnet(h,[xyStator(kk,1) xyStator(kk,2)],eleName,l,BLKLABELS.materials(xyStator(kk,3),:),'None',[0, 0, 1],-1);
        case codMatFeSta    % stator iron
            eleName = cell2mat(namesStator.FeYoke{indexFe});
            indexFe = indexFe+1;
            h = MakeComponentMagnet(h,[xyStator(kk,1) xyStator(kk,2)],eleName,l,BLKLABELS.materials(xyStator(kk,3),:),'None',[0, 0, 1],xyStator(kk,4));
    end
end


% Rotore
xyRotor = BLKLABELS.rotore.xy;
%namesRotor = BLKLABELS.names;
indexAir = 1;
indexPM  = 1;

for kk=1:length(xyRotor(:,1))
    switch xyRotor(kk,3)
        case codMatAirRot % rotor air
            eleName = BLKLABELS.rotore.BarName{kk};
            indexAir = indexAir+1;
            h = MakeComponentMagnet(h,[xyRotor(kk,1), xyRotor(kk,2)],eleName,l,BLKLABELS.materials(xyRotor(kk,3),:),'Uniform',[xyRotor(kk,6),xyRotor(kk,7), xyRotor(kk,8)],xyRotor(kk,4));
        case codMatFeRot % rotor iron
            h = MakeComponentMagnet(h,[xyRotor(kk,1), xyRotor(kk,2)],'rotor',l,BLKLABELS.materials(xyRotor(kk,3),:),'None',[0, 0, 1],xyRotor(kk,4));
        case codMatBar   % PM
            eleName = BLKLABELS.rotore.BarName{kk};
            indexPM = indexPM+1;
            if Br==0 % || strcmp(geo.RotType,'SPM')
                h = MakeComponentMagnet(h,[xyRotor(kk,1), xyRotor(kk,2)],eleName,l,BLKLABELS.materials(2,:),'None',[0, 0, 1],xyRotor(kk,4));
            else
                %MatName=BLKLABELS.materials(6,:);
                h = MakeComponentMagnet(h,[xyRotor(kk,1), xyRotor(kk,2)],eleName,l,BLKLABELS.materials(xyRotor(kk,3),:),'Uniform',[xyRotor(kk,6),xyRotor(kk,7), xyRotor(kk,8)],xyRotor(kk,4));
            end
        case codMatShaft % shaft
            h = MakeComponentMagnet(h,[xyRotor(kk,1), xyRotor(kk,2)],'shaft',l,BLKLABELS.materials(xyRotor(kk,3),:),'None',[0, 0, 1],-1);
    end
end
% Traferro
for kk=1:size(BLKLABELStraf.names,1)
    h = MakeComponentMagnet(h,[BLKLABELStraf.xy(kk,1), BLKLABELStraf.xy(kk,2)],BLKLABELStraf.names{kk,1},l,BLKLABELS.materials(BLKLABELStraf.xy(kk,3),:),'None',[0, 0, 1],-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Crea le bobine di statore
n3ph = geo.win.n3phase;
avv = geo.win.avv;


for ii=1:n3ph
    % COIL U
    coil_name = ['U' int2str(ii)];
    coil_number = 1+3*(ii-1);
    MakeSimpleCoilMagnet_UVW;
    
    %COSTRUZIONE COIL V
    coil_name = ['V' int2str(ii)];
    coil_number = 2+3*(ii-1);
    MakeSimpleCoilMagnet_UVW;
    
    %COSTRUZIONE COIL W
    coil_name = ['W' int2str(ii)];
    coil_number = 3+3*(ii-1);
    MakeSimpleCoilMagnet_UVW
    
end

% Cancella tutte le linee di costruzione
DeleteAllLinesMagnet(h, geo.R+50, -geo.R-50, -geo.R-50, geo.R+50);

%%%%%%%%%%%%%%%%%%
% Airbox
%%%%%%%%%%%%%%%%%%
[xArcTraf1,yArcTraf1] = rot_point(geo.r+1/3*geo.g,0,0*pi/180);
[xArcTraf2,yArcTraf2] = rot_point(geo.r+1/3*geo.g,0,(2*geo.Qs)*pc*pi/180);
% [xArcTraf1,yArcTraf1] = rot_point(geo.r+1/3*geo.g,0,-pc*pi/180);
% [xArcTraf2,yArcTraf2] = rot_point(geo.r+1/3*geo.g,0,(2*geo.Qs-1)*pc*pi/180);

[THETA1,RHO1] = cart2pol(xArcTraf1,yArcTraf1);
THETA1=THETA1*180/pi;
[THETA2,RHO2] = cart2pol(xArcTraf2,yArcTraf2);  %RHO2 non usato
THETA2=THETA2*180/pi;

DrawArcSectorMagnet(h,[0 0],RHO1,geo.R,THETA1,THETA2);
% DrawArcMagnet(h,[0 0],RHO1,THETA2,geo.ps*180/geo.p);
DrawArcMagnet(h,[0 0],RHO1,THETA1,THETA2);

DrawLineMagnet(h,[0 0],[RHO1*cos(geo.ps*pi/geo.p) RHO1*sin(geo.ps*pi/geo.p)]);
DrawLineMagnet(h,[0 0],[RHO1 0]);

h = MakeComponentMagnet(h,[0, 0],'Rotor_Airbox',l,BLKLABELS.materials(1,:),'None',[0, 0, 1],-1);
h = MakeComponentMagnet(h,[RHO1*1.001, RHO1*0.0001],'Stator_Airbox',l,BLKLABELS.materials(1,:),'None',[0, 0, 1],-1);

% Suddivisione mesh al traferro
clear edgeList;
edgeList{1}{1} = 'air_gap_rotor';
edgeList{1}{2} = 2;
edgeList{2}{1} = 'air_gap_stator1';
edgeList{2}{2} = 2;
edgeList{3}{1} = 'air_gap_stator2';
edgeList{3}{2} = 2;

% (MG) Machine periodicity selection
% t = gcd(round(Q),geo.p);  % number of periods
if rem(geo.ps,2)==0
    bdryType='Even';
else
    bdryType='Odd';
end

if (geo.ps*180/geo.p<=90)
    edgeList{1}{3} = 3;
    edgeList{2}{3} = 3;
    edgeList{3}{3} = 3;
else
    edgeList{1}{3} = 4;
    edgeList{2}{3} = 4;
    edgeList{3}{3} = 4;
end

DimMesh=round(2*pc*geo.Qs/mesh_res.res_traf);
SetUniformMeshOnEdgesMagnet(h,edgeList,DimMesh);

%================ Imposta condizioni al contorno ==========================

% Flux Tangential su superficie esterna di Statore
clear bdryFaces
bdryFaces{1}{1} = 'Stator_Airbox';
if (geo.ps*180/geo.p<=90)
    bdryFaces{1}{2} = 5;
else
    bdryFaces{1}{2} = 6;
end
bdryName = 'Flux Tangential';
SetBdryTangentialMagnet(h,bdryName,bdryFaces);

% Odd or even boundary aplication
clear bdryFaces
if (geo.ps*180/geo.p<=90)
    bdryFaces{1}{1} = 'Stator_Airbox';
    bdryFaces{1}{2} = 6;
    bdryFaces{2}{1} = 'Rotor_Airbox';
    bdryFaces{2}{2} = 5;
    bdryFaces{3}{1} = 'Rotor_Airbox';
    bdryFaces{3}{2} = 4;
elseif (geo.ps*180/geo.p==180)
    bdryFaces{1}{1} = 'Stator_Airbox';
    bdryFaces{1}{2} = 3;
    bdryFaces{2}{1} = 'Rotor_Airbox';
    bdryFaces{2}{2} = 4;
    bdryFaces{3}{1} = 'Rotor_Airbox';
    bdryFaces{3}{2} = 3;
else
    bdryFaces{1}{1} = 'Stator_Airbox';
    bdryFaces{1}{2} = 3;
    bdryFaces{2}{1} = 'Rotor_Airbox';
    bdryFaces{2}{2} = 5;
    bdryFaces{3}{1} = 'Rotor_Airbox';
    bdryFaces{3}{2} = 3;
end

bdryName = 'Periodic';
rotAxis = [0, 0, 1];
rotAngle = -(geo.ps*180/geo.p);
SetBdryRoundPeriodicMagnet(h,bdryName,bdryType,bdryFaces,rotAxis,[0, 0, 0],rotAngle)

% Imposta parte rotante

[var1 var2]=size(BLKLABELS.rotore.BarName);     %var2
clear Motion_Components
for kk=1:var1
    Motion_Components{kk}=(BLKLABELS.rotore.BarName{kk,1});
end
if strcmp(geo.RotType,'SPM')
    %     for jj = kk+1:kk+geo.ps*2
    %         Motion_Components{jj}=['Rotor_Air_Zone_',num2str(jj-kk)];
    %     end
    
    jj = kk;
    Motion_Components{jj+1}='rotor';
    Motion_Components{jj+2}='shaft';
    Motion_Components{jj+3}= 'air_gap_rotor';
    Motion_Components{jj+4}= 'Rotor_Airbox';
else
    Motion_Components{kk+1}='rotor';
    Motion_Components{kk+2}='shaft';
    Motion_Components{kk+3}= 'air_gap_rotor';
    Motion_Components{kk+4}= 'Rotor_Airbox';
end
h = MakeMotionComponentMagnet(h,Motion_Components,'Motion#1',6,0);

%% Imposta Parametri per simulazione Transient with Motion

% SetTransientOptions (h,num2str(tStart),num2str(tStep),num2str(tStop),'Yes');
SetTransientOptions (h,0,0.1,6,'Yes',0)

%% ===================== Salva e chiudi MagNet ============================
% pathname_DXF=[pathname,filename(1:end-4),'\DXF\'];
% pathname_DXF=[pathname,filename(1:end-4)];

% [h,f] = SaveDocumentMagnet(h,[pathname,filename(1:end-4),'.mn']);
% CloseMagnet(h)


