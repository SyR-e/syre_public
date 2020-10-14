%%%%%%%%%%%%%
% 2013/07/29 MG Script di costruzione geometria in magnet a partire da un file dxf.
%%%%%%%%%%%%%
%   Preparazione modello MagNet per simulazioni TWM con Arco di rotore
%   interno
%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%

addpath .\Boundaries
addpath .\Coils
addpath .\Draw
addpath .\Interface
addpath .\Meshing
addpath .\Selection
addpath .\Solution
addpath .\Utilities
addpath .\m1

%% Commonly use words
CGD='Call getDocument()';
CGDGV='Call getDocument().getView()';
IMCUS='infoMakeComponentUnionSurfaces';
IMCRV='infoMakeComponentRemoveVertices';
ITIS=' infoToggleInSelection';
IATS='infoAddToSelection';

%% === Apertura di Magnet e creazione modello =============================
load ultimo;
[docName,docPath]=uigetfile([docPath,'\*.mat'],'Select the Machine Mat-file');
save([cd,'\ultimo'],'docPath')
load([docPath,docName]);
dxfName=[docPath,docName(1:end-4),'.dxf'];
% CENTRI=BLKLABELS;

CENTRIstat=geo.BLKLABELS.statore;
CENTRIrot=geo.BLKLABELS.rotore;
CodificaMateriali;
CENTRI.materials=char(Air1, Air2, Stat_Cu,Stat_Fe,Rot_Fe,Stat_Magn,Rot_Magn,Shaft_Steel,Rot_Cu);
Mac=geo;
%% variable not present in geo
Mac.Q=6*Mac.p*Mac.q;
Mac.RtS=geo.r+geo.g;
%% ========================
L_assiale=Mac.l;
t=gcd(Mac.Q,Mac.p); % machine periodicity

h = OpenMagnet(1);  % 1 significa visibile, 0 invisibile

% Nuovo Documento 
h = NewDocumentMagnet(h);

% h = SaveDocumentMagnet(h,docName);

% Importa file dxf e crea geometria 
ImportDXFMagnet(h,dxfName);

pc=360/(Mac.Q)/2; % mezzo passo cava
xr=Mac.RtS-Mac.g;

[traferro,CENTRItraf]=traferroMatr(pc,xr,Mac.g,Mac.Q,Mac.Qs,Mac.ps,Mac.p,Mac.fem.res_traf);
[rig_traf, col_traf] = size(traferro);

% Se geometria MOGA
for ii=1:rig_traf
    if(traferro(ii,col_traf)==0)
        DrawLineMagnet(h,[traferro(ii,1) traferro(ii,2)],[traferro(ii,3) traferro(ii,4)]);
%         keyboard
    else
        DrawArcMagnetperPunti(h,traferro(ii,1),traferro(ii,2),traferro(ii,3),traferro(ii,4),traferro(ii,5),traferro(ii,6));
%         keyboard
    end
end
% % Se geometria MOGA
% DrawLineMagnet(h,[r_rot*cosd(360/Ksym) r_rot*sind(360/Ksym)],[r1*cosd(360/Ksym) r1*sind(360/Ksym)]);
% DrawLineMagnet(h,[r_rot 0],[r1 0]);
% DrawArcMagnet(h,[0 0],r1,0,360/Ksym);
% DrawArcSectorMagnet(h,[0 0],r2,r3,0-th_Stat_Rot,360/Ksym-th_Stat_Rot);
% DrawArcSectorMagnet(h,[0 0],r1,r2,0-th_Stat_Rot,360/Ksym-th_Stat_Rot);
% % DrawArcMagnet(h,[0 0],r3,0-th_Stat_Rot,360/Ksym-th_Stat_Rot);
% % DrawArcMagnet(h,[0 0],r2,0-th_Stat_Rot,360/Ksym-th_Stat_Rot);

%% Crea componenti nel modello 

% Assegna materiali alle diverse parti della macchina

%%%%%%%%%%%%%
%% Statore
%%%%%%%%%%%%%

[var1,var2]=size(CENTRIstat.names.air_slot);
for kk=1:var1
    h = MakeComponentMagnet(h,[CENTRIstat.xy(kk,1), CENTRIstat.xy(kk,2)],cell2mat(CENTRIstat.names.air_slot{kk,1}),L_assiale,CENTRI.materials(2,:),'None',[0, 0, 1],CENTRIstat.xy(kk,4));
    indice_air_slot=kk;
end

[var1,var2]=size(CENTRIstat.names.Cu_slot);
for kk=1:var1
    h = MakeComponentMagnet(h,[CENTRIstat.xy(kk+indice_air_slot,1), CENTRIstat.xy(kk+indice_air_slot,2)],cell2mat(CENTRIstat.names.Cu_slot{kk,1}),L_assiale,CENTRI.materials(3,:),'None',[0, 0, 1],-1);
    indice_Cu_slot=kk+indice_air_slot;
end
[var1,var2]=size(CENTRIstat.names.FeYoke);
for kk=1:var1
    h = MakeComponentMagnet(h,[CENTRIstat.xy(kk+indice_Cu_slot,1), CENTRIstat.xy(kk+indice_Cu_slot,2)],cell2mat(CENTRIstat.names.FeYoke{kk,1}),L_assiale,CENTRI.materials(4,:),'None',[0, 0, 1],CENTRIstat.xy(kk+indice_Cu_slot,4));
    
end

%%%%%%%%%%%
%% Rotore
%%%%%%%%%%%

% Barrier
% Mac.Br=0.4;
% CENTRI.materials(7,:)='NMF3F';
if (Mac.Br==0)
    [BarSize,var2]=size(CENTRIrot.BarName);
    for kk=1:BarSize
        h = MakeComponentMagnet(h,[CENTRIrot.xy(kk,1), CENTRIrot.xy(kk,2)],cell2mat(CENTRIrot.BarName{kk,1}),L_assiale,CENTRI.materials(2,:),'None',[0, 0, 1],CENTRIrot.xy(kk,4));
    end
else
    % Barrier Filled with magnet
    
    [BarSize,var2]=size(CENTRIrot.BarName);
    for kk=1:BarSize
        h = MakeComponentMagnet(h,[CENTRIrot.xy(kk,1), CENTRIrot.xy(kk,2)],cell2mat(CENTRIrot.BarName{kk,1}),L_assiale,CENTRI.materials(7,:),'Uniform',[CENTRIrot.xy(kk,6),CENTRIrot.xy(kk,7), CENTRIrot.xy(kk,8)],CENTRIrot.xy(kk,4));
    end
end

% rotor Iron
 kk=(BarSize+1);
h = MakeComponentMagnet(h,[CENTRIrot.xy(kk,1), CENTRIrot.xy(kk,2)],'rotor',L_assiale,CENTRI.materials(5,:),'None',[0, 0, 1],CENTRIrot.xy(kk,4));

% Shaft
 kk=BarSize+2;
h = MakeComponentMagnet(h,[CENTRIrot.xy(kk,1), CENTRIrot.xy(kk,2)],'shaft',L_assiale,CENTRI.materials(8,:),'None',[0, 0, 1],-1);


%%%%%%%%%%%%
%% Traferro
%%%%%%%%%%%%
for kk=1:size(CENTRItraf.names,1)
h = MakeComponentMagnet(h,[CENTRItraf.xy(kk,1), CENTRItraf.xy(kk,2)],CENTRItraf.names{kk,1},L_assiale,CENTRI.materials(CENTRItraf.xy(kk,3),:),'None',[0, 0, 1],-1);
end


%%
% % Crea le bobine di statore
% stringhe di uso comune

avv=Mac.avv(:,1:Mac.Qs);
  
%COSTRUZIONE COIL U
coil_name = 'U';
coil_number = 1;
MakeSimpleCoilMagnet_UVW;

%COSTRUZIONE COIL V
coil_name = 'V';
coil_number = 2;
MakeSimpleCoilMagnet_UVW;

%COSTRUZIONE COIL W
coil_name = 'W';
coil_number = 3;
MakeSimpleCoilMagnet_UVW

%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%
%% Cancella tutte le linee di costruzione

DeleteAllLinesMagnet(h, Mac.R+50, -Mac.R-50, -Mac.R-50, Mac.R+50);

%%%%%%%%%%%%%%%%%%
% Costruzione delle airbox
%%%%%%%%%%%%%%%%%%
%% 2013/07/29 MG in questa prima bozza si esegue una rivalutazione dei punti principali al traferro
[xArcTraf1,yArcTraf1] = rot_point(xr+1/3*Mac.g,0,-pc*pi/180);
[xArcTraf2,yArcTraf2] = rot_point(xr+1/3*Mac.g,0,(2*Mac.Qs-1)*pc*pi/180);

[THETA1,RHO1] = cart2pol(xArcTraf1,yArcTraf1);
THETA1=THETA1*180/pi;
[THETA2,RHO2] = cart2pol(xArcTraf2,yArcTraf2);
THETA2=THETA2*180/pi;

DrawArcSectorMagnet(h,[0 0],RHO1,Mac.R,THETA1,THETA2);
DrawArcMagnet(h,[0 0],RHO1,THETA2,Mac.ps*180/Mac.p);

DrawLineMagnet(h,[0 0],[RHO1*cos(Mac.ps*pi/Mac.p) RHO1*sin(Mac.ps*pi/Mac.p)]);
DrawLineMagnet(h,[0 0],[RHO1 0]);

h = MakeComponentMagnet(h,[0, 0],'Rotor_Airbox',L_assiale,CENTRI.materials(1,:),'None',[0, 0, 1],-1);
h = MakeComponentMagnet(h,[RHO1*1.001, RHO1*0.0001],'Stator_Airbox',L_assiale,CENTRI.materials(1,:),'None',[0, 0, 1],-1);

% Suddivisione mesh al traferro
clear edgeList;
edgeList{1}{1} = 'air_gap_rotor';
edgeList{1}{2} = 2;
edgeList{2}{1} = 'air_gap_stator1';
edgeList{2}{2} = 2;
edgeList{3}{1} = 'air_gap_stator2';
edgeList{3}{2} = 2;

% (MG) Machine periodicity selection
Q = Mac.Q;                    % number of slots 
t = gcd(round(Mac.Q),Mac.p);  % number of periods

if ((6*t/Q)>1)
    bdryType='Even';
else
    bdryType = 'Odd';
end

if (Mac.ps*180/Mac.p<=90)
    
    edgeList{1}{3} = 3;
    edgeList{2}{3} = 3;
    edgeList{3}{3} = 3;
else
    edgeList{1}{3} = 4;
    edgeList{2}{3} = 4;
    edgeList{3}{3} = 4;
end

DimMesh=round(2*pc*Mac.Qs/Mac.fem.res_traf);
SetUniformMeshOnEdgesMagnet(h,edgeList,DimMesh);


%================ Imposta condizioni al contorno ==========================

% Flux Tangential su superficie esterna di Statore
clear bdryFaces
bdryFaces{1}{1} = 'Stator_Airbox';
if (Mac.ps*180/Mac.p<=90)
    bdryFaces{1}{2} = 5;
else
    bdryFaces{1}{2} = 6;
end
bdryName = 'Flux Tangential';
SetBdryTangentialMagnet(h,bdryName,bdryFaces);

% Odd or even boundary aplication
clear bdryFaces
if (Mac.ps*180/Mac.p<=90)
    bdryFaces{1}{1} = 'Stator_Airbox';
    bdryFaces{1}{2} = 6;
    bdryFaces{2}{1} = 'Rotor_Airbox';
    bdryFaces{2}{2} = 5;
    bdryFaces{3}{1} = 'Rotor_Airbox';
    bdryFaces{3}{2} = 4;
elseif (Mac.ps*180/Mac.p==180)
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
rotAngle = -(Mac.ps*180/Mac.p);
SetBdryRoundPeriodicMagnet(h,bdryName,bdryType,bdryFaces,rotAxis,[0, 0, 0],rotAngle)

% Imposta parte rotante 

[var1 var2]=size(CENTRIrot.BarName);
clear Motion_Components
for kk=1:var1
   Motion_Components{kk}=cell2mat(CENTRIrot.BarName{kk,1});
end
Motion_Components{kk+1}='rotor';
Motion_Components{kk+2}='shaft';
Motion_Components{kk+3}= 'air_gap_rotor';
Motion_Components{kk+4}= 'Rotor_Airbox';


h = MakeMotionComponentMagnet(h,Motion_Components,'Motion#1',6,0);


%% Imposta Parametri per simulazione Transient with Motion

% SetTransientOptions (h,num2str(tStart),num2str(tStep),num2str(tStop),'Yes');
SetTransientOptions (h,0,0.1,6,'Yes',0)

%% ===================== Salva e chiudi MagNet ============================
[h,f] = SaveDocumentMagnet(h,[docPath,docName(1:end-4),'.mn']);
CloseMagnet(h)













