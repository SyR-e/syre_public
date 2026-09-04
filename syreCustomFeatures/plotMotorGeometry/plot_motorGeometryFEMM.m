function [hax,pShape] = plot_motorGeometryFEMM(hax,pathname,filename,fullRotor)

if nargin==0
    [filename,pathname] = uigetfile([cd '\*.fem'],'Select a FEMM model to plot');
    figure();
    figSetting(12,12)
    hax = gca;
end

if nargin()<4
    fullRotor = 1;
end


filename = [filename(1:end-4) '.fem'];

% if exist('openfemm.m','file')
%     flag_xfemm = 0;
% else
%     flag_xfemm = 1;
% end

flag_xfemm = 1;

if ~flag_xfemm
    disp('Importing machine geometry from FEMM...')
else
    disp('Importing machine geometry from XFEMM...')
end

if ~flag_xfemm
    openfemm(1)
    opendocument([pathname filename]);
    mi_createmesh;
    mi_analyze(1);
    closefemm;
    fileans = [pathname filename(1:end-4) '.ans'];
else
    syreDirectory = fileparts(which('GUI_Syre.mlapp'));
    newFile = checkPathSyntax([syreDirectory '\tmp\tmpPlotFEMM.fem']);
    FemmProblem = loadfemmfile([pathname filename]);
    FemmProblem.ProbInfo.SmartMesh = 0;
    FemmProblem = setcircuitcurrent(FemmProblem,'fase1',1);
    writefemmfile(newFile,FemmProblem)
    newFile = fmesher(newFile);
    fileans = fsolver(newFile,0,1);
end

% fileans = [filename(1:end-4) '.ans'];

keyWord = '[NumBlockLabels]';
fid=fopen([fileans],'r');
l=fgetl(fid);
while  ~contains(l,keyWord)
    l=fgetl(fid);
end
labels=cell2mat(textscan(fid,'','collectoutput',1));


keyWord = '[Solution]';
l=fgetl(fid);
while  ~contains(l,keyWord)
    l=fgetl(fid);
end
solution=cell2mat(textscan(fid,'','collectoutput',1));
fid=fclose(fid);

labels = labels(:,7);

numNodes    = solution(1,1);
numElements = solution(numNodes+2 ,1);

nodes    =  solution(2:numNodes+1 ,1:2);
elements =  solution(numNodes+3:numNodes+2+numElements ,1:4);
elements = elements +1;
elements(:,4) = labels(elements(:,4));

eleRotor    = elements(elements(:,4)==22,:);
eleStator   = elements(elements(:,4)==12,:);
eleSlot     = elements(elements(:,4)==1,:);
eleMagnet   = elements(elements(:,4)>199,:);
eleSleeve   = elements(elements(:,4)==199,:);
eleStruct   = [eleRotor ; eleMagnet];

vertRotor   = [nodes(eleRotor(:,1),1:2)  nodes(eleRotor(:,2),1:2)  nodes(eleRotor(:,3),1:2)];
vertStator  = [nodes(eleStator(:,1),1:2)  nodes(eleStator(:,2),1:2)  nodes(eleStator(:,3),1:2)];
vertSlot    = [nodes(eleSlot(:,1),1:2)  nodes(eleSlot(:,2),1:2)  nodes(eleSlot(:,3),1:2)];
vertMagnet  = [nodes(eleMagnet(:,1),1:2)  nodes(eleMagnet(:,2),1:2)  nodes(eleMagnet(:,3),1:2)];
vertSleeve  = [nodes(eleSleeve(:,1),1:2)  nodes(eleSleeve(:,2),1:2)  nodes(eleSleeve(:,3),1:2)];

psRotor  = polyshape;
psMagnet = polyshape;
psStator = polyshape;
psSlot   = polyshape;
psSleeve = polyshape;

for ii=1:length(vertRotor)
    psRotor(ii)     = polyshape([vertRotor(ii,1) vertRotor(ii,3) vertRotor(ii,5)] , [vertRotor(ii,2) vertRotor(ii,4) vertRotor(ii,6)]);
end

for ii=1:length(vertStator)
    psStator(ii)     = polyshape([vertStator(ii,1) vertStator(ii,3) vertStator(ii,5)] , [vertStator(ii,2) vertStator(ii,4) vertStator(ii,6)]);
end

for ii=1:length(vertSlot)
    psSlot(ii)     = polyshape([vertSlot(ii,1) vertSlot(ii,3) vertSlot(ii,5)] , [vertSlot(ii,2) vertSlot(ii,4) vertSlot(ii,6)]);
end

for ii=1:length(vertMagnet)
    psMagnet(ii)     = polyshape([vertMagnet(ii,1) vertMagnet(ii,3) vertMagnet(ii,5)] , [vertMagnet(ii,2) vertMagnet(ii,4) vertMagnet(ii,6)]);
end

for ii=1:length(vertSleeve)
    psSleeve(ii)     = polyshape([vertSleeve(ii,1) vertSleeve(ii,3) vertSleeve(ii,5)] , [vertSleeve(ii,2) vertSleeve(ii,4) vertSleeve(ii,6)]);
end

psRotor     = union(psRotor);
psStator    = union(psStator);
psSlot      = union(psSlot);
psMagnet    = union(psMagnet);
psSleeve    = union(psSleeve);

disp(['Geometry imported from FEMM!'])


colors.lamination = [0.5 0.5 0.5];
colors.conductors = [1.0 0.5 0.0];
colors.magnets    = [0.0 0.0 1.0];
colors.sleeve     = [0.0 0.6 0.0];

plot(hax,psMagnet,'FaceColor',colors.magnets,'EdgeColor','none','FaceAlpha',0.8,'DisplayName','PM')
plot(hax,psRotor,'FaceColor',colors.lamination,'EdgeColor','none','FaceAlpha',0.8,'DisplayName','Rotor core')
plot(hax,psStator,'FaceColor',colors.lamination,'EdgeColor','none','FaceAlpha',0.8,'DisplayName','Stator core')
plot(hax,psSlot,'FaceColor',colors.conductors,'EdgeColor','none','FaceAlpha',0.8,'DisplayName','Stator windings')
plot(hax,psSleeve,'FaceColor',colors.sleeve,'EdgeColor','none','FaceAlpha',0.8,'DisplayName','Rotor sleeve')

pShape.rotor  = psRotor;
pShape.stator = psStator;
pShape.magnet = psMagnet;
pShape.slot   = psSlot;
pShape.sleeve = psSleeve;

Rmax = max(abs(pShape.stator.Vertices(:,1)+j*pShape.stator.Vertices(:,2)));
tol = Rmax*0.05;

set(hax,...
    'XLim',[-tol Rmax+tol],...
    'YLim',[-tol Rmax+tol],...
    'DataAspectRatio',[1 1 1]);

if fullRotor
    angleMax = max(angle(pShape.rotor.Vertices(:,1)+j*pShape.rotor.Vertices(:,2)))*180/pi;
    numRep = round(360/angleMax)-1;
    for ii=1:numRep
        angleRot = angleMax*ii*pi/180;
        psTmp.stator = rotate_pShape(pShape.stator,angleRot);
        psTmp.rotor  = rotate_pShape(pShape.rotor,angleRot);
        psTmp.slot   = rotate_pShape(pShape.slot,angleRot);
        psTmp.magnet = rotate_pShape(pShape.magnet,angleRot);
        psTmp.sleeve = rotate_pShape(pShape.sleeve,angleRot);
        plot(hax,psTmp.magnet,'FaceColor',colors.magnets,'EdgeColor','none','FaceAlpha',0.8)
        plot(hax,psTmp.rotor,'FaceColor',colors.lamination,'EdgeColor','none','FaceAlpha',0.8)
        plot(hax,psTmp.stator,'FaceColor',colors.lamination,'EdgeColor','none','FaceAlpha',0.8)
        plot(hax,psTmp.slot,'FaceColor',colors.conductors,'EdgeColor','none','FaceAlpha',0.8)
        plot(hax,psTmp.sleeve,'FaceColor',colors.sleeve,'EdgeColor','none','FaceAlpha',0.8)
    end
    set(hax,...
    'XLim',(Rmax+tol)*[-1 1],...
    'YLim',(Rmax+tol)*[-1 1]);
end




