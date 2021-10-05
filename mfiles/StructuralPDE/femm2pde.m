% Copyright 2021
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

function [structModel] = femm2pde(geo,mat,simSetup)

meshSize  = simSetup.meshSize;
shaftBC   = simSetup.shaftBC;
evalSpeed = simSetup.evalSpeed;
filename  = simSetup.filename;
pathname  = simSetup.pathname;


% load data from input structures
% rotor = geo.rotor;
% nmax  = geo.nmax;
% meshSize = 'coarse';

switch meshSize
    case 'coarse'
        Hmax  = 5*geo.pont0/1e3;%2/1e3;
        Hmin  = geo.pont0/1e3/2;
        Hgrad = 1.1;
        warning('Coarse mesh selected!!!')
    case 'fine'
        Hmax = geo.pont0/1e3/2;
        Hmin = geo.pont0/1e3/10;
        Hgrad = 1.1;
    otherwise
        Hmax = geo.pont0/1e3/2;
        Hmin = geo.pont0/1e3/10;
        Hgrad = 1.1;
end

% additional data
if ~isfield(mat.Rotor,'E')
    E = [];
else
    E = mat.Rotor.E*1e9;
end

if isempty(E)
    disp('Young module not set for the selected iron. Set to 210 GPa')
    E = 210e9;
end

nu = 0.3;

structModel = createpde(2);

filename  = strrep(filename,'.mat','.fem');
fileans   = strrep(filename,'.fem','.ans');

openfemm(1);
opendocument([pathname filename]);
if isfile([pathname fileans])
     mi_loadsolution;
else
     mi_createmesh;
     mi_analyze(1);
     mi_loadsolution;
end

numNodes    = mo_numnodes;
numElements = mo_numelements;

nodes       = zeros(2,numNodes);
elements    = zeros(3,numElements);
eleGroup    = zeros(1,numElements);

for ii=1:numNodes
tmp=mo_getnode(ii);
nodes(1,ii)=tmp(1)/1000;
nodes(2,ii)=tmp(2)/1000;
end


k=1; %zz=1; jj=1; 
psMagnet = polyshape;

for ii=1:numElements
    tmp=mo_getelement(ii);
    elements(1,ii)  = tmp(1);
    elements(2,ii)  = tmp(2);
    elements(3,ii)  = tmp(3);
%     elements(4,ii)  = tmp(7);
    eleGroup(ii)    = tmp(7);
    if (eleGroup(ii)~=22 && eleGroup(ii)<200)
        filt(k)=ii;
        k=k+1;
    end
    
%     if (eleGroup(ii)~=22)
%         filtRotor(jj)=ii;
%         jj=jj+1;
%     end
%     
%     if (eleGroup(ii)<200)
%         filtPM(zz)=ii;
%         zz=zz+1;
%     end
    
   if (eleGroup(ii)>199)
      index         = elements(:,ii);
      vertex        = nodes(:,index);
      vertex        = vertex';
      psMagnetTemp  = polyshape(vertex);
      psMagnet      = union(psMagnet,psMagnetTemp);
    end
end
% eleRotor = elements;
% elePM    = elements;
% eleRotor(:,filtRotor) = [];
% elePM(:,filtPM)    = [];

elements(:,filt) = [];

closefemm

geometryFromMesh(structModel,nodes,elements);
generateMesh(structModel,'Hmax',Hmax,'Hmin',Hmin,'Hgrad',Hgrad,'GeometricOrder','quadratic');

% pmModel     = createpde();
% geometryFromMesh(pmModel,nodes,elePM);
% pmGeo       = pmModel.Geometry;
% vertices    = pmGeo.vertexCoordinates(1:pmGeo.NumVertices);

% jj=1;
% psMagnet = polyshape;
% for ii=1:(length(vertices)/4)
%       vertex        = vertices(jj:(jj+3),:);
%       vertex(:,3)=[];
%       kconv         = convhull(vertex);
%       pmVertex      = [vertex(kconv,1) vertex(kconv,2)];
%       psMagnetTemp  = polyshape(pmVertex);
%       psMagnet      = union(psMagnet,psMagnetTemp);
%       jj = jj+4;
% end

% for ii=1:(length(elePM))
%       index         = elePM(:,ii);
%       vertex        = nodes(:,index);
%       vertex        = vertex';
%       kconv         = convhull(vertex);
%       pmVertex      = [vertex(kconv,1) vertex(kconv,2)];
%       psMagnetTemp  = polyshape(pmVertex);
%       psMagnet      = union(psMagnet,psMagnetTemp);
% end

% psRotor = polyshape;
% for ii=1:(length(eleRotor))
%       index         = eleRotor(:,ii);
%       vertex        = nodes(:,index);
%       vertex        = vertex'*1000;
%       kconv         = convhull(vertex);
%       rotVertex     = [vertex(kconv,1) vertex(kconv,2)];
%       psRotorTemp   = polyshape(rotVertex);
%       psRotor       = union(psRotor,psRotorTemp);
% end

% figure
% figSetting
% plot(psMagnet)
% axis equal
% title('Magnets Area')

% Boundary condition: shaft is fixed

hfig = figure();
figSetting();
hax = axes('OuterPosition',[0 0 1 1]);
set(hax,'DataAspectRatio',[1 1 1]);
pdegplot(structModel,'EdgeLabels','on');
title('Mesh geometry')

hchild = get(hax,'Children');
xy = hchild(1).VertexData(1:2,:);
r  = abs(xy(1,:)+j*xy(2,:));
th = angle(xy(1,:)+j*xy(2,:));
index = 1:1:length(r);
indexShaft = index(r==min(r));
indexSpider1 = index(th==min(th));
indexSpider2 = index(th==max(th));

close(hfig);

if geo.ps==2*geo.p
    if shaftBC
        applyBoundaryCondition(structModel,...
            'dirichlet',...
            'Edge',indexShaft,...
            'u',[0 0]);
    end
else
    %'Edge',[indexShaft indexSpider1 indexSpider2],...
    if shaftBC
        applyBoundaryCondition(structModel,...
            'dirichlet',...
            'Edge',[indexShaft],...
            'u',[0,0]);
    end
%     applyBoundaryCondition(structModel,...
%         'dirichlet',...
%         'Edge',indexSpider1,...
%         'u',@(location,state)slidingBC(location,state,[cos(min(th));sin(min(th))]));
%     applyBoundaryCondition(structModel,...
%         'dirichlet',...
%         'Edge',indexSpider2,...
%         'u',@(location,state)slidingBC(location,state,[cos(max(th));sin(max(th))]));
    applyBoundaryCondition(structModel,...
        'dirichlet',...
        'Edge',[indexSpider1 indexSpider2],...
        'u',[0,0]);
end
% Create c matrix for PDE (reference from pdeModeler for plane stress)
cMat = [
    E/(1-nu^2)
    0
    E/(2*(1+nu))
    0
    E/(2*(1+nu))
    E*nu/(1-nu^2)
    0
    E/(2*(1+nu))
    0
    E/(1-nu^2)];

% Preparation of the data for centrifugal stress computation

meshData = structModel.Mesh;
xyNodes  = meshData.Nodes(1,:)+j*meshData.Nodes(2,:);
eleData  = [xyNodes(meshData.Elements(1,:));xyNodes(meshData.Elements(2,:));xyNodes(meshData.Elements(3,:))];
% eleG     = mean(eleData,1);



dataForCF.kgm3_Fe  = mat.Rotor.kgm3;
dataForCF.kgm3_PM  = mat.LayerMag.kgm3;
dataForCF.w        = evalSpeed*pi/30;
dataForCF.psMagnet = psMagnet;
% dataForCF.eleData = eleData;


specifyCoefficients(structModel,...
    'm',0,...
    'd',0,...
    'c',cMat,...
    'a',0,...
    'f',@(x,y)centrifugalForce(x,y,dataForCF));


% keyboard
