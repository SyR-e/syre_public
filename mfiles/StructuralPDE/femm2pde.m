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

function [structModel,data4GeoMat] = femm2pde(geo,mat,simSetup)

meshSize  = simSetup.meshSize;
shaftBC   = simSetup.shaftBC; 
flagFull  = simSetup.flagFull;
evalSpeed = simSetup.evalSpeed;
filename  = simSetup.filename;
pathname  = simSetup.pathname;


% load data from input structures
% rotor = geo.rotor;
% nmax  = geo.nmax;
% meshSize = 'coarse';

switch meshSize
    case 'coarse'
        Hmax  = geo.pont0/1e3*4;
        Hmin  = geo.pont0/1e3*0.8;
        Hgrad = 2;
        Hedge = Hmin;
        warning('Coarse mesh selected!!!')
    case 'fine'
        Hmax  = geo.pont0/1e3*2;
        Hmin  = geo.pont0/1e3*0.5;
        Hgrad = 2;
        Hedge = geo.pont0/1e3*0.5;
        %         Hmax  = geo.pont0/1e3*2;
        %         Hmin  = geo.pont0/1e3*1;
        %         Hgrad = 1;
        %         Hedge = geo.pont0/1e3*1.5;
    otherwise
        Hmax  = geo.pont0/1e3/2;
        Hmin  = geo.pont0/1e3/5;
        Hgrad = 2;
        Hedge = Hmin;
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
% E  = 200e9;
% warning('Young module imposed at 200 GPa')

structModel = createpde(2);

filename  = strrep(filename,'.mat','.fem');
% fileans   = strrep(filename,'.fem','.ans');

% syreDirectory = fileparts(which('GUI_Syre.mlapp'));

[~,tmpath]=createTempDir();
copyfile([pathname filename],[tmpath filename])

newFile = [tmpath filename];

openfemm(1);
opendocument(newFile);
% reduce mesh size (if not custom)
xy = geo.BLKLABELS.rotore.xy;
for ii=1:size(xy,1)
    if xy(ii,3)==1 % Air
        if ~geo.custom
            mi_selectlabel(xy(ii,1),xy(ii,2));
            mi_setblockprop('Air',0,0,'None',0,2,0);
            mi_clearselected;
        end
    elseif xy(ii,3)==7 % Shaft
        if simSetup.meshShaft
            mi_selectlabel(xy(ii,1),xy(ii,2));
            mi_setblockprop('Air',1,0,'None',0,400,0);
            mi_clearselected;
        end
    elseif xy(ii,3)==6 % PM
        if ~geo.custom
            mi_selectlabel(xy(ii,1),xy(ii,2));
            magdir=atan2(xy(ii,7),xy(ii,6))*180/pi;
            mi_setblockprop(mat.LayerMag.MatName,1,0,'None',magdir,201,0);
            mi_clearselected;
        end
    elseif xy(ii,3)==5 % rotor iron
        if ~geo.custom
            mi_selectlabel(xy(ii,1),xy(ii,2));
            mi_setblockprop(mat.Rotor.MatName,1,0,'None',0,22,0);
            mi_clearselected;
        end
    elseif xy(ii,3)==8 % rotor bar
        if ~geo.custom
            mi_selectlabel(xy(ii,1),xy(ii,2));
            mi_setblockprop(mat.BarCond.MatName,1,0,'None',0,201,0);
            mi_clearselected;
        end
    elseif xy(ii,3)==9 % sleeve
        if ~geo.custom
            mi_selectlabel(xy(ii,1),xy(ii,2));
            mi_setblockprop(mat.Sleeve.MatName,1,0,'None',0,199,0);
            mi_clearselected;
        end
    end
end

% % add shaft ring spring
if simSetup.meshShaft
    if shaftBC==2
        ArIn = geo.Ar/10;
        if geo.ps==2*geo.p
            mi_addnode(+ArIn,0);
            mi_addnode(-ArIn,0);
            mi_addarc(+ArIn,0,-ArIn,0,180,2);
            mi_addarc(-ArIn,0,+ArIn,0,180,2);
            mi_addblocklabel(ArIn/2,0);
            mi_selectlabel(ArIn/2,0);
            mi_setblockprop('Air',1,0,'None',0,2,0);
            mi_clearselected();
        else
            x1 = ArIn;
            y1 = 0;
            [x2,y2] = rot_point(x1,y1,pi/geo.p*geo.ps);
            mi_addnode(x1,y1);
            mi_addnode(x2,y2);
            mi_addarc(x1,y1,x2,y2,180/geo.p*geo.ps,2);
            [x0,y0] = rot_point(x1/2,y1/2,pi/geo.p*geo.ps/2);
            mi_addblocklabel(x0,y0);
            mi_selectlabel(x0,y0);
            mi_setblockprop('Air',1,0,'None',0,2,0);
            mi_clearselected();
            mi_selectsegment(x1/2,y1/2);
            mi_selectsegment(x2/2,y2/2);
            mi_setsegmentprop('None',1,0,0,2);
            mi_clearselected();
        end
    end
end

mi_createmesh;
mi_analyze(1);
mi_loadsolution;

numNodes    = mo_numnodes;
numElements = mo_numelements;

nodes       = zeros(2,numNodes);
elements    = zeros(3,numElements);
eleGroup    = zeros(1,numElements);
elementID   = ones(1,numElements); % 1-->Fe, 2-->PM/Al, 3-->sleeve

for ii=1:numNodes
    tmp=mo_getnode(ii);
    nodes(1,ii)=tmp(1)/1000;
    nodes(2,ii)=tmp(2)/1000;
end


k=1; %zz=1; jj=1; 
psMagnet = polyshape;
psSleeve = polyshape;
% if simSetup.meshShaft
    psShaft  = polyshape;
% end

for ii=1:numElements
    tmp=mo_getelement(ii);
    elements(1,ii)  = tmp(1);
    elements(2,ii)  = tmp(2);
    elements(3,ii)  = tmp(3);
    eleGroup(ii)    = tmp(7);
    if (eleGroup(ii)~=22 && eleGroup(ii)~=199 && eleGroup(ii)<200)
        filt(k)=ii;
        k=k+1;
    end

    if eleGroup(ii)==199
        index    = elements(:,ii);
        vertex   = nodes(:,index);
        vertex   = vertex';
        psTemp   = polyshape(vertex);
        psSleeve = union(psSleeve,psTemp);
        elementID(ii) = 3;
    elseif (eleGroup(ii)==400)
%         if simSetup.meshShaft
            index = elements(:,ii);
            vertex   = nodes(:,index);
            vertex   = vertex';
            psTemp   = polyshape(vertex);
            psShaft = union(psShaft,psTemp);
            elementID(ii) = 4;
%         end
    elseif (eleGroup(ii)>199)
        index    = elements(:,ii);
        vertex   = nodes(:,index);
        vertex   = vertex';
        psTemp   = polyshape(vertex);
        psMagnet = union(psMagnet,psTemp);
        elementID(ii) = 2;
    end
end

if flagFull
    nRep = 2*geo.p/geo.ps-1;
    base.nodes    = nodes;
    base.elements = elements;
    base.group    = eleGroup;
    for ii=1:nRep
        [xTmp,yTmp] = rot_point(base.nodes(1,:),base.nodes(2,:),2*pi*2*geo.p/geo.ps*ii);
        nodes = [nodes, [xTmp;yTmp]];
        elements = [elements, elements+length(base.nodes)];
        eleGroup = [eleGroup, base.group];
    end
    elementID = ones(1,numElements*(nRep+1));
    elementID(eleGroup==199)=3; % sleeve
    eleGroup(eleGroup==199)=0;
    elementID(eleGroup>199)=2; % PM
end


% eleRotor = elements;
% elePM    = elements;
% eleRotor(:,filtRotor) = [];
% elePM(:,filtPM)    = [];

elements(:,filt) = [];
elementID(:,filt) = [];

closefemm
delete(newFile)

% geometryFromMesh(structModel,nodes,elements,elementID);
geometryFromMesh(structModel,nodes,elements);
generateMesh(structModel,'Hmax',Hmax,'Hmin',Hmin,'Hgrad',Hgrad,'GeometricOrder','quadratic','Hedge',{1:1:structModel.Geometry.NumEdges, Hedge});

% Boundary condition

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
% indexSpider1 = index(th==min(th));
% indexSpider2 = index(th==max(th));
indexAirgap  = index(r==max(r));

thMin = min(th);
thMax = max(th);
thTol = 1*pi/180;
indexSpider1 = index(th<thMin+thTol);
indexSpider2 = index(th>thMax-thTol);


close(hfig);

if geo.ps==2*geo.p
    if shaftBC>0
        applyBoundaryCondition(structModel,...
            'dirichlet',...
            'Edge',indexShaft,...
            'u',[0 0]);
    end
else
    %'Edge',[indexShaft indexSpider1 indexSpider2],...
    if shaftBC>0
        applyBoundaryCondition(structModel,...
            'dirichlet',...
            'Edge',[indexShaft],...
            'u',[0,0]);
    end
%     applyBoundaryCondition(structModel,...
%         'mixed',...
%         'Edge',indexSpider1,...
%         'u',0,...
%         'EquationIndex',1);
%     applyBoundaryCondition(structModel,...
%         'mixed',...
%         'Edge',indexSpider2,...
%         'u',0,...
%         'EquationIndex',2);
%     applyBoundaryCondition(structModel,...
%         'dirichlet',...
%         'Edge',[indexSpider1 indexSpider2],...
%         'u',@(location,state)slidingBC2(location,state));
    applyBoundaryCondition(structModel,...
        'dirichlet',...
        'Edge',[indexSpider1 indexSpider2],...
        'u',[0,0]);
    applyBoundaryCondition(structModel,...
        'dirichlet',...
        'Edge',[indexSpider1 indexSpider2],...
        'r',@(location,state)slidingBCr(location,state),...
        'h',@(location,state)slidingBCh(location,state));
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

dataForCF.kgm3_Fe     = mat.Rotor.kgm3;
dataForCF.kgm3_PM     = mat.LayerMag.kgm3;
dataForCF.w           = evalSpeed*pi/30;
dataForCF.psMagnet    = psMagnet;
dataForCF.psSleeve    = psSleeve;
dataForCF.kgm3_sleeve = mat.Sleeve.kgm3;

if strcmp(geo.RotType,'IM')
    dataForCF.kgm3_PM = mat.BarCond.kgm3;
end

% dataForCF.eleData = eleData;

% warning('Mass density set by default!!!')

dataForCmatrix.E_Fe     = mat.Rotor.E*1e9;
dataForCmatrix.E_PM     = mat.Rotor.E*1e9*1e-3;
dataForCmatrix.nu       = 0.3;
dataForCmatrix.psMagnet = psMagnet;
dataForCmatrix.E_sleeve = mat.Sleeve.E*1e9;
dataForCmatrix.psSleeve = psSleeve;

data4GeoMat.kgm3_Fe        = mat.Rotor.kgm3;
data4GeoMat.kgm3_PM        = mat.LayerMag.kgm3;
data4GeoMat.kgm3_sleeve    = mat.Sleeve.kgm3;
data4GeoMat.w              = evalSpeed*pi/30;
data4GeoMat.psMagnet       = psMagnet;
data4GeoMat.psSleeve       = psSleeve;
data4GeoMat.E_Fe           = mat.Rotor.E*1e9;
data4GeoMat.E_PM           = mat.Rotor.E*1e9*1e-3;
data4GeoMat.E_sleeve       = mat.Sleeve.E*1e9;
data4GeoMat.nu             = 0.3;
data4GeoMat.sigmaMaxFe     = mat.Rotor.sigma_max*1e6;
data4GeoMat.sigmaMaxSleeve = mat.Sleeve.sigma_max*1e6;
% if simSetup.meshShaft
    data4GeoMat.psShaft        = psShaft;
% end
data4GeoMat.kgm3_shaft     = 0;
data4GeoMat.E_shaft        = mat.Rotor.E*1e9*1e-6;

if strcmp(geo.RotType,'IM')
    data4GeoMat.kgm3_PM = mat.BarCond.kgm3;
end


% if (strcmp(geo.RotType,'IM')||geo.hs>0||~strcmp(mat.LayerMag.MatName,'Air'))
%     specifyCoefficients(structModel,...
%         'm',0,...
%         'd',0,...
%         'c',@(x,y)defineCmatrix(x,y,data4GeoMat),...
%         'a',0,...
%         'f',@(x,y)centrifugalForce(x,y,data4GeoMat));

 specifyCoefficients(structModel,...
        'm',0,...
        'd',0,...
        'c',cMat,...
        'a',0,...
        'f',@(x,y)centrifugalForce(x,y,data4GeoMat));

%     specifyCoefficients(structModel,...
%     'm',0,...
%     'd',0,...
%     'c',cMat,...
%     'a',0,...
%     'f',@(x,y)centrifugalForce(x,y,dataForCF));

% else
%     specifyCoefficients(structModel,...
%         'm',0,...
%         'd',0,...
%         'c',cMat,...
%         'a',0,...
%         'f',@(x,y)centrifugalForce(x,y,data4GeoMat));
% end

% if strcmp(geo.RotType,'IM')
%     specifyCoefficients(structModel,...
%         'm',0,...
%         'd',0,...
%         'c',@(x,y)defineCmatrix(x,y,dataForCF),...
%         'a',0,...
%         'f',@(x,y)centrifugalForce(x,y,dataForCmatrix));
% else
%     specifyCoefficients(structModel,...
%         'm',0,...
%         'd',0,...
%         'c',cMat,...
%         'a',0,...
%         'f',@(x,y)centrifugalForce(x,y,dataForCF));
% end

% keyboard
