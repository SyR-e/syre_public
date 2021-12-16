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

function [structModel] = syre2pde(geo,mat,simSetup)

flagFull  = simSetup.flagFull;
meshSize  = simSetup.meshSize;
shaftBC   = simSetup.shaftBC;
evalSpeed = simSetup.evalSpeed;


% load data from input structures
rotor = geo.rotor;
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

% create the geometry of the problem

nEle = max(rotor(:,9));

density = zeros(1,nEle);
matCode = zeros(1,nEle);
materialCodes;


for ii=1:nEle
    rotTmp = rotor(rotor(:,9)==ii,:);
    rotTmp(:,1:6) = rotTmp(:,1:6)/1e3;
    [~,~,ps(ii)] = calcAreaShape(rotTmp);
    switch rotTmp(1,8)
        case codMatAirRot
            density(ii) = 0;
            matCode(ii) = codMatAirRot;
        case codMatBar
            density(ii) = mat.LayerMag.kgm3;
            matCode(ii) = codMatBar;
        case codMatFeRot
            density(ii) = mat.Rotor.kgm3;
            matCode(ii) = codMatFeRot;
        case codMatShaft
            density(ii) = mat.Shaft.kgm3;
            matCode(ii) = codMatShaft;
    end
end


psRotor = ps(nEle-1);

for ii=1:nEle
    if matCode(ii)==codMatAirRot
        if ~exist('psHole','var')
            psHole = ps(ii);
        else
            psHole = union(psHole,ps(ii));
        end
    elseif matCode(ii)==codMatBar
        psHole = subtract(psHole,ps(ii));
        if ~exist('psMagnet','var')
            psMagnet = ps(ii);
        else
            psMagnet = union(psMagnet,ps(ii));
        end
    end
end

psRotor = subtract(ps(nEle-1),psHole);
psRotor = subtract(psRotor,ps(end));

if ~exist('psMagnet','var')
    psMagnet = [];
end


if flagFull
    if geo.ps~=(2*geo.p)
        nRep = 2*geo.p/geo.ps-1;
        psBaseR = psRotor;
        psBaseM = psMagnet;
        for ii=1:nRep
            th = (pi/geo.p*geo.ps)*ii*180/pi;

            psNew = rotate(psBaseR,th);
            psRotor = union(psRotor,psNew);
            if ~isempty(psMagnet)
                psNew = rotate(psBaseM,th);
                psMagnet = union(psMagnet,psNew);
            end
        end

        geo.ps = 2*geo.p;

    end
    warning('Full motor selected!!!')
end


tr = triangulation(psRotor);
geometryFromMesh(structModel,tr.Points',tr.ConnectivityList');
generateMesh(structModel,'Hmax',Hmax,'Hmin',Hmin,'Hgrad',Hgrad,'GeometricOrder','quadratic');


% Boundary conditions: shaft is fixed

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
%         'Edge',[indexSpider1 indexSpider2],...
%         'u',@(x,y)slidingBC(x,y),...
%         'Vectorized','on');

%     applyBoundaryCondition(structModel,...
%         'dirichlet',...
%         'Edge',indexSpider1,...
%         'h',[1,0;0,1],...
%         'r',[cos(min(th));sin(min(th))],...
%         'Vectorized','on');
%     applyBoundaryCondition(structModel,...
%         'dirichlet',...
%         'Edge',indexSpider2,...
%         'h',[1,0;0,1],...
%         'r',[cos(max(th));sin(max(th))],...
%         'Vectorized','on');
%     applyBoundaryCondition(structModel,...
%         'dirichlet',...
%         'Edge',[indexSpider1 indexSpider2],...
%         'h',[1,0;0,1],...
%         'r',@(x,y)slidingBC_r(x,y),...
%         'Vectorized','off');
    applyBoundaryCondition(structModel,...
        'dirichlet',...
        'Edge',[indexSpider1],...
        'u',[0,0]);
    applyBoundaryCondition(structModel,...
        'dirichlet',...
        'Edge',[indexSpider2],...
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


setInitialConditions(structModel,[0;0]);
structModel.SolverOptions.ReportStatistics = 'on';

% keyboard


















