% Copyright 2019
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expres_trafs or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dataSet] = skinEffect_draw(dataSet)


clc
% read slot conductor data from dataSet

kcu           = dataSet.SlotFillFactor;
condType      = dataSet.SlotConductorType;
condIns       = dataSet.SlotConductorInsulation;
condShape     = dataSet.SlotConductorShape;
condRadius    = dataSet.SlotConductorRadius;
condWidth     = dataSet.SlotConductorWidth;
condHeight    = dataSet.SlotConductorHeight;
condNumber    = dataSet.SlotConductorNumber;
condBottomGap = dataSet.SlotConductorBottomGap;

% draw stator model (from GUI_APP_DrawMachine)
[~, ~, geo,per,mat] = data0(dataSet);
[geo,~,mat] = interpretRQ(geo.RQ,geo,mat);
geo.x0 = geo.r/cos(pi/2/geo.p);
fem.res = 0;
fem.res_traf = 0;
% nodes
[~,~,geo] = ROTmatr(geo,fem,mat);
[geo,~,~] = STATmatr(geo,fem);

%load([dataSet.currentpathname dataSet.currentfilename]);

% update slot conductor data
dataSet.SlotFillFactor          = kcu;
dataSet.SlotConductorType       = condType;
dataSet.SlotConductorInsulation = condIns;
dataSet.SlotConductorShape      = condShape;
dataSet.SlotConductorRadius     = condRadius;
dataSet.SlotConductorWidth      = condWidth;
dataSet.SlotConductorHeight     = condHeight;
dataSet.SlotConductorNumber     = condNumber;
dataSet.SlotConductorBottomGap  = condBottomGap;

geo.win.kcu        = kcu;
geo.win.condType   = condType;
geo.win.condIns    = condIns;
geo.win.condShape  = condShape;
geo.win.rCond      = condRadius;
geo.win.wCond      = condWidth;
geo.win.hCond      = condHeight;
geo.win.nCond      = condNumber;
geo.win.GapBotCond = condBottomGap;

materialCodes;  % load material codes

% load data from geo
condType  = geo.win.condType;   % conductor type (Round/Square)
tol       = geo.win.condIns;    % conductor insulation thickness [mm]
condHB    = geo.win.condShape; % conductor shape factor (h/b)

kcu          = geo.win.kcu;
nCondIn      = geo.win.nCond;

Aslot  = geo.Aslot;
r      = geo.r;
g      = geo.g;
l      = geo.l;
p      = geo.p;
q      = geo.q;
Qs     = geo.Qs;
stator = geo.stator;
ttd    = geo.ttd;
n3ph   = geo.win.n3phase;

alphaSlot = 2*pi/(6*p*q*n3ph);   % slot pitch [rad]

if strcmp(condType,'Round')
    slotOut = slotModel_eval_roundCond(geo);
else
    geo.win.kcu = NaN;
    slotOut = slotModel_eval_rectCond(geo);
end

rCond     = slotOut.rCond;
wCond     = slotOut.wCond;
hCond     = slotOut.hCond;
nCond     = slotOut.nCond;
% nCondMax  = slotOut.nCondMax;
xyCenters = slotOut.centers;

% Acond = pi*rCond^2+(wCond+2*rCond)*(hCond-2*rCond)+2*(wCond-2*rCond)*rCond+2*(hCond-2*rCond)*rCond;
Acond = wCond*hCond-(4*rCond^2-pi*rCond^2);

% define slot geometry
indCu  = 1;
indAir = indCu+Qs;

slotMat = stator(stator(:,9)==indCu,:);
airMat  = stator(stator(:,9)==indAir,:);
slotMat = slotMat(2:end-1,:);   % remove the first and the last line (boundary copper-air)
slotMat = [airMat(1:end/2,:);slotMat; airMat(end/2:end,:)];
slotMat(:,8) = codMatFeSta;  % change the material code of the slot side from copper to iron
slotMat(:,9) = 0;            % group zero for the iron

slotMat = rotateMatrix(slotMat,-alphaSlot/2);

% add the missing lines of the matrix (rotor airgap, stator airgap, tooth
% bottom

x1Fe = r*cos(alphaSlot/2);
y1Fe = r*sin(alphaSlot/2);
x2Fe = (r+g)*cos(alphaSlot/2);
y2Fe = (r+g)*sin(alphaSlot/2);

mesh_res = dimMesh(geo,'singt');
res_traf = mesh_res.res_traf;


slotMat = [slotMat
    +x1Fe +y1Fe +x2Fe +y2Fe NaN   NaN    0 codMatFeSta 0
    0     0     +x2Fe +y2Fe +x2Fe -y2Fe -1 codMatFeSta 0
    +x2Fe -y2Fe +x1Fe -y1Fe NaN   NaN    0 codMatFeSta 0
    0     0     +x1Fe -y1Fe +x1Fe +y1Fe +1 codMatFeSta 0
    ];

slotMat = rotateMatrix(slotMat,pi/2);


xySlot = [
    0 r+g/2     codMatAirSta res_traf 1
    0 r+g+ttd/2 codMatAirSta res_traf 1
    ];


matCond = [];
xyCond  = [];
indexCond = 1;

dh = wCond/2;
dw = hCond/2;
rc = rCond;

for ii=1:length(xyCenters)
    x0 = real(xyCenters(ii));
    y0 = imag(xyCenters(ii));
    xyCond = [xyCond; x0 y0 codMatCu res_traf 1];
    if rCond>0
        matCond = [matCond
            x0-dh    y0-dw+rc x0-dh    y0+dw-rc NaN      NaN       0 codMatCu indexCond
            x0-dh+rc y0+dw-rc x0-dh    y0+dw-rc x0-dh+rc y0+dw    -1 codMatCu indexCond
            x0-dh+rc y0+dw    x0+dh-rc y0+dw    NaN      NaN       0 codMatCu indexCond
            x0+dh-rc y0+dw-rc x0+dh-rc y0+dw    x0+dh    y0+dw-rc -1 codMatCu indexCond
            x0+dh    y0+dw-rc x0+dh    y0-dw+rc NaN      NaN       0 codMatCu indexCond
            x0+dh-rc y0-dw+rc x0+dh    y0-dw+rc x0+dh-rc y0-dw    -1 codMatCu indexCond
            x0+dh-rc y0-dw    x0-dh+rc y0-dw    NaN      NaN       0 codMatCu indexCond
            x0-dh+rc y0-dw+rc x0-dh+rc y0-dw    x0-dh    y0-dw+rc -1 codMatCu indexCond
            ];
    else
        matCond = [matCond
            x0-dh y0-dw x0-dh y0+dw NaN NaN 0 codMatCu indexCond
            x0-dh y0+dw x0+dh y0+dw NaN NaN 0 codMatCu indexCond
            x0+dh y0+dw x0+dh y0-dw NaN NaN 0 codMatCu indexCond
            x0+dh y0-dw x0-dh y0-dw NaN NaN 0 codMatCu indexCond
            ];
    end
    indexCond = indexCond+1;
end


disp(['Number of input conductor = ' int2str(nCondIn)])
disp(['Max number of conductor = ' int2str(nCond)])
if nCond<nCondIn
    warning('Actual conductor number lower than the input number');
else
    nCond = nCondIn;
    xyCond = xyCond(1:nCond,:);
    tmp = matCond(:,9);
    index = 1:1:numel(tmp);
    index = index(tmp<=nCond);
    matCond = matCond(index,:);
end

kcu = nCond*Acond/Aslot;

% update dataSet and geo
geo.win.kcu       = kcu;
geo.win.nCond     = nCond;
geo.win.rCond     = rCond;
geo.win.wCond     = wCond;
geo.win.hCond     = hCond;
geo.win.condShape = hCond/wCond;

dataSet.SlotFillFactor          = geo.win.kcu;
dataSet.SlotConductorType       = geo.win.condType;
dataSet.SlotConductorShape      = geo.win.condShape;
dataSet.SlotConductorRadius     = geo.win.rCond;
dataSet.SlotConductorWidth      = geo.win.wCond;
dataSet.SlotConductorHeight     = geo.win.hCond;
dataSet.SlotConductorNumber     = geo.win.nCond;

% save model
pathname = dataSet.currentpathname;
filename = dataSet.currentfilename;
outFolder = checkPathSyntax([filename(1:end-4) '_results\FEA results\']);
if ~exist([pathname outFolder],'dir')
    mkdir([pathname outFolder]);
end

resFolder = checkPathSyntax([pathname outFolder filename(1:end-4) '_slotModel\']);

mkdir(resFolder);

slot.xy      = xyCond;
slot.slotMat = slotMat;
slot.matCond = matCond;
slot.nCond   = nCond;

save([resFolder 'slotModel.mat'],'dataSet','geo','per','mat','slot');

copyfile([pathname filename(1:end-4) '.mat'],[resFolder filename(1:end-4) '.mat'])
copyfile([pathname filename(1:end-4) '.fem'],[resFolder filename(1:end-4) '.fem'])

if ~geo.XFEMMcreation
    openfemm(1);
    FEMM_initialize(geo,mat);
    mi_probdef(0,'millimeters','planar',1e-8,l,15);
    
    % draw slot (air)
    draw_lines_arcs(slotMat,1,mesh_res.res_traf);
    for ii=1:size(xySlot,1)
        mi_addblocklabel(xySlot(ii,1),xySlot(ii,2));
        mi_selectlabel(xySlot(ii,1),xySlot(ii,2));
        mi_setblockprop('Air',0,res_traf,'None',0,1,0);
        mi_clearselected;
    end
    
    mi_addboundprop('A=0', 0, 0, 0, 0, 0, 0, 0, 0, 0);  % inner and outer circles
    xBou = (x1Fe+x2Fe)/2;
    yBou = (y1Fe+y2Fe)/2;
    mi_selectarcsegment(xBou,yBou);
    mi_setarcsegmentprop(res_traf,'A=0',0,1);
    mi_clearselected;
    mi_selectarcsegment(xBou,-yBou);
    mi_setarcsegmentprop(res_traf,'A=0',0,1);
    mi_clearselected;
    
    % draw conductors
    draw_lines_arcs(matCond,3,res_traf);
    for ii=1:size(xyCond,1)
        circuitName = ['conductor_' int2str(ii)];
        mi_addcircprop(circuitName,0,1);
        mi_addblocklabel(xyCond(ii,1),xyCond(ii,2));
        mi_selectlabel(xyCond(ii,1),xyCond(ii,2));
        mi_setblockprop(mat.SlotCond.MatName,0,res_traf,circuitName,0,3,0);
        mi_clearselected;
    end
    
    % save model
    
    mi_saveas([resFolder 'slotModel.fem'])
    closefemm();
else
    error('Slot model creation with XFEMM not yet implemented!')
end
disp(['Slot model saved in:'])
disp([resFolder 'slotModel.fem'])

figure()
figSetting()
set(gca,'DataAspectRatio',[1 1 1]);
GUI_Plot_Machine(gca,slotMat);
GUI_Plot_Machine(gca,matCond);
saveas(gcf,[resFolder filename(1:end-4) '_slotModel.fig']);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['kcu = ' num2str(kcu)])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
