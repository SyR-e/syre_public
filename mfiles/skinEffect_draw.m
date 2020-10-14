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

kcu        = dataSet.SlotFillFactor;
condType   = dataSet.SlotConductorType;
condIns    = dataSet.SlotConductorInsulation;
condShape  = dataSet.SlotConductorShape;
condRadius = dataSet.SlotConductorRadius;
condWidth   = dataSet.SlotConductorWidth;
condHeight = dataSet.SlotConductorHeight;
condNumber = dataSet.SlotConductorNumber;

% draw stator model (from GUI_APP_DrawMachine)
[~, ~, geo,per,mat] = data0(dataSet);
[geo,~,mat] = interpretRQ(dataSet.RQ,geo,mat);
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

geo.win.kcu       = kcu;
geo.win.condType  = condType;
geo.win.condIns   = condIns;
geo.win.condShape = condShape;
geo.win.rCond     = condRadius;
geo.win.wCond     = condWidth;
geo.win.hCond     = condHeight;
geo.win.nCond     = condNumber;


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

if ~isnan(kcu) % conductor size computed from kCu
    Acond = Aslot*kcu/nCondIn;
    if strcmp(condType,'Round')
        rCond = (Acond/pi)^0.5;
        wCond = 2*rCond;
        hCond = 2*rCond;
    else
        wCond = (Acond/condHB)^0.5;
        hCond = condHB*wCond;
        rCond = geo.win.rCond;
        if rCond>wCond/2
            rCond=wCond/2;
        end
        if rCond>hCond/2
            rCond=hCond/2;
        end
    end
else          % conductor size imposed, kCu computed at the end
    if strcmp(condType,'Round')
        rCond = geo.win.rCond;
        wCond = 2*rCond;
        hCond = 2*rCond;
        kcu = (nCondIn*pi*rCond^2)/Aslot;
        Acond = pi*rCond^2;
    else
        wCond = geo.win.wCond;
        hCond = geo.win.hCond;
        rCond = geo.win.rCond;
        if rCond>wCond/2
            rCond=wCond/2;
        end
        if rCond>hCond/2
            rCond=hCond/2;
        end
        kcu = (nCondIn*wCond*hCond)/Aslot;
        Acond = wCond*hCond;
    end
end

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


xySlot = [
    r+g/2     0 codMatAirSta res_traf 1
    r+g+ttd/2 0 codMatAirSta res_traf 1
    ];

% search the main slot size

x1 = slotMat(3,1);
y1 = slotMat(3,2);
x2 = slotMat(3,3);
y2 = slotMat(3,4);

% computation of the slot profile
hSlot    = x2-x1;
wSlotMin = 2*y1;
wSlotMax = 2*y2;

if wSlotMin==wSlotMax
    m = 0;
    q = wSlotMin/2;
else
    m = (y2-y1)/(x2-x1);
    q = y1-m*x1;
end

% if strcmp(condType,'Round')
%     dh = 2*(rCond+tol)*sin(pi/3);
%     db = 2*(rCond+tol);
% else
%     dh = hCond+2*tol;
%     db = wCond+2*tol;
% end

if strcmp(condType,'Round')
    nh = floor(hSlot/(hCond)); % max number of height divisions
else
    nh = floor(hSlot/(hCond+2*tol)); % max number of height divisions
end

% nb = floor(bSlot/(wCond+2*tol)); % max number of base divisions


% slot filling, from the bottom slot
matCond = [];
xyCond  = [];
indexCond = 1;

dh = hCond/2;
dw = wCond/2;
rc = rCond;

for xx=1:1:nh
    if xx==1
        x0 = x2-tol-hCond/2;
        xMin = x0-hCond;
        % xMax = x0+hCond;
        yLim = m*xMin+q;
        nw = floor(2*yLim/(wCond+2*tol));
        if rem(nw,2)==0
            flagEven = 1;
        else
            flagEven = 0;
        end
    else
        if strcmp(condType,'Round')
            x0 = x0-2*(rCond+tol)*sin(pi/3);
            xMin = x0-hCond;
            % xMax = x0+hCond;
            yLim = m*xMin+q;
            nw = floor(2*yLim/(wCond+2*tol));
            if flagEven
                if rem(nw,2)==0
                    nw = nw-1;
                end
                flagEven = 0;
            else
                if rem(nw,2)>0
                    nw = nw-1;
                end
                flagEven = 1;
            end
        else
            x0 = x0-hCond-2*tol;
            xMin = x0-hCond;
            % xMax = x0+hCond;
            yLim = m*xMin+q;
            nw = floor(2*yLim/(wCond+2*tol));
            if nw/2==0
                flagEven = 1;
            else
                flagEven = 0;
            end
        end
    end
    
    % y0 = (nb-1)*(wCond/2+tol);
    
    for yy=1:1:nw
        y0 = (nw-1)*(wCond/2+tol)-(wCond+2*tol)*(yy-1); % conductor center
        xyCond = [xyCond; x0 y0 codMatCu res_traf 1];
        % conductor definition
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
    
end

nCond = indexCond-1;
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
outFolder = [filename(1:end-4) '_results\FEA results\'];
if ~exist([pathname outFolder],'dir')
    mkdir([pathname outFolder]);
end

resFolder = [pathname outFolder filename(1:end-4) '_slotModel\'];

mkdir(resFolder);

slot.xy      = xyCond;
slot.slotMat = slotMat;
slot.matCond = matCond;
slot.nCond   = nCond;

save([resFolder filename(1:end-4) '_slotModel.mat'],'dataSet','geo','per','mat','slot');

copyfile([pathname filename(1:end-4) '.mat'],[resFolder filename(1:end-4) '.mat'])
copyfile([pathname filename(1:end-4) '.fem'],[resFolder filename(1:end-4) '.fem'])

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

mi_saveas([resFolder filename(1:end-4) '_slotModel.fem'])
closefemm();
disp(['Slot model saved in:'])
disp([resFolder filename(1:end-4) '_slotModel.fem'])

figure()
figSetting()
set(gca,'DataAspectRatio',[1 1 1]);
GUI_Plot_Machine(gca,slotMat);
GUI_Plot_Machine(gca,matCond);
saveas(gcf,[resFolder filename(1:end-4) '_slotModel.fig']);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(['kcu = ' num2str(kcu)])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
