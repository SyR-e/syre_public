% Copyright 2026
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

function [slotOut] = slotModel_eval_rectCond(geo)

kcu      = geo.win.kcu;
nCondIn  = geo.win.nCond;
tol      = geo.win.condIns;  % conductor insulation thickness [mm]
Aslot    = geo.Aslot;

r        = geo.r;
g        = geo.g;
l        = geo.l;
p        = geo.p;
q        = geo.q;
Qs       = geo.Qs;
ws       = geo.st;
stator   = geo.stator;
ttd      = geo.ttd;
n3ph     = geo.win.n3phase;
liner    = geo.win.liner;
offset   = geo.win.gapBotCond;

condDist = 0;
grid_pad = 1;

alphaSlot = 2*pi/(6*p*q*n3ph);   % slot pitch [rad]

index = stator(:,end);
slotMat = stator(index==1,:);

slotMat = rotateMatrix(slotMat,-alphaSlot/2);
slotMat = rotateMatrix(slotMat,pi/2);

hfig(1) = figure();
figSetting()
GUI_Plot_Machine(gca,slotMat);
hchild = get(gca,'Children');
X = [];
Y = [];
for ii=1:length(hchild)
    X = [X hchild(ii).XData];
    Y = [Y hchild(ii).YData];
end
close(hfig);

P_slot = polyshape(X, Y);
if liner > 0
    P_usable = polybuffer(P_slot,-liner);
else
    P_usable = P_slot;
end
A_usable = area(P_usable);

if isnan(kcu)
    rCond = geo.win.rCond;
    wCond = geo.win.wCond;
    hCond = geo.win.hCond;
else
    A_cu_total_req = kcu * A_usable;
    A_cu_per_wire  = A_cu_total_req / nCondIn;

    wCond = ws-2*liner-2*tol;
    hCond = A_cu_per_wire/wCond;
    rCond = 0;
end

wOut = wCond + 2*tol;
hOut = hCond + 2*tol;

erode_dist  = hOut/2 + condDist;

P_fit = polybuffer(P_usable,-erode_dist);

pitch_x = wOut + condDist;
pitch_y = hOut + condDist;

V = P_fit.Vertices;
xmin = min(V(:,1));
xmax = max(V(:,1));
ymin = min(V(:,2));
ymax = max(V(:,2));

% Expand bbox a bit (optional)
dx = (xmax-xmin);
dy = (ymax-ymin);
xmin = xmin - (grid_pad-1)*0.5*dx;
xmax = xmax + (grid_pad-1)*0.5*dx;
ymin = ymin - (grid_pad-1)*0.5*dy+offset;
ymax = ymax + (grid_pad-1)*0.5*dy+offset;

xs = (xmin:pitch_x:xmax).';
xs = 0;
ys = (ymin:pitch_y:ymax).';

[XX,YY] = meshgrid(xs,ys);
cand = [XX(:), YY(:)];

% Keep candidates inside feasible center region
in = isinterior(P_fit,cand(:,1),cand(:,2));
cand = cand(in,:);

% Choose closest to centroid (compactness)
c0 = centroid(P_fit);
if size(cand,1) >= nCondIn
    dist2 = sum((cand - c0).^2,2);
    [~, idx] = sort(dist2,'ascend');
    centers = cand(idx(1:nCondIn),:);
    N_fit = nCondIn;
else
    centers = cand;
    N_fit = size(cand,1);
end


slotOut.nCondMax = N_fit;
slotOut.nCond    = min([N_fit nCondIn]);
slotOut.nCondIn  = nCondIn;
slotOut.centers  = centers(:,1)+j*centers(:,2);
slotOut.rCond    = rCond;
slotOut.wCond    = wCond;
slotOut.hCond    = hCond;
