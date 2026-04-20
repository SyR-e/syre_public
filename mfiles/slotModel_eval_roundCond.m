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

function [slotOut] = slotModel_eval_roundCond(geo)

kcu      = geo.win.kcu;
nCondIn  = geo.win.nCond;
condType = geo.win.condType;   % conductor type (Round/Square)
tol      = geo.win.condIns;    % conductor insulation thickness [mm]
Aslot    = geo.Aslot;
r        = geo.r;
g        = geo.g;
l        = geo.l;
p        = geo.p;
q        = geo.q;
Qs       = geo.Qs;
stator   = geo.stator;
ttd      = geo.ttd;
n3ph     = geo.win.n3phase;
liner    = geo.win.liner;
condDist = 0;               % distance between conductors
offset   = geo.win.gapBotCond;

grid_pad = 1; % oversampling factor (1-3)

alphaSlot = 2*pi/(6*p*q*n3ph);   % slot pitch [rad]

index = geo.stator(:,end);
slotMat = geo.stator(index==1,:);

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

% build slot polyshape
P_slot = polyshape(X, Y);
if liner > 0
    P_usable = polybuffer(P_slot,-liner);
else
    P_usable = P_slot;
end
A_usable = area(P_usable);

if isnan(kcu)
    rCond = geo.win.rCond;
    wCond = 2*rCond;
    hCond = 2*rCond;
    d_cu  = 2*rCond;
    d_out = d_cu + 2*tol;
    r_out = d_out/2;
else
    A_cu_total_req = kcu * A_usable;
    A_cu_per_wire  = A_cu_total_req / nCondIn;

    d_cu  = 2*sqrt(A_cu_per_wire/pi);       % [mm] bare copper diameter
    d_out = d_cu + 2*tol;              % [mm] outer diameter including enamel

    r_out = d_out/2;
    rCond = d_cu/2;
    wCond = 2*rCond;
    hCond = 2*rCond;
end

erode_dist = r_out+condDist;

P_fit = polybuffer(P_usable,-erode_dist);

pitch_x = d_out+condDist;
pitch_y = (sqrt(3)/2)*(d_out+condDist);

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

% Create staggered rows
y_vals = (ymin:pitch_y:ymax).';
cand = [];

for r = 1:numel(y_vals)
    yy = y_vals(r);
    if mod(r,2)==1
        x_start = xmin;
    else
        x_start = xmin + 0.5*pitch_x;
    end
    x_vals = (x_start : pitch_x : xmax).';
    cand = [cand; [x_vals, yy*ones(size(x_vals))]];
end

in = isinterior(P_fit,cand(:,1),cand(:,2));
cand = cand(in,:);

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