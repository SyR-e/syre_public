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

function [out]=eval_maxStress(structModel,sVonMises,R,geo,mat)

x         = structModel.Mesh.Nodes(1,:);
y         = structModel.Mesh.Nodes(2,:);
sigma_max = mat.Rotor.sigma_max*1e6;
MaxStress = max(sVonMises);
pos_max = find(sVonMises==MaxStress);

x_max = x(pos_max)*10^3;
y_max = y(pos_max)*10^3;

nodesOver = sum (sVonMises>sigma_max);
pos_over = find(sVonMises>sigma_max);
x_over = x(pos_over)*10^3;
y_over = y(pos_over)*10^3;

% Deformation calculations
% Displacements
Ux = R.NodalSolution(:,1);
Uy = R.NodalSolution(:,2);
% Original nodal positions
x0 = structModel.Mesh.Nodes(1,:)';
y0 = structModel.Mesh.Nodes(2,:)';
% Deformation magnitude
DefMagn = sqrt(Ux.^2 + Uy.^2);
MaxDef  = max(DefMagn);

% Deformed positions
De_x = x0 + Ux;
De_y = y0 + Uy;
% Radial deformed position
De_r = sqrt(De_x.^2 + De_y.^2);
max_r = (max(De_r))*1e3;
% Air-gap Clearance
r0 = geo.r;
ag = geo.g;
agclear = ag - (max_r - r0);

prc = 99; % percentile for GalFer Contest structural index

% Search for the stress in each ribs
sigmaRadMax = zeros(1,geo.nlay);
sigmaTanMax = zeros(1,geo.nlay);
nPointOverRad = zeros(1,geo.nlay);
nPointOverTan = zeros(1,geo.nlay);
sTmpTot = [];
sTmpTan = [];
sTmpRad = [];

for ii=1:geo.nlay
    if (strcmp(geo.RotType,'Circular')||strcmp(geo.RotType,'Seg'))
        % check tangential ribs
        x1 = geo.xxD1k(ii);
        y1 = geo.yyD1k(ii); % First edge of the rib
        x2 = geo.xxD2k(ii);
        y2 = geo.yyD2k(ii); % Second edge of the rib
        x0 = geo.xpont(ii);
        y0 = geo.ypont(ii); % Rib center
        r1 = ((x0-x1)^2+(y0-y1)^2)^0.5; % Distance center-first edge
        r2 = ((x0-x2)^2+(y0-y2)^2)^0.5; % Distance center-second edge
        r = max([r1,r2,geo.pontT(ii)])+geo.pont0;
        fi = linspace(0,2*pi,51);
        X = x0+r*cos(fi);
        Y = y0+r*sin(fi); % Circle based on ribs size
        [X,Y] = rot_point(X,Y,90/geo.p*pi/180);
        index_T = inpolygon(x*1000,y*1000,X,Y); % Points inside circle
        %xTmp = x(index)/1000;
        %yTmp = y(index)/1000;
        sTmp_T = sVonMises(index_T); % Nodes inside circle
        sigmaTanMax(ii) = max(sTmp_T);
        % sigmaTanAvg(ii) = mean(sTmp);
        nPointOverTan(ii) = sum(sTmp_T>sigma_max);
        sTmpTot = [sTmpTot;sTmp_T];
        sTmpTan = [sTmpTan;sTmp_T];

        % check radial ribs
        if geo.pontR(ii)>0
            if geo.radial_ribs_split(ii) % Quadrilateral
                X = [geo.XpontSplitBarSx(1,ii) geo.XpontSplitBarDx(1,ii) geo.XpontSplitBarDx(2,ii) geo.XpontSplitBarSx(2,ii)];
                Y = [geo.YpontSplitBarSx(1,ii) geo.YpontSplitBarDx(1,ii) geo.YpontSplitBarDx(2,ii) geo.YpontSplitBarSx(2,ii)];
            else % Trapezoidal
                X = [geo.XpontRadBarSx(ii) geo.XpontRadBarDx(ii) geo.XpontRadBarDx(ii) geo.XpontRadBarSx(ii)];
                Y = [geo.YpontRadBarSx(ii) geo.YpontRadBarDx(ii) 0 0];
            end
            [X,Y] = rot_point(X,Y,90/geo.p*pi/180);
            index_R = inpolygon(x*1000,y*1000,X,Y); % Points inside circle
            %xTmp = x(index)/1000;
            %yTmp = y(index)/1000;
            sTmp_R = sVonMises(index_R); % Nodes inside circle
            if ~isempty(sTmp_R)
                sigmaRadMax(ii) = max(sTmp_R);
                % sigmaRadAvg(ii) = mean(sTmp);
                nPointOverRad(ii) = sum(sTmp_R>sigma_max);
            end
        else
            sTmp = [];
            sTmp_R = [];
        end
        sTmpTot = [sTmpTot;sTmp_R];
        sTmpRad = [sTmpRad;sTmp_R];
    end
end

out.MaxStress       = MaxStress;
out.x_max           = x_max;
out.y_max           = y_max;
out.nodesOver       = nodesOver;
out.x_over          = x_over;
out.y_over          = y_over;
out.MaxDef          = MaxDef*1e3;
out.agclear         = agclear;
out.sigmaRadMax     = sigmaRadMax;
out.sigmaTanMax     = sigmaTanMax;
out.nPointOverRad   = nPointOverRad;
out.nPointOverTan   = nPointOverTan;
out.sigmaTotPrc     = prctile(sTmpTot,prc);
out.Tan_sigmaTotPrc = prctile(sTmpTan,prc);
out.Rad_sigmaTotPrc = prctile(sTmpRad,prc);


% out.stress_T  = stress_T;
% out.stress_R  = stress_R;

