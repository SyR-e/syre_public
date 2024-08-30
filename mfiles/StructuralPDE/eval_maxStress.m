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

function [out]=eval_maxStress(structModel,sVonMises,geo,mat)
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

prc = 99; % percentile for GalFer Contest structural index

% Search for the stress in each ribs
sigmaRadMax = zeros(1,geo.nlay);
sigmaTanMax = zeros(1,geo.nlay);
sigmaRadAvg = zeros(1,geo.nlay);
sigmaTanAvg = zeros(1,geo.nlay);
nPointOverRad = zeros(1,geo.nlay);
nPointOverTan = zeros(1,geo.nlay);
sTmpTot = [];

for ii=1:geo.nlay
    if (strcmp(geo.RotType,'Circular')||strcmp(geo.RotType,'Seg'))
        % chech tangential ribs
        x1 = geo.xxD1k(ii);
        y1 = geo.yyD1k(ii);
        x2 = geo.xxD2k(ii);
        y2 = geo.yyD2k(ii);
        x0 = geo.xpont(ii);
        y0 = geo.ypont(ii);
        r1 = ((x0-x1)^2+(y0-y1)^2)^0.5;
        r2 = ((x0-x2)^2+(y0-y2)^2)^0.5;
        r = max([r1,r2,geo.pontT(ii)])+geo.pont0;
        fi = linspace(0,2*pi,51);
        X = x0+r*cos(fi);
        Y = y0+r*sin(fi);
        [X,Y] = rot_point(X,Y,90/geo.p*pi/180);
        index = inpolygon(x*1000,y*1000,X,Y);
        %xTmp = x(index)/1000;
        %yTmp = y(index)/1000;
        sTmp = sVonMises(index);
        sigmaTanMax(ii) = max(sTmp);
        sigmaTanAvg(ii) = mean(sTmp);
        nPointOverTan(ii) = sum(sTmp>sigma_max);
        sTmpTot = [sTmpTot;sTmp];

        % check radial ribs
        if geo.pontR(ii)>0
            if geo.radial_ribs_split(ii)
                X = [geo.XpontSplitBarSx(1,ii) geo.XpontSplitBarDx(1,ii) geo.XpontSplitBarDx(2,ii) geo.XpontSplitBarSx(2,ii)];
                Y = [geo.YpontSplitBarSx(1,ii) geo.YpontSplitBarDx(1,ii) geo.YpontSplitBarDx(2,ii) geo.YpontSplitBarSx(2,ii)];
            else
                X = [geo.XpontRadBarSx(ii) geo.XpontRadBarDx(ii) geo.XpontRadBarDx(ii) geo.XpontRadBarSx(ii)];
                Y = [geo.YpontRadBarSx(ii) geo.YpontRadBarDx(ii) 0 0];
            end
            [X,Y] = rot_point(X,Y,90/geo.p*pi/180);
            index = inpolygon(x*1000,y*1000,X,Y);
            %xTmp = x(index)/1000;
            %yTmp = y(index)/1000;
            sTmp = sVonMises(index);
            if ~isempty(sTmp)
                sigmaRadMax(ii) = max(sTmp);
                sigmaRadAvg(ii) = mean(sTmp);
                nPointOverRad(ii) = sum(sTmp>sigma_max);
            end
        else
            sTmp = [];
        end
        sTmpTot = [sTmpTot;sTmp];
    end
end



% % stress_R = nan(1,geo.nlay);
% % stress_T = nan(1,geo.nlay);
%
% if any(geo.pontR & strcmp(geo.RotType,'Seg'))
%     ii = (geo.pontR>0 & geo.radial_ribs_split==0);
%     xR_xo(ii) = (geo.B2k(ii)-geo.hc(ii)/2);
%     yR_xo(ii) = 0;
%
%     jj = (geo.pontR>0 & geo.radial_ribs_split==1);
%     xR_xo(jj) = (geo.B2k(jj)-geo.hc(jj)/2);
%     yR_xo(jj) = geo.YpBar1(jj)-geo.pontR(jj)/4;
%
%     [xR,yR] = rot_point(xR_xo,yR_xo,pi/2/geo.p);
%
%     for kk=1:geo.nlay
%         if xR(kk)>0
%             [~,pos_mid]  = min(abs((x-xR(kk)*10^-3))+abs((y-yR(kk)*10^-3)));
%             stress_R(kk) = sVonMises(pos_mid);
%         end
%     end
% end
%
% if strcmp(geo.RotType,'Seg')
%     xT_xo = geo.xpont;
%     yT_xo = geo.ypont;
%
%     [xT,yT] = rot_point(xT_xo,yT_xo,pi/2/geo.p);
%
%     for kk=1:geo.nlay
%         if xT(kk)>0
%             [~,pos_mid]  = min(abs((x-xT(kk)*10^-3))+abs((y-yT(kk)*10^-3)));
%             stress_T(kk) = sVonMises(pos_mid);
%         end
%     end
% end

out.MaxStress     = MaxStress;
out.x_max         = x_max;
out.y_max         = y_max;
out.x_over        = x_over;
out.y_over        = y_over;
out.nodesOver     = nodesOver;
out.sigmaRadMax   = sigmaRadMax;
out.sigmaTanMax   = sigmaTanMax;
out.sigmaRadAvg   = sigmaRadAvg;
out.sigmaTanAvg   = sigmaTanAvg;
out.nPointOverRad = nPointOverRad;
out.nPointOverTan = nPointOverTan;
out.sigmaTotPrc   = prctile(sTmpTot,prc);

% out.stress_T  = stress_T;
% out.stress_R  = stress_R;

