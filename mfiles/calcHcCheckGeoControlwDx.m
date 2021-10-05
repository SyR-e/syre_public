% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.


% Input: dalpha [deg], hc_pu [p.u.], dx [p.u.]
% Output: hc [mm]

function geo = calcHcCheckGeoControlwDx(geo)

r = geo.r;              % Raggio del rotore al traferro
p = geo.p;              % Paia poli
nlay = geo.nlay;        % N° layers
R = geo.R;              % Raggio ext
g = geo.g;              % Traferro
lt = geo.lt;            % Lunghezza denti
pont0 = geo.pont0;      % Ponticelli al traferro
Ar    = geo.Ar;         % Shaft radius
hfe_min = geo.hfe_min;        % min tickness of each steel flux guide
Bx0=geo.Bx0;
pontT = geo.pontT;
x0 = geo.x0;

dalpha = geo.dalpha;
hc_pu = geo.hc_pu;
dx=geo.dx;

alpha = cumsum(dalpha);

% max shaft radius
ArLim = x0 - r * tan(pi/2/p);
if Ar>ArLim
    Ar=ArLim;
    disp('#1Ar Shaft radius is bigger than the maximum allowed')
    geo.Ar=ArLim;
end

% rotor space available radialwise (air + steel)
% htot = r - ArLim;
ly = R - r - g - lt;    % stator yoke
lyr = 1.0 * ly;         % lower limit of the total steel tickness
la = r - ArLim -lyr;    % upper limit of the total air insulation

% min hc
hc_half_min = la/nlay/8;% occhio che nn deve essere troppo piccolo se no le barriere verranno sempre eccessivamente piccole, ma?! :-|

beta = 180/pi * calc_apertura_cerchio(pi/180*alpha,r,x0);
rbeta = (x0 - r * cos(alpha*pi/180))./(cos(beta*pi/180));

hc_pu(hc_pu>1) = 1;
hc_pu(hc_pu<-1) = -1;

% Per il momento la taratura è a mano:
hfeqMax=rbeta(end)-rbeta(1)-(nlay-1)*2*hc_half_min; % max iron (carriers)
hfeqMin=(nlay-1)*hfe_min;                           % min iron (carriers)
hfeq=hfeqMax-(hfeqMax-hfeqMin)*(mean(hc_pu));
la=rbeta(end)-rbeta(1)-hfeq;
% coeff_gamma=hc_pu(1)/2+hc_pu(nlay)/2+sum(hc_pu(2:nlay-1));

%% preliminary evaluation of hc
hc = [];
% error_type1=1;
% conta=1;
if (nlay==1)
    %% max hc according to alpha min
    hc_half_max1 = (alpha*pi/180/(1+alpha*pi/180)*(r-pont0));
    % (needs division by 2 .. don't know why but it works)
    hc_half_max1 = hc_half_max1 * 2;
    %% max hc according to alpha max (27 Jan 2011)
    temp_alpha_hfemin = hfe_min/r; % rad
    temp_alpha_hc_2 = pi/(2*p) - alpha*pi/180 - temp_alpha_hfemin;
    hc_half_max2 = (temp_alpha_hc_2/(1+temp_alpha_hc_2)*(r-pont0));
    hc_half_max = min(hc_half_max1,hc_half_max2);
    hc(1) = hc_pu(1) * hc_half_max * 2;
    if hc(1)<2*hc_half_min
        hc(1)=2*hc_half_min;
    end
    if hc(1)>2*hc_half_max
        hc(1)=2*hc_half_max;
    end
else
    % layer one is the outermost, layer nlay is the innermost
    for jj = nlay:-1:1
        switch jj
            case nlay
                hc_half_max       = 0.5*la/(0.5+(sum(hc_pu(2:nlay-1)))/(hc_pu(nlay))+hc_pu(1)/hc_pu(nlay)/2);
                hc_nlay_temp_half = 0.5*la/(0.5+(sum(hc_pu(2:nlay-1)))/(hc_pu(nlay))+hc_pu(1)/hc_pu(nlay)/2);
                hc(jj) = 2*hc_half_max(1);
                if hc(jj)<2*hc_half_min
                    hc(jj)=2*hc_half_min;
                end
            case 1
                hc_half_max = min((hc_pu(jj)/hc_pu(nlay))*hc_nlay_temp_half,(r-pont0-x0+rbeta(jj)));
                hc(jj) = 2*hc_half_max(1);
                if hc(jj)<2*hc_half_min && hc_half_min<=(r-pont0-x0+rbeta(jj))
                    hc(jj)=2*hc_half_min;
                end
            otherwise
                hc_half_max =(hc_pu(jj)/hc_pu(nlay))*hc_nlay_temp_half;
                hc(jj) = 2*hc_half_max(1);
                if hc(jj)<2*hc_half_min
                    hc(jj)=2*hc_half_min;
                end
        end
    end
end

%% evaluation of B1k and B2k
dx(dx>0.99) = 0.99;
dx(dx<-0.99) = -0.99;

B1k=Bx0-hc/2+dx.*hc/2;
B2k=Bx0+hc/2+dx.*hc/2;

for ii=1:nlay-1
    if B2k(ii+1)>Bx0(ii)
        B2k(ii+1) = Bx0(ii)-hfe_min;
    end
    if B1k(ii)<Bx0(ii+1)
        B1k(ii) = Bx0(ii+1)+hfe_min;
    end
end


beta = 180/pi * calc_apertura_cerchio(pi/180*alpha,r,x0);
rbeta = (x0 - r * cos(alpha*pi/180))./(cos(beta*pi/180));
[xpont,ypont] = calc_intersezione_cerchi(r-pontT, rbeta, x0);

% sort the B1k, B2k nodes
all_nodes = [ArLim];
for j = 1:nlay
    all_nodes = [all_nodes B1k(end-j+1) B2k(end-j+1)];
end

all_nodes(all_nodes>xpont(1)) = xpont(1); % check outer barrier interference 
all_nodes(all_nodes<ArLim) = ArLim; % check shaft interference

% sort nodes twist (inherent interference)
all_nodes = [all_nodes(1) sort(all_nodes(2:end))];
% check min dimension (hfe_min for Fe and Air)
delta_nodes = diff(all_nodes);
delta_nodes(delta_nodes<hfe_min)=hfe_min;

% reconstruct B1k and B2k vectors
all_nodes = ArLim + cumsum(delta_nodes);
if (all_nodes(end)>xpont(1))
    delta_nodes = delta_nodes * ((xpont(1)-ArLim)/(all_nodes(end)-ArLim));
    all_nodes = ArLim + cumsum(delta_nodes);
end

for j = 1:nlay
    B2k(j) = all_nodes(end-2*j+2);
    B1k(j) = all_nodes(end-2*j+1);
end

%% output
if B1k(B1k>=Bx0)
    B1k(B1k>=Bx0) = Bx0(B1k>=Bx0) - 0.1;
end
if B2k(B2k<=Bx0)
    B2k(B2k<=Bx0) = Bx0(B2k<=Bx0) + 0.1;
end

hc = B2k-B1k;
dx = (mean([B1k;B2k])-Bx0)./(hc/2);

geo.hc_pu = hc_pu;
geo.hc=hc;
geo.dx=dx;

geo.Ar = Ar;
geo.B1k=B1k;
geo.B2k=B2k;

