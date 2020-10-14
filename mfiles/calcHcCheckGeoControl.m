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


% Input: dalpha [deg], hc_pu [p.u.]
% Output: alpha [deg], hc [mm]

function geo = calcHcCheckGeoControl(geo)

r = geo.r;              % Raggio del rotore al traferro
p = geo.p;              % Paia poli
nlay = geo.nlay;        % N° layers
R = geo.R;              % Raggio ext
g = geo.g;              % Traferro
lt = geo.lt;            % Lunghezza denti
pont0 = geo.pont0;      % Ponticelli al traferro
% dalpha = geo.dalpha;
alpha = cumsum(geo.dalpha);
hc_pu = geo.hc_pu;
hfe_min = geo.hfe_min;        % min tickness of each steel flux guide

x0 = geo.x0;                % center of the circular barriers profiles
Ar = x0 - r * tan(pi/2/p);  % max allowed shaft radius
% geo.ArMaxAdmis = Ar;
htot = r - Ar;              % rotor space available radialwise (air + steel)
ly = R - r - g - lt;        % stator yoke
lyr = 1.0 * ly;             % rotor yoke ( = ly)
la = r - Ar -lyr;           % maximum allowed insulation (total air, base value)
hc_half_min = la/nlay/8;    % min hc/2
% hc_half_min = pont0;

beta = 180/pi * calc_apertura_cerchio(pi/180*alpha,r,x0);   % arc angles respect to center x0
rbeta = (x0 - r * cos(alpha*pi/180))./(cos(beta*pi/180));   % radius of circular barriers profiles

delta=(1/(nlay)*sum(hc_pu));
hfeqMax=rbeta(end)-rbeta(1)-(nlay-1)*2*hc_half_min;         % total Fe, max (accounts for nlay-1 flux carriers)
hfeqMin=(nlay-1)*hfe_min;                                   % total Fe, min

hfeq=hfeqMax-(hfeqMax-hfeqMin)*(delta);

la=rbeta(end)-rbeta(1)-hfeq;
% coeff_gamma=hc_pu(1)/2+hc_pu(nlay)/2+sum(hc_pu(2:nlay-1));
% la=rbeta(end)-rbeta(1)-(nlay-1)*hfe_min;
%
% laprimo=la/(nlay-1);

hc = [];
error_type1=1;
% conta=1;
if (nlay==1)

    % max hc according to alpha min
    hc_half_max1 = (alpha*pi/180/(1+alpha*pi/180)*(r-pont0));
    % (needs division by 2 .. don't know why but it works)
    hc_half_max1 = hc_half_max1 * 2;

    % max hc according to alpha max (27 Jan 2011)
    temp_alpha_hfemin = hfe_min/r; % rad
    temp_alpha_hc_2 = pi/(2*p) - alpha*pi/180 - temp_alpha_hfemin;
    hc_half_max2 = (temp_alpha_hc_2/(1+temp_alpha_hc_2)*(r-pont0));
    hc_half_max = min(hc_half_max1,hc_half_max2);
    %     hc_pu(1) = 1;
    hc(1) = hc_pu(1) * hc_half_max * 2;
    if hc(1)<2*hc_half_min
        hc(1)=2*hc_half_min;
    end

    if hc(1)>2*hc_half_max
        hc(1)=2*hc_half_max;
    end

end
% layer one is the outermost, layer nlay is the innermost
for jj = nlay:-1:1
    if (jj == nlay)

        hc_half_max = 0.5*la/(0.5+(sum(hc_pu(2:nlay-1)))/(hc_pu(nlay))+hc_pu(1)/hc_pu(nlay)/2);
        hc_nlay_temp_half=0.5*la/(0.5+(sum(hc_pu(2:nlay-1)))/(hc_pu(nlay))+hc_pu(1)/hc_pu(nlay)/2);
        hc(jj) = 2*hc_half_max(1);
        if hc(jj)<2*hc_half_min
            hc(jj)=2*hc_half_min;
        end        
    else
        if jj == 1
            hc_half_max = min((hc_pu(jj)/hc_pu(nlay))*hc_nlay_temp_half,(r-pont0-x0+rbeta(jj)));
            hc(jj) = 2*hc_half_max(1);
            if hc(jj)<2*hc_half_min && hc_half_min<=(r-pont0(jj)-x0+rbeta(jj))
                hc(jj)=2*hc_half_min;
            end
        else
            hc_half_max =(hc_pu(jj)/hc_pu(nlay))*hc_nlay_temp_half;
            hc(jj) = 2*hc_half_max(1);
            if hc(jj)<2*hc_half_min
                hc(jj)=2*hc_half_min;
            end
        end
    end
end
% end

% 2014/02/25 MG Determinazione dei punti di barriera sull'asse q:
hcIni=hc;
hc=abs(hcIni);
Bx0=x0-rbeta;
B1k=Bx0-hc./2;
B2k=Bx0+hc./2;

for k=1:nlay-1
    %% #4 vincolo 1-n° layer overlap);
    if (B2k(k+1)>=B1k(k))   % questa condizione vale invece per tutte le barriere
        Dint=B2k(k+1)-B1k(k);
        B2p=B1k(k+1)+(hc(k)+hc(k+1)-Dint-hfe_min(k))/(1+hc_pu(k)/hc_pu(k+1));
        B1p=B2k(k)-(hc(k)+hc(k+1)-Dint-hfe_min(k))/(1+hc_pu(k+1)/hc_pu(k));
        error_type1=1;
        disp('#1 conflict between two barrier')
        %         keyboard
        %% #5 vincolo 1-n° intersezione arie, spessore lato barriera<pont0/2 --> equa ripartizione aria ferro');
        % condizione vale nel caso in cui muovendosi non c'è più spazio per l'aria condizione critica,
        % scelte casuali erronee, ma si privilegia la fattibilità della
        % macchina per proseguire nell'ottimizzazione...
        if ((Bx0(k)<B1p)||(Bx0(k+1)>B2p)|| ((B2p-Bx0(k+1))<hc_half_min) || ((Bx0(k)-B1p)<hc_half_min))
            B1p=Bx0(k)-(Bx0(k)-Bx0(k+1))/3;
            B2p=Bx0(k+1)+(Bx0(k)-Bx0(k+1))/3;
            disp('#2 low space for iron and air')
        end % end #5
        B2k(k+1)=B2p;
        B1k(k)=B1p;
    else
        disp('');
        error_type1(k)=0;
    end % end #4
end % end for nlay
%% Safety control for last flux barrier check feseability space between air and iron of the spyder at the air-gap and along the q axis
[xc_temp,yc_temp]=calc_intersezione_cerchi(r,rbeta(nlay),x0);
dPointEndBar=calc_distanza_punti([xc_temp,yc_temp],[r*cos(pi/2/p), r*sin(pi/2/p)]);

if (dPointEndBar<(Bx0(nlay)-B1k(nlay)))
    B1k(nlay)=Bx0(nlay)-dPointEndBar+hfe_min(nlay)/2;
end

if (B1k(nlay)<geo.Ar+hfe_min(end))    % questa condizione vale per l'ultima barriera di flux
    geo.Ar=max(B1k(nlay)-hfe_min(end),Bx0(nlay)-dPointEndBar);
end

hc=B2k-B1k;

%% 2014/02/25 MG Condition fo the first flux barrier if hc is to high and nlay bigger the length of the barrier too smal to be drawn:
if ((r-pont0(1)-hc(1)/2)<=x0-rbeta(1))
    temp_hc1=2*(r+rbeta(1)-x0-(1.5)*pont0(1));
    if (temp_hc1>0)
        hc(1)=temp_hc1;
        B2k(1)=B1k(1)+temp_hc1;
    end
end
geo.hcIni=hcIni;    % For control the intial hc value are saved

%% end seurity control 1th flux barrier
%%
% SALVO hc NELLA STRUTTURA 'geo'.
geo.hc = hc;
% geo = rmfield(geo,'hc_pu');

% CALCOLO 'nlay_effettivo', 'delta' E 'nr' E LI SALVO IN GEO.
% geo.delta = [alpha(1) diff(alpha) (90/p)-alpha(end)];
geo.B1k=B1k;
geo.B2k=B2k;
% geo.nr = ceil(90 ./ geo.delta ) * 2;
