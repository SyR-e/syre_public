function kMsat = compute_kMsat_EESM(geo, mat, per)
%compute_kMsat_EESM  (refactor veloce, no parfor)
% - riduce drasticamente le chiamate a interp1 usando griddedInterpolant
% - vettorializza il loop su alpha in STEP2
% - fattorizza la valutazione AT_tot(Bg) in una funzione locale

% --- Impostazioni manuali (per ora)
Nr        = geo.win.Nf;
Ir        = per.if0;
Ns        = geo.win.Ns;
if isfield(per,'id')
    Id        = per.id;
else
    Id = 0;
end
AT_target = Nr*Ir+Ns/2/geo.p*Id;
Bg_guess  = 1.2;

% =========================================================
% INPUT (da geo/mat/per) + conversioni mm -> m dove serve
p    = geo.p;
q    = geo.q;

n3ph = geo.win.n3phase;
m    = 3*n3ph;                 % coerente con Nslot = 6*p*q*n3ph

L    = geo.l * 1e-3;           % [m]

Dso  = (geo.R*2) * 1e-3;       % [m]
Dsi  = ((geo.r+geo.g)*2) * 1e-3; % [m]
g0   = geo.g * 1e-3;           % [m]
Dri  = (geo.Ar*2) * 1e-3;      % [m]

kb   = geo.dalpha_pu;          % [-]

hst  = geo.lt * 1e-3;          % [m]
h1   = geo.ttd * 1e-3;         % [m]

wst  = geo.wt * 1e-3;          % [m]
acs  = geo.acs;                % [-]
wp   = geo.wp * 1e-3;          % [m]

wst0 = (pi*(geo.r+geo.g)/p/q/3/n3ph)*(1-acs) * 1e-3; % [m]

h2_mm = (pi*(geo.r+geo.g)/geo.p/geo.q/3/geo.win.n3phase - geo.wt)/2 * tand(geo.tta);
h2    = h2_mm * 1e-3;          % [m]

hp1  = 0*geo.hph * 1e-3;       % [m]
hp2  = geo.hpb * 1e-3;         % [m]

BH   = mat.Stator.BH;          % Nx2 [B(T)  H(A/m)]
b    = BH(:,1);
h    = BH(:,2);

% =========================================================
% Preliminari
Nslot = m*q*2*p;
Dro   = Dsi - 2*g0;
tau_s = pi*Dsi/Nslot;
RY    = geo.lyr;
SY    = geo.ly;
mu_0  = 4*pi*1e-7;

dmax  = g0/cosd(kb*90);
davg  = (dmax+g0)/2;
ee    = (tau_s-wst0)/2/davg;
kc    = tau_s/(tau_s-2/pi*davg*(2*ee*atan(ee)-log(1+ee^2)));

% =========================================================
% Interpolanti (molto più veloci di interp1 ripetuto)
H_of_B = griddedInterpolant(b, h, 'linear', 'linear');  % H(B)

% =========================================================
% STEP 1: costruisci AT_td(Bg) e AT_d(Bg) su griglia
Bg_max_v = linspace(0, 5, 100);
AT_t     = zeros(size(Bg_max_v));

% costanti per STEP1
npt = 20;
dx  = (hst - h2)/npt;
% x campionato al centro di ciascun segmento
xVec = ( (1:npt) - 0.5 ) * dx;

% ws/gamma dipendono solo da x (non da Bg)
wsVec    = pi*(Dsi + 2*xVec)/Nslot - wst;
gammaVec = (wsVec./wst) * mu_0;     % semplificato

H0 = H_of_B(0); % utile in bisezione

for k = 1:numel(Bg_max_v)
    Bd = Bg_max_v(k);
    if Bd == 0
        AT_t(k) = 0;
        continue
    end

    % ---- "AT_teeth"
    Phi_t = Bd * tau_s * L;

    % Piece 1 (h1, wst0)
    B1  = Phi_t/(wst0*L);
    H1  = H_of_B(B1);
    AT1 = H1*h1;

    % Piece 2 (integrale su hst-h2)
    B0  = Phi_t/(wst*L);
    AT2 = 0;

    for ii = 1:npt
        gamma = gammaVec(ii);

        % Risolvi g(B)=B + gamma*H(B) - B0 = 0 con bisezione robusta
        Blo = 0;
        Bhi = max(B0, 1e-6);

        glo = Blo + gamma*H0 - B0;
        ghi = Bhi + gamma*H_of_B(Bhi) - B0;

        % espandi Bhi finché bracketta
        kexp = 0;
        while glo*ghi > 0 && kexp < 30
            Bhi = 2*Bhi;
            ghi = Bhi + gamma*H_of_B(Bhi) - B0;
            kexp = kexp + 1;
        end

        % bisezione (iterazioni fisse)
        Bt = Bhi; %#ok<NASGU>
        Ht = H_of_B(Bhi);
        for itB = 1:40
            Bmid = 0.5*(Blo + Bhi);
            Hmid = H_of_B(Bmid);
            gmid = Bmid + gamma*Hmid - B0;

            if abs(gmid) < 1e-6 || (Bhi - Blo) < 1e-6
                Ht = Hmid;
                break
            end

            if glo*gmid <= 0
                Bhi = Bmid;
                ghi = gmid; %#ok<NASGU>
            else
                Blo = Bmid;
                glo = gmid;
            end
        end

        AT2 = AT2 + Ht*dx;
    end

    AT_t(k) = AT1 + AT2;
end

AT_d  = (Bg_max_v/mu_0) * (g0*kc);
AT_td = AT_t + AT_d;

% interpolanti per STEP2 (AT_td(Bg) e inversa Bg(AT))
ATtd_of_Bg = griddedInterpolant(Bg_max_v, AT_td, 'linear', 'linear');
ATd_of_Bg  = griddedInterpolant(Bg_max_v, AT_d,  'linear', 'linear');

% Inversa: serve x monotono -> ordina per AT_td
[AT_td_s, idx] = sort(AT_td(:));
Bg_s = Bg_max_v(idx);
% rimuovi duplicati per sicurezza
[AT_td_u, iu] = unique(AT_td_s, 'stable');
Bg_u = Bg_s(iu);
Bg_of_AT = griddedInterpolant(AT_td_u, Bg_u, 'linear', 'linear');

% =========================================================
% STEP 2: inversione AT_tot(Bg) = AT_target con bisezione

Nalpha = 100;
Dalpha = (pi/2)/Nalpha;
alpha  = ((0:Nalpha-1) + 0.5) * Dalpha;
cosa   = cos(alpha);
maskKb = alpha < kb*pi/2;

% costanti per vettorializzazione
const_phi = (Dro + g0) * Dalpha / (2*p);
const_sy  = (Dso - SY) * Dalpha / 2;

% funzione locale: AT_tot(Bg)
    function [AT_tot, AT_d_val] = eval_ATtot(Bg)
        AT_d_val = ATd_of_Bg(Bg);

        if Bg == 0
            AT_tot = 0;
            return
        end

        AT_max = ATtd_of_Bg(Bg);

        % Ba(alpha)
        ATa = AT_max * cosa;         % 1xNalpha
        Ba  = zeros(size(ATa));
        if any(maskKb)
            Ba(maskKb) = Bg_of_AT(ATa(maskKb));
        end

        % cumulata e flusso
        Stot = cumsum(Ba);
        Phi_alpha_vec = Stot * const_phi;              % 1xNalpha
        Bsy = Phi_alpha_vec / SY;                      % 1xNalpha

        AT_sy1 = const_sy * sum(H_of_B(Bsy));

        Phi_alpha_end = Phi_alpha_vec(end);
        Phi_p = 2 * Phi_alpha_end * L;

        Bp    = Phi_p/(wp*L);
        AT_pv = H_of_B(Bp) * (hp1 + hp2);

        Bry   = (Phi_p/2)/(RY*L);
        AT_ryv = H_of_B(Bry) * (Dri + RY) * pi/(2*p)/2;

        AT_tot = AT_max + AT_sy1 + AT_pv + AT_ryv;
    end

% --- Bracketing
Bg_lo = 0;
Bg_hi = max(1e-6, Bg_guess);

F_lo = -AT_target;

[AT_tot_hi, ~] = eval_ATtot(Bg_hi);
F_hi = AT_tot_hi - AT_target;

Bg_hi_max = 10;
while F_lo*F_hi > 0 && Bg_hi < Bg_hi_max
    Bg_hi = 2*Bg_hi;
    [AT_tot_hi, ~] = eval_ATtot(Bg_hi);
    F_hi = AT_tot_hi - AT_target;
end

if F_lo*F_hi > 0
    kMsat = NaN;
    return
end

tolBg = 1e-4;
tolF  = 1e-2;
maxIt = 60;

AT_tot_found = NaN;
AT_d_found   = NaN;

for it = 1:maxIt
    Bg_mid = 0.5*(Bg_lo + Bg_hi);

    [AT_tot_mid, AT_d_mid] = eval_ATtot(Bg_mid);
    F_mid = AT_tot_mid - AT_target;

    if abs(F_mid) < tolF || (Bg_hi - Bg_lo) < tolBg
        AT_tot_found = AT_tot_mid;
        AT_d_found   = AT_d_mid;
        break
    end

    if F_lo*F_mid <= 0
        Bg_hi = Bg_mid;
        F_hi  = F_mid;
    else
        Bg_lo = Bg_mid;
        F_lo  = F_mid;
    end
end

kMsat = AT_tot_found / AT_d_found;

end
