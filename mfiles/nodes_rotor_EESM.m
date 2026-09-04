function [geo,mat,temp] = nodes_rotor_EESM(geo,mat)

ri         = geo.Ar;                % internal rotor radius     (mm)
p          = geo.p;                 % pole pairs                (-)
beta_pu    = geo.dalpha_pu;         % pole arc extension        (pu)
hry        = geo.lyr;               % rotor yoke height         (mm)
hp1        = geo.hpb;               % rotor pole body length    (mm)
hp2        = geo.hph;               % rotor pole head length    (mm)
wp         = geo.wp;                % rotor pole body width     (mm)
wb         = geo.wb;                % coil width                (mm)
hb         = geo.hb;                % coil height               (mm)
MinTol     = geo.pont0;             % minimun tolerance         (mm)
g0         = geo.g;                 % air gap                   (mm)
re_lim     = geo.r;                 % air gap radius            (mm)
thHead_deg = geo.thHead_deg;        % Head slope                (°)
r_fillet   = geo.r_fillet;          % Pole head fillet          (mm)
thYoke_deg = geo.thYoke_deg;        % Yoke slope                (°)
r_bfillet  = geo.r_bfillet;         % Pole yoke fillet          (mm)
PoleShape  = geo.headShape;

narcs = 15;    % # of arcs for approximation of the pole head
beta_pu = min(max(beta_pu, 0.5), 0.95); % beta saturation
beta = pi/p*beta_pu;                    % pole arc extension                    (rad)
sat = 1/(1+p)+0.5;                      % variable pole body width saturation   (pu)

theta = pi/p/2;                 % bisector                              (rad)

% parameters' check 
if MinTol     <  0.01    MinTol     = 0.01;   end
if p          <= 0       p          = 1;      end
if hry        <= MinTol  hry        = MinTol; end
if hp1        <= MinTol  hp1        = MinTol; end
if hp2        <= 2       hp2        = 2;      end
ry = ri + hry;                                      % rotor's yoke external radius          (mm)
rb = ry + hp1;                                      % rotor's body's tooth external radius  (mm)
re = rb + hp2;                                      % rotors external radius                (mm)
wp = max(MinTol,min(wp,2*ry*sin(theta)));
if wb         <= MinTol  wb         = MinTol; end
if hb         <= MinTol  hb         = MinTol; end
if thHead_deg <= 0       thHead_deg = 0;      end
if r_fillet   <= 0       r_fillet   = 0;      end
if thYoke_deg <= 0       thYoke_deg = 0;      end
if r_bfillet  <= 0       r_bfillet  = 0;      end


if re > re_lim
    re = re_lim;
    if re*cos(beta/2) < rb
        rb = re*cos(beta/2);
        [xp4, yp4] = intersezione_retta_circonferenza(0,0,ry,0,wp/2);
        if rb < xp4 + MinTol
            xp4 = rb - MinTol;
            ry = sqrt(xp4^2+yp4^2);
            if ry < ri
                warning('Change beta')
            end
        end
    end
else
    hp2 = hp2+ re_lim -re;
    re = re_lim;
end

% Calcolo testa del polo
xM = re;
yM = 0;
tau = beta/2/narcs;
gamma = beta/2-(0:narcs-1)*tau;

switch PoleShape
    case 0      % CRF
        K = ones(1,narcs);

    case 1      % 1/cos
        K = 1./cos(gamma*p);

    case 2      % Third Harmonic
        k3 = 0.23;
        den = cos(gamma*p)+k3*cos(3*gamma*p);
        K = max(den)./den;
end

xp7 = zeros(1,narcs);
yp7 = zeros(1,narcs);

% for ii = 1:narcs
%     [xp7(ii),yp7(ii)] = rot_point(xM,yM,gamma(ii));
%     xp7(ii) = xp7(ii)-g0*(K(ii)-1)*cos(gamma(ii));
%     yp7(ii) = yp7(ii)-g0*(K(ii)-1)*sin(gamma(ii));
% end

xp7 = (xM - g0 .* (K - 1)) .* cos(gamma);% - yM .* sin(gamma);
yp7 = (xM - g0 .* (K - 1)) .* sin(gamma);% + yM .* cos(gamma);

xp1 = ri;
yp1 = 0;

xp2 = ri*cos(theta);
yp2 = ri*sin(theta);

xp3 = ry*cos(theta);
yp3 = ry*sin(theta);

[xp4, yp4] = intersezione_retta_circonferenza(0,0,ry,0,wp/2);

xp4o = xp4;
yp4o = yp4;

xp5 = rb;
yp5 = yp4;

xp6 =rb;
headCutParallel = 1;
if headCutParallel == 1
    yp6 = yp7(1);
else
    yp6 = xp6*tan(beta/2);
end


% xM = re;
% yM = 0;
% tau = beta/2/narcs;
% kmax = 0;
% xp7 = zeros(1,narcs);
% yp7 = zeros(1,narcs);
% PoleShape = geo.headShape;
% 
% if PoleShape == 0       % CRF Pole Shape
%     for ii=1:narcs
%         gamma = beta/2 - (ii-1)*tau;
%         [xp7(ii),yp7(ii)] = rot_point(xM,yM,gamma);
%     end
% 
% elseif PoleShape == 1   % 1/cos Pole Shape
%     for ii=1:narcs
%         gamma = beta/2 - (ii-1)*tau;
%         [xp7(ii),yp7(ii)] = rot_point(xM,yM,gamma);
%         xp7(ii) = xp7(ii) - g0*(1/cos(gamma*p)-1)*cos(gamma);
%         yp7(ii) = yp7(ii) - g0*(1/cos(gamma*p)-1)*sin(gamma);
%     end
% 
% elseif PoleShape == 2   % Third Harmonic Pole Shape
%     for ii=1:narcs
%         gamma = beta/2 - (ii-1)*tau;
%         k3 = 0.23;
%         if cos(gamma*p)+k3*cos(3*gamma*p) > kmax
%             kmax = cos(gamma*p)+k3*cos(3*gamma*p);
%         end
%     end
%     for ii=1:narcs
%         gamma = beta/2 - (ii-1)*tau;
%         [xp7(ii),yp7(ii)] = rot_point(xM,yM,gamma);
%         xp7(ii) = xp7(ii) - g0*(kmax/(cos(gamma*p)+k3*cos(3*gamma*p))-1)*cos(gamma);
%         yp7(ii) = yp7(ii) - g0*(kmax/(cos(gamma*p)+k3*cos(3*gamma*p))-1)*sin(gamma);
%     end
% end

if xp7(1) < xp5     %rotor pole beta saturation (extreme situation for p = 1)
    if xp7(1) < ri+3*MinTol
        xTemp = ri+3*MinTol;
        yTemp = sqrt(re^2-xTemp^2);
        beta = 2*atan(yTemp/xTemp);
        beta_pu = beta*p/pi;
        tau = beta/2/narcs;
        if PoleShape == 0       % CRF Pole Shape
            for ii=1:narcs
                gamma = beta/2 - (ii-1)*tau;
                [xp7(ii),yp7(ii)] = rot_point(xM,yM,gamma);
            end

        elseif PoleShape == 1   % 1/cos Pole Shape
            for ii=1:narcs
                gamma = beta/2 - (ii-1)*tau;
                [xp7(ii),yp7(ii)] = rot_point(xM,yM,gamma);
                xp7(ii) = xp7(ii) - g0*(1/cos(gamma*p)-1)*cos(gamma);
                yp7(ii) = yp7(ii) - g0*(1/cos(gamma*p)-1)*sin(gamma);
            end

        elseif PoleShape == 2   % Third Harmonic Pole Shape
            for ii=1:narcs
                gamma = beta/2 - (ii-1)*tau;
                [xp7(ii),yp7(ii)] = rot_point(xM,yM,gamma);
                xp7(ii) = xp7(ii) - g0*(kmax/(cos(gamma*p)+k3*cos(3*gamma*p))-1)*cos(gamma);
                yp7(ii) = yp7(ii) - g0*(kmax/(cos(gamma*p)+k3*cos(3*gamma*p))-1)*sin(gamma);
            end
        end
        xp4 = xp7(1) - MinTol;
    end
    xp5 = xp7(1);
    xp6 = xp7(1);
    yp6 = yp7(1);
    if xp5 < xp4 + 2*MinTol
        xp4 = xp5 - 2*MinTol;
        dist = sqrt(xp4^2+yp4^2);
        xp3 = dist*cos(theta);
        yp3 = dist*sin(theta);
    end
end


% if xp7(1) < xp5     %beta saturation per 1/cos()
%     test = 50;
%     tau = beta/2/test;
%     flag = 0;
%     ii = 1;
%     while flag == 0
%         gamma = (ii-1)*tau;
%         [xp7_t(ii),yp7_t(ii)] = rot_point(xM,yM,gamma);
%         xp7_t(ii) = xp7_t(ii) - g0*(1/cos(gamma*p)-1)*cos(gamma);
%         yp7_t(ii) = yp7_t(ii) - g0*(1/cos(gamma*p)-1)*sin(gamma);
%         if xp7_t(ii) < xp6
%             [xc,yc] = calc_center_given_3pts(xp7_t(ii),yp7_t(ii),xp7_t(ii-1),yp7_t(ii-1),xp7_t(ii-2),yp7_t(ii-2));
%             d = calc_distanza_punti_altern(xc,yc,xp7_t(ii),yp7_t(ii));
%             [xc, yc] = rot_point(xc,yc,pi/2);
%             a = -1;
%             b = 0;
%             c = rb;
%             a1 = b;
%             b1 = a;
%             [xTemp,yTemp] = calc_int_retta_circ_gen(xc,yc,d,a1,b1,c);
%             [xTemp,yTemp] = rot_point(xTemp,yTemp,-pi/2);
%             xp6 = xTemp(1);
%             yp6 = yTemp(1);
%             beta =2*(atan(yTemp(1)/xTemp(1)));
%             beta_pu = beta/pi*p;
%             tau = beta/2/narcs;
%             for ii=1:narcs
%                 gamma = beta/2 - (ii-1)*tau;
%                 [xp7(ii),yp7(ii)] = rot_point(xM,yM,gamma);
%                 xp7(ii) = xp7(ii) - g0*(1/cos(gamma*p)-1)*cos(gamma);
%                 yp7(ii) = yp7(ii) - g0*(1/cos(gamma*p)-1)*sin(gamma);
%             end
%             flag = 1;
%         end
%         ii = ii + 1;
%     end
% end

if yp5/yp6 > sat
    yp5 = sat*yp6;
    [xp4, yp4] = intersezione_retta_circonferenza(0,0,ry,0,yp5);
    xp4o = xp4;
    yp4o = yp4;
    wp = 2*yp5;
end

[xp3,yp3] = intersezione_tra_rette(1,0,-xp4,-yp2/xp2,1,0); % Yoke non CRF ma dritto (modificato anche build_matrix_EESM p3-p4)

% Pole yoke angle
if thYoke_deg > 89
    thYoke_deg = 89;
end

if thYoke_deg ~= 0
    [a1, b1, c1] = retta_ruotata_attorno_punto(1, 0, -xp4, xp4, yp4, thYoke_deg);
    [xp3,yp3] = intersezione_tra_rette(a1,b1,c1,-yp2/xp2, 1, 0);
    if thYoke_deg > 0 && xp3 < xp2
        xp3 = xp2;
        yp3 = yp2;
        thYoke_deg = -atand((xp4-xp3)/(yp4-yp3));
    end
end

% Pole head angle
if thHead_deg > 89
    thHead_deg = 89;
end
if thHead_deg ~= 0
    x_off = (yp6-yp5)*tand(thHead_deg);
    if x_off > xp5-xp4
        x_off = xp5-xp4;
        thHead_deg = atand(x_off/(yp6-yp5));
    end
    xp5 = xp5 - x_off;
end

temp = struct();
% Pole Yoke Fillet
if r_bfillet ~= 0
    dist_p3p4 = hypot(xp3-xp4o,yp3-yp4o);
    if r_bfillet > min(abs(xp5-xp4),dist_p3p4)
        r_bfillet = min(abs(xp5-xp4),dist_p3p4);
    end
    if r_bfillet ~= 0
        [a1,b1,c1] = retta_per_2pti(xp3,yp3,xp4,yp4);
        [a2,b2,c2] = retta_per_2pti(xp4,yp4,xp5,yp5);
        [a1o,b1o,c1o] = calc_retta_offset(a1,b1,c1,-r_bfillet);
        [a2o,b2o,c2o] = calc_retta_offset(a2,b2,c2,-r_bfillet);
        [xc_bfillet,yc_bfillet] = intersezione_tra_rette(a1o,b1o,c1o,a2o,b2o,c2o);
        [xp4(1),yp4(1)] = proiezione_punto_retta(a1,b1,c1,xc_bfillet,yc_bfillet);
        [xp4(2),yp4(2)] = proiezione_punto_retta(a2,b2,c2,xc_bfillet,yc_bfillet);
        temp.xc_bfillet = xc_bfillet;
        temp.yc_bfillet = yc_bfillet;
    end
end

% Pole Head Fillet
if r_fillet ~= 0
    if r_fillet > min(abs(xp5-xp4(end)),abs(yp6-yp5))
        r_fillet = min(abs(xp5-xp4(end)),abs(yp6-yp5));
    end
        [a1,b1,c1] = retta_per_2pti(xp6,yp6,xp5,yp5);
        [a2,b2,c2] = retta_per_2pti(xp4(end),yp4(end),xp5,yp5);
        [a1o,b1o,c1o] = calc_retta_offset(a1,b1,c1,+r_fillet);
        [a2o,b2o,c2o] = calc_retta_offset(a2,b2,c2,-r_fillet);
        [xc_fillet,yc_fillet] = intersezione_tra_rette(a1o,b1o,c1o,a2o,b2o,c2o);
        [xp5(2),yp5(2)] = proiezione_punto_retta(a1,b1,c1,xc_fillet,yc_fillet);
        if a2 <= 1e-13
            a2 = 0;
        end
        [xp5(1),yp5(1)] = proiezione_punto_retta(a2,b2,c2,xc_fillet,yc_fillet); % ERRORE SE a2 è circa 0 ma non proprio zero tipo e-16
        temp.xc_fillet = xc_fillet;
        temp.yc_fillet = yc_fillet;
    
end

% Coil design
wb_min = max(yp5(end) - yp5(1),yp4(1) - yp4(end));
wb_max = max(yp6-yp5(1));
wb = max(min(wb,wb_max),wb_min);

xc2 = xp5;
yc2 = yp5;

[a2,b2,c2] = retta_per_2pti(xp5(end),yp5(end),xp6,yp6);
[xc3,yc3] = intersezione_tra_rette(0,1,-(wb+wp/2),a2,b2,c2);

% if round(hb,2) >= round(xp6-xp3,2)
if hb >= (xp6-xp4(end))
    hb = round(xp6-xp3,2);
    yc4 = min(wb+wp/2,yp3);
    [a1,b1,c1] = retta_per_2pti(xp3,yp3,xp4(1),yp4(1));
    xc4 = (-b1*yc4-c1)/a1;

    xc1 = xp4;
    yc1 = yp4;
else%if hb < (xp5(1)-xp4(end))
    % if hb > (xp6-xp4(end))
    %     hb = (xp6-xp4(end));
    % end
    % hb = max(hb,xp6-xp4(end));
    xc4 = xc3 - hb;
    yc4 = min(wb+wp/2,xc4*tan(theta)-MinTol);
    
    xc1 = xc4;
    yc1 = yp4(end);
end

for ii = 1:1:7
    x_var = sprintf('xp%d',ii);
    y_var = sprintf('yp%d',ii);
    temp.(x_var) = eval(x_var);
    temp.(y_var) = eval(y_var);
end
for ii = 1:1:4
    x_var = sprintf('xc%d',ii);
    y_var = sprintf('yc%d',ii);
    temp.(x_var) = eval(x_var);
    temp.(y_var) = eval(y_var);
end
temp.xM = xM;
temp.yM = yM;
temp.narcs = narcs;
temp.rcirc = re+g0;
temp.theta = theta;

temp.xc = mean([xc1,xc2,xc3,mean(xc4)]);
temp.yc = mean([yc1,yc2,yc3,mean(yc4)]);

temp.xair = rb*cos(0.98*theta);
temp.yair = rb*sin(0.98*theta);

temp.xmag = NaN;
temp.ymag = NaN;
temp.zmag = NaN;

temp.mirrorFlag = 1;
temp.mirrorFlagAir = 1;

ry = sqrt(xp4o^2+yp4o^2);

% computation of coil cross-section (circa slot area)
if r_bfillet ~= 0 && length(xc1) > 1 % Se la bobina è più bassa non interseca il raccordo allo yoke -> xc1 è un punto
    ang = linspace(atan2(yc1(1)-yc_bfillet,xc1(1)-xc_bfillet),atan2(yc1(2)-yc_bfillet,xc1(2)-xc_bfillet), 21);
    xc1 = [xc1(1) xc_bfillet+r_bfillet*cos(ang) xc1(2)];
    yc1 = [yc1(1) yc_bfillet+r_bfillet*sin(ang) yc1(2)];
end
if r_fillet ~= 0
    ang = linspace(atan2(yc2(1)-yc_fillet,xc2(1)-xc_fillet), atan2(yc2(2)-yc_fillet,xc2(2)-xc_fillet), 21);
    xc2 = [xc2(1) xc_fillet+r_fillet*cos(ang) xc2(2)];
    yc2 = [yc2(1) yc_fillet+r_fillet*sin(ang) yc2(2)];
end
X = [xc1 xc2 xc3 xc4 xc1(1)];
Y = [yc1 yc2 yc3 yc4 yc1(1)];

geo.Acoilf = polyarea(X,Y);


geo.p  = p;
geo.dalpha_pu  = round(beta_pu,3);
geo.dalpha     = round(beta,3);
geo.lyr        = round(ry-ri,2);
geo.hpb        = round(xp6-ry,2);
geo.hph        = round(xM-xp6,2);
geo.wp         = round(wp,2);
geo.wb         = wb;
geo.hb         = round(hb,2);
geo.r_fillet   = round(r_fillet,2);
geo.thHead_deg = round(thHead_deg,1);
geo.r_bfillet  = round(r_bfillet,2);
geo.thYoke_deg = round(thYoke_deg,1);
% geo.pont0      = MinTol;
% geo.g          = g0;

% Additional parameters for compatibility
geo.AreaC = 0;
geo.AreaE = 0;

% For COMSOL sweep parameters - COMSOL Parametrization with Mechs
geo.hre = calc_distanza_punti([xp6,yp6],[xp7,yp7]);
geo.wcu = calc_distanza_punti([xp6,yp6],[xc3,yc3]);

% For COMSOL Fillet numeric stability
geo.filsub = ry-xp3;

%geo.wd = wd;
%geo.beta_pu = beta_pu;
%geo.beta = beta;
% = ri;
% = p;
% = hry;
% = hp1;
% = hp2;
% = wp;
% = wb;
% = MinTol

% figure
% figSetting
% axis equal
% plot([xp1,xp2,xp3,xp4,xp5,xp6,xp7,xM],[yp1,yp2,yp3,yp4,yp5,yp6,yp7,yM],'b')
% plot([xc1,xc2,xc3,xc4,xc1],[yc1,yc2,yc3,yc4,yc1],'r')
% plot([temp.xc, temp.xair],[temp.yc, temp.yair],'g*')
