% Copyright 2024
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

narcs = 15;                     % # of arcs for approximation of the pole head

% parameters' check 
if MinTol     <  0.01    MinTol     = 0.01;   end
if p          <= 0       p          = 1;      end
if hry        <= MinTol  hry        = MinTol; end
if hp1        <= MinTol  hp1        = MinTol; end
if hp2        <= 2       hp2        = 2;      end
if wp         <= MinTol  wp         = MinTol; end
if wb         <= MinTol  wb         = MinTol; end
if hb         <= MinTol  hb         = MinTol; end
if thHead_deg <= 0       thHead_deg = 0;      end
if r_fillet   <= 0       r_fillet   = 0;      end

% beta saturation
if beta_pu > 0.95 
    beta_pu = 0.95;
elseif beta_pu < 0.5
    beta_pu = 0.5;
end

% parameters' calculation
theta = pi/p/2;                 % bisector                              (rad)
beta = pi/p*beta_pu;            % pole arc extension                    (rad)
sat = 1/(1+p)+0.5;              % variable pole body width saturation   (pu)
ry = ri + hry;                  % rotor's yoke external radius          (mm)
rb = ry + hp1;                  % rotor's body's tooth external radius  (mm)
re = rb + hp2;                  % rotors external radius                (mm)

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


xp1 = ri;
yp1 = 0;

xp2 = ri*cos(theta);
yp2 = ri*sin(theta);

xp3 = ry*cos(theta);
yp3 = ry*sin(theta);

if yp3 > wp/2
    [xp4, yp4] = intersezione_retta_circonferenza(0,0,ry,0,wp/2);
else
    xp4 = xp3;
    yp4 = yp3;
    wp=2*yp3;
end

xp5 = rb;
yp5 = yp4;

xp6 =rb;
yp6 = xp6*tan(beta/2);

xM = re;
yM = 0;
tau = beta/2/narcs;
kmax = 0;
xp7 = zeros(1,narcs);
yp7 = zeros(1,narcs);
PoleShape = geo.headShape;

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
        k3 = 0.23;
        if cos(gamma*p)+k3*cos(3*gamma*p) > kmax
            kmax = cos(gamma*p)+k3*cos(3*gamma*p);
        end
    end
    for ii=1:narcs
        gamma = beta/2 - (ii-1)*tau;
        [xp7(ii),yp7(ii)] = rot_point(xM,yM,gamma);
        xp7(ii) = xp7(ii) - g0*(kmax/(cos(gamma*p)+k3*cos(3*gamma*p))-1)*cos(gamma);
        yp7(ii) = yp7(ii) - g0*(kmax/(cos(gamma*p)+k3*cos(3*gamma*p))-1)*sin(gamma);
    end
end

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
    wp = 2*yp5;
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
% Fillet
if r_fillet ~= 0
    if r_fillet > min(abs(xp5-xp4),abs(yp6-yp5))
        r_fillet = min(abs(xp5-xp4),abs(yp6-yp5));
    end
    if r_fillet ~= 0
        [a1,b1,c1] = retta_per_2pti(xp6,yp6,xp5,yp5);
        [a2,b2,c2] = retta_per_2pti(xp4,yp4,xp5,yp5);
        [a1o,b1o,c1o] = calc_retta_offset(a1,b1,c1,+r_fillet);
        [a2o,b2o,c2o] = calc_retta_offset(a2,b2,c2,-r_fillet);
        [xc_fillet,yc_fillet] = intersezione_tra_rette(a1o,b1o,c1o,a2o,b2o,c2o);
        [xp5(2),yp5(2)] = proiezione_punto_retta(a1,b1,c1,xc_fillet,yc_fillet);
        [xp5(1),yp5(1)] = proiezione_punto_retta(a2,b2,c2,xc_fillet,yc_fillet);
        temp.xc_fillet = xc_fillet;
        temp.yc_fillet = yc_fillet;
    end
end

% Coil design
if hb > calc_distanza_punti([xp5(1) yp5(1)],[xp4 yp4])
    hb = calc_distanza_punti([xp5(1) yp5(1)],[xp4 yp4]);
end
if wb > yp6-yp5(1)
    wb = yp6-yp5(1);
elseif wb < yp5(end)-yp5(1)
    wb = yp5(end)-yp5(1);
end

xc1 = xp5(1) - hb;
yc1 = yp4;

xc2 = xp5;
yc2 = yp5;

% xc3 = xp5(end);
yc3 = yp5(1) + wb;
xc3 = xp6 - (yp6-yc3)*tand(thHead_deg); % In case of pole head with angle

xc4 = xc1;
yc4 = yc3;

% if yc4 > yp6
%     yc4 = yp6;
%     yc3 = yp6;
% end

if yc4 > xc4*tan(0.98*theta)
    yc4 = xc4*tan(0.98*theta);%-MinTol;
end

if yc4 < yc1
    if yp4 > yp6
        [a1,b1,c1] = retta_per_2pti(0,0,xp6,yp6);
        [a2,b2,c2] = retta_per_2pti(xp4,yp4,xp5(1),yp5(1));
        [xc1,yc1] = intersezione_tra_rette(a1,b1,c1,a2,b2,c2); 
    else
        xc1 = yp4/tan(beta/2);
    end
    xc4 = xc1;
    yc4 = yc1;
    if yc2 < yc1
        [a,b,c] = retta_per_2pti(xp4,yp4,xp5(1),yp5(1));
        yc1 = -(a*xc1+c)/b;
    end
end
if xc4 < xp4
    xc4 = xp4;
end
if xc1 < xp4
    xc1 = xp4;
    yc1 = yp4;
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

ry = sqrt(xp4^2+yp4^2);


% computation of coil cross-section (circa slot area)
if r_fillet==0
    X = [xc1 xc2 xc3 xc4 xc1];
    Y = [yc1 yc2 yc3 yc4 yc1];
else
    ang1 = angle(xc2(1)-xc_fillet+j*yc2(1)-j*yc_fillet);
    ang2 = angle(xc2(2)-xc_fillet+j*yc2(2)-j*yc_fillet);
    xRacc = xc_fillet+r_fillet*cos(linspace(ang1,ang2,101));
    yRacc = yc_fillet+r_fillet*sin(linspace(ang1,ang2,101));
    X = [xc1 xc2(1) xRacc xc2(2) xc3 xc4 xc1];
    Y = [yc1 yc2(1) yRacc yc2(2) yc3 yc4 yc1];
end

geo.Acoilf = polyarea(X,Y);


geo.p  = p;
geo.dalpha_pu  = round(beta_pu,3);
geo.dalpha     = round(beta,3);
geo.lyr        = round(ry-ri,2);
geo.hpb        = round(xp6-ry,2);
geo.hph        = round(xM-xp6,2);
geo.wp         = round(wp,2);
geo.wb         = round(yc3-yc2(1),2);
geo.hb         = round(hb,2);
geo.r_fillet   = round(r_fillet,2);
geo.thHead_deg = round(thHead_deg,1);
geo.pont0      = MinTol;
geo.g          = g0;

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
