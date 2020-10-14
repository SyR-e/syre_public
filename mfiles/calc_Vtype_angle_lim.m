function [betaLim] = calc_Vtype_angle_lim(geo,alpha,r,Ar,hc)
% 
% [betaLim] = calc_Vtype_angle_lim(geo,alpha,r,Ar,hc)
% 
% Compute the minimum value for the beta angle for Vtype geometry.
% NB: the hc value MUST be correct BEFORE this function
% SF - 20181012

% Load data from geo
if nargin==1
    tmp     = cumsum(geo.dalpha);
    alpha   = tmp(end);     % angle alpha
    r       = geo.r;        % rotor radius
    Ar      = geo.Ar;       % shaft radius
    hc      = geo.hc;       % barrier width
end
pontT   = geo.pontT;    % minimum mechanical tolerance and tangential ribs thickness
hfe_min = geo.hfe_min;  % minimum iron thickness


% min value computation
xc = (r-pontT-hc/2).*cos(alpha*pi/180);
yc = (r-pontT-hc/2).*sin(alpha*pi/180);
R  = hc/2;

x1 = Ar+hfe_min;
y1 = 0;

% the minimum angle is the slope of the line tangent to the circle at the
% barrier end and that cross the inner possible point (B1k limit).
% The circumference has center in (xc,yc) and radius hc/2, while the
% point is defined as (x1,y1)=(B1kLim,0)

a = (x1-xc).^2-R.^2;
b = -2*(x1-xc).*(y1-yc);
c = (y1-yc).^2-R.^2;

% m = roots([a b c]);
% m = real(m);

m1 = (-b+(b.^2-4*a.*c).^0.5)./(2*a);
m2 = (-b-(b.^2-4*a.*c).^0.5)./(2*a);

beta1 = atan(m1);
beta2 = atan(m2);

betaLim = beta1;
betaLim(beta2>beta1)=beta2(beta2>beta1);













