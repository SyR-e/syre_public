
%   staircaseFunction.m
%   Function for MMF staircase evaluation
% 
%   called twice, both for rotor and stator MMF discretization dictated by
%   respective slots
% 
%   This script was used to produce the figures presented at the
%   "PM Machine Design Boot Camp", UW Madison, November 11 -14, 2014
% 
%   Copyright 2014 Gianmario Pellegrino - gianmario.pellegrino@polito.it


function [df,da] = staircaseAnyAlpha(alphaElt)

% nlay = floor(nr/4);    % regular stairs only

xmax = pi; ph = 0;
x = linspace(0,xmax,1000);
Fs = sin(x - ph);

nlay = length(alphaElt);

angles = fliplr(alphaElt);
angles = pi/2 - angles;
angles = [0 angles];
if angles(end) < pi/2
    angles = [angles pi/2];
end
% angles = [angles(1:nlay+1) angles(end)];
levels = (cos(angles(1:end-1)) - cos(angles(2:end)))./(angles(2:end) - angles(1:end-1));
levels(1) = 0;
temp = fliplr(levels);
temp = [levels temp(2:end)];
levels = temp;
angles = angles(2:end);
angles_all = [angles(1:end-1) pi-(fliplr(angles(1:end-1)))];

y = zeros(1,length(x));
for jj = 2:length(angles_all)
    y((x < angles_all(jj)) & (x > angles_all(jj-1))) = levels(jj);
end

if (0)
    figure(1)
    plot(x*180/pi,abs(Fs),'LineWidth',2), grid on, axis([ph 180 0 1]);
    hold on
    plot(x*180/pi,abs(y),'k','LineWidth',2)
    set(gca,'PlotBoxAspectRatio',[1 0.5 1])
    xlabel('rotor angular coordinate - elt degrees')
    ylabel('MMF - Am')
end

f = levels(1:nlay+1);
df = diff(f);
da = diff(angles);