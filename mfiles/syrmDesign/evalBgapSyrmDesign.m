% Copyright 2016
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

function [x,y,kf1,kfm] = evalBgapSyrmDesign(q,kt,debug)
% 
% [x,y,kf1,kfm] = evalBgapSyrmDesign(q,kt)
% 
% Evaluate the airgap induction with stair case.
% Main function take from staircaseRegular.m

if nargin==2
    debug=0;
end

x = linspace(-pi,pi,1001);

if debug
    figure()
    hold on
    grid on
    box on
    xlabel('\theta_{elt} [°]')
    ylabel('B_{gap} [pu]')
    set(gca,'XLim',[-180 180],'XTick',-180:30:180);
    set(gca,'YLim',[-1.6 1.6],'YTick',-1.6:0.2:1.6);
    set(gca,'PlotBoxAspectRatio',[1 0.8 1]);
    plot([-180 180],[0 0],'-k','LineWidth',0.5,'HandleVisibility','off')
    plot([0 0],[-2 2],'-k','LineWidth',0.5,'HandleVisibility','off')
    plot(x*180/pi,1/kt*cos(x),'-b','LineWidth',2,'DisplayName','generator: y=1/kt*cos(x)')
end


%% From staircaseRegular.m
nlay = floor(6*q/4);    % regular stairs only

xmax = pi; ph = 0;
xtmp = linspace(0,xmax,1000);
Fs = sin(xtmp - ph);

delta = 360/(6*q);
angles = (delta/2:delta:90) * pi/180;
angles = [0 angles];
if angles(end) < pi/2
    angles = [angles pi/2];
end
angles = [angles(1:nlay+1) angles(end)];
levels = (1/kt*cos(angles(1:end-1)) - 1/kt*cos(angles(2:end)))./(angles(2:end) - angles(1:end-1));
levels(1) = 0;
temp = fliplr(levels);
temp = [levels temp(2:end)];
levels = temp;
angles = angles(2:end);
angles_all = [angles(1:end-1) pi-(fliplr(angles(1:end-1)))];

ytmp = zeros(1,length(xtmp));
for jj = 2:length(angles_all)
    ytmp((xtmp < angles_all(jj)) & (xtmp > angles_all(jj-1))) = levels(jj);
end

f = levels(1:nlay+1);
df = diff(f);
da = diff(angles);
%%

% y=zeros(1,length(x));
xtmp=xtmp-pi/2;
y=interp1(xtmp,ytmp,x);
y(x>+pi/2)=-interp1(xtmp,ytmp,fliplr(abs(x(x>+pi/2)-pi/2)));
y(x<-pi/2)=-interp1(xtmp,ytmp,fliplr(abs(x(x<-pi/2)+pi/2)));

if debug
    plot(x*180/pi,y,'-.g','LineWidth',1.5,'DisplayName','unsaturated staircase')
end

% saturation
y(y>+1) = +1;
y(y<-1) = -1;

% mean and first harmonic factor
kfm=mean(abs(y));
kf1=1/pi*trapz(x,y.*cos(x));

if (debug)
    plot(x*180/pi,y,'-g','LineWidth',2,'DisplayName','saturated staircase')
    plot([-180 180],+[1 1],'--k','LineWidth',1,'HandleVisibility','off')
    plot([-180 180],-[1 1],'--k','LineWidth',1,'HandleVisibility','off')
    plot(x*180/pi,kfm*ones(1,length(x)),'--r','LineWidth',1,'DisplayName','mean: y=b*kfm')
    plot(x*180/pi,kf1*cos(x),'-r','LineWidth',1,'DisplayName','1^{st} harmonics: y=b*kf1*cos(x)')
    legend('show','Location','South')
end
