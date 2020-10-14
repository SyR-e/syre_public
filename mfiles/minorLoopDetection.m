% Copyright 2020
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

function [major,minor] = minorLoopDetection(x,y,debug)

% Detect major and minor loops for a given waveform.
% Useful for iron loss calcularion according to the iGSE model
% If y is complex, the component on its main axis is considered.
% The input vectors must be column vectors

if nargin==2
    debug=0;
end

if max(abs(imag(y)))>0
    ang = angle(y);
    ang(ang<0) = ang(ang<0)+pi;
    ang = mean(ang);
    y = y*exp(-j*ang);
    y = real(y);
end

x = reshape(x,[1,numel(x)]);
y = reshape(y,[1,numel(y)]);


% tic

x = x-x(1); % remove x offset
dx = diff(x(1:2));
index = 1:1:length(x);

% y = y-mean(y);

iMin = index(y==min((y)));
iMin = iMin(1);
y = [y(iMin:end) y(1:iMin-1)]; % sorting of the waveform from the positive peak


if debug
    figure();
    figSetting(16,16);
    hax = axes('OuterPosition',[0 0 1 1]);
    set(gca,'XLim',[0 360],'XTick',0:60:360);
    set(gca,'GridAlpha',1);
    plot(hax,x,y,'-o','DisplayName','$y = f(x)$');
    legend('show','Location','northeastoutside');
end

%xTmp = [-dx x x(end)+dx];
yTmp = [y(end) y y(1)];

index = 1:1:length(x);

dS = sign(diff(sign(diff(yTmp))));

iMax = index(y==max(y));
iMin = index(y==min(y));

iMax = iMax(1);
iMin = iMin(1);

% iMax = 1;
% iMin = max(index)/2+1;

dSM = zeros(size(x));
dSM([iMax iMin]) = dS([iMax iMin]);

if debug
    stem(hax,x,dS,'-o','LineWidth',1.5,'DisplayName','$d\Delta$')
    stem(hax,x,dSM,'--x','LineWidth',1.5,'DisplayName','$d\Delta_{major}$')
end

iD = index(dS~=0);

iD = [iD index(end)];

flagMin = 0;
mm = 1;

y = [y y(1)];
dy = diff(y);
index = [index NaN];

if debug
    keyboard
end

major.y = [];

for ii = 1:length(iD)-1
    if((dS(iD(ii))==dSM(iD(ii))))
        % sono sicuramente all'inizio della zona major loop
        major.y = [major.y y(iD(ii):iD(ii+1))];
        if debug
            plot(hax,x(iD(ii):iD(ii+1)),y(iD(ii):iD(ii+1)),'go','DisplayName','major loop')
        end
        flagMinor = 0;
        dM = (major.y(end)-major.y(end-1))/abs(major.y(end)-major.y(end-1));
    else
        if ~flagMin
            % sicuramente all'inizio di un minor loop
            minor(mm).y = y(iD(ii):iD(ii+1)-1);
            if debug
                plot(hax,x(iD(ii):iD(ii+1)-1),y(iD(ii):iD(ii+1)-1),'rx','DisplayName','minor loops')
            end
            flagMin = 1;
        else
            yTmp = y;
            dyTmp = dy;
            dyTmp(1:iD(ii)) = NaN;
            iFin = index(sign(dyTmp)~=dM);
            yTmp(iFin) = NaN;
            yTmp([iMax iMin]) = y([iMax iMin]);
            yTmp(1:iD(ii))  = NaN;
            % detect when the major loop restart
            if dM>0
                iFin = index(yTmp>=major.y(end));
            else
                iFin = index(yTmp<=major.y(end));
            end
            iFin = iFin(1);
            iFin(isnan(iFin)) = max(index(1:end-1));
            
            if iFin<iD(ii+1)
                % il major loop ricomincia prima del prossimo flesso
                minor(mm).y = [minor(mm).y y(iD(ii):iFin)];
                mm = mm+1;
                flagMin = 0;
                if debug
                    plot(hax,x(iD(ii):iFin),y(iD(ii):iFin),'r+','DisplayName','minor loops')
                end
                % dopo ho un major
                major.y = [major.y y(iFin:iD(ii+1))];
                if debug
                    plot(hax,x(iFin:iD(ii+1)),y(iFin:iD(ii+1)),'go','DisplayName','major loop')
                end
                dM = (major.y(end)-major.y(end-1))/abs(major.y(end)-major.y(end-1));
            else
                % il major loop ricomincia dopo il prossimo flesso
                 minor(mm).y = [minor(mm).y y(iD(ii):iD(ii+1))];
                 if debug
                     plot(hax,x(iD(ii):iD(ii+1)),y(iD(ii):iD(ii+1)),'rs','DisplayName','minor loops')
                 end
            end
        end
    end
end



major.x = (0:1:length(major.y)-1)*dx;

if ~exist('minor','var')
    minor.x = [];
    minor.y = [];
end
if length(minor)>1
    for ii=1:length(minor)
        minor(ii).x = (0:1:length(minor(ii).y)-1)*dx;
    end
else
    if ~isempty(minor.y)
        for ii=1:length(minor)
            minor(ii).x = (0:1:length(minor(ii).y)-1)*dx;
        end
    else
        minor.x = [];
    end
end




% toc
