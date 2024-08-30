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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_flxdn_fig(geo,out,newDir,filemot)

% 
% plot_flxdn_fig(geo,out,newDir,filemot)
% 


%% setup

Bg = out.SOL.Bg;
Bt = out.SOL.Bt;
By = out.SOL.By;

Q = 6*geo.p*geo.q*geo.win.n3phase;
alpha_th = geo.wt/(geo.r+geo.g)*180/pi;
alpha_slot = 360/Q;
Qs = geo.Qs;

names{1} = '$B_{g}$ (T)';
names{2} = '$B_{t}$ (T)';
names{3} = '$B_{y}$ (T)';

filenames{1} = 'Bgap';
filenames{2} = 'Btooth';
filenames{3} = 'Byoke';

%% create figures

for ii=1:3
    hfig(ii)=figure();
    figSetting()
    hax(ii)=gca;
    xlabel('$\xi$ ($^\circ$ mech)')
    ylabel(names{ii})
    
    for jj=1:Qs
        X1=[0
            alpha_th/2
            alpha_th/2
            0];
        Y= [-0.1
            -0.1
            +0.1
            +0.1];
        X2=[alpha_slot-alpha_th/2
            alpha_slot
            alpha_slot
            alpha_slot-alpha_th/2];
        X1=X1+alpha_slot*(jj-1)*ones(length(Y),1);
        X2=X2+alpha_slot*(jj-1)*ones(length(Y),1);
        
        fill(hax(ii),X1,Y,'b','HandleVisibility','off');
        fill(hax(ii),X2,Y,'b','HandleVisibility','off');
    end
end

%% plot figures and formatting

plot(hax(1),Bg(:,1),Bg(:,2:end));
plot(hax(2),Bt(:,1),Bt(:,2:end));
plot(hax(3),By(:,1),By(:,2:end));

for ii=1:3
    hchild=get(hax(ii),'Children');
    for jj=1:length(hchild)
        set(hchild(end-jj+1),'DisplayName',['$\theta= ' num2str(out.SOL.th(jj)-geo.th0(1)) '^\circ$']);
    end
    legend(hax(ii),'show');
    set(hax(ii),'XLim',[min(Bg(:,1)) max(Bg(:,1))]);
end

%% save

for ii=1:3
    saveas(hfig(ii),[newDir filemot(1:end-4) '_' filenames{ii} '.fig']);
end