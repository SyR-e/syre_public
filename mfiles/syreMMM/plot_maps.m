
% Copyright 2018
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

function plot_maps(pathname,filename)

% plot_maps

% input file contains:
% - Id Iq Fd Fq
% - optionally: loss maps

% output:
% - mesh of all maps in the input file (Id Iq Fd Fq and maybe loss maps)
% - mesh of T, amp(F), amp(I), IPF
% - saves all figures


if (nargin<1)
    % pathname = cd;
    % filename = 'fdfq_idiq_n256.mat'
    [FILENAME, pathname, FILTERINDEX] = uigetfile(['.mat'], 'LOAD DATA');
    load([pathname FILENAME]);
else
    load([pathname '\' filename])
end

% Fd Fq curves
figure, figSetting
[value, index] = min(abs(Iq(:,1)));
plot(Id(index,:),Fd(index,:)),
plot(Id(1,:),Fd(1,:)),
plot(Id(end,:),Fd(end,:)), hold on
[value, index] = min(abs(Id(1,:)));
plot(Iq(:,index),Fq(:,index)), hold all
plot(Iq(:,1),Fq(:,1)),
plot(Iq(:,end),Fq(:,end)),
xlabel('[A]'), ylabel('[Vs]')
h=gcf(); %AS
% if isoctave()
%     fig_name=strcat(pathname, 'fdfq curves');
%     hgsave(h,[fig_name]);
% else
%     saveas(h,[pathname 'fdfq curves'])
% end

% Fd [Vs]
figure, figSetting
mesh(Id,Iq,Fd), grid on, xlabel('$$ i_d [A] $$'), ylabel('$$ i_q [A] $$'), zlabel('$$ \lambda_d [Vs] $$')
view(3)
plot3(Id(1,:),Iq(1,:),Fd(1,:),'k','LineWidth',2),
plot3(Id(end,:),Iq(end,:),Fd(end,:),'k--','LineWidth',2),
[value, index] = min(abs(Iq(:,1)));
plot3(Id(index,:),Iq(index,:),Fd(index,:),'k','LineWidth',2),
h=gcf(); %AS
% if isoctave()
%     fig_name=strcat(pathname, 'Fd mesh');
%     hgsave(h,[fig_name]);
% else
%     saveas(h,[pathname 'Fd mesh'])
% end

% Fq [Vs]
figure, figSetting  
mesh(Id,Iq,Fq), grid on, xlabel('$$ i_d [A] $$'), ylabel('$$ i_q [A] $$'), zlabel('$$ \lambda_q [Vs] $$')
view(3)
plot3(Id(:,1),Iq(:,1),Fq(:,1),'k','LineWidth',2),
plot3(Id(:,end),Iq(:,end),Fq(:,end),'k--','LineWidth',2),
[value, index] = min(abs(Id(1,:)));
plot3(Id(:,index),Iq(:,index),Fq(:,index),'k','LineWidth',2),
h=gcf(); %AS
% if isoctave()
%     fig_name=strcat(pathname, 'Fq mesh');
%     hgsave(h,[fig_name]);
% else
%     saveas(h,[pathname 'Fq mesh'])
% end

% TORQUE
if exist('T','var')
    figure, figSetting 
    mesh(Id,Iq,abs(T)), grid on, xlabel('$$ i_d [A] $$'), ylabel('$$ i_q [A] $$'), zlabel('Torque [Nm]')
    view(3)
    h=gcf(); %AS
%     if isoctave()
%         fig_name=strcat(pathname, 'Torque mesh');
%         hgsave(h,[fig_name]);
%     else
%         saveas(h,[pathname 'Torque mesh'])
%     end
    
end
if exist('dT','var')
    figure, figSetting
    mesh(Id,Iq,dT), grid on, xlabel('$$ i_d [A] $$'), ylabel('$$ i_q [A] $$'), zlabel('Torque ripple (std) [Nm]')
    view(3)
    h=gcf(); %AS
%     if isoctave()
%         fig_name=strcat(pathname, 'Torque Ripple mesh');
%         hgsave(h,[fig_name]);
%     else
%         saveas(h,[pathname 'Torque Ripple mesh'])
%     end
    
end

if exist('dTpp','var')
    figure, figSetting
    mesh(Id,Iq,dTpp), grid on, xlabel('$$ i_d [A] $$'), ylabel('$$ i_q [A] $$'), zlabel('Torque ripple (peak-to-peak) [Nm]')
    view(3)
    h=gcf(); %AS
%     if isoctave()
%         fig_name=strcat(pathname, 'Torque Ripple pp mesh');
%         hgsave(h,[fig_name]);
%     else
%         saveas(h,[pathname 'Torque Ripple pp mesh'])
%     end
    
end

% core loss
if exist('Pfes_h','var')
    % Loss Map
    figure, figSetting
    subplot(2,2,1)
    view(3)
    mesh(Id,Iq,Pfes_c); grid on, hold on
    xlabel('$$ i_d [A] $$'), ylabel('$$ i_q [A] $$'),zlabel('Pc-stat [W]')
    subplot(2,2,3)
    view(3)
    mesh(Id,Iq,Pfes_h); grid on, hold on
    xlabel('$$ i_d [A] $$'), ylabel('$$ i_q [A] $$'),zlabel('Ph-stat [W]'),
    subplot(2,2,2)
    view(3)
    mesh(Id,Iq,Pfer_c); grid on, hold on
    xlabel('$$ i_d [A] $$'), ylabel('$$ i_q [A] $$'),zlabel('Pc-rot [W]')
    subplot(2,2,4)
    view(3)
    mesh(Id,Iq,Pfer_h); grid on, hold on
    xlabel('$$ i_d [A] $$'), ylabel('$$ i_q [A] $$'),zlabel('Ph-rot [W]')
    h=gcf(); %AS
%     if isoctave()
%         fig_name=strcat(pathname, 'Loss mesh');
%         hgsave(h,[fig_name]);
%     else
%         saveas(h,[pathname 'Loss mesh'])
%     end
    
end

% pm loss
if exist('Ppm','var')
    figure, figSetting  
    mesh(Id,Iq,Ppm), grid on, xlabel('i_d [A]'), ylabel('i_q [A]'), zlabel('pm loss [W]')
    view(3)
    h=gcf(); %AS
%     if isoctave()
%         fig_name=strcat(pathname, 'PM loss mesh');
%         hgsave(h,[fig_name]);
%     else
%         saveas(h,[pathname 'PM loss mesh'])
%     end
    clear  fig_name
end





