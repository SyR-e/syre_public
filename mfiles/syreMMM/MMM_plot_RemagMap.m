% Copyright 2025
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

function MMM_plot_RemagMap(motorModel)

% load data
Remag = motorModel.VFMdata.data.Remag;
MS = Remag.MS;                                        ;
Id = Remag.id;
Iq = Remag.iq;

% Check for MS_ph field and convert to degrees if present
if isfield(Remag, 'MS_ph')
    MS_ph = Remag.MS_ph;
    MS_ph_deg = rad2deg(MS_ph); % Convert radians to degrees
else
    MS_ph = [];
    MS_ph_deg = [];
end

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = checkPathSyntax([motName '_results\MMM results\' 'Remag Map\']);

% Prepare figure
hfig(1) = figure();
figSetting()
hax(1,1) = axes('OuterPosition',[0 0 1 1]);

% Create contour plot of MS values
[C, h] = contourf(hax(1,1), Id, Iq, MS, 'LineStyle', '-', 'LineWidth', 1, 'DisplayName', 'MS Contour');
colorbar(hax(1,1));

% Add MS_ph iso-lines in red if available
if ~isempty(MS_ph_deg)
    hold(hax(1,1), 'on');
    
    % Create contour lines for MS_ph values
    [C_ph, h_ph] = contour(hax(1,1), Id, Iq, MS_ph_deg, 'LineColor', 'red', 'LineWidth', 2, ...
                          'LineStyle', '-', 'DisplayName', 'MS Phase (deg)');
    
    
    clabel(C_ph, h_ph, 'FontSize', 8, 'Color', 'red', 'LabelSpacing', 300, ...
           'BackgroundColor', 'white', 'EdgeColor', 'red');
    
    hold(hax(1,1), 'off');
end

xlabel(hax(1,1), '$I_d$ (A)')
ylabel(hax(1,1), '$I_q$ (A)')
title(hax(1,1), 'MS Remag')
if ~isempty(MS_ph_deg)
    title(hax(1,1), 'MS Amplitude and Phase Remag')
else
    title(hax(1,1), 'MS Amplitude Remag')
end

set(hfig(1), 'FileName', [pathname resFolder 'RemagMap.fig']);
set(hfig(1), 'Name', 'RemagMap');

grid(hax(1,1), 'on');
grid(hax(1,1), 'minor');

legend(hax(1,1), 'show', 'Location', 'best');

set(hax(1,1), 'XLim', [min(Id(:)) max(Id(:))]);
set(hax(1,1), 'YLim', [min(Iq(:)) max(Iq(:))]);

%% Save figure
answer = 'No';
answer = questdlg('Save figure?', 'Save', 'Yes', 'No', answer);
if strcmp(answer, 'Yes')
    if ~exist([pathname resFolder], 'dir')
        mkdir([pathname resFolder]);
    end
    savePrintFigure(hfig(1));
end

end