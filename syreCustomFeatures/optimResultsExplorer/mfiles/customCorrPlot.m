% Copyright 2026
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

function customCorrPlot(hfig,front,designIdx,objIdx,designBounds,feasibilityLimits,varNames)
% customCorrPlot Plot correlation grid with design bounds and feasibility limits.
%
% Inputs:
%   front               - data matrix (rows = samples, cols = variables+objectives)
%   designIdx           - indices of design variables (e.g. 1:6)
%   objIdx              - indices of objectives  (e.g. 7:10)
%   designBounds        - Mx2 array of [min max] for each design variable (M = num design)
%   feasibilityLimits   - Kx2 array [min max] for each objective (K = num objectives)
%   varNames (optional) - cell or string array of length numVars with names
%
% Example:
%   customCorrPlot(front,1:6,7:10,designBounds,feasLimits,varNames)

% Validate sizes
vars = [designIdx(:); objIdx(:)];
numVars = numel(vars);
if size(front,2) < max(vars)
    error('front does not have enough columns for given indices.');
end

if nargin < 6 || isempty(varNames)
    varNames = arrayfun(@(v) sprintf('Var %d', v), vars, 'UniformOutput', false);
else
    if numel(varNames) ~= numVars
        error('varNames must have length equal to number of variables (design+obj).');
    end
    varNames = cellstr(varNames);
end

% Create figure with minimal spacing
if isempty(hfig)
    hfig = figure();
    figSetting(14,14)
end
tiledlayout(numVars, numVars, 'Padding','none', 'TileSpacing','none');

% helper to compute data limits with small padding
padLimits = @(v) deal( ...
    (min(v(~isnan(v))) - 0.02*range(v(~isnan(v)))), ...
    (max(v(~isnan(v))) + 0.02*range(v(~isnan(v)))) );

for i = 1:numVars
    yi = vars(i);
    for j = 1:numVars
        xi = vars(j);
        ax = nexttile((i-1)*numVars + j);
        xdata = front(:, xi);
        ydata = front(:, yi);
        hold(ax, 'on');

        % Remove numeric tick marks and labels from all axes
        ax.XTick = [];
        ax.YTick = [];
        ax.XTickLabel = {};
        ax.YTickLabel = {};

        if i == j
            % Diagonal: histogram + boundary shading + objective target lines
            nbins = max(10, round(sqrt(size(front,1))));
            histogram(ax, xdata, nbins, 'EdgeColor','none','FaceColor',[0 0.8 0]);

            if ismember(xi, designIdx)
                idxDesign = find(designIdx==xi,1);
                if ~isempty(idxDesign) && ~any(isnan(designBounds(idxDesign,:)))
                    lb = designBounds(idxDesign,1); ub = designBounds(idxDesign,2);
                    yl = get(ax,'YLim');
                    if isfinite(lb) && isfinite(ub)
                        xlim(ax, [lb ub]);
                    else
                        axis(ax,'tight');
                    end
                else
                    axis(ax,'tight');
                end
            elseif ismember(xi, objIdx)
                idxObj = find(objIdx==xi,1);
                if ~isempty(idxObj) && ~any(isnan(feasibilityLimits(idxObj,:)))
                    fl = feasibilityLimits(idxObj,:);
                    yl = get(ax,'YLim');
                    if isfinite(fl(1))
                        plot(ax, [fl(1) fl(1)], yl, 'r-');
                    end
                    if isfinite(fl(2))
                        plot(ax, [fl(2) fl(2)], yl, 'r-');
                    end
                    if any(isfinite(xdata))
                        [xmin, xmax] = padLimits(xdata);
                        if xmin==xmax, xmin = xmin-1; xmax = xmax+1; end
                        xlim(ax, [xmin xmax]);
                    else
                        axis(ax,'tight');
                    end
                else
                    axis(ax,'tight');
                end
            else
                axis(ax,'tight');
            end

        else
            % Off-diagonal: scatter + boundary shading + objective target lines
            plot(ax, xdata, ydata,'.');

            % X-axis bounds (design or objective)
            if ismember(xi, designIdx)
                idxDesign = find(designIdx==xi,1);
                if ~isempty(idxDesign) && ~any(isnan(designBounds(idxDesign,:)))
                    lb = designBounds(idxDesign,1); ub = designBounds(idxDesign,2);
                    yl = get(ax,'YLim'); patch(ax, [lb lb ub ub], [yl(1) yl(2) yl(2) yl(1)], [0.95 1 0.95], 'FaceAlpha',0.22, 'EdgeColor','none');
                    if isfinite(lb) && isfinite(ub)
                        xlim(ax, [lb ub]);
                    end
                end
            end
            if ismember(xi, objIdx)
                idxObj = find(objIdx==xi,1);
                if ~isempty(idxObj) && ~any(isnan(feasibilityLimits(idxObj,:)))
                    fl = feasibilityLimits(idxObj,:);
                    yl = get(ax,'YLim');
                    if isfinite(fl(1)) && isfinite(fl(2))
                        patch(ax, [fl(1) fl(1) fl(2) fl(2)], [yl(1) yl(2) yl(2) yl(1)], [0.95 1 0.95], 'FaceAlpha',0.18, 'EdgeColor','none');
                    else
                        if isfinite(fl(2))
                            patch(ax, [min(xdata) min(xdata) fl(2) fl(2)], [yl(1) yl(2) yl(2) yl(1)], [0.95 1 0.95], 'FaceAlpha',0.18, 'EdgeColor','none');
                        end
                        if isfinite(fl(1))
                            patch(ax, [fl(1) fl(1) max(xdata) max(xdata)], [yl(1) yl(2) yl(2) yl(1)], [0.95 1 0.95], 'FaceAlpha',0.18, 'EdgeColor','none');
                        end
                    end
                    if isfinite(fl(1)), plot(ax, [fl(1) fl(1)], yl, 'r-', 'LineWidth',1.2); end
                    if isfinite(fl(2)), plot(ax, [fl(2) fl(2)], yl, 'r-', 'LineWidth',1.2); end
                    if any(isfinite(xdata))
                        [xmin, xmax] = padLimits(xdata);
                        if xmin==xmax, xmin = xmin-1; xmax = xmax+1; end
                        xlim(ax, [xmin xmax]);
                    end
                end
            end

            % Y-axis bounds (design or objective)
            if ismember(yi, designIdx)
                idxDesignY = find(designIdx==yi,1);
                if ~isempty(idxDesignY) && ~any(isnan(designBounds(idxDesignY,:)))
                    lb = designBounds(idxDesignY,1); ub = designBounds(idxDesignY,2);
                    if isfinite(lb) && isfinite(ub)
                        ylim(ax, [lb ub]);
                    end
                end
            end
            if ismember(yi, objIdx)
                idxObjY = find(objIdx==yi,1);
                if ~isempty(idxObjY) && ~any(isnan(feasibilityLimits(idxObjY,:)))
                    fl = feasibilityLimits(idxObjY,:);
                    xl = get(ax,'XLim');
                    if isfinite(fl(1)) && isfinite(fl(2))
                        patch(ax, [xl(1) xl(2) xl(2) xl(1)], [fl(1) fl(1) fl(2) fl(2)], [1 0.95 0.95], 'FaceAlpha',0.18, 'EdgeColor','none');
                    else
                        if isfinite(fl(2))
                            patch(ax, [xl(1) xl(2) xl(2) xl(1)], [ -Inf fl(2) fl(2) -Inf ], [1 0.95 0.95], 'FaceAlpha',0.18, 'EdgeColor','none');
                        end
                        if isfinite(fl(1))
                            patch(ax, [xl(1) xl(2) xl(2) xl(1)], [ fl(1) Inf Inf fl(1) ], [1 0.95 0.95], 'FaceAlpha',0.18, 'EdgeColor','none');
                        end
                    end
                    if isfinite(fl(1)), plot(ax, xl, [fl(1) fl(1)], 'r-'); end
                    if isfinite(fl(2)), plot(ax, xl, [fl(2) fl(2)], 'r-'); end
                    if any(isfinite(ydata))
                        [ymin, ymax] = padLimits(ydata);
                        if ymin==ymax, ymin = ymin-1; ymax = ymax+1; end
                        ylim(ax, [ymin ymax]);
                    end
                end
            end

            xlabel(''); ylabel('');
        end

        % Show variable names only on left-most column and bottom row (no numeric ticks)
        if j == 1
            ylabel(ax, varNames{i}, 'Interpreter','none', 'FontSize',9);
        end
        if i == numVars
            xlabel(ax, varNames{j}, 'Interpreter','none', 'FontSize',9);
        end

        title(ax,'');
        box(ax,'on');
        hold(ax,'off');
    end
end

% sgtitle('Custom Correlation Matrix');

end
