% Copyright 2025
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

function f = surfplot_slider(X,Y,Z,Ir,figureFullName,labels) % X,Y: [N x N]; Z: [N x N x M]; Ir [L x 1] L-> Number of shown layers
    
    % --- Input Handling for Optional Labels ---
    if nargin < 6 || isempty(labels)
        % If 'labels' is not provided or is empty, use default labels
        xlabel_str = '$i_d$ (A)';
        ylabel_str = '$i_q$ (A)';
        zlabel_str = '$T$ (Nm)';
    else
        % If 'labels' is provided, use them
        % Ensure 'labels' is a cell array with 3 elements (for x, y, z)
        if ~iscell(labels) || numel(labels) ~= 3
            warning('surfplot_slider:InvalidLabels', 'Labels input must be a 1x3 cell array. Using default labels.');
            xlabel_str = '$i_d$ (A)';
            ylabel_str = '$i_q$ (A)';
            zlabel_str = '$T$ (Nm)';
        else
            xlabel_str = labels{1};
            ylabel_str = labels{2};
            zlabel_str = labels{3};
        end
    end

    % --- Figure and Axes Setup ---
    f = figure;
    figSetting(); % Assumed to be defined elsewhere

    ax = axes('OuterPosition',[0 0 1 1],...
         'XLim',[min(X(:)) max(X(:))],...
         'YLim',[min(Y(:)) max(Y(:))],...
         'PlotBoxAspectRatio',[1 1 0.8]);
    xlabel(xlabel_str, 'Interpreter', 'latex')
    ylabel(ylabel_str, 'Interpreter', 'latex')
    zlabel(zlabel_str, 'Interpreter', 'latex')
    set(f,'FileName',figureFullName);
    [~,figureName] = fileparts(figureFullName);
    set(f,'Name',figureName)
    view(3)
    hold(ax, 'on'); % Keep hold on for subsequent plots

    % --- Initialize GUI Data Structure (for UserData) ---
    % This struct holds all shared variables and UI handles.
    guiData.X = X;
    guiData.Y = Y;
    guiData.Z = Z;
    guiData.Ir = Ir;
    guiData.current_Ir_index = numel(Ir); % Start at max Ir value (last index)

    % Preallocate handles for surfaces and contours
    num_ir_values = numel(Ir);
    guiData.hSurf = gobjects(num_ir_values, 1); % Use gobjects for graphics handles
    guiData.hCont = gobjects(num_ir_values, 1);

    % Create all surfaces and contours, initially hidden
    for ii = 1:num_ir_values
        guiData.hSurf(ii) = surf(ax, X, Y, Z(:,:,ii), 'FaceColor', 'interp', 'EdgeColor', 'none', 'Visible', 'off');
        [~,guiData.hCont(ii)] = contour3(ax, X, Y, Z(:,:,ii),'EdgeColor','k','ShowText','off', 'Visible', 'off');
    end

    guiData.hIrTextBox = uicontrol('Parent', f, 'Style', 'edit', ...
                                   'String', num2str(round(Ir(guiData.current_Ir_index),2)), ...
                                   'Position', [55 45 50 20], ...
                                   'HorizontalAlignment', 'center', ...
                                   'FontSize', 14, 'BackgroundColor', 'white', ...
                                   'Visible', 'off'); 

    guiData.hButtonUp = uicontrol('Parent', f, 'Style', 'pushbutton', 'String', '▲', ...
                                  'Position', [20 55 30 20], ...
                                  'FontSize', 16, 'FontWeight', 'bold');
    
    guiData.hButtonDown = uicontrol('Parent', f, 'Style', 'pushbutton', 'String', '▼', ...
                                    'Position', [20 20 30 20], ...
                                    'FontSize', 16, 'FontWeight', 'bold');

    set(f, 'UserData', guiData); % Store guiData in figure's UserData

    % --- Assign Callbacks ---
    % Callbacks only receive the figure handle 'f' and the direction/src
    set(guiData.hIrTextBox, 'Callback', @(src, event) updateValueFromTextBox(f)); % Pass figure handle
    set(guiData.hButtonUp, 'Callback', @(src, event) incrementDecrementValue(1, f));
    set(guiData.hButtonDown, 'Callback', @(src, event) incrementDecrementValue(-1, f));

    % --- Nested Callback Functions ---

    function incrementDecrementValue(direction, fHandle)
        % Retrieve current GUI state
        current_guiData = get(fHandle, 'UserData'); 
        
        new_index = current_guiData.current_Ir_index + direction;
        
        % Clamp index to valid range
        if new_index < 1
            new_index = 1;
            warning('MatlabPlot:MinIrReached', 'Minimum Ir reached.');
        elseif new_index > numel(current_guiData.Ir) % Use numel for vector length
            new_index = numel(current_guiData.Ir);
            warning('MatlabPlot:MaxIrReached', 'Maximum Ir reached.');
        end
        
        current_guiData.current_Ir_index = new_index; % Update index
        
        % Update text box string
        set(current_guiData.hIrTextBox, 'String', num2str(round(current_guiData.Ir(new_index),2))); 
        
        % Save updated GUI state
        set(fHandle, 'UserData', current_guiData); 
        
        % Refresh plot visibility
        updateSurfaceVisibility(fHandle);
    end

    function updateValueFromTextBox(fHandle)
        % Retrieve current GUI state
        current_guiData = get(fHandle, 'UserData'); 
        
        inputString = get(current_guiData.hIrTextBox, 'String'); % Get string directly from guiData handle
        inputVal = str2double(inputString);
        
        if ~isnan(inputVal) && isreal(inputVal)
            [~, closest_index] = min(abs(current_guiData.Ir - inputVal)); 
            
            % Validate input against predefined Ir values
            if ~isempty(closest_index) && abs(current_guiData.Ir(closest_index) - inputVal) < 1e-3
                current_guiData.current_Ir_index = closest_index; 
                set(current_guiData.hIrTextBox, 'String', num2str(round(current_guiData.Ir(closest_index),2))); 
                
                % Save updated GUI state
                set(fHandle, 'UserData', current_guiData); 
                
                updateSurfaceVisibility(fHandle);
            else
                % Revert to previous valid value if input is not close to any Ir
                set(current_guiData.hIrTextBox, 'String', num2str(round(current_guiData.Ir(current_guiData.current_Ir_index),2)));
                warning('MatlabPlot:InvalidInput', 'Ir value must match one of the predefined values.');
            end
        else
            % Revert if input is not a number
            set(current_guiData.hIrTextBox, 'String', num2str(round(current_guiData.Ir(current_guiData.current_Ir_index),2)));
            warning('MatlabPlot:InvalidInput', 'Invalid input. Please enter a number.');
        end
    end

    function updateSurfaceVisibility(fHandle)
        % Retrieve current GUI state
        current_guiData = get(fHandle, 'UserData'); 
        
        irValue = current_guiData.Ir(current_guiData.current_Ir_index);
        title(['Map @ Ir = ', num2str(round(irValue,2)), ' A']);
        
        % Update visibility of all surfaces and contours
        for jj = 1:numel(current_guiData.hSurf) % Use numel for number of elements
            if jj == current_guiData.current_Ir_index
                set(current_guiData.hSurf(jj), 'Visible', 'on');
                set(current_guiData.hCont(jj), 'Visible', 'on');
            else
                set(current_guiData.hSurf(jj), 'Visible', 'off');
                set(current_guiData.hCont(jj), 'Visible', 'off');
            end
        end
        
        % Enable/disable buttons based on current index
        if current_guiData.current_Ir_index == 1
            set(current_guiData.hButtonDown, 'Enable', 'off'); 
        else
            set(current_guiData.hButtonDown, 'Enable', 'on');  
        end
        
        if current_guiData.current_Ir_index == numel(current_guiData.Ir)
            set(current_guiData.hButtonUp, 'Enable', 'off');   
        else
            set(current_guiData.hButtonUp, 'Enable', 'on');    
        end
    end

    % Initial call to set correct visibility and UI state
    updateSurfaceVisibility(f); 
end