function f = slider_inductance(X, Y, Z, MS_values, figureFullName, labels, plotType)
    % Interactive 3D plot with magnetization state (MS) slider control
    % Supports both incremental and apparent inductance maps
    
    % Set default plot type if not provided
    if nargin < 7
        plotType = 'inductance';
    end
    
    % Set default axis labels if not provided
    if nargin < 6 || isempty(labels)
        xlabel_str = '$i_d$ (A)';
        ylabel_str = '$i_q$ (A)';
        zlabel_str = 'Value';
    else
        xlabel_str = labels{1};
        ylabel_str = labels{2};
        zlabel_str = labels{3};
    end

    % Create figure with standard settings
    f = figure();
    figSetting();

    % Set up main axes with consistent 3D view
    hax = axes('OuterPosition',[0 0 0.85 1],...  % Leave space for colorbar
        'XLim',[min(X,[],'all') max(X,[],'all')],...
        'YLim',[min(Y,[],'all') max(Y,[],'all')],...
        'PlotBoxAspectRatio',[1 1 0.8]);
    xlabel(xlabel_str)
    ylabel(ylabel_str)
    zlabel(zlabel_str)
    view(3)
    
    % Store figure name for saving
    [~, figureName] = fileparts(figureFullName);
    set(f, 'FileName', figureFullName, 'Name', figureName);

    % Initialize data structure for GUI state management
    guiData.X = X;
    guiData.Y = Y;
    guiData.Z = Z;
    guiData.MS_values = MS_values;
    guiData.current_MS_index = length(MS_values);  % Start at maximum MS
    guiData.hax = hax;
    guiData.plotType = plotType;

    % Create appropriate visualization based on plot type
    if strcmp(plotType, 'app') && contains(figureName, 'FluxM')
        % Special handling for flux maps - show edges for better visibility
        guiData.hSurf = surf(X(:,:,guiData.current_MS_index), Y(:,:,guiData.current_MS_index), Z(:,:,guiData.current_MS_index), 'FaceColor', 'interp', 'EdgeColor', 'k', 'LineWidth', 1.5);
        guiData.hContour = []; % No contour lines for flux maps
    else
        % Standard inductance maps - surface with contour lines
        guiData.hSurf = surf(X(:,:,guiData.current_MS_index), Y(:,:,guiData.current_MS_index), Z(:,:,guiData.current_MS_index), 'FaceColor', 'interp', 'EdgeColor', 'none');
        hold on;
        [~, guiData.hContour] = contour3(X(:,:,guiData.current_MS_index), Y(:,:,guiData.current_MS_index), Z(:,:,guiData.current_MS_index), 'EdgeColor', 'k', 'ShowText', 'off');
        hold off;
    end

    % Add colorbar with consistent scaling across all MS values - FIXED CLIM
    guiData.hColorbar = colorbar('eastoutside');
    ylabel(guiData.hColorbar, zlabel_str, 'Interpreter', 'latex');
    clim_min = min(Z(:));
    clim_max = max(Z(:));
    
    % Ensure valid color limits (min < max)
    if clim_min == clim_max
        if clim_min == 0
            clim_lim = [-0.1 0.1]; % Default range for zero data
        else
            clim_lim = [clim_min-0.1*abs(clim_min) clim_max+0.1*abs(clim_max)];
        end
    else
        clim_lim = [clim_min clim_max];
    end
    set(hax, 'CLim', clim_lim);

    % Compact UI controls positioned at bottom of figure
    controlY = 10;
    controlHeight = 20;
    controlWidth = 25;

    % Navigation buttons and MS value display
    guiData.hButtonPrev = uicontrol('Parent', f, 'Style', 'pushbutton', 'String', '◀', ...
                                'Position', [50, controlY, controlWidth, controlHeight], ...
                                'FontSize', 10, 'FontWeight', 'bold');

    guiData.hMSTextBox = uicontrol('Parent', f, 'Style', 'edit', ...
                                   'String', num2str(MS_values(guiData.current_MS_index)), ...
                                   'Position', [80, controlY, 50, controlHeight], ...
                                   'HorizontalAlignment', 'center', ...
                                   'FontSize', 10, 'BackgroundColor', 'white');
    
    uicontrol('Parent', f, 'Style', 'text', 'String', 'MS:', ...
              'Position', [135, controlY, 25, controlHeight], 'FontSize', 9, ...
              'HorizontalAlignment', 'left', 'BackgroundColor', get(f, 'Color'));
    
    % Main slider for MS navigation
    guiData.hSlider = uicontrol('Parent', f, 'Style', 'slider', ...
                                'Min', 1, 'Max', length(MS_values), 'Value', guiData.current_MS_index, ...
                                'Position', [165, controlY, 200, controlHeight], ...
                                'SliderStep', [1/(length(MS_values)-1) 1/(length(MS_values)-1)]);
    
    guiData.hButtonNext = uicontrol('Parent', f, 'Style', 'pushbutton', 'String', '▶', ...
                                    'Position', [370, controlY, controlWidth, controlHeight], ...
                                    'FontSize', 10, 'FontWeight', 'bold');
    
    % Display MS value range for reference
    uicontrol('Parent', f, 'Style', 'text', 'String', sprintf('[%.2f-%.2f]', min(MS_values), max(MS_values)), ...
              'Position', [400, controlY, 80, controlHeight], 'FontSize', 8, ...
              'HorizontalAlignment', 'left', 'BackgroundColor', get(f, 'Color'));
    
    set(f, 'UserData', guiData);

    % --- Nested callback functions (must be defined before callbacks are set) ---
    
    function changeMS(direction, fHandle)
        % Change MS index by specified direction (-1 or +1)
        guiData = get(fHandle, 'UserData');
        new_index = guiData.current_MS_index + direction;
        
        % Clamp index to valid range
        if new_index < 1
            new_index = 1;
        elseif new_index > length(guiData.MS_values)
            new_index = length(guiData.MS_values);
        end
        
        guiData.current_MS_index = new_index;
        set(fHandle, 'UserData', guiData);
        set(guiData.hSlider, 'Value', new_index);
        updateMSDisplay(fHandle);
    end

    function updateMSFromSlider(fHandle)
        % Update display when slider is moved
        guiData = get(fHandle, 'UserData');
        new_index = round(get(guiData.hSlider, 'Value'));
        
        % Ensure index is within bounds
        if new_index < 1, new_index = 1; end
        if new_index > length(guiData.MS_values), new_index = length(guiData.MS_values); end
        
        guiData.current_MS_index = new_index;
        set(fHandle, 'UserData', guiData);
        updateMSDisplay(fHandle);
    end

    function updateMSFromTextBox(fHandle)
        % Update display when MS value is entered in text box
        guiData = get(fHandle, 'UserData');
        inputVal = str2double(get(guiData.hMSTextBox, 'String'));
        
        if ~isnan(inputVal)
            % Find closest matching MS value
            [~, closest_index] = min(abs(guiData.MS_values - inputVal));
            
            if abs(guiData.MS_values(closest_index) - inputVal) < 1e-3
                guiData.current_MS_index = closest_index;
                set(fHandle, 'UserData', guiData);
                set(guiData.hSlider, 'Value', closest_index);
                updateMSDisplay(fHandle);
            else
                % Revert to current value if input doesn't match predefined values
                set(guiData.hMSTextBox, 'String', num2str(guiData.MS_values(guiData.current_MS_index)));
                warndlg('MS value must match one of the predefined values.', 'Invalid MS');
            end
        else
            % Revert if input is not a valid number
            set(guiData.hMSTextBox, 'String', num2str(guiData.MS_values(guiData.current_MS_index)));
            warndlg('Please enter a valid number.', 'Invalid Input');
        end
    end

    function updateMSDisplay(fHandle)
        % Update all visual elements for current MS value
        guiData = get(fHandle, 'UserData');
        current_index = guiData.current_MS_index;
        current_MS = guiData.MS_values(current_index);
        
        % Update z-axis limits with safe range calculation - FIXED HERE
        current_Z = guiData.Z(:,:,current_index);
        z_min = min(current_Z(:));
        z_max = max(current_Z(:));
        
        % Ensure valid z-limits (min < max)
        if z_min == z_max
            z_lim = [z_min-0.001 z_max+0.001];
        else
            z_lim = [z_min z_max];
        end
        
        set(guiData.hax, 'ZLim', z_lim);
        
        % Update surface data
        set(guiData.hSurf, 'ZData', guiData.Z(:,:,current_index));
        set(guiData.hSurf, 'CData', guiData.Z(:,:,current_index));
        
        % Update contour lines (except for flux maps)
        if ~isempty(guiData.hContour) && isvalid(guiData.hContour)
            delete(guiData.hContour);
        end
        
        if ~strcmp(guiData.plotType, 'app') || ~contains(get(fHandle, 'Name'), 'FluxM')
            hold(guiData.hax, 'on');
            [~, guiData.hContour] = contour3(guiData.X(:,:,current_index), guiData.Y(:,:,current_index), guiData.Z(:,:,current_index), 'EdgeColor', 'k', 'ShowText', 'off');
            hold(guiData.hax, 'off');
        end
        
        % Update title and UI elements
        title(sprintf('%s - MS = %.3f', get(fHandle, 'Name'), current_MS));
        set(guiData.hMSTextBox, 'String', num2str(current_MS));
        set(guiData.hButtonPrev, 'Enable', ifelse(current_index > 1, 'on', 'off'));
        set(guiData.hButtonNext, 'Enable', ifelse(current_index < length(guiData.MS_values), 'on', 'off'));
        
        set(fHandle, 'UserData', guiData);
        drawnow;
    end

    % --- Set callbacks after nested functions are defined ---
    set(guiData.hMSTextBox, 'Callback', @(src, event) updateMSFromTextBox(f));
    set(guiData.hSlider, 'Callback', @(src, event) updateMSFromSlider(f));
    set(guiData.hButtonPrev, 'Callback', @(src, event) changeMS(-1, f));
    set(guiData.hButtonNext, 'Callback', @(src, event) changeMS(1, f));
    
    % Initialize display with current MS value
    updateMSDisplay(f);
end

% Helper function for conditional assignment
function out = ifelse(condition, true_val, false_val)
    if condition
        out = true_val;
    else
        out = false_val;
    end
end