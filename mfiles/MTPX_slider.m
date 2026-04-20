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

function f = MTPX_slider(X,Y,Z1,Z2,Ir,MTPA,MTPV,figureFullName,labels) % X,Y: [N x N]; Z: [N x N x M]; Ir [L x 1]; ... ; label 0 for current, 1 for flux
    
    % --- Input Handling for Optional Labels ---
    if nargin < 9 || isempty(labels) || labels == 0
        % If 'labels' is not provided or is empty, use default labels
        labels = 0;
        xlabel_str = '$i_d$ (A)';
        ylabel_str = '$i_q$ (A)';
    else
        xlabel_str = '$\lambda_d$ (Vs)';
        ylabel_str = '$\lambda_q$ (Vs)';
    end
% Plot in function of currents/fluxes instead of more generic x,y,z -> lables now is a flag!!!
    %     % If 'labels' is provided, use them
    %     % Ensure 'labels' is a cell array with 3 elements (for x, y, z)
    %     if ~iscell(labels) || numel(labels) ~= 3
    %         warning('surfplot_slider:InvalidLabels', 'Labels input must be a 1x3 cell array. Using default labels.');
    %         xlabel_str = '$i_d$ (A)';
    %         ylabel_str = '$i_q$ (A)';
    %     else
    %         xlabel_str = labels{1};
    %         ylabel_str = labels{2};
    %     end
    % end

    % --- Figure and Axes Setup ---
    f = figure;
    figSetting(); % Assumed to be defined elsewhere

    ax = axes( ...
            'XLim',[min(X(:),[],'all') max(X(:),[],'all')],...
            'YLim',[min(Y(:),[],'all') max(Y(:),[],'all')],...
            'DataAspectRatio',[1 1 1]);
    if (~sum(isnan(Z2(:)))&&(max(Z2(:))~=min(Z2(:))))
        set(gca,'Clim',[min(Z2(:),[],'all') max(Z2(:),[],'all')]);
    else
        set(gca,'Clim',[0 1]);
    end
    xlabel(xlabel_str, 'Interpreter', 'latex')
    ylabel(ylabel_str, 'Interpreter', 'latex')
    iLevels = get(gca,'XTick');
    set(ax,'XTick',iLevels,'YTick',iLevels);
    set(f,'FileName',figureFullName);
    [~,figureName] = fileparts(figureFullName);
    set(f,'Name',figureName)
    hold(ax, 'on'); % Keep hold on for subsequent plots

    % --- Initialize GUI Data Structure (for UserData) ---
    % This struct holds all shared variables and UI handles.
    guiData.X = X;
    guiData.Y = Y;
    guiData.Z1 = Z1;
    guiData.Z2 = Z2;
    guiData.Ir = Ir;
    guiData.current_Ir_index = numel(Ir); % Start at max Ir value (last index)

    % Preallocate handles for surfaces and contours
    num_Ir_values = numel(Ir);
    guiData.hCont_dT = gobjects(num_Ir_values, 1); % Use gobjects for graphics handles
    guiData.hCont_T = gobjects(num_Ir_values, 1);
    guiData.hCont_I = gobjects(num_Ir_values, 1);
    guiData.hPlot_MTPA = gobjects(num_Ir_values, 1);
    guiData.hPlot_MTPV = gobjects(num_Ir_values, 1);
    
    MTPA_id = interp1(MTPA.Ir_layers,MTPA.id,Ir(:));
    MTPA_iq = interp1(MTPA.Ir_layers,MTPA.iq,Ir(:));
    MTPV_id = interp1(MTPV.Ir_layers,MTPV.id,Ir(:));
    MTPV_iq = interp1(MTPV.Ir_layers,MTPV.iq,Ir(:));

    MTPA_fd = interp1(MTPA.Ir_layers,MTPA.fd,Ir(:));
    MTPA_fq = interp1(MTPA.Ir_layers,MTPA.fq,Ir(:));
    MTPV_fd = interp1(MTPV.Ir_layers,MTPV.fd,Ir(:));
    MTPV_fq = interp1(MTPV.Ir_layers,MTPV.fq,Ir(:));

    % Create all surfaces and contours, initially hidden
    for ii = 1:num_Ir_values
        if ~sum(isnan(Z2(:,:,ii)))
            [~,guiData.hCont_dT(ii)] = contourf(X(:,:,ii),Y(:,:,ii),Z2(:,:,ii),'LineWidth',1,'DisplayName','$\Delta T_{pp}$ (Nm)','ShowText','on', 'Visible', 'off');
        else
            [~,guiData.hCont_dT(ii)] = contourf(X(:,:,ii),Y(:,:,ii),zeros(size(X(:,:,ii))),'LineWidth',1,'DisplayName','$\Delta T_{pp}$ (Nm)','ShowText','on', 'Visible', 'off');
        end
        
        [~,guiData.hCont_T(ii)] = contour(X(:,:,ii),Y(:,:,ii),Z1(:,:,ii),'-','LineColor',[0.0 0.0 0.0],'LineWidth',1,'DisplayName','$T$ (Nm)','showText','on', 'Visible', 'off');
        if labels == 0
            legLabel = '$I$ (A)';
        else
            legLabel = '$\lambda$ (Vs)';
        end
        [~,guiData.hCont_I(ii)] = contour(X(:,:,ii),Y(:,:,ii),abs(X(:,:,ii)+j*Y(:,:,ii)),abs(unique(iLevels)),'-','LineColor',[1.0 0.0 0.0],'LineWidth',0.5,'DisplayName',legLabel,'ShowText','on', 'Visible', 'off');
        
        if labels == 0
            guiData.hPlot_MTPA(ii) = plot(MTPA_id(ii,:),MTPA_iq(ii,:),'-r','DisplayName','MTPA', 'Visible', 'off');
            guiData.hPlot_MTPV(ii) = plot(MTPV_id(ii,:),MTPV_iq(ii,:),'-b','DisplayName','MTPV', 'Visible', 'off');
        else
            guiData.hPlot_MTPA(ii) = plot(MTPA_fd(ii,:),MTPA_fq(ii,:),'-r','DisplayName','MTPA', 'Visible', 'off');
            guiData.hPlot_MTPV(ii) = plot(MTPV_fd(ii,:),MTPV_fq(ii,:),'-b','DisplayName','MTPV', 'Visible', 'off');
        end
        

    end
    legend('show','Location','northeast');
    

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
        
        IrValue = current_guiData.Ir(current_guiData.current_Ir_index);
        title(['Torque map @ Ir = ', num2str(round(IrValue,2)), ' A']);
        
        % Update visibility of all surfaces and contours
        for jj = 1:numel(current_guiData.hCont_dT) % Use numel for number of elements
            if jj == current_guiData.current_Ir_index
                set(current_guiData.hCont_dT(jj), 'Visible', 'on');
                set(current_guiData.hCont_T(jj), 'Visible', 'on');
                set(current_guiData.hCont_I(jj), 'Visible', 'on');
                set(current_guiData.hPlot_MTPA(jj), 'Visible', 'on');
                set(current_guiData.hPlot_MTPV(jj), 'Visible', 'on');
                set(current_guiData.hPlot_MTPV(jj), 'Visible', 'on');

                visibleHandles = [ ...
                    current_guiData.hCont_dT(current_guiData.current_Ir_index), ...
                    current_guiData.hCont_T(current_guiData.current_Ir_index), ...
                    current_guiData.hCont_I(current_guiData.current_Ir_index), ...
                    current_guiData.hPlot_MTPA(current_guiData.current_Ir_index), ...
                    current_guiData.hPlot_MTPV(current_guiData.current_Ir_index)];
                legend(visibleHandles, 'Location', 'northeast');
                %legend(legend().String(5*(jj-1)+1:5*jj));
            else
                set(current_guiData.hCont_dT(jj), 'Visible', 'off');
                set(current_guiData.hCont_T(jj), 'Visible', 'off');
                set(current_guiData.hCont_I(jj), 'Visible', 'off');
                set(current_guiData.hPlot_MTPA(jj), 'Visible', 'off');
                set(current_guiData.hPlot_MTPV(jj), 'Visible', 'off');
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