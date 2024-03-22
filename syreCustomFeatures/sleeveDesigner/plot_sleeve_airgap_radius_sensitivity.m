function plot_sleeve_airgap_radius_sensitivity(in,app)
%PLOT_SLEEVE_AIRGAP_RADIUS_SENSITIVITY 

h_1_vector   = in.h_1_vector;
du_12_vector = in.du_12_vector;
r_1i_vector  = in.r_1i_vector;

in.du_12_min = min(du_12_vector);
in.du_12_max = max(du_12_vector);

in.h_1_min = min(h_1_vector);
in.h_1_max = max(h_1_vector);

du_12_vector_length = double((in.du_12_min:0.001:in.du_12_max)*10^-3);
h_1_vector_length   = (in.h_1_min:0.1:in.h_1_max)*10^-3;

speed_over_radius_and_prestress = NaN(length(r_1i_vector),length(du_12_vector_length));
speed_over_radius_and_thickness = NaN(length(r_1i_vector),length(h_1_vector_length));

dataSet = app.dataSet;

preValueStatorOuterRadius = dataSet.StatorOuterRadius; %Stator outer radius [mm]
preValueAirgapRadius = dataSet.AirGapRadius;
splitRatio = preValueAirgapRadius/preValueStatorOuterRadius;

if length(h_1_vector) == 1

    for i = 1:length(r_1i_vector)
        
        dataSet.StatorOuterRadius = r_1i_vector(i)/splitRatio;
        % %Load dataSet
        
        % 
        % %Assign the loop parameters
        % dataSet.AirGapRadius    = r_1i_vector(i);
        % dataSet.SleeveThickness = h_1_vector;
        
        if i == 1
            [dataSet] = scaleMachine(dataSet,preValueStatorOuterRadius);
        else
            [dataSet] = scaleMachine(dataSet,r_1i_vector(i-1)/splitRatio);
        end 
        %Update dataSet
        
        % app.dataSet = dataSet;
        % %draw the machine to get the lamination parametes
        % app = GUI_APP_DrawMachine(app);
         
        %Derive equivalent model parameters
        [~, ~, geo, ~, mat] = data0(dataSet);
        [par,in] = eval_sleeveThickness(dataSet, geo, mat);
    
        %Assign the vector boundaries
        in.h_1_min   = min(h_1_vector);
        in.h_1_max   = max(h_1_vector);
        in.du_12_min = min(du_12_vector);
        in.du_12_max = max(du_12_vector);
        
        %Derive the useful maps
        [map] = plot_sleeve_design(par,in);
 
        speed_over_radius_and_prestress(i,:) = map.n_max_sl(1,:);
        breakage_area(i) = map.du_12_break;
        danger_area(i) = map.du_12_danger_limit;

    end

    hfig=figure('Name','Plot speed over radius','NumberTitle','off');
    figSetting(15,12)
    ax = gca;
    ax.XAxis.Exponent = 0;
    ax.YAxis.Exponent = 0;
    xlim([min(r_1i_vector) max(r_1i_vector)]);
    ylim([in.du_12_min in.du_12_max])
    xlabel("Airgap radius $r_{1i}$ in $mm$");
    ylabel("Prestress $du_{12}$ in $mm$");
    titletext = sprintf('Sleeve thickness = %.2g [mm]', (h_1_vector));
    title(titletext);
    
    if length(du_12_vector) == 1
        plot(r_1i_vector,speed_over_radius_and_prestress);
    else
        contour(r_1i_vector,map.du_12*10^3,speed_over_radius_and_prestress','DisplayName','$n_{max}$','ShowText','on','LevelStep',2000);
    end

    if min(map.du_12_danger_limit) < max(breakage_area)
    %Plot the danger zone:
    %the sleeve stress reach the maximum stress before the lift-off speed
    plot(r_1i_vector,danger_area*10^3,'color','k','HandleVisibility','off','LineWidth',0.5) %contour line
        if in.du_12_max > min(danger_area*10^3)
            patch([r_1i_vector fliplr(r_1i_vector)], [danger_area*10^3 max(ylim)*ones(size(danger_area*10^3))], [0.8500 0.3250 0.0980],'DisplayName','Danger area','HandleVisibility','on','facealpha',0.6);
        end
    end

    plot(r_1i_vector,breakage_area*10^3,'color','k','HandleVisibility','off','LineWidth',0.5); %contour line
    if in.du_12_max > min(breakage_area*10^3)
        patch([r_1i_vector fliplr(r_1i_vector)], [breakage_area*10^3 max(ylim)*ones(size(breakage_area*10^3))], [0.6350 0.0780 0.1840],'DisplayName','Breakage area','HandleVisibility','on','facealpha',0.6);
    end

    legend

    dataSet = app.dataSet;
    dataSet.AirGapRadius            = eval(app.AirGapRadiusEdit.Value);
    dataSet.SleeveThickness         = eval(app.RotorSleeveThicknessEditField.Value);
    dataSet.ShaftRadius             = eval(app.ShaftRadEdit.Value);
    app.dataSet = dataSet;
    app = GUI_APP_DrawMachine(app);
end

if length(du_12_vector) == 1

   for i = 1:length(r_1i_vector)

        dataSet = app.dataSet;
        dataSet.AirGapRadius    = r_1i_vector(i);
        dataSet.SleeveThickness = h_1_vector(1);
        app.dataSet = dataSet;
        app = GUI_APP_DrawMachine(app);
    
        [~, ~, geo, ~, mat] = data0(app.dataSet);
        [par,in] = eval_sleeveThickness(dataSet, geo, mat);
    
        in.h_1_min   = min(h_1_vector);
        in.h_1_max   = max(h_1_vector);
        in.du_12_min = min(du_12_vector);
        in.du_12_max = max(du_12_vector);
    
        [map] = plot_sleeve_design(par,in);
    
        speed_over_radius_and_thickness(i,:) = map.n_max_sl(:,1)';

    end

    hfig=figure('Name','Plot speed over radius','NumberTitle','off');
    figSetting(15,12)
    ax = gca;
    ax.XAxis.Exponent = 0;
    ax.YAxis.Exponent = 0;
    xlim([min(r_1i_vector) max(r_1i_vector)]);
    ylim([in.h_1_min in.h_1_max])
    xlabel("Airgap radius $r_{1i}$ in $mm$");
    ylabel("Sleeve thickness $h_1$ in $mm$");
    titletext = sprintf('Prestress = %.2g [mm]', du_12_vector);
    title(titletext);
    
    if length(h_1_vector) == 1
        plot(r_1i_vector,speed_over_radius_and_prestress);
    else
        contour(r_1i_vector,map.h_1*10^3,speed_over_radius_and_thickness','DisplayName','$n_{max}$','ShowText','on','LevelStep',2000);
    end
    legend

    dataSet = app.dataSet;
    dataSet.AirGapRadius            = eval(app.AirGapRadiusEdit.Value);
    dataSet.SleeveThickness         = eval(app.RotorSleeveThicknessEditField.Value);
    dataSet.ShaftRadius             = eval(app.ShaftRadEdit.Value);
    app.dataSet = dataSet;
    app = GUI_APP_DrawMachine(app);
 
% elseif length(in.n_max_vector) == 1

end


