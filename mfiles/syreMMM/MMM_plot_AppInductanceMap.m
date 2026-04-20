function MMM_plot_AppInductanceMap(motorModel)

% load data
Id = motorModel.AppInductanceMap_dq.Id;
Iq = motorModel.AppInductanceMap_dq.Iq;
Ld = motorModel.AppInductanceMap_dq.Ld;
Lq = motorModel.AppInductanceMap_dq.Lq;
Tr = motorModel.AppInductanceMap_dq.Tr;
Tm = motorModel.AppInductanceMap_dq.Tm;

if strcmp(motorModel.data.motorType, 'VFM')
    Fdm = motorModel.AppInductanceMap_dq.Fdm;
    Fqm = motorModel.AppInductanceMap_dq.Fqm;
else
    Fm = motorModel.AppInductanceMap_dq.Fm;
end

axisType = motorModel.data.axisType;

pathname = motorModel.data.pathname;
motName  = motorModel.data.motorName;
resFolder = checkPathSyntax([motName '_results\MMM results\' 'Apparent Inductance Maps - ' int2str(motorModel.data.tempPM) 'deg\']);

%% Check if motor is VFM and has MS dimension
if strcmp(motorModel.data.motorType, 'VFM')

    MS_values = motorModel.FluxMap_dq.MS;

    %% Create figures with MS slider for VFM
    figNames{1} = 'InducD';
    figNames{2} = 'InducQ';
    figNames{3} = 'FluxD';
    figNames{4} = 'FluxQ';
    figNames{5} = 'Anisotropy';
    figNames{6} = 'TorqueReluctance';
    figNames{7} = 'TorquePM';

    switch axisType
        case 'SR'
            csi = Ld./Lq;
            csiName = '$L_d/L_q$';
        case 'PM'
            csi = Lq./Ld;
            csiName = '$L_q/L_d$';
    end

    % Create all figures with sliders
    for ii = 1:length(figNames)
        switch ii
            case 1
                Z_data = Ld;
                zlabel_str = '$L_{d}$ (H)';
            case 2
                Z_data = Lq;
                zlabel_str = '$L_{q}$ (H)';
            case 3
                Z_data = Fdm;
                zlabel_str = '$\lambda_{dm}$ (Vs)';
            case 4
                Z_data = Fqm;
                zlabel_str = '$\lambda_{qm}$ (Vs)';
            case 5
                Z_data = csi;
                zlabel_str = csiName;
            case 6
                Z_data = Tr;
                zlabel_str = '$T_{rel}$ (Nm)';
            case 7
                Z_data = Tm;
                zlabel_str = '$T_{PM}$ (Nm)';
        end
        
        figureFullName = [pathname resFolder figNames{ii} '.fig'];
        labels = {'$i_d$ (A)', '$i_q$ (A)', zlabel_str};
        

        hfig(ii) = slider_inductance(Id, Iq, Z_data, reshape(MS_values(1,1,:),[],1), figureFullName, labels, 'app');
        
    end

else
    %% Original code for non-VFM motors
    figNames{1} = 'InducD';
    figNames{2} = 'InducQ';
    figNames{3} = 'FluxM';
    figNames{4} = 'Anisotropy';
    figNames{5} = 'TorqueReluctance';
    figNames{6} = 'TorquePM';
    
    switch axisType
        case 'SR'
            csi = Ld./Lq;
            csiName = '$L_d/L_q$';
        case 'PM'
            csi = Lq./Ld;
            csiName = '$L_q/L_d$';
    end
    
    for ii=1:length(figNames)
        hfig(ii) = figure();
        figSetting();
        hax(ii) = axes('OuterPosition',[0 0 1 1],...
            'XLim',[min(Id,[],'all') max(Id,[],'all')],...
            'YLim',[min(Iq,[],'all') max(Iq,[],'all')],...
            'PlotBoxAspectRatio',[1 1 0.8]);
        xlabel('$i_d$ (A)')
        ylabel('$i_q$ (A)')
        view(3)
        set(hfig(ii),'FileName',[pathname resFolder figNames{ii} '.fig'])
        set(hfig(ii),'Name',figNames{ii});
        switch ii
            case 1
                zlabel('$L_{d}$ (H)')
                set(gca,'ZLim',[min(Ld,[],'all') max(Ld,[],'all')]);
            case 2
                zlabel('$L_{q}$ (H)')
                set(gca,'ZLim',[min(Lq,[],'all') max(Lq,[],'all')]);
            case 3
                zlabel('$\lambda_m$ (Vs)')
                set(gca,'ZLim',[min(Fm,[],'all') max(Fm,[],'all')]);
                switch axisType
                    case 'SR'
                        view(0,0)
                    case 'PM'
                        view(90,0)
                end
            case 4
                zlabel(csiName)
                set(gca,'ZLim',[min(csi,[],'all') max(csi,[],'all')]);
            case 5
                zlabel('$T_{rel}$ (Nm)')
            case 6
                zlabel('$T_{PM}$ (Nm)')
        end
    end

    surf(hax(1),Id,Iq,Ld,'FaceColor','interp','EdgeColor','none')
    contour3(hax(1),Id,Iq,Ld,'EdgeColor','k','ShowText','off')
    surf(hax(2),Id,Iq,Lq,'FaceColor','interp','EdgeColor','none')
    contour3(hax(2),Id,Iq,Lq,'EdgeColor','k','ShowText','off')
    surf(hax(3),Id,Iq,Fm,'FaceColor','interp','EdgeColor','k','LineWidth',1.5)
    surf(hax(4),Id,Iq,csi,'FaceColor','interp','EdgeColor','none')
    contour3(hax(4),Id,Iq,csi,'EdgeColor','k','ShowText','off')
    surf(hax(5),Id,Iq,Tr,'FaceColor','interp','EdgeColor','none')
    contour3(hax(5),Id,Iq,Tr,'EdgeColor','k','ShowText','off')
    surf(hax(6),Id,Iq,Tm,'FaceColor','interp','EdgeColor','none')
    contour3(hax(6),Id,Iq,Tm,'EdgeColor','k','ShowText','off')
end

%% Save figures
answer = 'No';
answer = questdlg('Save figures?','Save','Yes','No',answer);
if strcmp(answer,'Yes')
    if ~exist([pathname resFolder],'dir')
        mkdir([pathname resFolder]);
    end
    
    for ii=1:length(hfig)
        savePrintFigure(hfig(ii));
    end
end
end