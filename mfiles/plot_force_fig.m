function plot_force_fig(geo,out,newDir,filemot)

% 
% plot_flxdn_fig(geo,out,newDir,filemot)
% 


%% setup

Fr = out.SOL.Fr;
Ft = out.SOL.Ft;

names{1} = '$F_{r}$ [$Nm$]';
names{2} = '$F_{t}$ [$Nm$]';

filenames{1} = 'Fr';
filenames{2} = 'Ft';

%% create figures

for ii=1:2
    hfig(ii)=figure();
    figSetting()
    hax(ii)=gca;
    xlabel('$\alpha$ [$^\circ mech$]')
    ylabel(names{ii})
end

%% plot figures and formatting

plot(hax(1),Fr(:,1),Fr(:,2:end));
plot(hax(2),Ft(:,1),Ft(:,2:end));

for ii=1:2
    hchild=get(hax(ii),'Children');
    format long
    for jj=1:length(hchild)
%         set(hchild(end-jj+1),'DisplayName',['$\theta= ' num2str(out.SOL.th(jj)-geo.th0(1),2) '^\circ$']);
        set(hchild(end-jj+1),'DisplayName',['$\theta= ' sprintf('%1.0f',out.SOL.th(jj)-geo.th0(1)) '^\circ$']);
    end
    legend(hax(ii),'show');
    set(hax(ii),'XLim',[min(Fr(:,1)) max(Fr(:,1))]);
end

%% save

for ii=1:2
    saveas(hfig(ii),[newDir filemot(1:end-4) '_' filenames{ii} '.fig']);
end

%% create figures

for ii=1:2
    hfig(ii)=figure();
    figSetting()
    hax(ii)=gca;
    GUI_Plot_Machine(gca,geo.rotor);
    hchild=get(gca,'Children');
    for jj=1:length(hchild)
        set(hchild(jj),'LineWidth',1,'HandleVisibility','off');
    end
    set(gca,'DataAspectRatio',[1 1 1]);
    set(gca,'XLimMode','auto','YLimMode','auto');
    set(gca,'XGrid','off','XMinorGrid','off');
    set(gca,'YGrid','off','YMinorGrid','off');
end

%% plot figures and formatting

am = (geo.r+geo.g-Fr(:,2:end));
th = Fr(:,1)*pi/180;
plot(hax(1),am.*cos(th),am.*sin(th));
rMax = geo.r+geo.g+max(max(-Fr(:,2:end)));
set(hax(1),'XLim',rMax*[-1 1],'XTick',linspace(-rMax,rMax,11));
set(hax(1),'YLim',rMax*[-1 1],'YTick',linspace(-rMax,rMax,11));

am = (geo.r+geo.g-Ft(:,2:end));
th = Ft(:,1)*pi/180;
plot(hax(2),am.*cos(th),am.*sin(th));
rMax = geo.r+geo.g+max(max(Ft(:,2:end)));
set(hax(2),'XLim',rMax*[-1 1],'XTick',linspace(-rMax,rMax,11));
set(hax(2),'YLim',rMax*[-1 1],'YTick',linspace(-rMax,rMax,11));

for ii=1:2
    hchild=get(hax(ii),'Children');
    for jj=1:length(hchild)
        %         set(hchild(end-jj+1),'DisplayName',['$\theta= ' num2str(out.SOL.th(jj)-geo.th0(1),2) '^\circ$']);
        set(hchild(end-jj+1),'DisplayName',['$\theta= ' sprintf('%1.0f',out.SOL.th(jj)-geo.th0(1)) '^\circ$']);        
    end
    legend(hax(ii),'show');
    %set(hax(ii),'XLim',[min(Fr(:,1)) max(Fr(:,1))]);
end

%% save

for ii=1:2
    saveas(hfig(ii),[newDir filemot(1:end-4) '_' filenames{ii} 'Geometry.fig']);
end


