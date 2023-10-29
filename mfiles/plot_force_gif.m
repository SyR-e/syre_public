function plot_force_gif(geo,out,newDir,filemot)

% 
% plot_force_gif(geo,out,newDir,filemot)
% 


%% setup

Fr = out.SOL.Fr;
Ft = out.SOL.Ft;

names{1} = '$F_{r}$ [Nm/m$^2$]';
names{2} = '$F_{t}$ [Nm/m$^2$]';

filenames{1} = [newDir filemot(1:end-4) '_Fr.gif'];
filenames{2} = [newDir filemot(1:end-4) '_Ft.gif'];

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

am = (geo.r+geo.g-Fr(:,2));
th = Fr(:,1)*pi/180;
hplot(1) = plot(hax(1),am.*cos(th),am.*sin(th),'-b');
rMax = geo.r+geo.g+max(max(-Fr(:,2:end)));
set(hax(1),'XLim',rMax*[-1 1],'XTick',[]);
set(hax(1),'YLim',rMax*[-1 1],'YTick',[]);

am = (geo.r+geo.g-Ft(:,2));
th = Ft(:,1)*pi/180;
hplot(2) = plot(hax(2),am.*cos(th),am.*sin(th),'-b');
rMax = geo.r+geo.g+max(max(Ft(:,2:end)));
set(hax(2),'XLim',rMax*[-1 1],'XTick',[]);
set(hax(2),'YLim',rMax*[-1 1],'YTick',[]);

for ii=1:length(Fr(1,2:end))
    am = (geo.r+geo.g-Fr(:,ii+1));
    th = Fr(:,1)*pi/180;
    set(hplot(1),'XData',am.*cos(th),'YData',am.*sin(th));
    %     title(hax(1),['$F_r$ @ $\theta=' num2str(out.SOL.th(ii)-geo.th0(1),2) '^\circ$'])
    title(hax(1),['$p_r$ @ $\theta=' sprintf('%1.0f',out.SOL.th(ii)-geo.th0(1)) '^\circ$'])
    
    
    am = (geo.r+geo.g-Ft(:,ii+1));
    th = Fr(:,1)*pi/180;
    set(hplot(2),'XData',am.*cos(th),'YData',am.*sin(th));
%     title(hax(2),['$F_t$ @ $\theta=' num2str(out.SOL.th(ii)-geo.th0(1),2) '^\circ$'])
    title(hax(2),['$p_t$ @ $\theta=' sprintf('%1.0f',out.SOL.th(ii)-geo.th0(1)) '^\circ$'])
    
    for jj=1:2
        frame = getframe(hfig(jj));
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if ii==1
            imwrite(imind,cm,filenames{jj},'gif','Loopcount',inf);
        else
            imwrite(imind,cm,filenames{jj},'gif','WriteMode','append');
        end
    end
end