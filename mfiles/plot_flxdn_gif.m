function plot_flxdn_gif(geo,out,newDir,filemot)

% 
% plot_flxdn_gif(geo,out,newDir,filemot)
% 


%% setup

Bg = out.SOL.Bg;
Bt = out.SOL.Bt;
By = out.SOL.By;

Q = 6*geo.p*geo.q;
alpha_th = geo.wt/(geo.r+geo.g)*180/pi;
alpha_slot = 360/Q;
Qs = geo.Qs;

names{1} = '$B_{g}$ [$T$]';
names{2} = '$B_{t}$ [$T$]';
names{3} = '$B_{y}$ [$T$]';

filenames{1} = [newDir filemot(1:end-4) 'Bgap.gif'];
filenames{2} = [newDir filemot(1:end-4) 'Btooth.gif'];
filenames{3} = [newDir filemot(1:end-4) 'Byoke.gif'];

%% create figures

for ii=1:3
    hfig(ii)=figure();
    figSetting()
    hax(ii)=gca;
    xlabel('$\alpha$ [$^\circ mech$]')
    ylabel(names{ii})
    set(hax(ii),'XLim',[min(Bg(:,1)) max(Bg(:,1))]);
    
    for jj=1:Qs
        X1=[0
            alpha_th/2
            alpha_th/2
            0];
        Y= [-0.1
            -0.1
            +0.1
            +0.1];
        X2=[alpha_slot-alpha_th/2
            alpha_slot
            alpha_slot
            alpha_slot-alpha_th/2];
        X1=X1+alpha_slot*(jj-1)*ones(length(Y),1);
        X2=X2+alpha_slot*(jj-1)*ones(length(Y),1);
        
        fill(hax(ii),X1,Y,'b','HandleVisibility','off');
        fill(hax(ii),X2,Y,'b','HandleVisibility','off');
    end
end

%% axis limit
tmp=max(max(abs(Bg(:,2:end))));
set(hax(1),'Ylim',tmp*[-1 1]);
tmp=max(max(abs(Bt(:,2:end))));
set(hax(2),'Ylim',tmp*[-1 1]);
tmp=max(max(abs(By(:,2:end))));
set(hax(3),'Ylim',tmp*[-1 1]);


%% plot

plot(hax(1),Bg(:,1),Bg(:,2:end));
plot(hax(2),Bt(:,1),Bt(:,2:end));
plot(hax(3),By(:,1),By(:,2:end));

for ii=1:length(Bg(1,2:end))
    
    cla(hax(1))
    plot(hax(1),Bg(:,1),Bg(:,ii+1),'-b')
    title(hax(1),['$\theta=' num2str(out.SOL.th(ii)-geo.th0,2) ' ^\circ$'])
    drawnow
    
    cla(hax(2))
    plot(hax(2),Bt(:,1),Bt(:,ii+1),'-b')
    title(hax(2),['$\theta=' num2str(out.SOL.th(ii)-geo.th0,2) ' ^\circ$'])
    drawnow
    
    cla(hax(3))
    plot(hax(3),By(:,1),By(:,ii+1),'-b')
    title(hax(3),['$\theta=' num2str(out.SOL.th(ii)-geo.th0,2) ' ^\circ$'])
    drawnow
    
    for jj=1:3
        % GIF save
        frame=getframe(hfig(jj));
        im=frame2im(frame);
        [imind,cm]=rgb2ind(im,256);
        if ii==1
            imwrite(imind,cm,filenames{jj},'gif','Loopcount',inf);
        else
            imwrite(imind,cm,filenames{jj},'gif','WriteMode','append');
        end
    end
end