function [dataSet,flagS,hfig] = sleeveDesign(dataSet)
%SLEEVEDESIGN

[~, ~, geo, ~, mat] = data0(dataSet);

%Design equations
[par,in] = eval_sleeveThickness(dataSet, geo, mat);

%Chose the area of interest
prompt = {'Min Sleeve thickness:','Max Sleeve thickness:','Min Prestress:','Max Prestress:'};
dlgtitle = 'Input';
dims = [1 35; 1 35; 1 35; 1 35];
definput = {'0.5','5','0','1'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(answer)
    flagS=0;
    return
end
in.h_1_min   = str2double(answer{1});
in.h_1_max   = str2double(answer{2});
in.du_12_min = str2double(answer{3});
in.du_12_max = str2double(answer{4});

%Obtain the results for the plot
[map] = plot_sleeve_design(par,in);

%Figure settings
hfig=figure('Name','Sleeve Designer','NumberTitle','off');
figSetting(15,12)
ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
xlim([in.h_1_min in.h_1_max]);
ylim([in.du_12_min in.du_12_max])
xlabel("Sleeve thickness $h_1$ in $mm$");
ylabel("Prestress $du_{12}$ in $mm$");
titletext = sprintf('Operating temperature = %d °C', (in.dT+20));
title(titletext);

%Plot the max speed
contour(map.h_1*10^3,map.du_12*10^3,map.n_max_sl','DisplayName','$n_{max}$','ShowText','on','LevelStep',1000);

if min(map.du_12_danger_limit) < max(map.du_12_break)
%Plot the danger zone:
%the sleeve stress reach the maximum stress before the lift-off speed
plot(map.h_1*10^3,map.du_12_danger_limit*10^3,'color','k','HandleVisibility','off','LineWidth',0.5) %contour line
    if in.du_12_max > min(map.du_12_danger_limit*10^3)
        patch([map.h_1*10^3 fliplr(map.h_1*10^3)], [map.du_12_danger_limit*10^3 max(ylim)*ones(size(map.du_12_danger_limit*10^3))], [0.8500 0.3250 0.0980],'DisplayName','Danger area','HandleVisibility','on','facealpha',0.6);
    end
end

%Plot the breakage zone:
%the sleeve prestress-stress reach the maximum stress
plot(map.h_1*10^3,map.du_12_break*10^3,'color','k','HandleVisibility','off','LineWidth',0.5); %contour line
if in.du_12_max > min(map.du_12_break*10^3)
    patch([map.h_1*10^3 fliplr(map.h_1*10^3)], [map.du_12_break*10^3 max(ylim)*ones(size(map.du_12_break*10^3))], [0.6350 0.0780 0.1840],'DisplayName','Breakage area','HandleVisibility','on','facealpha',0.6);
end

%Plot the actual airgap
% plot(ones(1,10)*dataSet.AirGapThickness,linspace(in.du_12_min,in.du_12_max,10),'k','DisplayName','Actual airgap','HandleVisibility','on')

legend
set(hfig,'UserData',map);

button = questdlg('pick up a machine?','SELECT','Yes','No','Yes');

    while isequal(button,'Yes')
        figure(hfig)
        [dataSet.SleeveThickness,dataSet.SleeveInterference] = ginput(1);
        setup=inputdlg({'h_1','du_12'},'(h_1,du_12) values',1,{num2str(dataSet.SleeveThickness,3),num2str(dataSet.SleeveInterference,3)});
    

        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp(['h_1 = ' num2str(dataSet.SleeveThickness) ';']);
        disp(['du_12 = ' num2str(dataSet.SleeveInterference) ';']);
        
        button = questdlg('pick up another machine?','SELECT','Yes','No','Yes');
    
        if isequal(button,'No')
            buttonS = questdlg('save the last machine?','SELECT','Yes','No','Yes');
            figure(hfig)
        end
    end

    if ~exist('buttonS')
        buttonS='No';
    end
    
    flagS=0;

    if isequal(buttonS,'Yes')
        % save new machine
        flagS=1;
        newnamestring = ['h_1' num2str(dataSet.SleeveThickness,3) 'du_12' num2str(dataSet.SleeveInterference,3)];
        newnamestring(newnamestring=='.') = '';
        dataSet.currentfilename = strrep(dataSet.currentfilename,'.mat',[newnamestring '.mat']);
    end

 figure(hfig)