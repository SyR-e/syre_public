function plotSleeveDesignStats(dataSet,app)
%PLOTSLEEVEDESIGNSTATS Summary of this function goes here
%   Detailed explanation goes here

clc,
[~, ~, geo, ~, mat] = data0(dataSet);
[geo, ~, mat] = interpretRQ(geo.RQ,geo,mat);

button = questdlg('Choose a plot','SELECT','Sleeve Design Explict','Plots versus speed','Airgap radius sensitivity','Sleeve Design Explict');

if isequal(button,'Sleeve Design Explict')

    [par,in] = eval_sleeveThickness(dataSet, geo, mat);

    %Chose the area of interest
    prompt = {'Sleeve thickness vector [mm]:','Interference vector [mm]:'};
    dlgtitle = 'Input';
    dims = [1 35; 1 35];
    definput = {'1:0.5:5','0.1:0.1:1'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    in.h_1_vector   = str2num(answer{1});
    in.du_12_vector = str2num(answer{2});

    plot_sleeve_design_explicit(par,in)

elseif isequal(button,'Plots versus speed')

    [par,in] = eval_sleeveThickness(dataSet, geo, mat);
    plot_sleeve_versus_speed(par,in)

elseif isequal(button,'Plots versus radius')

    [par,in] = eval_sleeveThickness(dataSet, geo, mat);
    plot_sleeve_versus_radius(par,in)

elseif isequal(button,'Airgap radius sensitivity')

    %Chose the area of interest
    prompt = {'Airgap radius vector [mm]:','Sleeve thickness vector [mm]:','Interference vector [mm]:'};
    dlgtitle = 'Input';
    dims = [1 35;1 35; 1 35];
    definput = {'50:2:80','2','0.1:0.1:1'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    in.r_1i_vector  = str2num(answer{1});
    in.h_1_vector   = str2num(answer{2});
    in.du_12_vector = str2num(answer{3});

    plot_sleeve_airgap_radius_sensitivity(in,app)
end