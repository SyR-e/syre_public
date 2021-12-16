% Copyright 2019
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

function [ichOut,ich]=eval_ich(setup)
%
% [ichOut,ich]=eval_ich(setup)
%
% Compute the characteristic current of the selected machine using FEMM
%
% There are three ways to use this function:
% 1) as a stand-alone function (no inputs):
%    The motor and the temperature vector will be set with proper windows.
%    FEA simulations are done with 6 position over 60 elt deg and results
%    are automatically saved.
% 2) from SyR-e GUI:
%    The input is the dataSet structure, with the settings from the GUI.
%    Results are automatically saved.
% 3) as function from script:
%    in this case the input is the structure "setup" and it is possible to
%    avoid the results save. FEA simulations are done with 6 positions
%    over 60 elt deg. The "setup" structure is composed as follows:
%     setup.pathname --> name of the path
%          .filename --> name of the .mat file
%          .tempVect --> vector of the PM temperatures
%          .flagSave --> =1-->save the data
%                        =0--> do not save the data but return the results
%

clc

if nargin()==0
    load LastPath.mat
    [filename,pathname,] = uigetfile([pathname '\.mat'],'Select a machine');
    
    load([pathname filename])
    
    % inputs
    prompt={
        'Temperature vector'
        };
    name = 'Characteristic current';
    numlines = 1;
    answer = {mat2str(mat.LayerMag.temp.temp)};
    
    answer = inputdlg(prompt,name,numlines,answer);
    
    tempVect = eval(answer{1});
    flagSave = 1;
    
    setup.filename = filename;
    setup.pathname = pathname;
    setup.tempVect = tempVect;
    setup.flagSave = flagSave;
    
    dataSet.tempPP        = tempVect;
    dataSet.BrPP          = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,tempVect);
    dataSet.NumOfRotPosPP = 6;
    dataSet.AngularSpanPP = 60;
    
    
else
    if isfield(setup,'flagSave') % input con setup
        pathname = setup.pathname;
        filename = setup.filename;
        tempVect = setup.tempVect;
        flagSave = setup.flagSave;
        
        load([pathname filename]);
        
        dataSet.tempPP        = tempVect;
        dataSet.BrPP          = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,tempVect);
        dataSet.NumOfRotPosPP = 6;
        dataSet.AngularSpanPP = 60;
        
        
    else % input da GUI, con dataSet
        
        load([setup.currentpathname setup.currentfilename]);
        dataSet.tempPP          = setup.tempPP;
        dataSet.BrPP            = setup.BrPP;
        dataSet.NumOfRotPosPP   = setup.NumOfRotPosPP;
        dataSet.AngularSpanPP   = setup.AngularSpanPP;
        dataSet.currentfilename = setup.currentfilename;
        dataSet.currentpathname = setup.currentpathname;
        
        if ~isfield(dataSet,'axisType')
            if strcmp(dataSet.TypeOfRotor,'SPM') || strcmp(dataSet.TypeOfRotor,'Vtype')
                dataSet.axisType = 'PM';
            else
                dataSet.axisType = 'SR';
            end
        end
        
        if ~strcmp(dataSet.axisType,dataIn.axisType)
            if strcmp(dataSet.axisType,'PM')
                geo.th0 = geo.th0 + 90;
            else
                geo.th0 = geo.th0 - 90;
            end
        end
        
        tempVect = dataSet.tempPP;
        pathname = dataSet.currentpathname;
        filename = dataSet.currentfilename;
        flagSave = 1;
    end
end


per.nsim_singt      = dataSet.NumOfRotPosPP;
per.delta_sim_singt = dataSet.AngularSpanPP;
per.tempPP          = dataSet.tempPP;
per.BrPP            = dataSet.BrPP;

dataSet.currentfilename = filename;
dataSet.currentpathname = pathname;

if (strcmp(dataSet.TypeOfRotor,'SPM')||strcmp(dataSet.TypeOfRotor,'Vtype'))
    per.gamma=-180;
    axes_type='PM';
else
    per.gamma=90;
    axes_type='SR';
end

motname=[pathname filename(1:end-4) '.fem'];

MaxIter = 10;
i0 = per.i0;

ichVect = nan(size(tempVect));
FmVect  = nan(size(tempVect));

disp('Starting FEMM simulations...')

for tt=1:length(tempVect)
    ichTest{tt} = nan(1,MaxIter);
    FmTest{tt}  = nan(1,MaxIter);
    
    
    disp(['Temperature ' int2str(tt) ' of ' int2str(length(tempVect)) ' - ' int2str(tempVect(tt)) ' Celsius degree'])
    
    done=0;
    ii=1;
    
    while ~done
        if ii==1
            ichTest{tt}(ii)=0;
            tol=-inf;
        elseif ii==2
            if tt==1
                ichTest{tt}(ii) = 1;
            else
                ichTest{tt}(ii)=ichVect(tt-1)/i0;
            end
            tol=FmTest{tt}(1)/100;
        else
            ichTest{tt}(ii)=ichTest{tt}(ii-2)-FmTest{tt}(ii-2)/(FmTest{tt}(ii-1)-FmTest{tt}(ii-2))*(ichTest{tt}(ii-1)-ichTest{tt}(ii-2));
            tol=FmTest{tt}(1)/100;
        end
        
        per.overload = ichTest{tt}(ii);
        per.BrPP     = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,tempVect(tt));
        
        [~,~,~,out,~] = FEMMfitness([],geo,per,mat,'singt',motname);
        
        if strcmp(axes_type,'PM')
            FmTest{tt}(ii) = out.fd;
        else
            FmTest{tt}(ii) = -out.fq;
        end
        
        disp([' Iteration ' int2str(ii) ': ' num2str(ichTest{tt}(ii)*i0,4) ' A --> ' num2str(FmTest{tt}(ii),4) ' Vs'])
        
        if abs(FmTest{tt}(ii))>tol
            ii=ii+1;
            done=0;
        else
            done=1;
            ichVect(tt) = ichTest{tt}(ii)*i0;
            FmVect(tt)  = FmTest{tt}(1);
        end
        
        if ii>MaxIter
            done=1;
        end
    end
    
end

disp('FEMM simulations done!!!')

ich=ichVect(1);

ichOut.ichVect  = ichVect;
ichOut.tempVect = tempVect;
ichOut.FmVect   = FmVect;
ichOut.i0       = i0;
ichOut.ichTest  = ichTest;
ichOut.fMTest   = FmTest;
ichOut.pathname = pathname;
ichOut.filename = filename;

if flagSave
    
    outFolder = [filename(1:end-4) '_results\FEA results\'];
    if ~exist([pathname outFolder],'dir')
        mkdir([pathname outFolder]);
    end
    
    resFolder = [pathname outFolder 'charCurr - ' datestr(now,30) '\'];
    mkdir(resFolder)
    
    save([resFolder 'ichOut.mat'],'ich','ichOut');
    
    figure()
    figSetting()
    xlabel('$\Theta_{PM}$ [$^\circ C$]')
    ylabel('$i_{ch}$ [$A$]')
    plot(tempVect,ichVect,'-bo');
    saveas(gcf,[resFolder 'ichVStemp.fig'])
    
    figure()
    figSetting()
    xlabel('$\Theta_{PM}$ [$^\circ C$]')
    ylabel('$\lambda_{m}$ [$Vs$]')
    plot(tempVect,FmVect,'-bo');
    saveas(gcf,[resFolder 'FmVStemp.fig'])
    
    figure()
    figSetting()
    xlabel('iteration')
    ylabel('$i_{ch}$ [$A$]')
    for ii=1:length(ichTest)
        plot(ichTest{ii}*i0,'-o','DisplayName',['$\Theta_{PM}=' num2str(tempVect(ii)) '^\circ C$']);
    end
    legend('show')
    saveas(gcf,[resFolder 'ichIterations.fig'])
    
    figure()
    figSetting()
    xlabel('iteration')
    ylabel('$\lambda_{m}$ [$Vs$]')
    for ii=1:length(FmTest)
        plot(FmTest{ii},'-o','DisplayName',['$\Theta_{PM}=' num2str(tempVect(ii)) '^\circ C$']);
    end
    legend('show')
    saveas(gcf,[resFolder 'FmIterations.fig'])
end












