% Copyright 2019
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expres_trafs or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results] = skinEffect_eval(dataSet)

clc

% defaultCurrent   = [0.5 1 2 4 8 12 25 50 100 150 200 250 300];

if nargin==0
    defaultFrequency = [1 10 50 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500];
    load LastPath
    [filename,pathname,success] = uigetfile([pathname '_slotModel.fem'],'Select FEMM slot model');
    if success
        save LastPath pathname
    else
        error('No file selected')
    end
    
    prompt = {
        'Stator frequency vector [Hz]',...
%         'Stator current vector [Apk]',...
        };
    name = 'Skin Effect Evaluation';
    numlines = 1;
    answers = {
        'defaultFrequency'
%         '1'
        };
    
    answers = inputdlg(prompt,name,numlines,answers);
    fRef = eval(answers{1});
%     iRef = eval(answers{2});
    
    load([pathname filename],'dataSet')
    dataSet.SlotConductorFrequency = fRef;
    
else
    pathname = [dataSet.currentpathname dataSet.currentfilename(1:end-4) '_results\FEA results\' dataSet.currentfilename(1:end-4) '_slotModel\'];
    filename = ['slotModel.fem'];
    fRef     = dataSet.SlotConductorFrequency;
    tRef     = dataSet.SlotConductorTemperature;
end

[~,~,~,per,mat] = data0(dataSet);


% resFolder
date  = fix(clock);
year  = date(1);
month = date(2);
day   = date(3);
hour  = date(4);
mins  = date(5);

resFolder=[pathname 'evaluation' int2str(year) int2str(month) int2str(day) '-' int2str(hour) int2str(mins) ' - ' int2str(per.tempcu) 'deg\'];
mkdir(resFolder);

copyfile([pathname filename(1:end-4) '.fem'],[resFolder filename(1:end-4) '.fem'])
copyfile([pathname filename(1:end-4) '.mat'],[resFolder filename(1:end-4) '.mat'])

% initialization
[fRef,tRef] = meshgrid(fRef,tRef);
[nR,nC] = size(fRef);
fRef = fRef(:);
tRef = tRef(:);

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Slot model evaluation...')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

% DC simulation
disp('DC FEMM simulation...')
tic
outDC = skinEffect_point(20,0,per,mat,[pathname filename]);
time = toc;
outDC.R = outDC.P/outDC.I^2;
disp(['DC FEMM simulations done in ' num2str(time) ' s'])

% AC simulations
ppState = parallelComputingCheck(1);

disp('AC FEMM simulations...')
if ppState==0
    estTime = time*length(fRef);
else
    estTime = time*length(fRef)/ppState;
end
disp(['Estimated simulation time: ' num2str(estTime) ' s']);
tic
if ppState==0
    for ii=1:length(fRef)
        [OUT{ii}]=skinEffect_point(tRef(ii),fRef(ii),per,mat,[pathname filename]);
    end
else
    parfor ii=1:length(fRef)
        [OUT{ii}]=skinEffect_point(tRef(ii),fRef(ii),per,mat,[pathname filename]);
    end
end
time=toc;
disp(['AC FEMM simulations done in ' num2str(time) ' s'])

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Slot model evaluated!')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

[~,dirName] = skinEffect_point(max(tRef),0,per,mat,[pathname filename]);
copyfile([dirName 'slot_0.fem'],[resFolder filename(1:end-4) '_0Hz.fem'])
copyfile([dirName 'slot_0.ans'],[resFolder filename(1:end-4) '_0Hz.ans'])

[~,dirName] = skinEffect_point(max(tRef),max(fRef),per,mat,[pathname filename]);
copyfile([dirName 'slot_0.fem'],[resFolder filename(1:end-4) '_' int2str(max(fRef)) 'Hz.fem'])
copyfile([dirName 'slot_0.ans'],[resFolder filename(1:end-4) '_' int2str(max(fRef)) 'Hz.ans'])


% post-processing
results.I = zeros(size(OUT));
results.f = zeros(size(OUT));
results.P = zeros(size(OUT));
results.Q = zeros(size(OUT));
results.T = zeros(size(OUT));
results.T0 = 20;

for ii=1:length(OUT)
    results.I(ii) = OUT{ii}.I;
    results.f(ii) = OUT{ii}.f;
    results.P(ii) = OUT{ii}.P;
    results.Q(ii) = OUT{ii}.Q;
    results.T(ii) = OUT{ii}.T;
end

results.R = results.P./results.I.^2;
results.L = results.Q./results.I.^2./(2*pi*results.f);
results.k = results.R./(outDC.R.*(1+mat.SlotCond.alpha.*(results.T-results.T0)));

results.I = reshape(results.I,[nR,nC]);
results.f = reshape(results.f,[nR,nC]);
results.P = reshape(results.P,[nR,nC]);
results.Q = reshape(results.Q,[nR,nC]);
results.R = reshape(results.R,[nR,nC]);
results.L = reshape(results.L,[nR,nC]);
results.k = reshape(results.k,[nR,nC]);
results.T = reshape(results.T,[nR,nC]);

save([resFolder 'skinEffectResults.mat'],'results','OUT','outDC','dataSet');



if nR==1
    % single current level
    figure()
    figSetting()
    xlabel('$f$ [Hz]')
    ylabel('$k_{AC}=\frac{R_{AC}}{R_{DC}}$');
    plot(results.f,results.k,'-bo')
    saveas(gcf,[resFolder 'skinEffectFactor.fig']);
else
    figure()
    figSetting()
    view(3)
    surf(results.f,results.T,results.k)
    xlabel('$f$ [Hz]')
    ylabel('$\Theta_{Cu}$ [$^\circ$C]')
    zlabel('$k_{AC}=\frac{R_{AC}}{R_{DC}}$');
    saveas(gcf,[resFolder 'skinEffectFactor3D.fig']);

    figure()
    figSetting()
    xlabel('$f$ [Hz]')
    ylabel('$k_{AC}=\frac{R_{AC}}{R_{DC}}$');
    for ii=1:size(results.T,1)
        plot(results.f(ii,:),results.k(ii,:),'-o','DisplayName',['$\Theta_{PM}=' num2str(results.T(ii,1)) '^\circ$C'])
    end
    legend('show','Location','northwest')
    saveas(gcf,[resFolder 'skinEffectFactor.fig']);
end
