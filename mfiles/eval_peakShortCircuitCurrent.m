% Copyright 2021
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

function [pkSCout] = eval_peakShortCircuitCurrent(dataIn)

if isfield(dataIn,'flagSave')
    flagSave = dataIn.flagSave;
else
    flagSave = 1;
end

if isfield(dataIn,'flagFEAfix')
    clc
    flagFEAfix = dataIn.flagFEAfix;
    geo        = dataIn.geo;
    mat        = dataIn.mat;
    per        = dataIn.per;
    axis_type  = geo.axisType;
    flagSave   = 0;
    RQ         = dataIn.RQ;
else
    flagFEAfix = 0;
end

if nargin()~=1
    error('Wrong number of function input!')
end

if ~flagFEAfix
    load([dataIn.currentpathname dataIn.currentfilename])


    if ~isfield(geo,'axisType')
        if strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype')
            geo.axisType = 'PM';
        else
            geo.axisType = 'SR';
        end
    end

    if ~strcmp(geo.axisType,dataIn.axisType)
        %geo.axisType = dataIn.axisType;
        if strcmp(dataIn.axisType,'PM')
            geo.th0 = geo.th0 - 90;
        else
            geo.th0 = geo.th0 + 90;
        end
    end

    % update loaded dataSet fields with GUI dataSet (dataIn) fields
    dataSet.RatedCurrent     = dataIn.RatedCurrent;
    dataSet.CurrLoPP         = dataIn.CurrLoPP;
    % dataSet.SimulatedCurrent = dataIn.SimulatedCurrent;
    dataSet.SimulatedCurrent = dataSet.RatedCurrent*dataSet.CurrLoPP;
    dataSet.GammaPP          = dataIn.GammaPP;
    dataSet.BrPP             = dataIn.BrPP;
    dataSet.tempPP           = dataIn.tempPP;
    dataSet.NumOfRotPosPP    = dataIn.NumOfRotPosPP;
    dataSet.AngularSpanPP    = dataIn.AngularSpanPP;
    dataSet.EvalSpeed        = dataIn.EvalSpeed;
    dataSet.axisType         = dataIn.axisType;

    dataSet.currentpathname = dataIn.currentpathname;
    dataSet.currentfilename = dataIn.currentfilename;

    filename = dataSet.currentfilename;
    pathname = dataSet.currentpathname;

    axis_type  = dataSet.axisType;
end

% create result folder
if flagSave
    outFolder = [filename(1:end-4) '_results\FEA results\'];
    if ~exist([pathname outFolder],'dir')
        mkdir([pathname outFolder]);
    end

    resFolder = ['peakShortCircuitCurrent_' datestr(now,30) '\'];
    mkdir([pathname outFolder],resFolder);
    resFolder = [pathname outFolder resFolder];


    save([resFolder 'pkScOut.mat'],'dataSet');
end

if ~flagFEAfix
    % update per structure
    per.nsim_singt      = dataSet.NumOfRotPosPP;
    per.delta_sim_singt = dataSet.AngularSpanPP;
    per.tempPP          = dataSet.tempPP;
    per.BrPP            = dataSet.BrPP;
    per.i0              = dataSet.RatedCurrent;

    % initialize variables
    idq0 = dataSet.SimulatedCurrent.*cosd(dataSet.GammaPP)+j*dataSet.SimulatedCurrent.*sind(dataSet.GammaPP);
    fdq0 = nan(size(idq0));
    T0   = nan(size(idq0));
    idq  = nan(size(idq0));
    fdq  = nan(size(idq0));


    i0 = per.i0;

    motname = [dataSet.currentpathname dataSet.currentfilename(1:end-4) '.fem'];
else
    motname = dataIn.filemot;
    idq0 = dataIn.idq0;
    i0 = per.i0;
end

tol = 100;
maxIter = 10;



for ii=1:length(idq0)
        
    if ~flagFEAfix
        disp('Starting FEMM simulations...')
        disp(' Healthy point simulation...')
        disp(['Operating point ' int2str(ii) ' of ' int2str(length(idq0)) ' --> ' int2str(real(idq0(ii))) '+j*' int2str(imag(idq0(ii))) ' A'])
        per.overload = abs(idq0(ii))/per.i0;
        per.gamma = angle(idq0(ii))*180/pi;
        [~,~,~,out,~] = FEMMfitness([],geo,per,mat,'singt',motname);
        fdq0(ii) = out.fd+j*out.fq;
        T0(ii)   = out.T;
        disp([' Healthy point simulation done: |fdq| = ' num2str(abs(fdq0(ii)),4) ' Vs'])
    else
        per.overload = abs(idq0(ii))/per.i0;
        fdq0(ii) = dataIn.fdq0;
        T0(ii)   = dataIn.T0;
    end

   

    tol = abs(fdq0(ii))/100;

    done = 0;
    jj=1;

    iTest{ii} = nan(1,maxIter);
    fTest{ii} = nan(1,maxIter);

    while ~done
        if ii==1 && jj<3
            if jj==1
                iTest{ii}(jj) = 1;
            elseif jj==2
                iTest{ii}(jj) = 2;
            end
        else
            fVect = [];
            iVect = [];
            for zz=1:ii
                if zz==1
                    fVect = fTest{zz};
                    iVect = iTest{zz};
                else
                    fVect = [fVect (fTest{zz})];
                    iVect = [iVect (iTest{zz})];
                end
            end
            fVect = fVect(~isnan(fVect));
            iVect = iVect(~isnan(iVect));
            [fVect,index] = sort(fVect);
            iVect = iVect(index);

            if strcmp(axis_type,'SR')
                iTest{ii}(jj) = interp1(fVect,iVect,abs(fdq0(ii)),'linear','extrap');
            else
                iTest{ii}(jj) = interp1(fVect,iVect,-abs(fdq0(ii)),'linear','extrap');
            end
        end

        per.overload = abs(iTest{ii}(jj));
        if strcmp(axis_type,'SR')
            per.gamma = 90;
        else
            per.gamma = -180;
        end
        if ~flagFEAfix
            [~,~,~,out,~] = FEMMfitness([],geo,per,mat,'singt',motname);
        else
            [~,~,~,out,~] = FEMMfitness(RQ,geo,per,mat,'singt',motname);
        end

        if strcmp(axis_type,'SR')
            fTest{ii}(jj) = out.fq;
        else
            fTest{ii}(jj) = -out.fd;
        end

        if ~flagFEAfix
            disp(['  - Iteration ' int2str(jj) ': ' num2str(abs(iTest{ii}(jj)*i0),4) ' A --> ' num2str(abs(fTest{ii}(jj)),4) ' Vs'])
        end

        if (abs(fTest{ii}(jj))>abs(fdq0(ii))+tol)||(abs(fTest{ii}(jj))<abs(fdq0(ii))-tol)
            done = 0;
            jj = jj+1;
        else
            done = 1;
            if strcmp(axis_type,'SR')
                idq(ii) = j*iTest{ii}(jj);
                fdq(ii) = j*fTest{ii}(jj);
            else
                idq(ii) = -iTest{ii}(jj);
                fdq(ii) = fTest{ii}(jj);
            end
        end

        if ii>maxIter
            done=1;
        end
    end
end

if ~flagFEAfix
    disp('FEMM simulations done!!!')
end

idq = idq*per.i0;


pkSCout.idq0    = idq0;
pkSCout.fdq0    = fdq0;
pkSCout.T0      = T0;
pkSCout.idq     = idq;
pkSCout.fdq     = fdq;
pkSCout.iTest   = iTest;
pkSCout.fTest   = fTest;

if flagSave
    save([resFolder 'pkScOut.mat'],'pkSCout','-append');

    hfig(1) = figure();
    figSetting(16,10);
    hax(1) = axes(...
        'OuterPosition',[0 0 1 1],...
        'XLim',[0.5 length(idq)+0.5],'XTick',1:1:length(idq));
    xlabel('working point')
    ylabel('[A]')
    set(hfig(1),'FileName',[resFolder 'peakShortCircuitCurrent.fig'])
    hleg(1) = legend(hax(1),'show','Location','southoutside','Orientation','horizontal');


    hfig(2) = figure();
    figSetting(16,10);
    hax(2) = axes('OuterPosition',[0 0 1 1]);
    xlabel('[A]')
    ylabel('[Vs]')
    set(hfig(2),'FileName',[resFolder 'iterations.fig'])
    hleg(2) = legend(hax(2),'show','Location','northeastoutside');

    bar(hax(1),abs(idq),'BarWidth',0.8,'EdgeColor',0.5*[1 0 0],'FaceColor',[1 0 0],'DisplayName','peak short circuit current')
    bar(hax(1),abs(idq0),'BarWidth',0.5,'EdgeColor',0.5*[0 0 1],'FaceColor',[0 0 1],'DisplayName','initial current')


    fVect = [];
    iVect = [];
    for ii=1:length(iTest)
        fVect = [fVect (fTest{ii})];
        iVect = [iVect (iTest{ii})];
        iPlot = iTest{ii}(~isnan(iTest{ii}));
        fPlot = fTest{ii}(~isnan(fTest{ii}));
        plot(hax(2),iPlot*per.i0,fPlot,'-o','DisplayName',['working point ' int2str(ii)]);
    end

    iVect = iVect(~isnan(iVect));
    fVect = fVect(~isnan(fVect));
    [iVect,index] = sort(iVect);
    fVect = fVect(index);
    plot(hax(2),iVect*per.i0,fVect,':k','DisplayName','all simulations')

    for ii=1:length(hfig)
        savePrintFigure(hfig(ii));
    end
end