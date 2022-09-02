% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

%% MODE
% Multi-objective Evolutionary Algorithm (MOEA) based on Differential
% Evolution (DE).
% It implements a greedy selection based on pure dominance.
% DE algorithm has been introduced in:
%
% Storn, R., Price, K., 1997. Differential evolution: A simple and
% efficient heuristic for global optimization over continuous spaces.
% Journal of Global Optimization 11, 341 – 359.
%%


function OUT=MODE2(options, dataSet)

%% Reading parameters from options
generations     = options.MAXGEN;    % Maximum number of generations.
populationSize  = options.XPOP;      % Population size.
numVariables    = options.NVAR;      % Number of decision variables.
numObjectives   = options.NOBJ;      % Number of objectives.
Bounds          = options.FieldD;    % Optimization bounds.
Initial         = options.Initial;   % Initialization bounds.
scaleFactor     = options.Esc;       % Scaling fator in DE algorithm.
crossOverRate   = options.Pm;        % Crossover probability in DE algorithm.
costFunction    = options.mop;       % Cost function.
OUT.MatrixPop   = [];
OUT.MatrixFitness = [];
OUT.MatrixPFront  = [];
OUT.MatrixPset    = [];
tau1            = options.tau1;

[success,message,messageid] = mkdir('partial_optimization');

if ~isempty(message)
    rmdir('partial_optimization','s');
    mkdir('partial_optimization');
end

%% Initial random population
Parent = zeros(populationSize,numVariables);  % Parent population.
Mutant = zeros(populationSize,numVariables);  % Mutant population.
Child  = zeros(populationSize,numVariables);  % Child population.
functionEvaluations    = 0;                 % Function Evaluation.


if(generations==1)
    rng default
    Xlhc = lhsdesign(numVariables,populationSize);
    for i=1:populationSize
        for nvar=1:numVariables
            Parent(i,nvar) = Bounds(nvar,1)+Xlhc(nvar,i)*(Bounds(nvar,2)-Bounds(nvar,1));
        end
    end

    [JxParent,GeoxParent] = costFunction(Parent,options);

    geo0=options.geo0;
    per=options.per;
    mat=options.mat;
    thisfilepath = fileparts(which('syre.png'));
    
    if strcmp(geo0.RQnames(end),'gamma')
        numVariables = numVariables-1;
    end
    XpopReal = zeros(populationSize,numVariables);
    for i=1:populationSize
        for nvar=1:numVariables
            tmp = strcat('GeoxParent{i}.', geo0.RQnames(nvar));
            XpopReal(i,nvar) = eval(tmp{1});
        end
    end

    OUT.Xpop           = Parent;   % Population
    OUT.Jpop           = JxParent; % Poopulation's Objective Vector
    OUT.Gpop           = GeoxParent;
    OUT.XpopReal       = XpopReal;
    if strcmp(geo0.RQnames(end),'gamma')
        OUT.XpopReal = [OUT.XpopReal OUT.Xpop(:,end)];
    end
    

    filename=fullfile(thisfilepath,'results',['OUT_LHC_' datestr(now,30)]);
    save(filename,'OUT','per','geo0','dataSet','mat');
    close

else
    for i=1:populationSize
        for nvar=1:numVariables
            Parent(i,nvar) = Initial(nvar,1)+(Initial(nvar,2) - Initial(nvar,1))*rand();
        end
    end

    if size(options.InitialPop,1)>=1
        Parent(1:size(options.InitialPop,1),:)=options.InitialPop;
    end
    disp('------------------------------------------------')
    disp(['Initialization process']); %#ok<*NBRAK>
    disp(['Evaluating initial population']);
    disp('------------------------------------------------')

    JxParent = costFunction(Parent,options);
    functionEvaluations = functionEvaluations+populationSize;



    %% Evolution process

    for n=1:generations
        Tstart=now;
        for i=1:populationSize
            rev=randperm(populationSize);

            %% Mutant vector calculation
            if rand<tau1
                sf=rand;
            else
                sf=scaleFactor;
            end
            Mutant(i,:)= Parent(rev(1),:)+sf*(Parent(rev(2),:)-Parent(rev(3),:));

            for nvar=1:numVariables %Bounds are always verified
                if Mutant(i,nvar)<Bounds(nvar,1)
                    Mutant(i,nvar) = Bounds(nvar,1);
                elseif Mutant(i,nvar)>Bounds(nvar,2)
                    Mutant(i,nvar)=Bounds(nvar,1);
                end
            end

            %% Crossover operator
            for nvar=1:numVariables
                if rand > crossOverRate
                    Child(i,nvar) = Parent(i,nvar);
                else
                    Child(i,nvar) = Mutant(i,nvar);
                end
            end

        end
        disp('------------------------------------------------')
        disp(['Evaluating ' num2str(n) '-th generation...']);
        disp('------------------------------------------------')

        JxChild = costFunction(Child,options);
        functionEvaluations=functionEvaluations+populationSize;

        %% Selection
        for k=1:populationSize
            if JxChild(k,:) <= JxParent(k,:)
                Parent(k,:) = Child(k,:);
                JxParent(k,:) = JxChild(k,:);
            end
        end
        tempPopulation = [Parent JxParent; Child JxChild];
        tempPopulation = nonDominationSort(tempPopulation, numObjectives, numVariables);
        tempPopulation = cutPopulation(tempPopulation, numObjectives, numVariables,populationSize);

        Parent     = tempPopulation(:,1:numVariables);
        Objectives = tempPopulation(:,numVariables+1:numVariables+numObjectives);
        Ranks      = tempPopulation(:,numVariables+1+numObjectives);
        %CrowdDist  = tempPopulation(:,end);

        Tend = now;
        EvalTime = Tend-Tstart;
        EndTime = datevec(EvalTime*(generations-n)+now);

        PFront = Objectives(Ranks==1,:);
        PSet   = Parent(Ranks==1,:);
        OUT.MatrixPop      = cat(3, OUT.MatrixPop, Parent);
        OUT.MatrixFitness  = cat(3, OUT.MatrixFitness, Objectives);
        OUT.Xpop           = Parent;   % Population
        OUT.Jpop           = Objectives; % Poopulation's Objective Vector
        OUT.PSet           = PSet;     % Pareto Set
        OUT.PFront         = PFront;   % Pareto Front
        OUT.MatrixPset{n}  = PSet;
        OUT.MatrixPFront{n}= PFront;
        OUT.Param          = options;  % MODE Parameters
        OUT.eval_type      = options.eval_type;
        options.CounterGEN = n;
        options.CounterFES = functionEvaluations;
        options.EvalTime   = datevec(EvalTime);
        options.EndTime    = EndTime;

        [OUT, options]=PrinterDisplay(OUT,options); % To print results on screen

        %        if isoctave()  %OCT
        %            current_dir = pwd();
        %            file_name2 = strcat(current_dir,'\partial_optimization\generation_', int2str(n),'.mat');
        %            save ('-mat-binary', file_name2);
        %            clear file_name2 current_dir
        %        else
        save(['partial_optimization\generation_' int2str(n)]);
        %        end

        if functionEvaluations>options.MAXFUNEVALS || n>options.MAXGEN
            disp('Termination criteria reached.')
            break;
        end

        % Delete tmp files
        if isoctave() %OCT
            %     if strcmp(dataSet.RMVTmp,'ON')
            %        disp('Deleting tmp files...')
            %        ttt1=pwd();
            %        tempdirname = strcat(ttt1,'\tmp');
            %        if exist(tempdirname,'dir')
            %            rmdir(tempdirname,'s');
            %            clear tempdirname ttt1
            %        end
            %        if exist([cd,'\tmp'],'dir') == 0
            %            mkdir([cd,'\tmp']);
            %        end
            %        disp('tmp files deleted')
            %      end
        else
            if strcmp(dataSet.RMVTmp,'ON')
                disp('Deleting tmp files...')
                if exist([cd,'\tmp'],'dir')
                    [status,msg] = rmdir([cd,'\tmp'],'s');
                end
                if exist([cd,'\tmp'],'dir') == 0
                    [status,msg] = mkdir([cd,'\tmp']);
                end
                disp('tmp files deleted')
            end
        end
        %
    end

    OUT.Xpop=Parent;
    OUT.Jpop=Objectives;
    %[OUT.PFront, OUT.PSet]=DominanceFilter(PFront,PSet); %A Dominance Filter
    s=size(OUT.MatrixPop);

    M1=[];
    M2=[];
    if numel(s)<3
        M1=[M1;OUT.MatrixPset{1}]; % #ok<AGROW>
        M2=[M2;OUT.MatrixPFront{1}]; % #ok<AGROW>
    else
        for i=1:s(3)
            M1=[M1;OUT.MatrixPset{i}]; %#ok<AGROW>
            M2=[M2;OUT.MatrixPFront{i}]; %#ok<AGROW>
        end
    end
    if isoctave() %OCT
        I=paretoset(M2);
    else
        I=isparetosetMember(M2);
    end

    OUT.PFront=(M2(I==1,:));
    OUT.PSet=(M1(I==1,:));

    if strcmp(options.SaveResults,'yes')
        geo0=options.geo0;
        per=options.per;
        mat=options.mat;
        thisfilepath = fileparts(which('syre.png'));
        %    if isoctave() %OCT
        %    filename=fullfile(thisfilepath,'results',['OUT_' datestr(now,30) '.mat']);
        %    save('-mat7-binary', filename,'OUT','per','geo0','dataSet','mat'); %Results are saved
        %    else
        filename=fullfile(thisfilepath,'results',['OUT_' datestr(now,30)]);
        save(filename,'OUT','per','geo0','dataSet','mat');
        %     if per.MechStressOptCheck
        %         save(filename,'structModel','sVonMises','-append');
        %     end%Results are saved
        %    end

        clear geo0 per
    end

    disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    disp('Red  asterisks : Set Calculated.')
    disp('Black diamonds : Filtered Set.')
    if strcmp(options.SaveResults,'yes')
        disp(['Check out OUT_' datestr(now,30) ...
            ' variable on folder for results.'])
    end
    disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')

    F = OUT.PFront;

    if options.NOBJ==3
        figure
        figSetting
        plot3(F(:,1),F(:,2),F(:,3),'dk','MarkerFaceColor','k');
        grid on;
        hold on;
    elseif options.NOBJ==2
        figure
        figSetting
        plot(F(:,1),F(:,2),'dk','MarkerFaceColor','k');
        grid on;
        hold on;
    elseif options.NOBJ==1
        figure
        figSetting
        plot(options.CounterGEN,log(min(F(:,1))),'dk','MarkerFaceColor','k');
        grid on;
        hold on;
    end

    % plot(F(:,1),F(:,2),'dk','MarkerFaceColor','k');...
    %     grid on; hold on;
    % PostProcessing of current optimization result
    % postProcExecution = questdlg('Do you want to run post processing??','Postprocessing','Yes','No','Yes');
    % if strcmp(postProcExecution,'Yes')
    disp('PostProcessing of current optimization result...');
    [~,ff]=fileparts(filename);
    % if strcmp(func2str(options.CostProblem),
    % '@(x)FEMMfitness(x,geo,per,mat,eval_type)')
    evalParetoFront(ff,dataSet);
    % else
    %     evalParetoFrontX(ff,dataSet);
    % end

end

%% Print and Display information
% Modify at your convenience
%
function [OUT, Dat]=PrinterDisplay(OUT,Dat)

disp('------------------------------------------------')
disp(['Generation: ' num2str(Dat.CounterGEN)]);
disp(['FEs: ' num2str(Dat.CounterFES)]);
disp(['Pareto Front Size: ' mat2str(size(OUT.PFront,1))]);
disp(['Evaluation Time: ' int2str(Dat.EvalTime(4)) ' h ' num2str(Dat.EvalTime(5)) ' min ' num2str(round(Dat.EvalTime(6))) ' sec']);
disp(['Actual time             : ' datestr(now())]);
disp(['End of evolution process: ' datestr(Dat.EndTime)]);
disp('------------------------------------------------')

if mod(Dat.CounterGEN,1)==0
    if Dat.NOBJ==3
        stem3(OUT.PFront(:,1),OUT.PFront(:,2),OUT.PFront(:,3),'*r');
        grid on; hold on; drawnow
    elseif Dat.NOBJ==2
        figSetting
        plot(OUT.PFront(:,1),OUT.PFront(:,2),'*r');
        grid on; hold on; drawnow
    elseif Dat.NOBJ==1
        figSetting
        plot(Dat.CounterGEN,log(min(OUT.PFront(:,1))),'*r');
        grid on; hold on; drawnow
        %     else
        %         figSetting
        %         plot(OUT.PFront(:,1),OUT.PFront(:,2),'*r');
        %         grid on; hold on; drawnow
        %         title('Pareto front on the first two objs only')
    end
end


%% Release and bug report:
%
% November 2012: Initial release