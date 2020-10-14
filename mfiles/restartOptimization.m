function OUT=restartOptimization(flag)
% 
% OUT=restartOptimization(flag)
% 
% Useful when the optimization is interrupted and you must restart.
% (code copied and edited from MODE2.m)

clc


%% Load partial save of the optimization
if flag
    pathname=[cd '\partial_optimization\'];
    tmp=dir(pathname);
    FILENAME=tmp(end).name;
else
    [FILENAME, pathname, FILTERINDEX] = uigetfile('.mat', 'Load the last generation of partial save');
end

load([pathname FILENAME]);

%% Setting parameters of the optimization
% options = OUT.Param;
% 
% generations       = options.MAXGEN;    % Maximum number of generations.
% populationSize    = options.XPOP;      % Population size.
% numVariables      = options.NVAR;      % Number of decision variables.
% numObjectives     = options.NOBJ;      % Number of objectives.
% Bounds            = options.FieldD;    % Optimization bounds.
% Initial           = options.Initial;   % Initialization bounds.
% scaleFactor       = options.Esc;       % Scaling fator in DE algorithm.
% crossOverRate     = options.Pm;        % Crossover probability in DE algorithm.
% costFunction      = options.mop;       % Cost function.
% OUT.MatrixPop     = [];
% OUT.MatrixFitness = [];
% OUT.MatrixPFront  = [];
% OUT.MatrixPset    = [];
% tau1              = options.tau1;

save('dataSet','dataSet');

startGen=size(OUT.MatrixPop,3)+1;


figure()
if mod(options.CounterGEN,1)==0
    if options.NOBJ==3
        
        plot3(OUT.PFront(:,1),OUT.PFront(:,2),OUT.PFront(:,3),'*r');
        grid on;hold on,drawnow
    elseif options.NOBJ==2
        
        plot(OUT.PFront(:,1),OUT.PFront(:,2),'*r'); grid on;hold on,drawnow
    elseif options.NOBJ==1
        
        plot(options.CounterGEN,log(min(OUT.PFront(:,1))),'*r'); ...
            grid on; hold on,drawnow
    end
end



%% Evolution process (from MODE2.m)

for n=startGen:generations
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
    OUT.MatrixPFront{n}   = PFront;
    OUT.Param          = options;  % MODE Parameters
    OUT.eval_type      = options.eval_type;
    options.CounterGEN = n;
    options.CounterFES = functionEvaluations;
    options.EvalTime   = datevec(EvalTime);
    options.EndTime    = EndTime;
    
    [OUT, options]=PrinterDisplay(OUT,options); % To print results on screen
    
    save(['partial_optimization\generation_' int2str(n)]);
    
    if functionEvaluations>options.MAXFUNEVALS || n>options.MAXGEN
        disp('Termination criteria reached.')
        break;
    end
    
    
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
I=isparetosetMember(M2);
OUT.PFront=(M2(I==1,:));
OUT.PSet=(M1(I==1,:));

if strcmp(options.SaveResults,'yes')
    geo0=options.geo0;
    per=options.per;
    mat=options.mat;
    thisfilepath = fileparts(which('syre.png'));
    filename=fullfile(thisfilepath,'results',['OUT_' datestr(now,30)]);
    save(filename,'OUT','per','geo0','dataSet','mat'); %Results are saved
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

F=OUT.PFront;

plot(F(:,1),F(:,2),'dk','MarkerFaceColor','k');...
    grid on; hold on;
% PostProcessing of current optimization result
% postProcExecution = questdlg('Do you want to run post processing??','Postprocessing','Yes','No','Yes');
% if strcmp(postProcExecution,'Yes')
disp('PostProcessing of current optimization result...');
[~,ff]=fileparts(filename);
% if strcmp(func2str(options.CostProblem), '@(x)FEMMfitness(x,geo,per,mat,eval_type)')
evalParetoFront(ff,dataSet);
% else
%     evalParetoFrontX(ff,dataSet);
% end


delete('dataSet.mat');
% end

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
        
        plot3(OUT.PFront(:,1),OUT.PFront(:,2),OUT.PFront(:,3),'*r');
        grid on;hold on,drawnow
    elseif Dat.NOBJ==2
        
        plot(OUT.PFront(:,1),OUT.PFront(:,2),'*r'); grid on;hold on,drawnow
    elseif Dat.NOBJ==1
        
        plot(Dat.CounterGEN,log(min(OUT.PFront(:,1))),'*r'); ...
            grid on; hold on,drawnow
    end
end
