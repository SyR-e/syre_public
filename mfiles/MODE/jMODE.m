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


function OUT=jMODE(options)

%% Reading parameters from options
generations     = options.MAXGEN;    % Maximum number of generations.
populationSize  = options.XPOP;      % Population size.
numVariables    = options.NVAR;      % Number of decision variables.
numObjectives   = options.NOBJ;      % Number of objectives.
Bounds          = options.FieldD;    % Optimization bounds.
Initial         = options.Initial;   % Initialization bounds.
costFunction    = options.mop;       % Cost function.
tau1            = options.tau1;
tau2            = options.tau2;
fl              = options.fl;
fu              = options.fu;
OUT.MatrixPop   = [];
OUT.MatrixFitness = [];



%% Initial random population
Parent = zeros(populationSize,numVariables);  % Parent population.
Mutant = zeros(populationSize,numVariables);  % Mutant population.
Child  = zeros(populationSize,numVariables);  % Child population.
functionEvaluations    = 0;                 % Function Evaluation.

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
F   = rand(populationSize,1)*(fu-fl) + fl;
CR  = rand(populationSize,1); %Crossover rates

for n=1:generations
    for i=1:populationSize
        rev=randperm(populationSize);
        if rand < tau1
            F(i)=fl + fu * rand;
        end
        %% Mutant vector calculation
        Mutant(i,:)= Parent(rev(1),:)+F(i)*(Parent(rev(2),:)-Parent(rev(3),:));
        
        for nvar=1:numVariables %Bounds are always verified
            if Mutant(i,nvar)<Bounds(nvar,1)
                Mutant(i,nvar) = Bounds(nvar,1);
            elseif Mutant(i,nvar)>Bounds(nvar,2)
                Mutant(i,nvar)=Bounds(nvar,1);
            end
        end
        
        if rand < tau2
            CR(i) = rand;
        end
        %% Crossover operator
        for nvar=1:numVariables
            if rand > CR(i)
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
     for i=1:populationSize
         if JxChild(i,:) <= JxParent(i,:)
             Parent(i,:) = Child(i,:);
             JxParent(i,:) = JxChild(i,:);
         end
     end
    tempPopulation = [Parent F CR JxParent; Child F CR JxChild];
    tempPopulation = nonDominationSort(tempPopulation, numObjectives, numVariables+2);
    tempPopulation = cutPopulation(tempPopulation, numObjectives, numVariables + 2,populationSize);

    Parent     = tempPopulation(:,1:numVariables);
    F          = tempPopulation(:,numVariables+1);
    CR          = tempPopulation(:,numVariables+2);
    Objectives = tempPopulation(:,numVariables+3:numVariables+2+numObjectives);
    Ranks      = tempPopulation(:,numVariables+3+numObjectives);
    %CrowdDist  = tempPopulation(:,end);
    
    PFront = Objectives(Ranks==1,:);
    PSet   = Parent(Ranks==1,:);
    OUT.MatrixPop      = cat(3, OUT.MatrixPop, Parent);
    OUT.MatrixFitness  = cat(3, OUT.MatrixFitness, Objectives);
    OUT.Xpop           = Parent;   % Population
    OUT.Jpop           = Objectives; % Poopulation's Objective Vector
    OUT.PSet           = PSet;     % Pareto Set
    OUT.PFront         = PFront;   % Pareto Front
    OUT.Param          = options;  % MODE Parameters
    OUT.eval_type      = options.eval_type;
    options.CounterGEN = n;
    options.CounterFES = functionEvaluations;
    
    [OUT, options]=PrinterDisplay(OUT,options); % To print results on screen
    
    if functionEvaluations>options.MAXFUNEVALS || n>options.MAXGEN
        disp('Termination criteria reached.')
        break;
    end
end

OUT.Xpop=Parent;
OUT.Jpop=Objectives;
[OUT.PFront, OUT.PSet]=DominanceFilter(PFront,PSet); %A Dominance Filter


if strcmp(options.SaveResults,'yes')
    thisfilepath = fileparts(which('data0.m'));
    filename=fullfile(thisfilepath,'results',['OUT_' datestr(now,30)]);
    save(filename,'OUT'); %Results are saved
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
for i=1:size(F,1)
    if numObjectives==1
        hold on;
        plot(options.CounterGEN,log(min(F(:,1))),'dk', ...
            'MarkerFaceColor','k'); grid on; hold on;
    elseif numObjectives==2
        hold on;
        plot(F(i,1),F(i,2),'dk','MarkerFaceColor','k');...
            grid on; hold on;
    elseif numObjectives==3
        hold on;
        plot3(F(i,1),F(i,2),F(i,3),'dk','MarkerFaceColor','k');...
            grid on; hold on;
    end
end

%% Print and Display information
% Modify at your convenience
%
function [OUT, Dat]=PrinterDisplay(OUT,Dat)

disp('------------------------------------------------')
disp(['Generation: ' num2str(Dat.CounterGEN)]);
disp(['FEs: ' num2str(Dat.CounterFES)]);
disp(['Pareto Front Size: ' mat2str(size(OUT.PFront,1))]);
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

%% Dominance Filter
%
% A filter based on dominance criteria
%
function [F, C]=DominanceFilter(F,C)

Xpop=size(F,1);
Nobj=size(F,2);
Nvar=size(C,2);
k=0;

for xpop=1:Xpop
    dom=0;
    
    for i=1:Xpop
        if F(xpop,:)==F(i,:)
            if xpop > i
                dom=1;
                break;
            end
        else
            if F(xpop,:)>=F(i,:)
                dom=1;
                break;
            end
        end
    end
    
    if dom==0
        k=k+1;
        F(k,:)=F(xpop,:);
        C(k,:)=C(xpop,:);
    end
end
F=F(1:k,:);
C=C(1:k,:);

%% Release and bug report:
%
% November 2012: Initial release