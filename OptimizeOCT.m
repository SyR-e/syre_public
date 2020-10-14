function OptimizeOCT(filename,pathname)
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

%% Overall Description
% This code implements a basic multi-objective optimization algorithm based
% on Diferential Evolution (DE) algorithm:
%
% Storn, R., Price, K., 1997. Differential evolution: A simple and 
% efficient heuristic for global optimization over continuous spaces. 
% Journal of Global Optimization 11, 341 – 359.
%
% When one objective is optimized, the standard DE runs; if two or more
% objectives are optimized, the greedy selection step in DE algorithm is 
% performed using a dominance relation.
%%

[filename, pathname, fltidx]=uigetfile(' *.fem', 'Pick a motor');
load(strrep(filename,'.fem','.mat'));   
[bounds, objs, geo, per, mat]=data0(dataSet);
dataSet.RQ = buildDefaultRQ(bounds);
   
dat.geo0=geo;
dat.per=per;
dat.mat=mat;
filemot=strcat(pathname,filename(1:end-4),'.fem');
if isoctave()  %OCT
        save ('-mat7-binary', 'dataSet.mat');
    else
        save('dataSet.mat');
    end
    
if isoctave()
  more off;
end
%%%%%%%%%% FEMM fitness handle %%%%%%%%%%%%%%%%%%%%%%%%%%
eval_type = 'MO_OA';                         % you can choose between "MO_OA" and "MO_GA"
                                             % "MO_GA" use multi-objective algorithm from matlab ga toolbox
                                             % "MO_OA" use multi-objective de algorithm
FitnessFunction = @(x)FEMMfitness(x,geo,per,mat,eval_type,filename); %OCT
%FitnessFunction = @(x)FEMMfitness(x,geo,per,mat,eval_type,filemot); %OCT
% FitnessFunction = @(x)FEMMfitness(x,geo,per,eval_type);
%FitnessFunction = @(x)zdtTestFunctions(x,1);
%bounds = [zeros(10,1) ones(10,1)];
%bounds = [0 1;-5*ones(9,1) 5*ones(9,1)]; ZDT4
dat.CostProblem = FitnessFunction;           % Cost function instance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variables regarding the optimization algorithm
% For guidelines for the parameter tuning see:
%
% Storn, R., Price, K., 1997. Differential evolution: A simple and 
% efficient heuristic for global optimization over continuous spaces. 
% Journal of Global Optimization 11, 341 – 359.
%
% Das, S., Suganthan, P. N., 2010. Differential evolution: A survey of the 
% state-of-the-art. IEEE Transactions on Evolutionary Computation. Vol 15,
% 4 - 31.

NOBJ = sum(objs(:,2));                       % Number of objectives
XPOP = dataSet.XPop;                         % Population size
Esc = 0.75;                                  % Scaling factor
Pm= 0.2;                                     % Croosover Probability

NVAR = size(bounds,1);                       % Number of decision variables
MAXGEN = dataSet.MaxGen;                     % Generation bound
MAXFUNEVALS = 20000*NVAR*NOBJ;               % Function evaluations bound
    
%% Variables regarding the optimization problem

switch eval_type
    case 'MO_OA'
        dat.FieldD = bounds;                  % Initialization 
        dat.Initial = bounds;                 % Optimization bounds (see data0.m)
        dat.NOBJ = NOBJ;
        dat.NRES = 0;                         % Number of constraints
        dat.NVAR = NVAR;
        dat.mop = str2func('evaluateF');      % Cost function
        dat.eval_type = eval_type;

        dat.XPOP = XPOP;
        dat.Esc = Esc;
        dat.Pm  = Pm;
        dat.fl  = 0.1;
        dat.fu  = 0.9;
        dat.tau1= 0.1;
        dat.tau2= 0.1;
        
        dat.InitialPop=[];                    % Initial population (if any)
        dat.MAXGEN = MAXGEN;
        dat.MAXFUNEVALS = MAXFUNEVALS;                  
        dat.SaveResults='yes';                % Write 'yes' if you want to 
                                              % save your results after the
                                              % optimization process;
                                              % otherwise, write 'no';
        % Initialization (don't modify)
        dat.CounterGEN=0;
        dat.CounterFES=0;

        % Run the algorithm.
        OUT=MODE2(dat,dataSet);        
        
    case 'MO_GA'                              %'MO_GA' use multi-objective algorithm from matlab ga toolbox
        A = []; b = [];
        Aeq = []; beq = [];
        lb = bounds(:,1);
        ub = bounds(:,2);
        numberOfVariables = NVAR;

        % Adding Visualization
        options = gaoptimset('PlotFcns',@gaplotpareto);

        %
        options = gaoptimset(options,'DistanceMeasureFcn',{@distancecrowding,'genotype'});
        options = gaoptimset(options,'ParetoFraction',0.4);

        % popolazione
        options = gaoptimset(options,'PopulationSize',XPOP);%200

        %
        options = gaoptimset(options,'PopInitRange',[lb';ub']);
        options = gaoptimset(options,'TolFun',0.0001,'StallGenLimit',1000);

        % generazioni
        options = gaoptimset(options,'TimeLimit',MAXFUNEVALS,'Generations',MAXGEN); %100
        options = gaoptimset(options,'UseParallel','always'); 
        options = gaoptimset(options,'Display','diagnose'); 
        % Run the algorithm.
        [x,fval,exitFlag,output,population,scores] = gamultiobj(FitnessFunction,numberOfVariables,[],[],[],[],lb,ub,options);
        fprintf('The number of points on the Pareto front was: %d\n', size(x,1));

        OUT.eval_type = eval_type;
        OUT.PSet = x;
        OUT.PFront = fval;
        OUT.Xpop = population;
        OUT.Jpop = scores;
        OUT.Param = output;
        % salvataggi
        thisfilepath = fileparts(which('syre.png'));
        filename=fullfile(thisfilepath,'results',['OUT_' datestr(now,30)]);
        
%        if isoctave()  %OCT Octave can't save function handles in Matlab .mat file
%            file_name1 = strcat(thisfilepath,'\results','\OUT_',datestr(now,30),'.mat');
%            save ('-mat7-binary', file_name1,'OUT','geo0','per');
%            clear file_name1
%        else
            save(filename,'OUT','geo0','per');
%        end
        
        if isoctave() %OCT
            h=gcf();
            fig_name=strcat(thisfilepath,'\results','\OUT_',datestr(now,30));
            hgsave(h,[fig_name]);
        else
            saveas(gcf,filename);
        end
    
        Pareto_Fvals = OUT.PFront;
        Pareto_front = OUT.PSet;
        disp('PostProcessing of current optimization result...')      
%         evalParetoFront(filename,'data0.m')
end
        
%% matlabpool close
%
%% Release and bug report:
%
% November 2012: Initial release
% October 2016: Second release