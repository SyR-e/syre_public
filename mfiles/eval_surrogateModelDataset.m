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


function [OUT] = eval_surrogateModelDataset(options,dataSet)

%% Reading parameters from options
generations       = options.MAXGEN;    % Maximum number of generations. --> used for skip value for Sobol
populationSize    = options.XPOP;      % Population size.
numVariables      = options.NVAR;      % Number of decision variables.
%numObjectives     = options.NOBJ;      % Number of objectives.
Bounds            = options.FieldD;    % Optimization bounds.
%Initial           = options.Initial;   % Initialization bounds.
%scaleFactor       = options.Esc;       % Scaling fator in DE algorithm.
%crossOverRate     = options.Pm;        % Crossover probability in DE algorithm.
costFunction      = options.mop;       % Cost function.
OUT.MatrixPop     = [];
OUT.MatrixFitness = [];
OUT.MatrixPFront  = [];
OUT.MatrixPset    = [];
%tau1              = options.tau1;
optionsOK = options;
optionsOK.NOBJ = optionsOK.NOBJ+3;



%% Initial random population
Parent = zeros(populationSize,numVariables);  % Parent population.
Mutant = zeros(populationSize,numVariables);  % Mutant population.
Child  = zeros(populationSize,numVariables);  % Child population.
functionEvaluations    = 0;                 % Function Evaluation.



rng default
switch dataSet.optType
    case 'Surrogate model dataset (LHS)'
        X = lhsdesign(numVariables,populationSize);
        X = X';
    case 'Surrogate model dataset (Sobol)'
        pSobol = sobolset(numVariables,'Skip',generations);
        pSobol = scramble(pSobol,'MatousekAffineOwen');
        X = net(pSobol,populationSize);
end

for ii=1:populationSize
    for nvar=1:numVariables
        Parent(ii,nvar) = Bounds(nvar,1)+X(ii,nvar)*(Bounds(nvar,2)-Bounds(nvar,1));
    end
end

[JxParent,GeoxParent] = costFunction(Parent,optionsOK);

geo0=options.geo0;
per=options.per;
mat=options.mat;
thisfilepath = fileparts(which('syre.png'));

if strcmp(geo0.RQnames(end),'gamma')
    numVariables = numVariables-1;
end
XpopReal = zeros(populationSize,numVariables);
for ii=1:populationSize
    for nvar=1:numVariables
        tmp = strcat('GeoxParent{ii}.', geo0.RQnames(nvar));
        XpopReal(ii,nvar) = eval(tmp{1});
    end
end

OUT.Xpop           = Parent;   % Population
OUT.Jpop           = JxParent; % Poopulation's Objective Vector
OUT.Gpop           = GeoxParent;
OUT.XpopReal       = XpopReal;

if strcmp(geo0.RQnames(end),'gamma')
    OUT.XpopReal = [OUT.XpopReal OUT.Xpop(:,end)];
end


filename=fullfile(thisfilepath,'results',['OUT_SurogateModelDataSet_' datestr(now,30)]);
save(filename,'OUT','per','geo0','dataSet','mat');
close





