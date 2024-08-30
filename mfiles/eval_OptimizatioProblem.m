% Copyright 2024
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

function [OUT] = eval_OptimizatioProblem(dataSet)
%
% Launch SyR-e optimization or surrogate model dataset computation
% (previously implemented directly in the optimize button callback)

save('dataSet','dataSet');
figure()
[bounds, objs, geo, per, mat] = data0(dataSet);
filemot=([dataSet.currentpathname dataSet.currentfilename]);
dat.geo0 = geo;
per.objs = objs;
dat.per = per;
dat.mat = mat;
if (strcmp(dataSet.optType,'MODE Design')||strcmp(dataSet.optType,'MODE refine'))
    eval_type = 'MO_OA';
else
    eval_type = 'singt';
end

FitnessFunction = @(x)FEMMfitness(x,geo,per,mat,eval_type,strrep(filemot,'.mat','.fem'));

dat.CostProblem = FitnessFunction;           % Cost function instance
%% ========================================================================
NOBJ        = length(objs);
NOBJ        = sum(objs(:,2));
XPOP        = dataSet.XPop;
Esc         = 0.75;
Pm          = 0.2;
NVAR        = size(bounds,1);
MAXGEN      = dataSet.MaxGen;
MAXFUNEVALS = 20000*NVAR*NOBJ;

% Variables regarding the optimization problem
dat.FieldD      = bounds;
dat.Initial     = bounds;
dat.NOBJ        = NOBJ;
dat.NRES        = 0;
dat.NVAR        = NVAR;
dat.mop         = str2func('evaluateF');
dat.eval_type   = eval_type;
dat.XPOP        = XPOP;
dat.Esc         = Esc;
dat.Pm          = Pm;
dat.fl          = 0.1;
dat.fu          = 0.9;
dat.tau1        = 0.1;
dat.tau2        = 0.1;
dat.InitialPop  = [];
dat.MAXGEN      = MAXGEN;
dat.MAXFUNEVALS = MAXFUNEVALS;
dat.SaveResults = 'yes';
dat.CounterGEN  = 0;
dat.CounterFES  = 0;

% Run the algorithm.

switch dataSet.optType
    case {'MODE Design','MODE Refine'}
        OUT = MODE2(dat,dataSet);
        rmdir('partial_optimization','s');
    case {'Surrogate model dataset (LHS)','Surrogate model dataset (Sobol)'}
        OUT = eval_surrogateModelDataset(dat,dataSet);
end

delete('dataSet.mat');

if nargout()==0
    clear OUT;
end



