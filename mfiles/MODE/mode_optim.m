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

function varargout = mode_optim(funfcn, popsize, lb, ub, varargin)

% MODE               Optimization with Differential Evolution
%
% INPUT ARGUMENTS:
% ================
%
%   obj_fun     The objective function of which the global minimum
%               will be determined (function_handle). For multi-
%               objective optimization, several objective functions
%               may be provided as a cell array of function handles,
%               or alternatively, in a single function that returns
%               the different function values along the second
%               dimension.
%               Objective functions must accept either a [popsize x
%               dimensions] matrix argument, or a [1 x dimensions]
%               vector argument, and return a [popsize x number of
%               objectives] matrix or [1 x number of objective]
%               vector of associated function values (number of
%               objectives may be 1). With the first format, the
%               function is evaluated vectorized, in  the second
%               case CELLFUN() is used, which is a bit slower in
%               general.
%
%   popsize     positive integer. Indicates the TOTAL population
%               size, that is, the number of individuals of all
%               populations combined.
%
%   lb, ub      The lower and upper bounds of the problem's search
%               space, for each dimension. May be scalar in case all
%               bounds in all dimensions are equal. Note that at
%               least ONE of these must have a size of [1 x
%               dimensions], since the problem's dimensionality is
%               derived from it.
%
%   which_ones  DE.     %2013/11/04 MG some part of the sript hava also
%                       this formulation for indicationg the topology of the optimization
%                       algorithm used wich_ones is set to 'DE', it is the only possibility.
%
%   options/    Sets the options to be used by GODLIKE. Options may
%   'option',   be a structure created by set_options, or given as
%      value    individual ['option', value] pairs. See set_options
%               for a list of all the available options and their
%               defaults.
%
% OUTPUT ARGUMENTS:
% =================
%
%   sol         The solution that minizes the problem globally,
%               of size [1 x dimensions]. For multi-objective
%               optimization, this indicates the point with the
%               smallest distance to the (shifted) origin.
%
%   fval        The globally minimal function value
%
%   exitflag    Additional information to facilitate fully automated
%               optimization. Negative is `bad', positive `good'. A
%               value of '0' indicates GODLIKE did not perform any
%               operations and exited prematurely. A value of '1'
%               indicates normal exit conditions. A value of '-1'
%               indicates a premature exit due to exceeding the preset
%               maximum number of function evaluations. A value of
%               '-2' indicates that the amount of maximum GODLIKE
%               iterations has been exceeded, and a value of '-3'
%               indicates no optimum has been found (only for single-
%               objective optimization).
%
%   output      structure, containing much additional information
%               about the optimization as a whole; see the manual
%               for a more detailed description.
%
%   (For multi-objective optimization only)
%
%   Pareto_front, Pareto_Fvals
%               The full set of non-dominated solutions, and their
%               associated function values.
%
%   See also pop_single, pop_multi, set_options.


%      Author : Rody P.S. Oldenhuis
% Affiliation : Delft University of Technology
%               Faculty of Aerospace Engineering
%               Dep. of Astrodynamics & Satellite Systems
%     Contact : oldenhuis@dds.nl
%   Licensing/
%    (C) info : Frankly I don't care what you do with it,
%               as long as I get some credit when you copy
%               large portions of the code ^_^

%% Initialize

% basic check on input
error(nargchk(4, inf, nargin));

% more elaborate check on input (nested function)
check_input;

% resize and reshape boundaries and dimensions (nested function)
[lb, ub, sze, dimensions, options] = reformat_input(lb, ub, varargin{:});

% 2013/11/04 MG wich_ones is unecessary because refers to original
% GODLIKE algorithm and it defines the type of different optimization
% algorithm used, for as a precaution is set manually to DE but in
% future this line must be eliminated.
global which_ones;

which_ones={'DE'};

% test input objective function(s) to determine the problem's dimensions,
% number of objectives and proper input format (nested function)
[options, single, multi, test_evaluations] = test_funfcn(options);

% initialize more variables
algorithms = 1;                          % number of algorithms to use
generation = 1;                         % this is the first generation
pop = cell(algorithms,1);        % cell array of [population] objects
num_funevaluations  = 0;                         % number of function evaluations
[converged, output] = check_convergence;         % initial output structure

% Initially, [output] is a large structure used to move data to and from all the
% subfunctions. Later, it is turned into the output argument [output] by removing some
% obsolete entries from the structure.

% do an even more elaborate check (the behavior of this
% nested function is changed by passing the number of
% requested output arguments)
check_input(nargout);

%% 2013/10/16 MG initialization for saving population and fitness at each iteration
%%
output.MatrixPop=[];
output.MatrixFitness=[];

%% GODLIKE loop

% GODLIKE loop
while ~converged
    
    % randomize population sizes (minimum is 5 individuals)
    frac_popsize = break_value(popsize, 5);
    
    % randomize number of iterations per algorithm
    % ([options.GODLIKE.ItersUb] is the maximum TOTAL amount
    % of iterations that will be spent in all of the algorithms combined.)
    frac_iterations = break_value(options.MODE.ItersUb, options.MODE.ItersLb);
    
    % frac_iterations=5;
    % shuffle (or initialize) populations
    pop = interchange_populations(pop);
    %         keyboard
    
    % perform algorithm iterations
    if strcmpi(pop.algorithm, 'MS')
        % Multi-start behaves differently; its needs to
        % execute its iterations inside pop_single.
        
        % save previous value of number of function evaluations
        prev_FE = pop.funevals;
        
        % pass data via arguments
        pop.iterate(frac_iterations, num_funevaluations);
        
        % adjust number of function evaluations made
        num_funevaluations = num_funevaluations + pop.funevals - prev_FE;
        
    else % Perform single iterations for all other algorithms
        counter = 0; % used for single-objective optimization
%         frac_iterations
        for j = 1:frac_iterations
            % FC increase generation
            generation = j;
            %%
            %% 2013/10/16 MG save population at each iteration
            output.MatrixPop=cat(3,output.MatrixPop,pop.individuals);
            output.MatrixFitness=cat(3,output.MatrixFitness,pop.fitnesses);
            %%
            
            % do single iteration on this population
            pop.iterate;
            % keyboard
            % calculate total number of function evaluations
            % Appareantly, pop{:}.funevals doesn't work. So
            % we have to do a loop through all of them.
            funevaluations = 0;
            
            if ~isempty(pop),
                funevaluations=funevaluations+pop.funevals;
            end
            
            num_funevaluations = test_evaluations + funevaluations;
            
            % check for convergence of this iteration
            if multi
                % all are non-dominated, first display progress, then exit the loop
                if all(pop.pop_data.front_number == 0)
                    if ~isempty(options.display),
                        display_progress;
                    end,
                    break
                end
            elseif single
                % check algorithm convergence
                [alg_converged, output, counter] = check_convergence(false,output,counter);
                % if converged, first display progress, then exit the loop
                if alg_converged
                    if ~isempty(options.display), display_progress; end, break
                end
            end % if
            
            % check function evaluations, and exit if it
            % surpasses the preset maximum
            if (num_funevaluations >= options.MaxFunEvals)
                % also display last iteration
                if ~isempty(options.display), display_progress; end,
                converged = true; break,
            end
            
            % display progress at every iteration
            if ~isempty(options.display), display_progress; end
            
        end % algorithm loop
    end
    
    % if we have convergence inside the algorithm
    % loop, break the main loop
    % FC JAN 14: following line prevent creation of output
    % variables when max number of function calls is reached
    % if converged, break, end
    
    % % increase generation
    % generation = generation + 1;
    
    % check maximum iterations
    if (generation >= options.MaxIters), converged = true; end
    %         keyboard
    output.generation=generation;
    % check for convergence (and update output structure)
    
    [converged, output] = check_convergence(converged, output);
    
end % GODLIKE loop

% display final results
if ~isempty(options.display), display_progress; end

%% output values

% multi-objective optimization
if multi
    varargout{1} = output.most_efficient_point;
    varargout{2} = output.most_efficient_fitnesses;
    varargout{3} = output.pareto_front_individuals;
    varargout{4} = output.pareto_front_fitnesses;
    varargout{5} = output.exitflag;
    % remove some fields from output structure
    output = rmfield(output, {'pareto_front_individuals','pareto_front_fitnesses',...
        'exitflag','most_efficient_point','most_efficient_fitnesses'});
    % and output what's left
    varargout{6} = output;
    
    % single-objective optimization
elseif single
    
    % if all went normal
    if isfield(output, 'global_best_individual')
        varargout{1} = output.global_best_individual;
        varargout{2} = output.global_best_funval;
        varargout{3} = output.exitflag;
        % remove some fields from output structure
        outpt.algorithms = output.algorithms;   outpt.funcCount = output.funcCount;
        outpt.message    = output.message;      outpt.algorithm_info = output.algorithm_info;
        outpt.iterations = output.iterations;
        % and output
        varargout{4} = outpt;
        
        % but, no optimum might have been found
    else
        varargout{1} = NaN(1, dimensions);
        varargout{2} = NaN;
        varargout{3} = -3;
        % remove some fields from output structure
        output = rmfield(output, {'global_best_funval', 'exitflag','descent_counter',...
            'best_individuals','best_funcvalues','previous_global_best_funval',...
            'previous_best_funcvalues'});
        % adjust message
        output.message = sprintf('%s\n\n All function values encountered were INF or NaN.\n',...
            output.message);
        % output
        varargout{4} = output;
    end
end

%% nested functions

% =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =
% initialization shizzle
% =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

    function check_input(varargin)
        if (nargin == 0)
            if isempty(funfcn)
                error('MODE:function_not_defined',...
                    'MODE requires at least one objective function.');
            end
            if isempty(lb)||isempty(ub)||isempty(popsize)
                error('MODE:lbubpopsize_not_defined',...
                    'MODE requires arguments [lb], [ub] and [popsize].');
            end
            if ~isnumeric(lb)||~isnumeric(ub)||~isnumeric(popsize)
                error('MODE:lbubpopsize_not_numeric',...
                    'Arguments [lb], [ub] and [popsize] must be numeric.');
            end
            if any(~isfinite(lb)) || any(~isfinite(ub)) || ...
                    any(  ~isreal(lb)) || any(~isreal(ub))
                error('MODE:lbub_not_finite',...
                    'Values for [lb] and [ub] must be real and finite.');
            end
            if ~isvector(lb) || ~isvector(ub)
                error('MODE:lbub_mustbe_vector',...
                    'Arguments [lb] and [ub] must be given as vectors.');
            end
            if ~isa(funfcn, 'function_handle')
                % might be cell array
                if iscell(funfcn)
                    for ii = 1:numel(funfcn)
                        if ~isa(funfcn{ii}, 'function_handle')
                            error('MODE:funfcn_mustbe_function_handle',...
                                'All objective functions must be function handles.');
                        end
                    end
                    % otherwise, must be function handle
                else
                    error('MODE:funfcn_mustbe_function_handle',...
                        'Objective function must be given as a function handle.');
                end
            end
            if (nargin == 6) && ~isstruct(varargin{2})
                error('MODE:options_mustbe_structure',...
                    'Argument [options] must be a structure.')
            end
            if any(lb > ub)
                error('MODE:lb_larger_than_ub',...
                    'All entries in [lb] must be smaller than the corresponding entries in [ub].')
            end
            if ~isscalar(popsize) || ~isreal(popsize) || ~isfinite(popsize) || popsize < 0
                error('MODE:popsize_is_bad',...
                    'Argument [popsize] must be a real, positive and finite scalar.')
            end
        else
            if (options.MODE.ItersLb > options.MODE.ItersUb)
                warning('MODE:ItersLb_exceeds_ItersUb',...
                    ['Value of options.MODE.ItersLb is larger than value of\n',...
                    'options.GODLIKE.ItersUb. Values will simply be swapped.']);
                u_b = options.MODE.ItersUb;
                options.MODE.ItersUb = options.MODE.ItersLb;
                options.MODE.ItersLb = u_b;
            end
            if (options.MODE.ItersLb > options.MODE.ItersUb)
                warning('MODE:MaxIters_exceeds_MinIters',...
                    ['Value of options.MinIters is larger than value of\n',...
                    'options.MaxIters. Values will simply be swapped.']);
                u_b = options.MaxIters;
                options.MaxIters = options.MinIters;
                options.MinIters = u_b;
            end
            if single
                % single objective optimization has a maximum of 4 output arguments
                error(nargoutchk(0, 4, varargin{1}))
            elseif multi
                % multi-objective optimization has a maximum of 6 output arguments
                error(nargoutchk(0, 6, varargin{1}))
            end
            if strcmpi(options.display, 'plot') && single && dimensions > 2
                warning('MODE:Plotting_not_possible',...
                    ['Display type was set to ''Plot'', but the number of\n',...
                    'decision variables exceeds 2. The search space can note be\n',...
                    'displayed. Set options.display to ''off'' or ''on'' to \n',...
                    '''on'' to supress this message.'])
            end
            if strcmpi(options.display, 'plot') && multi && options.num_objectives > 3
                warning('MODE:Plotting_not_possible',...
                    ['Display type was set to ''Plot'', but the number of\n',...
                    'objective functions exceeds 3. The Pareto front can \n',...
                    'not be displayed. Set options.display to ''off'' or \n',...
                    '''on'' to supress this message.'])
            end
        end % if
    end % nested function

% reshape, resize and redefine input to predictable formats
    function [lb, ub, sze, dimensions, options] = ...
            reformat_input(lb, ub, varargin)
        %         keyboard
        
        % set options
        if nargin <= 1, options = set_options; end   % defaults
        %         if nargin == 2, options = varargin{2}; end   % structure provided
        if nargin >= 2 , options = set_options(varargin{1:end}); end
        % individually provided
        
        % save the original size of [lb] or [ub]
        max_elements = max(numel(lb),numel(ub));
        if (max_elements == numel(lb)), sze = size(lb); else  sze = size(ub); end
        
        % force [lb] and [ub] to be row vectors
        lb = lb(:).';  ub = ub(:).';
        
        % replicate one or the other when their sizes are not equal
        if ~all(size(lb) == size(ub))
            if     isscalar(lb)
                lb = repmat(lb, size(ub));
            elseif isscalar(ub)
                ub = repmat(ub, size(lb));
            else
                error('MODE:lbub_sizes_incorrect',...
                    ['If the size of either [lb] or [ub] is equal to the problem''s dimenions\n',...
                    'the size of the other must be 1x1.'])
            end
        end
        
        % define [dimensions]
        dimensions = numel(lb);
        
    end % nested function

% test the function, and determine the amount of objectives. Here
% it is decided whether the optimization is single-objective or
% multi-objective.
    function [options, single, multi, fevals] = test_funfcn(options)
        
        % initialize
        fevals = 0;
        
        % split multi/single objective
        fun = funfcn; % make a copy
        if iscell(funfcn) && (numel(funfcn) > 1)
            % no. of objectives is simply the amount of provided objective functions
            options.num_objectives = numel(funfcn);
            % single is definitely false
            single = false;
        elseif iscell(funfcn) && (numel(funfcn) == 1)
            % single it true but might still change to false
            single = true;
            % also convert function to function_handle in this case
            funfcn = funfcn{1};
        else
            % cast fun to cell
            fun = {funfcn};
            % single is true but might still change to false
            single = true;
        end
        
        % loop through all objective functions
        % (also works for single function)
        for ii = 1:numel(funfcn)
            
            % reshape to original size
            lb_original = reshape(lb, sze);
            
            % try to evaluate the function
            try
                %                 keyboard
                % simply evaluate the function with the lower bound
                sol = fun{ii}(lb_original);
                
                % keep track of the number of function evaluations
                fevals = fevals + 1;
                
                % see whether single must be changed to multi
                if single && (numel(sol) > 1), single = false; end
                %                 keyboard
                % it might happen that more than one function is provided,
                % but that one of the functions returns more than one function
                % value. GODLIKE cannot handle that case
                if (numel(sol) > 1) && (ii > 1)
                    error('MODE:multimulti_not_allowed',...
                        ['MODE cannot optimize multiple multi-objective problems ',...
                        'simultaneously.\nUse GODLIKE multiple times on each of your objective ',...
                        'functions separately.\n\nThis error is generated because the first of ',...
                        'your objective functions returned\nmultiple values, while ',...
                        'you provided multiple objective functions. Only one of\nthese formats ',...
                        'can be used for multi-objective optimization, not both.'])
                end
                
                % if evaluating the function fails, throw error
            catch userFcn_ME
                pop_ME = MException('MODE:function_doesnt_evaluate',...
                    'MODE cannot continue: failure during function evaluation.');
                userFcn_ME = addCause(userFcn_ME, pop_ME);
                rethrow(userFcn_ME);
            end % try/catch
            
        end % for
        
        % see if the optimization is multi-objective
        multi = ~single;
        
    end % nested function

% =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =
% functions used in the main loop
% =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =

% break up some [value] into a vector of random integers
% of length [algorithms], that sums up to [value]
    function frac_value = break_value(value, Lb)
        % NOTE: The case of these variables [Lb] and [Ub] is important.
        % The GODLIKE arguments [lb] or [ub] may get overwritten!
        %         keyboard
        % only one algorithm - just return value
        if algorithms == 1, frac_value = value; return; end
        
        % initially, the upper bound is the value minus
        % (algorithms-1) times the lower bound
        Ub = value - (algorithms-1)*Lb;
        
        % create array of length [algorithms] that
        % sums to [value]
        frac_value = zeros(algorithms, 1);
        %         for ii = 1:algorithms-1 % note the minus one
        %             % random value (make sure it's not zero)
        %             rnd = 0; while (rnd == 0), rnd = round(rand*(Ub-Lb) + Lb); end
        %             frac_value(ii) = rnd;
        %             % adjust max. value for next iteration
        %             Ub = round((value - sum(frac_value))/(algorithms-ii));
        %         end % for
        
        % last entry is the difference of the sum of all values and the original value
        frac_value(end) = value - sum(frac_value);
        
        % sort at random
        [dummy, inds] = sort(rand(size(frac_value,1),1));
        frac_value = frac_value(inds);
        
    end % nested function

% shuffle and (re)initialize the population objects
    function pop = interchange_populations(pop)
        %         keyboard
        % just initialize populations if this is the first iteration
        if (generation == 1)
            
            %% 2013/10/27 MG modifica definisco a priori DE e la riga
            %% di definizione dovrà essere successivamente eliminata
            % set algorithm for this iteration
            %                 options.algorithm = which_ones{ii};
            options.algorithm = which_ones;
            % initialize population
            %                 keyboard
            if single
                pop = pop_single(funfcn, frac_popsize, lb, ub, sze, dimensions, options);
            else
                pop = pop_multi(funfcn, frac_popsize, lb, ub, sze, dimensions, options);
            end
            
            % we're done
            return
        end
        
        % don't shuffle if there's only one algorithm
        if (algorithms == 1),
            return,
        end
        
    end % nested function

% update output values, and check for convergence
    function [converged, output, counter] = ...
            check_convergence(converged, output, varargin)
        %         keyboard
        % some algorithms might be doubly used.
        % save which ones they are
        %          persistent sames
        
        % no input - initialize
        if (nargin == 0)
            %   keyboard
            % initially, no convergence
            converged = false;
            
            % some algorithms might be doubly used. Find out
            % which ones, and create proper indices
            %              sames = ones(algorithms, 1);
            
            %% Commentate le seguenti linee possibile eliminazione successiva:
            %              for ii = 1:algorithms
            %                  same        = strcmpi(which_ones, which_ones{ii});
            %                  sames(same) = 1:nnz(same);
            %              end
            
            % general settings
            which_ones='DE';
            output.algorithms = upper(which_ones); % algorithms used
            output.exitflag   = 0;                 % neutral exitflag
            output.message    = sprintf('No iterations have been performed.');
            output.funcCount  = 0;
            
            output.algorithm_info.(upper(which_ones)).funcCount  = 0;
            output.algorithm_info.(upper(which_ones)).iterations = 0;
            
            
            % initialize [output] for single-objective optimization
            if single
                output.descent_counter               = 0;
                output.global_best_individual        = NaN(1,dimensions);
                output.previous_global_best_individual = NaN(1,dimensions);
                output.global_best_funval            = inf;
                output.previous_global_best_funval   = inf;
                output.best_funcvalues               = inf(1,algorithms);
                output.previous_best_funcvalues      = inf(1,algorithms);
                output.best_individuals              = NaN(algorithms,dimensions);
                output.previous_best_individuals     = NaN(algorithms,dimensions);
                
                output.algorithm_info.(upper(which_ones)).last_population = [];
                output.algorithm_info.(upper(which_ones)).last_fitnesses  = [];
                
            end
            
            % initialize [output] for multi-objective optimization
            if multi
                output.pareto_front_individuals = [];
                output.pareto_front_fitnesses   = [];
                output.most_efficient_point     = [];
                output.most_efficient_fitnesses = [];
                
            end
            
            % we're done
            return
            
            % otherwise, update according to the current status of [pops]
        else %if nargin
            %              keyboard
            % both per-algorithm and global check needs to be performed.
            % the mode of operation depends on the presence of a third
            % input argument. If given, only the current populations is
            % checked. If  omitted, all populations are checked.
            if (nargin == 3),
                alg_conv = true; algorithm = 1; counter = varargin{1};
            else alg_conv = false;
            end
            
            % general stuff
            output.funcCount  = num_funevaluations;
            output.iterations = generation;
            
            output.algorithm_info.(upper(which_ones)).iterations = pop.iterations;
            output.algorithm_info.(upper(which_ones)).funcCount  = pop.funevals;
            
            
            % convergence might already have occured. Determine the reason
            if converged
                % maximum function evaluations have been exceeded.
                if (num_funevaluations >= options.MaxFunEvals)
                    output.exitflag = -1;
                    output.message = sprintf(['Optimization terminated:\n',...
                        ' Maximum amount of function evaluations has been reached.\n',...
                        ' Increase ''MaxFunEvals'' option.']);
                end
                % maximum allowable iterations have been exceeded.
                if (generation >= options.MaxIters)
                    output.exitflag = -2;
                    output.message = sprintf(['Optimization terminated:\n',...
                        ' Maximum amount of iterations has been reached.\n',...
                        ' Increase ''MaxIters'' option.']);
                end % if
            end % if
            
            % stuff specific for single objective optimization
            if single
                
                % store previous global best function value
                output.previous_global_best_individual = output.global_best_individual;
                output.previous_global_best_funval = output.global_best_funval;
                output.previous_best_funcvalues    = output.best_funcvalues;
                output.previous_best_individuals   = output.best_individuals;
                
                % assign global best individuals and their function
                % values per algorithm
                
                [output.best_funcvalues, ind] = min(pop.fitnesses);
                output.best_individuals(1,:) = pop.individuals(ind, :);
                
                % save new global best individual and function value
                [min_func_val, index] = min(output.best_funcvalues);
                if (output.global_best_funval > min_func_val)
                    output.global_best_funval     = min_func_val;
                    output.global_best_individual = output.best_individuals(index, :);
                end
                
                % check convergence
                if ~converged
                    % per-algorithm convergence
                    if alg_conv
                        % update counter
                        if output.best_funcvalues(algorithm) < options.AchieveFunVal
                            if abs(output.previous_best_funcvalues(algorithm) - ...
                                    output.best_funcvalues(algorithm)) <= options.TolFun &&...
                                    all(abs(output.previous_best_individuals(algorithm) - ...
                                    output.best_individuals(algorithm))) <= options.TolX
                                counter = counter + 1;
                            else counter = 0;
                            end
                        end % if
                        
                        % if counter is larger than preset maximum,
                        % convergence has been achieved
                        if (counter > options.TolIters)
                            converged = true;
                        end
                        
                        % GODLIKE-convergence
                    else
                        % update counter
                        if output.global_best_funval < options.AchieveFunVal
                            if abs(output.previous_global_best_funval - ...
                                    output.global_best_funval) <= options.TolFun && ...
                                    all(abs(output.previous_global_best_individual - ...
                                    output.global_best_individual)) <= options.TolX
                                output.descent_counter = output.descent_counter + 1;
                            else output.descent_counter = 0;
                            end
                        end % if
                        
                        % if counter is larger than preset maximum, and the
                        % minimum amount of iterations has been performed,
                        % convergence has been achieved
                        if generation > options.MinIters && (output.descent_counter > 2)
                            converged = true;
                        end % if
                    end % if
                    
                    % finalize output
                    if converged && ~alg_conv
                        % insert the last population in the output
                        %                          for ii = 1:algorithms
                        %                              output.algorithm_info.(which_ones{ii})(sames(ii)).last_population = ...
                        %                                  pop{ii}.individuals;
                        %                              output.algorithm_info.(which_ones{ii})(sames(ii)).last_fitnesses = ...
                        %                                  pop{ii}.fitnesses;
                        %                          end
                        % insert the last population in the output
                        output.algorithm_info.(which_ones).last_population = ...
                            pop.individuals;
                        output.algorithm_info.(which_ones).last_fitnesses = ...
                            pop.fitnesses;
                        
                        
                        % finalize output structure
                        output.exitflag = 1;
                        output.message = sprintf(['Optimization terminated:\n\n',...
                            ' Coordinate differences were less than OPTIONS.TolX, and decrease\n',...
                            ' in function value was less than OPTIONS.TolFun for two consecutive\n',...
                            ' MODE-iterations. MODE algorithm converged without any problems.']);
                    end
                end % if
            end % if single
            
            % stuff specific for multi-objective optimization
            if multi
                
                % check convergence
                if ~converged
                    %                      keyboard
                    % see if the minimum amount of iterations has
                    % been performed yet
                    if generation > options.MinIters
                        % test if ALL populations are non-dominated
                        all_nd = false(algorithms, 1);
                        
                        all_nd = all(pop.pop_data.front_number == 0);
                        
                        % if we have not broken prematurely, all fronts are zero, and
                        % thus we have convergence
                        if all(all_nd), converged = true; end
                        
                        % finalize output structure
                        if converged
                            % complete output structure
                            output.exitflag = 1;
                            output.message = sprintf(['Optimization terminated:\n',...
                                'All trial solutions of all selected algorithms are non-dominated.\n',...
                                'MODE algorithm converged without any problems.']);
                        end % if (converged)
                    end % if (MinIters check)
                end % if ~converged
                
                % if converged, complete output structure
                if converged
                    % output complete Pareto front
                    
                    output.pareto_front_individuals = ...
                        [output.pareto_front_individuals; pop.individuals];
                    output.pareto_front_fitnesses = ...
                        [output.pareto_front_fitnesses; pop.fitnesses];
                    
                    % find most efficient point and its fitnesses
                    origin              = min(output.pareto_front_fitnesses);
                    shifted_fitnesses   = bsxfun(@minus, ...
                        output.pareto_front_fitnesses, origin);
                    distances_sq        = sum(shifted_fitnesses.^2,2);
                    [mindist_sq, index] = min(distances_sq);
                    output.most_efficient_point     = output.pareto_front_individuals(index, :);
                    output.most_efficient_fitnesses = output.pareto_front_fitnesses(index, :);
                    
                    %% 2013/10/16 MG modify to save geometry and fitness at each iterations
                    output.MatrixPop=cat(3,output.MatrixPop,pop.individuals);
                    output.MatrixFitness=cat(3,output.MatrixFitness,pop.fitnesses);
                    
                end % if converged
            end % if multi
        end % if
    end % nested function

% display the algorithm's progress
    function display_progress
%                 keyboard
        % if the algorithm is multistart, only print the header
        
        %TODO - display for MS
        if ~strcmpi(pop.algorithm, 'MS'), algorithm_index = j; end
        
        % Command window
        % ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·
        if strcmpi(options.display, 'on') || strcmpi(options.display, 'CommandWindow')
            
            % if not converged, display all relevant information
            % from the current iteration
            if ~converged
                % create counter string
                genstr = num2str(generation);
                if strcmp(genstr,'11')||strcmp(genstr,'12')||strcmp(genstr,'13')
                    counter_string = 'th';
                else
                    switch genstr(end)
                        case '1', counter_string = 'st';
                        case '2', counter_string = 'nd';
                        case '3', counter_string = 'rd';
                        otherwise, counter_string = 'th';
                    end
                end
                
                % display header if this is the first generation, first
                % algorithm and first algorithm iteration
                if (generation == 1) && (algorithm_index == 1)
                    
                    % 2013/11/04 MG assign to string for successive display
                    strings = '%s population.\n';
                    % output header
                    fprintf(1, ['\nMODE optimization ', strings], which_ones);
                    
                    % display single or multi-objective optimization, and
                    % population size, iterations low and high
                    if     single
                        fprintf(1,...
                            ['Performing single-objective optimization, with total population\n'...
                            'size of %d individuals. Lower bounds on algorithm iterations\n', ...
                            'is %d, upper bound is %d.\n'], popsize, options.MODE.ItersLb, ...
                            options.MODE.ItersUb);
                    elseif multi
                        fprintf(1,...
                            ['Performing multi-objective optimization, with %d objectives.\n',...
                            'Total population size is %d individuals. Lower bounds on\nalgorithm ',...
                            'iterations is %d, upper bound is %d.\n'], options.num_objectives,...
                            popsize, options.MODE.ItersLb, options.MODE.ItersUb);
                    end % if
                end % if
                
                % subsequent iterations
                % check if this is a new iteration
                
                %                         fprintf(1,...
                %                            ['\n==============================================================\n',...
                %                             '                         GENERATION %d\n'],generation);
                if single
                    fprintf(1, ...
                        '          Current global minimum: %+1.8e\n',...
                        output.global_best_funval);
                end
                %                         fprintf(1, ...
                %                             '==============================================================\n');
                
                % display new algorithm's header
                fprintf(1,...
                    ['· · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · ·\n',...
                    '                  %s algorithm, %d%s generation\n',...
                    '  popsize: %d, max.generations: %d, max.functioncalls: %d\n'],...
                    which_ones, generation, counter_string, ...
                    frac_popsize, frac_iterations, options.MaxFunEvals);
                if multi
                    fprintf(1, ...
                        '  #GEN  f.count      Pareto fronts       non-Pareto fronts\n');
                elseif single
                    fprintf(1, '  #  f.count       min(F)        std(F)         descent\n');
                end % if
                fprintf(1, '· · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · ·\n')
                
                
                if multi
                    fprintf(1, '%3d   %6d    %10d             %10d\n', ...
                        algorithm_index, pop.funevals, ...
                        nnz(pop.pop_data.front_number==0),...
                        nnz(pop.pop_data.front_number~=0))
                elseif single
                    fprintf(1, '%3d   %6d    %+1.5e  %+1.5e  %+1.5e\n',...
                        algorithm_index, pop.funevals, ...
                        min(pop.fitnesses),std(pop.fitnesses),...
                        output.previous_best_funcvalues -...
                        output.best_funcvalues)
                end % if
                
                % if we do have convergence, just display the output message
            else
                fprintf(1, '\n'), fprintf(1, output.message), fprintf(1, '\n\n'),
            end % if
        end % if (commandwindow)
        
        % Plot
        % ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·  ·
%                 keyboard
        if strcmpi(options.dispPLOT, 'on')
            
            % Check problem dimensionality (can not be larger than 2)
            if single && pop.dimensions > 2, return, end
            
            % Check number of objectives (can not be larger than 3)
            if ~single && pop.num_objectives > 3, return, end
            
            % initialize some stuff
            % (maximum of 16 algorithms can be displayed)
            clf, minfval = [];
            colors = {'r.';'b.';'g.';'k.';
                'ro';'bo';'go';'ko';
                'rx';'bx';'gx';'kx';
                'r^';'b^';'g^';'k^';};

            
            %% 2014/02/05 MG single optimization at the moment present no modification in the plot.
            % single-objective
            if single
                fvals = pop.fitnesses;
                % overall minimum and maximum function values
                if isempty(minfval)
                    minfval = min(fvals(:));  maxfval = max(fvals(:));
                    if ~isfinite(minfval), minfval = -1e-100; end
                    if ~isfinite(maxfval), maxfval = 1e100; end
                else
                    if min(fvals(:)) < minfval, minfval = min(fvals(:)); end
                    if max(fvals(:)) > maxfval, maxfval = max(fvals(:)); end
                end
                
                % also extract individuals
                inds = pop.individuals;
                
                % plot the variables versus their function value
                if (size(inds,2) == 1)      % one dimensional
                    plot(inds, fvals, colors{1});
                elseif (size(inds,2) == 2)  % two dimensional
                    plot3(inds(:, 1), inds(:, 2), fvals, colors{1});
                end % if
            end % if single
            
            % multi-objective
            if multi
                global ii kk pos;
                % extract function values
                fvals=output.MatrixFitness;
                output.last_case_to_plot=size(fvals,3);
                output.case_to_plot=round(linspace(1,output.last_case_to_plot,4));
                pos=find(output.case_to_plot==output.last_case_to_plot);
                output.case_to_plot(pos(2:end))=[];
                % plot the function values against each other
                if (size(fvals,2) == 2)     % two objectives
                    clf;
                    kk=1;
                    legend_entries={};
                    for ii=output.case_to_plot
                        plot(fvals(:, 1,ii), fvals(:, 2,ii), colors{kk});hold on;
                        
                        legend_entries{kk}=['Iteration ',num2str(ii),' of ',num2str(output.last_case_to_plot)];
                        kk=kk+1;
                    end
                    clear ii kk;
                    
                elseif (size(fvals,2) == 3) % three objectives
                    plot3(fvals(:, 1), fvals(:, 2), fvals(:, 3), colors{1});
                end % if
            end % if
            
            %% 2014/02/05 MG legend entries valid only for single optimization
            if single && ~converged
                legend_entries = ...
                    ['DE \bf{(evaluating)}'];
            end
            
            % make plot
            if single
                % plot & axes
                if     (size(inds,2) == 1) % one-dimensional
                    xlabel('Decision variable x'), ylabel('Function value F(x)')
                    axis([lb(1) ub(1) minfval-1e-100 maxfval+1e-100])
                    if converged
                        plot(output.global_best_individual, output.global_best_funval, 'ko',...
                            'MarkerFaceColor', 'g', 'MarkerSize', 10)
                    end
                elseif (size(inds,2) == 2) % two-dimensional
                    xlabel('Decision variable x_1'), ylabel('Decision variable x_2')
                    zlabel('Function value F(x)'), view(30,50)
                    axis([lb(1) ub(1) lb(2) ub(2) minfval-1e-16 maxfval+1e-16])
                    if converged
                        plot3(output.global_best_individual(1),output.global_best_individual(2),...
                            output.global_best_funval, 'ko','MarkerFaceColor', 'g', ...
                            'MarkerSize', 10), view(30,50)
                    end
                end
                % make nice title
                if ~converged
                    title({'Current population versus objective function';
                        ['Generation ', num2str(generation),...
                        ', Function evaluations = ', num2str(num_funevaluations)]});
                else
                    % adjust legend
                    legend_entries{end+1} = 'global optimum';
                    % create the title
                    % MG change in title black is the last cases analysed
                    title({'CONVERGED population versus objective function';
                        ['Generation ', num2str(generation),...
                        ', Function evaluations = ', num2str(num_funevaluations)];
                        '("Global" optimum is the black marker)';});
                end
            elseif multi
                % plot & axes
                if     (size(fvals,2) == 2)     % two objectives
                    xlabel('F_1(x)'), ylabel('F_2(x)'); grid on;
                    %% 2014/02/05 MG plot unecessary becouse the last solution is just take into account...
%                     if converged
%                         plot(output.most_efficient_fitnesses(1),...
%                             output.most_efficient_fitnesses(2), 'ko','MarkerFaceColor', 'g',...
%                             'MarkerSize', 10)
%                     end
                elseif (size(fvals,2) == 3)     % three objectives
                    xlabel('F_1(x)'), ylabel('F_2(x)'), zlabel('F_3(x)'), view(30,50)
                    if converged
                        plot3(output.most_efficient_fitnesses(1),...
                            output.most_efficient_fitnesses(2),output.most_efficient_fitnesses(3),...
                            'ko','MarkerFaceColor','g','MarkerSize',10), view(30,50)
                    end
                end
                % make nice title
                if ~ converged
                    title({'Current Pareto Front'; ['Generation ', num2str(generation),...
                        ', Function evaluations = ', num2str(num_funevaluations)]});
                else
                    % adjust legend
                    % FC legend_entries{end+1} = 'most efficient';
                    % make nice title
                    % MG change in title black is the last cases analysed
                    title({'Final Pareto Front'; ['Generation ', num2str(generation),...
                        ', Function evaluations = ', num2str(num_funevaluations)];
                        '(black marker is the most efficient point)'});
                end
            end
            
            % draw legend
            %keyboard
            legend(legend_entries);
            
            % do not delay plotting
            drawnow
            
        end % if
    end % nested function

end
