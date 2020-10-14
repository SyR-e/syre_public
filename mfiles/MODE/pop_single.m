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

classdef pop_single < handle
% =insert documentation here=

%      Author : Rody P.S. Oldenhuis
% Affiliation : Delft University of Technology
%               Faculty of Aerospace Engineering
%               Dep. of Astrodynamics & Satellite Systems 
%     Contact : oldnhuis@dds.nl
%   Licensing/
%    (C) info : Frankly I don't care what you do with it, 
%               as long as I get some credit when you copy 
%               large portions of the code ^_^
    
    % all properties are public
    properties 
        algorithm          % type of optimization algorithm used 
        funfcn             % objective function(s)
        individuals        % members of the population
        fitnesses          % corresponding fitnesses        
        size               % population size
        lb                 % lower bounds
        ub                 % upper bounds
        orig_size          % original size of the input
        dimensions         % dimensions
        funevals   = 0;    % number of function evaluations made        
        iterations = 0;    % iterations so far performed
        options            % options structure (see function [set_options] for info)            
        pop_data           % structure to store intermediate data                           
        % contents for single-objective optimization:
        %      pop_data.parent_population
        %      pop_data.offspring_population
        %      pop_data.function_values_parent
        %      pop_data.function_values_offspring
    end
        
    % public methods
    methods (Access = public)
        
        % constructor
        function pop = pop_single(varargin)
            
            % default check
            error(nargchk(2, 7, nargin));
            
            % input is ( new [pop_data] structure, previous [population] object, options )
            % (subsequent call from GODLIKE)
            % = = = = = = = = = = = = = = = = = = = = = = = = = = 
            
            if (nargin == 3) 
                
                % assign new pop_data structure
                pop.pop_data = varargin{1};
                
                % simply copy previous object
                pop.funfcn     = varargin{2}.funfcn;       pop.iterations = varargin{2}.iterations;
                pop.algorithm  = varargin{2}.algorithm;    pop.lb         = varargin{2}.lb;
                pop.funevals   = varargin{2}.funevals;     pop.ub         = varargin{2}.ub;
                pop.dimensions = varargin{2}.dimensions;   pop.orig_size  = varargin{2}.orig_size;  
                
                % copy individuals and fitnesses
                pop.individuals = pop.pop_data.parent_population;
                pop.fitnesses   = pop.pop_data.function_values_parent;
                
                % size and options might have changed
                pop.size = size(pop.individuals, 1);%#ok
                pop.options = varargin{3}; 
                % replicate [ub] and [lb]
                pop.lb = repmat(pop.lb(1, :), pop.size, 1);    
                pop.ub = repmat(pop.ub(1, :), pop.size, 1);
                
                return
            end
            
            % input is ( funfcn, popsize, lb, ub, dimensions, options )
            % (initialization call from GODLIKE)
            % = = = = = = = = = = = = = = = = = = = = = = = = = = 
             
            % parse input
            % · · · · · · · · · · · · · · · · · · · · · · 
                        
            % assign input
            pop.funfcn  = varargin{1};   pop.ub         = varargin{4};
            pop.size    = varargin{2};   pop.orig_size  = varargin{5};
            pop.lb      = varargin{3};   pop.dimensions = varargin{6};
            pop.options = varargin{7};
            
            % cast funfcn to cell if necessary
            if ~iscell(pop.funfcn), pop.funfcn = {pop.funfcn}; end
                        
            % replicate [lb] and [ub] to facilitate things a bit
            % (and speed it up some more)
            pop.lb = repmat(pop.lb, pop.size, 1);   pop.ub = repmat(pop.ub, pop.size, 1);
            
            % set optimization algorithm
            pop.algorithm = pop.options.algorithm;
                             
            % Initialize population             
            % · · · · · · · · · · · · · · · · · · · · · · 
            
            % initialize population 
            pop.individuals = pop.lb + rand(pop.size, pop.dimensions) .* (pop.ub-pop.lb);
                                      
            % insert copy into info structure
            pop.pop_data.parent_population = pop.individuals;
                        
            % temporarily copy parents to offspring positions 
            pop.pop_data.function_values_offspring = [];
            pop.pop_data.offspring_population      = pop.individuals;
            
            % evaluate function for initial population (parents only)             
            pop.evaluate_function;   
                        
            % copy function values into fitnesses properties
            pop.fitnesses = pop.pop_data.function_values_offspring;
            pop.pop_data.function_values_parent = pop.fitnesses;
                        
            % delete entry again
            pop.pop_data.function_values_offspring = [];
                                  
        end % function (constructor)  
        
        % single iteration
        function iterate(pop, times, FE)
            % [times] and [FE] are only used for the MultiStart algorithm
                        
            % select proper candiadates 
            
                pool = 1:pop.size;       % whole population otherwise
             
            % create offspring
            if nargin == 1
                pop.create_offspring(pool);
            else
                pop.create_offspring(pool, times, FE);  
            end
            
            %% 2013/10/28 MG probabilmente le seguenti righe possono essere eliminate
%             % if the algorithm is MS, this is the only step
%             if strcmpi(pop.algorithm, 'MS')
%                 % adjust iterations
%                 pop.iterations = pop.iterations + times;
%                 % then return                
%                 return            
%             end
            
            % carefully evaluate objective function(s)
            try
                pop.evaluate_function;   
            catch userFcn_ME                
                pop_ME = MException('pop_single:function_doesnt_evaluate',...
                    'MODE cannot continue: failure during function evaluation.');
                userFcn_ME = addCause(userFcn_ME, pop_ME);
                rethrow(userFcn_ME);
            end
            
            % replace the parents
            pop.replace_parents;         
                        
            % increase number of iterations made
            pop.iterations = pop.iterations + 1;
                        
        end % function (single iteration)              
                
    end % methods
    
    % % protected/hidden methods
    methods (Access = protected, Hidden)
        
        % tournament selection (only for GA)
        function pool = tournament_selection(pop, pool_size, tournament_size)
            
            % initialize mating pool
            pool = zeros(pool_size, 1);
            
            % total number of competitors
            rnd_inds = zeros(pool_size*tournament_size,1);
            
            % create random indices outside the loop (faster)
            for i = 1:floor(pool_size*tournament_size/pop.size)
                offset = pop.size*(i-1);
                [dummy, rnds] = sort(rand(pop.size,1));
                rnd_inds(offset+1:min(end,offset+pop.size), :) = rnds(1:min(end,nnz(~rnd_inds)));
            end
            
            % fill the mating pool
            for i = 1:pool_size    
                
                % select [tournament_size] individuals at random
                inds = rnd_inds(1:tournament_size);
                rnd_inds = rnd_inds(tournament_size+1:end);
                
                % let them compete according to
                % (xj < yj) if fit(xj) < fit(yj)
                [best, ind] = min(pop.fitnesses(inds));                
                
                % insert the index of the best one in the pool
                pool(i) = inds(ind);                
            
            end % for
            
        end % function (tournament selection)        
        
        %% generate new generation
%         function create_offspring(pop, pool, times, FE)
        function create_offspring(pop, pool)

            
            % get the size of the pool
            pool_size = length(pool);  
            
            % rename some stuff 
            parent_pop = pop.individuals(pool, :);
            parent_fit = pop.fitnesses(pool, :);
            
            % initialize
            newpop = zeros(pop.size, pop.dimensions);           % empty new population            
            newfit = NaN(pop.size, pop.options.num_objectives); % placeholder for the sites to 
                                                                % evaluate the function                        
            
            % generate offspring with selected algorithm
                %%
                %% Differential Evolution 
                %%
%                 case 'DE'  % Differential Evolution 
                                        
                    % I love DE for its elegance and simplicity, and 
                    % yet powerful optimization qualities
                    
                    % rename some stuff   
                    Flb = pop.options.DE.Flb;
                    Fub = pop.options.DE.Fub;
                    crossconst = pop.options.DE.CrossConst; 
%                     keyboard
                    % Neoteric Differential Evolution                    
                    for i = 1:pop.size   
                        % random indices
                        base = round(rand*(pool_size-1))+1; % RANDI is slower
                        d1   = round(rand*(pool_size-1))+1;
                        d2   = round(rand*(pool_size-1))+1;
                        % d2 may not be equal to d1
                        while (d1 == d2), d2 = round(rand*(pool_size-1))+1; end
                        % DE operator
                        if rand < crossconst || round(rand*(pool_size-1))+1 == i;
                            % DE operator when rnd < Cr
                            F = rand*(Fub-Flb) + Flb;
                            newpop(i, :) = parent_pop(base,:) +  ...
                                F*(parent_pop(d1,:) - parent_pop(d2,:));
                        else
                            % insert random parent otherwise
                            rnd_ind      = round(rand*(pool_size-1))+1;
                            newpop(i, :) = parent_pop(rnd_ind, :);
                            newfit(i, :) = parent_fit(rnd_ind, :);
                        end                        
                    end % for
                                    
                        
            % check constraints and boundaries after 
            % offspring generation
            [newpop, newfit] = pop.honor_bounds(newpop, newfit);
                         
            % insert result into pop
            pop.pop_data.offspring_population      = newpop;
            pop.pop_data.function_values_offspring = newfit;
            
        end % function (create offspring)
        
        % selectively replace the parent population with members 
        % from the offspring population (single-objective optimization)
        function replace_parents(pop)
            
            % rename for clarity
            new_fits = pop.pop_data.function_values_offspring;
            new_inds = pop.pop_data.offspring_population;
                
                % Differential Evolution 
                % Differential Evolution                     
                    % DE and GA both use simple greedy replacement 
                    better_inds = new_fits < pop.fitnesses;
                    pop.pop_data.parent_population(better_inds, :) = new_inds(better_inds, :);
                    pop.pop_data.function_values_parent(better_inds, :) = new_fits(better_inds, :);                   
            
            % copy individuals and fitnesses to respective properties            
            pop.fitnesses   = pop.pop_data.function_values_parent;       
            pop.individuals = pop.pop_data.parent_population;   
            
        end % function
        
        % evaluate the objective function(s) correctly
        function evaluate_function(pop)
                   
            % NOTE: suited for both single and multi-objective optimization
            
            % find evaluation sites
            if isempty(pop.pop_data.function_values_offspring)
                sites = 1:pop.size; % only in pop-initialization
            else
                sites = ~isfinite(pop.pop_data.function_values_offspring(:, 1));
            end
            
            % first convert population to cell
            true_pop = reshape(pop.pop_data.offspring_population(sites, :).', ...
                [pop.orig_size,nnz(sites)]); 
            % NOTE: for-loop is faster than MAT2CELL
            cell_pop = cell(1,1,nnz(sites));
            for ii = 1:size(true_pop,3), cell_pop{1,1,ii} = true_pop(:,:,ii); end %#ok

            % then evaluate all functions with cellfun
%             for ii = 1:numel(pop.funfcn)
%                 pop.pop_data.function_values_offspring(sites, ii) = ...
%                   cellfun(pop.funfcn{ii}, cell_pop);
%             end
% keyboard
%% CICLO MODIFICATO DA FC GEN 12
            for ii = 1:numel(pop.funfcn)
                pizzamat=cell2mat(cellfun(pop.funfcn{ii}, cell_pop,'UniformOutput',false));
                pizzasq=squeeze(pizzamat(1,:,:));
                
                pop.pop_data.function_values_offspring(sites,:) = pizzasq';
            end
%% FINE DELLA MODIFICA            
            % update number of function evaluations
            pop.funevals = pop.funevals + ...
                nnz(sites)*size(pop.pop_data.function_values_offspring, 2);% #ok

        end % function 
           
        % check boundaries
        function [newpop, newfit] = honor_bounds(pop, newpop, newfit)
            
            % find violation sites
            outsiders1 = false; outsiders2 = false;
            if ~isempty(newpop)
                outsiders1 = newpop < pop.lb;
                outsiders2 = newpop > pop.ub;
            end
            
            reinit = pop.lb + rand(pop.size, pop.dimensions).*(pop.ub-pop.lb);
            if any(outsiders1(:) | outsiders2(:))
                newpop(outsiders1) = reinit(outsiders1);
                newpop(outsiders2) = reinit(outsiders2);
                % also remove any function values
                newfit(any(outsiders1,2), :) = NaN;
                newfit(any(outsiders2,2), :) = NaN;
            end
            %             end % if
        end % function
        
                
    end % methods (private)
end % classdef
