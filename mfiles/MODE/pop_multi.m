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

classdef pop_multi < pop_single
% POP_MULTI         Class definition for a population to be
%                   used for multi-objective optimization 
%
% POP_MULTI is a SubClass of POP_SINGLE. The class 
% constructor works in the same way as that of POP_SINGLE,
% with the exception that an additional property is set:
%
%   pop.num_objectives      (number of objectives)
%
% All inputs and other properties are the same as for 
% POP_SINGLE -- type 'help pop_single' for more information. 
%
% The method ITERATE is now suited to optimize multi-
% objective problems. To that end, several other (hidden) 
% methods have been implemented: 
%
%   NON_DOMINATED_SORT()
%
%       a general implementation of NSGA-II's non-dominated 
%       sorting procedure. Sorts the current population
%       according to the domination level of the individuals.
%
%
%   pool = TOURNAMENT_SELECTION(pool_size, tournament_size) 
%
%       a general tournament selection procedure, that takes
%       [tournament_size] individuals randomly selected from
%       the offspring population and lets them compete with
%       the rankings and crowding distances as determining 
%       factors. The winning individual of each tournament
%       is inserted into [pool] until that [pool] contains
%       [pool_size] individuals. 
%
%   UPDATE_ALGORITHMS()
%    
%       Called from NON_DOMINATED_SORT(), updates some 
%       globaly changing values for the different algorithms. 
%       In pop_single, this is done in REPLACE_PARENTS(), but
%       as that step is not executed here, an extra method is
%       required. This updates for instance the [lbest],
%       [nbest] and [gbest] for PSO, and the temperature for
%       ASA. 
%
%
% See also pop_single, GODLIKE. 

%      Author : Rody P.S. Oldenhuis
% Affiliation : Delft University of Technology
%               Faculty of Aerospace Engineering
%               Dep. of Astrodynamics & Satellite Systems 
%     Contact : oldnhuis@dds.nl
%   Licensing/
%    (C) info : Frankly I don't care what you do with it, 
%               as long as I get some credit when you copy 
%               large portions of the code ^_^

% last edited 28/07/2009
    
    % additional properties are also public
    properties
        num_objectives     % number of objectives
        % contents of pop_data for multi-objective
        % optimization:
        %      pop_data.parent_population
        %      pop_data.offspring_population
        %      pop_data.function_values_parent
        %      pop_data.function_values_offspring
        %      pop_data.front_number
        %      pop_data.crowding_distance        
    end
    
    % public methods
    methods (Access = public)
                
        % simple constructor: create pop_single object, and 
        % just add the number of objectives
        function pop = pop_multi(varargin)
            % construct pop_single object
            pop = pop@pop_single(varargin{:});             
            % just set the amount of objectives
            pop.num_objectives = pop.options.num_objectives;              
        end % function (constructor)
        
        % perform one multi-objective iteration
        function iterate(pop)
            
            % one iteration is:
            pop.non_dominated_sort;      % non-dominated sort
            pool = ...               % binary tournament selection (half pop.size, not for PSO)
                  pop.tournament_selection(round(pop.size/2), 2);
                
%             keyboard
            pop.create_offspring(pool);  % create new offspring
            pop.evaluate_function;       % evaluate objective function(s)            
%             keyboard
            % increase number of iterations made
            pop.iterations = pop.iterations + 1; 
            
        end % function (iterate)
        
    end % public methods
    
    % protected/hidden methods
    methods (Access = protected, Hidden)
      
        % non-dominated sort, and crowding distance assignment
        function non_dominated_sort(pop)         
            
            % this function sorts the population according to non-domination.
            % At the same time, it will compute the crowding distance for every
            % individual. It will store these values for all individuals in the 
            % data structure [pop.pop_data]. 
            
            % NOTE: NON_DOMINATED_SORT is by far the most computationally costly
            % function of both POP_SINGLE and POP_MULTI. Any effort regarding
            % increasing efficiency should be focussed on this method. 
            
            % determine if this is the first call or a subsequent call, 
            % and create fitnesses and number of variables accordingly
            if ~isempty(pop.pop_data.function_values_offspring)  
                
                    inds = [pop.pop_data.parent_population;
                        pop.pop_data.offspring_population; ];
                    fits = [pop.pop_data.function_values_parent;
                        pop.pop_data.function_values_offspring];
                
                N    = 2*pop.size;
            else  
                inds = pop.individuals; 
                fits = pop.fitnesses;
                N    = pop.size;
            end
                             
            % pre-calculate some stuff for crowding distances
            crowding_dists = zeros(N, 1);              % initialize            
            [sorted_fitnesses, indices] = sort(fits);  % sort fitnesses            
            crowding_dists(indices(1, :))   = inf;     % always include boundaries
            crowding_dists(indices(end, :)) = inf;     % always include boundaries
                                         
            % count the number of solutions that dominate every other solution
            guy_is_dominated_by = inf(N, 1);            
            for guy = 1:N
                
                % extract its fitnesses
                this_guy = fits(guy, :);
                
                % This is where the largest portion of the computational cost is,
                %
                %                      -> !!BY FAR!! <-
                %
                % BSXFUN() has been used to compare everything in one go, which
                % is also much faster than using a nested for. Note that strictly 
                % speaking, TWO comparisons are required, since dominance is 
                % defined as
                %
                % (xi <= yi) for all i, AND (xi < yi) for at least one value of i
                %
                % I coded the two comparisons for completeness, but in practice 
                % two equal function values are almost never encountered, so I 
                % left it commented. To use the stricter version, just uncomment 
                % the second comparison (and the addition to the [dominated] array.
                less_eq   = bsxfun(@le, fits([1:guy-1, guy+1:N], :), this_guy);
                %less      = bsxfun(@lt, fits([1:guy-1, guy+1:N], :), this_guy);
                dominated = all(less_eq, 2);% & any(less, 2);
                
                % insert the domination count in dominatrix
                guy_is_dominated_by(guy) = nnz(dominated);
                
                % compute this guy's crowding distance
                for m = 1:pop.num_objectives                    
                    % current sorting indices
                    sort_inds = indices(:, m);                    
                    % find this guy's index
                    guys_ind = find(sort_inds == guy);                    
                    % compute crowding distance
                    if isfinite(crowding_dists(guy))
                        crowding_dists(guy) = crowding_dists(guy) + ...
                            (sorted_fitnesses(guys_ind+1, m)-sorted_fitnesses(guys_ind-1, m))/...
                            (sorted_fitnesses(end,m)-sorted_fitnesses(1,m));
                    else break
                    end % if
                end % for                
            end % for
%             keyboard
            % create new population
            new_pop = zeros(pop.size, pop.dimensions);
            new_fit = zeros(pop.size, pop.num_objectives);
            fronts  = zeros(pop.size, 1);
            not_filled = pop.size;
            for i = 0:max(guy_is_dominated_by)                
                % extract current front
                front = guy_is_dominated_by == i;                
                % number of entries
                entries = nnz(front);                
                % see if it still fits in the population 
                if entries <= not_filled                     
                    % if so, insert all individuals from this front in the population
                    % and keep track of their fitnesses 
                    
                    new_pop(pop.size-not_filled+1:pop.size-not_filled+entries, :) = inds(front, :);
                    new_fit(pop.size-not_filled+1:pop.size-not_filled+entries, :) = fits(front, :);
                    fronts (pop.size-not_filled+1:pop.size-not_filled+entries, 1) = i;
                    % adjust number of entries that have not yet been filled
                    not_filled = not_filled - entries;                  
                % if it does not fit, insert the remaining individuals based on 
                % their crowding distance
                else break 
                end % if                
            end % for
              
            % (for the first iteration, the WHOLE current population will fit 
            % in the new population. Skip that case)
            if (N ~= pop.size)          
%                 keyboard
                % sort last front encountered by crowding-distance
                front_inds = inds(front, :);       front_fits = fits(front, :);
                [sorted_front, indices] = sort(crowding_dists(front), 'descend');                
                % extract individuals & fitnesses in proper order
                sorted_inds = front_inds(indices, :); sorted_fits = front_fits(indices, :);                    
                % insert the remaining individuals in the new population 
                new_pop(pop.size-not_filled+1:end, :) = sorted_inds(1:not_filled, :);
                new_fit(pop.size-not_filled+1:end, :) = sorted_fits(1:not_filled, :);
                fronts (pop.size-not_filled+1:end, 1) = i; 
            end  

            % insert results in data structure                
            pop.pop_data.parent_population      = new_pop;            
            pop.pop_data.function_values_parent = new_fit;
            pop.pop_data.front_number           = fronts;
            pop.pop_data.crowding_distance      = crowding_dists;  
            
            % copy individuals and fitnesses to respective class properties
            pop.fitnesses   = pop.pop_data.function_values_parent;       
            pop.individuals = pop.pop_data.parent_population;  
            
            %  some algorithms need additional operations             
%             pop.update_algorithms;   % MG DE needs nothing so this line could be eliminated                    
            
        end % function (non-dominated sort)
        
        % tournament selection with crowding distances and rankings 
        % as competitive factors
        function pool = tournament_selection(pop, pool_size, tournament_size)
            
            % initialize mating pool
            pool = zeros(pool_size, 1);
%             keyboard
            % fill the mating pool
            for i = 1:pool_size    
                
                % select [tournament_size] individuals at random
                equal = true;
                while equal % make sure they're not equal
                    inds  = round( rand(tournament_size, 1)*(pop.size-1) + 1); % random indices
                    equal = numel(unique(inds)) ~= numel(inds);        % check if any are equal
                end
                
                % let them compete according to
                % (xj < yj) if (rank(xj) < rank(yj))
                %           or (rank(xj) == rank(yj) and distance(xj) > distance(yj)
                ranks     = pop.pop_data.front_number(inds);
                distances = pop.pop_data.crowding_distance(inds);
                for j = 1:tournament_size                    
                    % compare ranks
                    less_rank = ranks(j) < [ranks(1:j-1); ranks(j+1:end)];                                       
                    % rank is less than all others
                    if all(less_rank), best = inds(j); break; end                    
                    % compare distances and rank equality
                    more_dist = ranks(j) == [ranks(1:j-1); ranks(j+1:end)] &...
                                distances(j) >= [distances(1:j-1); distances(j+1:end)];                    
                    % rank is equal, but distance is less than all others
                    if any(~less_rank & more_dist), best = inds(j); break; end                        
                end % for
                
                % insert the index of the best one in the pool
                pool(i) = best;                
            
            end % for
%             keyboard
        end % function (tournament selection)        
         
    end % protected/hidden methods 
    
end % classdef
