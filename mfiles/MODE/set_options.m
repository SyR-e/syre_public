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

function options = set_options(varargin)
% SET_OPTIONS                 Set options for the various optimizers
%
% Usage:
%
%   options = set_options('option1', value1, 'option2', value2, ...)
%
%
%   SET_OPTIONS is an easy way to set all options for the global optimization
%   algorithms PSO, DE, GA, ASA in GODLIKE. All options, and their possible
%   values are:
%
%   ======================================================================
%   General Settings:
%   ======================================================================
%       Display : string, either 'off' (default), 'on' or 'CommandWindow',
%                 'Plot'. This option determines the type of display that
%                 is used to show the algorithm's progress. 'CommandWindow'
%                 (or simply 'on') will show relevant information in the
%                 command window, whereas 'Plot' will make a plot in every
%                 iteration of the current population. Note that 'Plot'
%                 will only work if the number of decision variables is 1
%                 or 2 in case of single-pbjective optimization, or between
%                 1 and 3 objectives for multi-objective optimization.
%                 Please note that using any other display setting than
%                 'off' can significantly slow down the optimization.
%   MaxFunEvals : positive scalar, defining the maximum number of
%                 allowable function evaluations. The default is 100,000.
%                 Note that every objective and constraint function
%                 evaluation will be counted as 1 function evaluation. For
%                 multi-objective optimization, each objective function
%                 will be counted.
%      MaxIters : positive scalar, defining the maximum number of
%                 iterations that can be performed. The default is 20.
%      MinIters : positive scalar. This option defines the minimum amount
%                 of iterations GODLIKE will perform. This is particularly
%                 useful in multi-objective problems with small population
%                 sizes, because this combination increases the probability 
%                 that GODLIKE reports convergence (all fronts are Pareto 
%                 fronts), while a Pareto front of much better quality is 
%                 obtained if some additional shuffles are performed. The
%                 default value is 2. 
%
%   ======================================================================
%   Options specific to the GODLIKE Algorithm:
%   ======================================================================
%       ItersLb : positive scalar. This sets the minimum number of
%                 iterations that will be spent in one of the selected
%                 heuristic optimizers, per GODLIKE iteration. The default
%                 value is 10.
%       ItersUb : positive scalar. This sets the maximum TOTAL amount of
%                 iterations that will be spent in all of the selected
%                 heuristic optimizers combined. The default value is 100.
%
%   ======================================================================
%   General Settings for Single-Objective Optimization:
%   ======================================================================  
%        TolIters: positive scalar. This option defines how many consecutive 
%                  iterations the convergence criteria must hold for each 
%                  individual algorithm, before that algorithm is said to 
%                  have converged. The default setting is 15 iterations. 
%           TolX : positive scalar. Convergence is assumed to be attained, 
%                  if the coordinate differences in all dimensions for a
%                  given amount of consecutive iterations is less than 
%                  [TolX]. This amount of iterations is [TolIters] for each 
%                  individual algorithm, and simply 2 for GODLIKE-iterations. 
%                  The default value is 1e-4.
%         TolFun : positive scalar. Convergence is said to have been 
%                  attained if the value of the objective function decreases 
%                  less than [TolFun] for a given amount of consecutive
%                  iterations. This amount of iterations is [TolIters] for 
%                  each individual algorithm, and simply 2 for the 
%                  GODLIKE-iterations. The default value is 1e-4.
%  AchieveFunVal : scalar. This value is used in conjunction with the
%                  [TolX] and [TolFun] settings. If set, the algorithm will 
%                  FIRST try to achieve this function value, BEFORE enabling
%                  the [TolX] and [TolFun] convergence criteria. By default, 
%                  it is switched off (equal to AchieveFunVal = inf).
%
%   ======================================================================
%   General Settings for Multi-Objective Optimization:
%   ======================================================================
%        SkipTest : If set to 'on', some initial tests that are performed on
%                   the objective and constraint functions. These tests
%                   automatically determine whether the function accepts
%                   vectorized input or not, and how many objectives the
%                   problem has. The default is 'on', but it may be switched
%                   'off'. In case it's switched 'off', the algorithm assumes
%                   all functions accept vectorized input, AND the number of
%                   objectives (the next option) has been given, AND the
%                   dimensionality of the problem is also given (two options
%                   down). The 'off'-switch will be ignored if either of these
%                   demands is not true.
%   NumObjectives : Positive scalar. Sets the number of objectives manually.
%                   When the objective function is a single function that
%                   returns multiple objectives, the algorithm has to first
%                   determine how many objectives there are. This takes some
%                   function evaluations, which may be skipped by setting this
%                   value manually.
%
%   ======================================================================
%   Options specific to the Differential Evolution algorithm:
%   ======================================================================
%          Flb : scalar. This value defines the lower bound for the range
%                from which the scaling parameter will be taken. The
%                default value is -1.5.
%          Fub : scalar. This value defines the upper bound for the range
%                from which the scaling parameter will be taken. The
%                default value is +1.5. These two values may be set equal
%                to each other, in which case the scaling parameter F is
%                simply a constant.
%   CrossConst : positive scalar. It defines the probability with which a
%                new trial individual will be inserted in the new
%                population. The default value is 0.95.
%

    
    % Last edited 30/Jul/2009
%     keyboard
    % create structure with default values if no input is given
    if (nargin == 0)
        
        % initialize
        options = struct;
        
        % general options
        options.display       = 'off';
        options.dispPLOT      = 'off';
        options.MaxFunEvals   = 1e5;
        options.MaxIters      = 20;
        options.MinIters      = 2;        
        options.TolIters      = 15;
        options.TolX          = 1e-4;
        options.TolFun        = 1e-4;
        options.AchieveFunVal = inf;
        
        % function evaluation
        options.num_objectives = 2;
        options.dimensions     = [];
        
        % Differential Evolution
        options.DE.Flb        = -1.5;
        options.DE.Fub        = 1.5;
        options.DE.CrossConst = 0.95;
                
        % GODLIKE
        options.MODE.ItersLb = options.MinIters;
        options.MODE.ItersUb = options.MaxIters;
        
        % finished
        return;
        
        % create structure with fields according to user input
    elseif (nargin > 0)
        
        % assign default values
        options = set_options;
        
%        2013/11/04 MG not necessary this control
        % errortrap
%         if (mod(nargin, 2) ~= 0)
%             error('Please provide values for all the options.')
%         end
%         keyboard
        % loop through all the inputs, and use an "if-else cancer" to
        % create the problem structure
        for i = 1:2:nargin
            option = varargin{i};
            value  = varargin{i+1};
            
            % if option is not recognized, continue to the next argument
            if ~isa(option, 'char')
                throwwarning(option, [], [], []);
                continue;
            end
            
            % General options
            if     strcmpi(option, 'Display')
                if ~ischar(value)
                    throwwarning('Display', 'char', value);
                    continue;
                end
                if     strcmpi(value, 'off')
                    options.display = [];
                elseif strcmpi(value, 'CommandWindow') || strcmpi(value, 'on')
                    options.display = 'CommandWindow';
                elseif strcmpi(value, 'Plot')
                    options.display = 'Plot';
                else
                    error('population:set_options:unknown_display_option',...
                        ['Unknown display type: ', '''', value, '''.'])
                end
            elseif strcmpi(option, 'MaxFunEvals')
                if ~isnumeric(value)
                    throwwarning('MaxFunEvals', 'double', value);
                    continue;
                end
                options.MaxFunEvals = value;
            elseif strcmpi(option, 'MaxIters')
                if ~isnumeric(value)
                    throwwarning('MaxIters', 'double', value);
                    continue;
                end
                options.MaxIters = value;
            elseif strcmpi(option, 'MinIters')
                if ~isnumeric(value)
                    throwwarning('MinIters', 'double', value);
                    continue;
                end
                options.MinIters = value;
            elseif strcmpi(option, 'TolIters')
                if ~isnumeric(value)
                    throwwarning('TolIters', 'double', value);
                    continue;
                end
                options.TolIters = value;
            elseif strcmpi(option, 'TolX')
                if ~isnumeric(value)
                    throwwarning('TolX', 'double', value);
                    continue;
                end
                options.TolX = value;
            elseif strcmpi(option, 'TolFun')
                if ~isnumeric(value)
                    throwwarning('TolFun', 'double', value);
                    continue;
                end
                options.TolFun = value;           
            elseif strcmpi(option, 'AchieveFunVal')
                if ~isnumeric(value)
                    throwwarning('AchieveFunVal', 'double', value);
                    continue;
                end
                options.AchieveFunVal = value;
            elseif strcmpi(option, 'dispPLOT')
                if ~ischar(value)
                    throwwarning('dispPLOT', 'char', value);
                    continue;
                end
                options.dispPLOT=value;
                
            % options specific to Differential Evolution
            elseif strcmpi(option, 'Flb')
                if ~isnumeric(value)
                    throwwarning('Flb', 'double', value);
                    continue;
                end
                options.DE.Flb = value;
            elseif strcmpi(option, 'Fub')
                if ~isnumeric(value)
                    throwwarning('Fub', 'double', value);
                    continue;
                end
                options.DE.Fub = value;
            elseif strcmpi(option, 'CrossConst')
                if ~isnumeric(value)
                    throwwarning('CrossConst', 'double', value);
                    continue;
                end
                options.DE.CrossConst = value;
                
            % options specific to MODE algorithm
            elseif strcmpi(option, 'ItersLb')
                if ~isnumeric(value)
                    throwwarning('algiters', 'double', value);
                    continue;
                end
                options.MODE.ItersLb = value;
            elseif strcmpi(option, 'ItersUb')
                if ~isnumeric(value)
                    throwwarning('ItersUb', 'double', value);
                    continue;
                end
                options.MODE.ItersUb = value;
                
                % General Settings
            elseif strcmpi(option, 'SkipTest')
                if ~isnumeric(value)
                    throwwarning('SkipTest', 'char', value);
                    continue;
                end
                if     strcmpi(value, 'on')
                    options.skip_function_test = true;
                elseif strcmpi(value, 'off')
                    options.skip_function_test = false;
                end
            elseif strcmpi(option, 'NumObjectives')
                if ~isnumeric(value)
                    throwwarning('NumObjectives', 'double', value);
                    continue;
                end
                options.num_objectives = value;
            else
                throwwarning(option);
            end % if
        end % for
    end % if
    
    % throw appropriate warning upon abuse of the function
    function throwwarning(option, required, given, varargin)%#ok
        
        % test type
        if nargin == 3
            provided = whos('given');
            provided = provided.class;
            warning('set_options:incorrectvalue', ...
                ['Incorrect class type given for option ''%s'';\n',...
                'required type is ''%s'', received ''%s''.\n',...
                'Using default value...'], option, required, provided);
            
        % unrecognized options will be ignored
        else
            warning('set_options:incorrectoption', ...
                'Unrecognized option, ''%s''. Ignoring it...', num2str(option))
        end % if
    end % nested function
end % function (set options)