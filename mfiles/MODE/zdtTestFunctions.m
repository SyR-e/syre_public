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

function [y, cons] = zdtTestFunctions(x,type)
switch type
    case 1
        [y,cons] = TP_ZDT1_objfun(x);
        return
    case  2
        [y,cons] = TP_ZDT2_objfun(x);
        return
    case 3
        [y,cons] = TP_ZDT3_objfun(x);
        return
    case 4
        [y,cons] = TP_ZDT4_objfun(x);
        return
    case 6
        [y,cons] = TP_ZDT6_objfun(x);
        return
    otherwise
        [y,cons] = TP_ZDT1_objfun(x);
        return
end

    function [y, cons] = TP_ZDT1_objfun(x)
        % Objective function : Test problem 'ZDT1'.
        %*************************************************************************
        y = [0, 0];
        cons = [];
        numVar = length(x);
        g = 1 + 9*sum(x(2:numVar))/(numVar-1);
        y(1) = x(1);
        y(2) = g*(1-sqrt(x(1)/g));
    end

    function [y, cons] = TP_ZDT2_objfun(x)
        % Objective function : Test problem 'ZDT2'.
        %*************************************************************************
        y = [0, 0];
        cons = [];
        numVar = length(x);
        g = 1 + 9*sum(x(2:numVar))/(numVar-1);
        y(1) = x(1);
        y(2) = g * ( 1-(x(1)/g)^2);
    end

    function [y, cons] = TP_ZDT3_objfun(x)
        % Objective function : Test problem 'ZDT3'.
        %*************************************************************************
        y = [0, 0];
        cons = [];
        numVar = length(x);
        g = 1 + 9*sum(x(2:numVar))/(numVar-1);
        y(1) = x(1);
        y(2) = g * ( 1-sqrt(x(1)/g) - x(1)/g*sin(10*pi*x(1)) );
    end
    function [y, cons] = TP_ZDT4_objfun(x)
        % Objective function : Test problem 'ZDT4'.
        %*************************************************************************
        y = [0, 0];
        cons = [];
        numVar = length(x);
        sum = 0;
        for i = 2:numVar
            sum = sum + x(i).^2 - 10 * cos(4*pi*x(i));
        end
        g = 1 + 10*(numVar-1)+sum;
        y(1) = x(1);
        y(2) = g * (1-sqrt(x(1)/g));
    end
    function [y, cons] = TP_ZDT6_objfun(x)
        % Objective function : Test problem 'ZDT6'.
        %*************************************************************************
        y = [0, 0];
        cons = [];
        numVar = length(x);
        g = 1 + 9 * (sum(x(2:numVar))/(numVar-1))^0.25;     
        y(1) = 1 - exp(-4*x(1)) * sin(6*pi*x(1))^6;
        y(2) = g * (1 - (y(1)/g)^2);
    end
end