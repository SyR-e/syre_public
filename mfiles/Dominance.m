% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function flag = Dominance(fitnessA,fitnessB,epsilon)

% Verify if a solution is Pareto dominant respect to the other.
% If epsilon is greater than zero, the epsilon-dominance is
% verified. For MINIMIZATION problems.
%
% Syntax:
%           flag = Dominance(fitnessA,fitnessB,epsilon)
%
% Inputs:
%           fitnessA = Nx1 vector containing fitness values for
%                      the solution A;
%
%           fitnessB = Nx1 vector containing fitness values for
%                      the solution B;
%
%           epsilon  = epsilon value. Must be greater or equal to zero.
%        
% Ouputs:
%           flag = equal to 0 if A and B are non-dominating, 1 if A
%                  dominates B and -1 if B dominates A.



dom_more = 0;
dom_eq   = 0;
dom_less = 0;

if (epsilon<0)
    error('Dominance: epsilon value must be greater or equal to zero')
end

[mA,nA] = size(fitnessA);
[mB,nB] = size(fitnessB);

if (mA>1)||(mB>1)
    error('Dominance: fitness values must be stored in a row vector.')
end

if (nA~=nA)
    error('Dominance: fitness values must be stored in a row vector of the same size.')
end

% Control accuracy using epsilon parameter
for k = 1 : nA
    if fitnessA(k)>0
        fitnessA(k) = 1/(1+epsilon)*fitnessA(k);
    else
        fitnessA(k) = (1+epsilon)*fitnessA(k);
    end
end

% Find how many components of fitnessA are improving components (dom_less), how many
% are depreciatory ones (dom_more) and how many are equivalents (dom_eq).
for k = 1 : nA
    if (fitnessA(k)<fitnessB(k))
        dom_less = dom_less+1;
    elseif (fitnessA(k)>fitnessB(k))
        dom_more = dom_more + 1;
    else
        dom_eq = dom_eq + 1;
    end
end

% Improved solution (there aren't depreciatory components and not all are equivalent)
if (dom_more==0)&&(dom_eq<nA) 
    flag = 1;
    
% Deprecated solution (there aren't depreciatory components and not all are equivalent)
elseif (dom_less==0)&&(dom_eq<nA) 
    flag = -1;
else
    flag = 0;
end