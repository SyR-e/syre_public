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

function f  = cutPopulation(intermediatePopulation, objNumber, varNumber, populationSize)

% This function replaces the chromosomes based on rank and crowding
% distance. Initially until the population size is reached each front is
% added one by one until addition of a complete front which results in
% exceeding the population size. At this point the chromosomes in that
% front is added subsequently to the population based on crowding distance.

[N] = size(intermediatePopulation);

% Get the index for the population sort based on the rank
[~,index] = sort(intermediatePopulation(:,objNumber + varNumber + 1));

% Now sort the individuals based on the index
for i = 1 : N
    sortedPopulation(i,:) = intermediatePopulation(index(i),:); %#ok<AGROW>
end

% Find the maximum rank in the current population
max_rank = max(intermediatePopulation(:,objNumber + varNumber + 1));

% Start adding each front based on rank and crowing distance until the
% whole population is filled.
previousIndex = 0;
for i = 1 : max_rank
    % Get the index for current rank i.e the last the last element in the
    % sorted_chromosome with rank i. 
    currentIndex = find(sortedPopulation(:,objNumber + varNumber + 1) == i, 1, 'last' );
    % Check to see if the population is filled if all the individuals with
    % rank i is added to the population. 
    if currentIndex > populationSize
        % If so then find the number of individuals with in with current
        % rank i.
        remaining = populationSize - previousIndex;
        % Get information about the individuals in the current rank i.
        tempPop = ...
            sortedPopulation(previousIndex + 1 : currentIndex, :);
        % Sort the individuals with rank i in the descending order based on
        % the crowding distance.
        [~,tempIndex] = ...
            sort(tempPop(:, objNumber + varNumber + 2),'descend');
        % Start filling individuals into the population in descending order
        % until the population is filled.
        for j = 1 : remaining
            f(previousIndex + j,:) = tempPop(tempIndex(j),:); %#ok<AGROW>
        end
        return;
    elseif currentIndex < populationSize
        % Add all the individuals with rank i into the population.
        f(previousIndex + 1 : currentIndex, :) = ...
            sortedPopulation(previousIndex + 1 : currentIndex, :);
    else
        % Add all the individuals with rank i into the population.
        f(previousIndex + 1 : currentIndex, :) = ...
            sortedPopulation(previousIndex + 1 : currentIndex, :);
        return;
    end
    % Get the index for the last added individual.
    previousIndex = currentIndex;
end
