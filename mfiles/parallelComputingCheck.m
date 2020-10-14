% Copyright 2019
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

function [ppState] = parallelComputingCheck(dispFlag)
% 
% [ppState] = parallelComputingCheck(dispFlag)
% 
% This function check the state of the parallel computing engine in Matlab.
% The output has the following meanings:
%   ppState == -1 --> parallel computing toolbox not available (not installed or not licenced)
%   ppState ==  0 --> parallel pool available and disabled
%   ppState >= +1 --> parallel pool available and enabled on ppState workers
% The input is a flag:
%   dispFlag == 0 --> no messages in the command window (default)
%   dispFlag == 1 --> display messages on the command window

if nargin()==0
    dispFlag=0;
end

ppState=-1;

test = ver();
for ii=1:length(test)
    if strcmp(test(ii).Name,'Parallel Computing Toolbox')
        ppState=0;
    end
end

if ppState==0
    if ~license('checkout','Distrib_Computing_Toolbox')
        ppState=-1;
    else
        if ~isempty(isprop(gcp('nocreate'),'NumWorkers'))
            tmp=gcp;
            ppState=tmp.NumWorkers;
        end
    end
end

if dispFlag
    switch ppState
        case -1
            disp('Parallel Computing Toolbox not available')
        case 0
            disp('Parallel pool not enabled')
        otherwise
            disp(['Parallel pool enabled on ' int2str(ppState) ' workers'])
    end
end

