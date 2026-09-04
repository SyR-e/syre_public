% Copyright 2026
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

function [optRes] = ORE_removeUnfeasible(optRes,flag)

if ~flag
    optRes.filt = ones(size(optRes.motNum));
else
    for ii=1:length(optRes.objLimits)
        optRes.filt(optRes.objParetoFront(:,ii)>optRes.objLimits(ii)) = NaN;
    end
end
