% Copyright 2021
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

function [out] = SDE_evalMotor(map,flagOut)

if nargin()==1
    flagOut = 0;
end

if ~flagOut
    out = [];
end

x = map.xSelect;
b = map.bSelect;

if isfield('dataAvailable',map)
    dataAvailable = map.dataAvailable;
else
    dataAvailable = fieldnames(map);
end

if ~flagOut
    disp(['-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-'])
end

for ii=1:length(dataAvailable)
    switch dataAvailable{ii}
        case {'hc_pu','NsI0_hc','hc','dx','geo','Br','xRaw','bRaw','dataSet','PMdim','sk','mechStressRad','mechStressTan','kmechrad','kmechtan','dataAvailable','dataSelect'}
            flagDisp=0;

        otherwise
            flagDisp=1;
    end
    if flagDisp
        % command = ['tmp = interp2(map.xx,map.bb,map.' dataAvailable{ii} '.*ones(size(map.xx)),x,b);'];
        % eval(command);
        tmp = interp2(map.xx,map.bb,map.(dataAvailable{ii}).*ones(size(map.xx)),x,b);
        if ~flagOut
            disp([dataAvailable{ii} ' = ' num2str(tmp)]);
        else
            out.(dataAvailable{ii}) = tmp;
        end
    end
end

if ~flagOut
    disp(['-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-'])
end


if nargout()==0
    clear out
end

