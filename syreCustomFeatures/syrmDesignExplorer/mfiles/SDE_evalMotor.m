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

function SDE_evalMotor(map)

x = map.xSelect;
b = map.bSelect;

dataAvailable = map.dataAvailable;

disp(['-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-'])
for ii=1:length(dataAvailable)
    switch dataAvailable{ii}
        case {'hc_pu','NsI0_hc','hc','dx','geo','Br','xRaw','bRaw','dataSet','PMdim','sk','mechStressRad','mechStressTan','kmechrad','kmechtan'}
            flagDisp=0;

        otherwise
            flagDisp=1;
    end
    if flagDisp
        command = ['tmp = interp2(map.xx,map.bb,map.' dataAvailable{ii} '.*ones(size(map.xx)),x,b);'];
        eval(command);
        disp([dataAvailable{ii} ' = ' num2str(tmp)]);
    end
end
disp(['-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-'])
