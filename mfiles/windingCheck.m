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


function [avv] = windingCheck(geo)

slot_layer_pos = geo.win.slot_layer_pos;    % side-by-side flag winding
avv = geo.win.avv;
[r,c] = size(avv);

if strcmp(slot_layer_pos,'side_by_side')
    for cc=1:c-1
        if ((abs(avv(1,cc))~=abs(avv(1,cc+1))) && (abs(avv(2,cc))==abs(avv(1,cc+1))))
            tmp = avv(:,cc);
            avv(1,cc) = tmp(2);
            avv(2,cc) = tmp(1);
        end
    end
end