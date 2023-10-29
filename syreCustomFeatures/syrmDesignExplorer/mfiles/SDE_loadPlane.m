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

function [map,filename] = SDE_loadPlane(pathname,filename)

if nargin()==1
    map = pathname;
    filename = 'Map not saved';
else
    if nargin()~=2
        [filename,pathname,~] = uigetfile([cd '\*.fig'],'Select xb design plane');
    end

    filename = erase(filename,'.fig');

    hfig = openfig([pathname filename],'invisible');
    map = get(hfig,'UserData');
    close(hfig);
end


if isfield(map,'dataAvailable')
    map = rmfield(map,'dataAvailable');
    map = rmfield(map,'dataSelect');
end

if isempty(map)
    error('User Data not present in the figure');
else

    filename = erase(filename,'.fig');
    
    map.xSelect = mean(unique(map.xx));
    map.bSelect = mean(unique(map.bb));
    map.Vdc = 565;
    if ~isfield(map,'Rs')
        map.Rs = map.kj.*(2*pi*map.geo.R/1000*map.geo.l/1000)./(map.geo.win.n3phase.*3/2.*map.i0.^2);
    end

    [map] = SDE_evalBaseSpeed(map);

    tmp = rmfield(map,'dataSet');
    tmp = rmfield(tmp,'xSelect');
    tmp = rmfield(tmp,'bSelect');
    tmp = rmfield(tmp,'Vdc');
    tmp = rmfield(tmp,'geo');
    
    tmp = orderfields(tmp);

    map.dataAvailable = fieldnames(tmp);
    map.dataSelect = [];
end
