% Copyright 2024
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

function [map] = mapsReInterpolation(map,xName,yName,nPoints,method)

if nargin()==4
    method = 'linear';
end

map0 = map;
names = fieldnames(map0);
xx0 = map.(xName);
yy0 = map.(yName);
xx = linspace(min(xx0(:)),max(xx0(:)),nPoints);
yy = linspace(min(yy0(:)),max(yy0(:)),nPoints);
[xx,yy] = meshgrid(xx,yy);
map.(xName) = xx;
map.(yName) = yy;

for ii=1:length(names)
    if isnumeric(map0.(names{ii}))
        if size(map0.(names{ii}))==size(map0.(xName))
            map.(names{ii}) = interp2(map0.(xName),map0.(yName),map0.(names{ii}),map.(xName),map.(yName),method);
        else
            map.(names{ii}) = map0.(names{ii});
        end
    end
end
