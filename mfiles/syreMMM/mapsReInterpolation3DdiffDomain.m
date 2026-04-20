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

function [map] = mapsReInterpolation3DdiffDomain(map,xName,yName,zName,XYPoints,ZPoints, Xnew, Ynew, Znew, method, extrapolation)


if nargin() < 10 || isempty(method)
    method = 'linear';
    extrapolation = 'none';
end


if nargin() < 9 || isempty(method)
    extrapolation = 'none';
end


map0 = map;
names = fieldnames(map0);

x_vec0 = reshape(map0.(xName)(1,:,1), [], 1);
y_vec0 = reshape(map0.(yName)(:,1,1), [], 1); 
z_vec0 = reshape(map0.(zName)(1,1,:), [], 1);

[X0_ndgrid, Y0_ndgrid, Z0_ndgrid] = meshgrid(x_vec0, y_vec0, z_vec0);

xx = reshape(linspace(min(Xnew, [], 'all'),max(Xnew, [], 'all'),XYPoints),[],1);
yy = reshape(linspace(min(Ynew, [], 'all'),max(Ynew, [], 'all'),XYPoints),[],1);
zz = reshape(linspace(min(Znew, [], 'all'),max(Znew, [], 'all'),ZPoints),[],1);

[xx_nd,yy_nd,zz_nd] = meshgrid(xx,yy,zz);

map.(xName) = xx_nd;
map.(yName) = yy_nd;
map.(zName) = zz_nd;

for ii=1:length(names)
    current_field_name = names{ii};
    if isnumeric(map0.(current_field_name))
        if size(map.(current_field_name))==size(map0.(xName))
            map.(current_field_name) = interp3(X0_ndgrid,Y0_ndgrid,Z0_ndgrid,map0.(current_field_name),map.(xName),map.(yName),map.(zName),method);
        end
    end
end
