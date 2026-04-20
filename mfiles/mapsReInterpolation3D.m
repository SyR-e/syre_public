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

function [map] = mapsReInterpolation3D(map,xName,yName,zName,nPoints,method)

if nargin() < 6 || isempty(method)
    method = 'linear';
end

map0 = map;
names = fieldnames(map0);

x_vec0 = map0.(xName)(1,:,1);
y_vec0 = map0.(yName)(:,1,1); 
z_vec0 = map0.(zName)(1,1,:);

x_vec0 = x_vec0(:);
y_vec0 = y_vec0(:);
z_vec0 = z_vec0(:);

[X0_ndgrid, Y0_ndgrid, Z0_ndgrid] = meshgrid(x_vec0, y_vec0, z_vec0);

xx0 = map.(xName);
yy0 = map.(yName);
zz0 = map.(zName);

xx = linspace(min(xx0(:)),max(xx0(:)),nPoints);
yy = linspace(min(yy0(:)),max(yy0(:)),nPoints);
zz = linspace(min(zz0(:)),max(zz0(:)),nPoints);

[yy_nd,xx_nd,zz_nd] = ndgrid(yy,xx,zz);
map.(xName) = xx_nd;
map.(yName) = yy_nd;
map.(zName) = zz_nd;

for ii=1:length(names)
    current_field_name = names{ii};
    if isnumeric(map0.(current_field_name))
        if size(map0.(current_field_name))==size(map0.(xName))
            % map.(current_field_name) = interp3(map0.(xName),map0.(yName),map0.(zName),map0.(current_field_name),map.(xName),map.(yName),map.(zName),method);
            map.(current_field_name) = interp3(X0_ndgrid,Y0_ndgrid,Z0_ndgrid,map0.(current_field_name),map.(xName),map.(yName),map.(zName),method);
        else
            map.(current_field_name) = map0.(current_field_name);
        end
    end
end
