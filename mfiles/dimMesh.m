% Copyright 2014
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

function fem = dimMesh(geo,eval_type)

% Define mesh size
% - thicker if in multi-objective optimization (according to K_mesh_MOOA)
% - thinner if in single machine evaluation (post-proc, according to K_mesh)

res_max=1/3*geo.g;

if strcmp(eval_type,'MO_GA')||strcmp(eval_type,'MO_OA')
    fem.res_traf=1/geo.p;
    if geo.q < 0.5
        fem.res_traf = fem.res_traf * 3;
    end
%     fem.res=geo.K_mesh_MOOA*fem.res_traf;
else
    fem.res_traf=1/geo.p * geo.mesh_K/geo.mesh_K_MOOA;
end

fem.res=geo.mesh_K_MOOA*fem.res_traf;

