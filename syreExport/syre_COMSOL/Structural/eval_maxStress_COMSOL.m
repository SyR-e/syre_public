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

function [out]=eval_maxStress_COMSOL(Stress,mat)

stress_rot =[];
x_stress = [];
y_stress = [];

for kk = 1:size(Stress, 1)
    if Stress(kk, 3) > mat.Rotor.sigma_max
       x_stress = Stress(kk, 1); 
       y_stress = Stress(kk, 2);    
    end
    stress_rot = [stress_rot; x_stress, y_stress];
end

x_max_stress = [];
y_max_stress = [];

[~, idx] = max(Stress(:, 3));

x_max_stress = Stress(idx, 1);
y_max_stress = Stress(idx, 2);

out.stressrot   = stress_rot;
out.x_over      = x_stress;
out.y_over      = y_stress;
out.x_max       = x_max_stress;
out.y_max       = y_max_stress;
