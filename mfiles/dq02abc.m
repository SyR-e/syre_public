% Copyright 2025
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

function i123 = dq02abc(idq0,theta)

theta_i = 0;    % was for multi 3 phase

T32 = 2/3 * [cos(theta_i)     cos(theta_i+2*pi/3)    cos(theta_i-2*pi/3)
             sin(theta_i)     sin(theta_i+2*pi/3)    sin(theta_i-2*pi/3);
             1/2              1/2                    1/2];
% Tt0 = [cos(theta) sin(theta)  0;
%        -sin(theta) cos(theta) 0;
%        0           0          1];

% idq0 to iab0
% iab0 = inv(Tt0)*idq0;
iab = (idq0(1,:) + j * idq0(2,:)) .* exp(j*theta);
iab0 = [real(iab);imag(iab);idq0(3,:)];

% iab0 to i123
i123 = inv(T32)*iab0;


