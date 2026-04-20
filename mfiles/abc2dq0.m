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

function idq0 = abc2dq0(i123,theta)

theta_i = 0;    % was for multi 3 phase

% abc -> alpha beta
T32 = 2/3 * [cos(theta_i)     cos(theta_i+2*pi/3)    cos(theta_i-2*pi/3)
             sin(theta_i)     sin(theta_i+2*pi/3)    sin(theta_i-2*pi/3);
             1/2              1/2                    1/2];
iab0 = T32 * i123;

% alpha beta --> dq
% Tt0 = [cos(theta) sin(theta)  0;
%        -sin(theta) cos(theta) 0;
%        0           0          1];
% idq0 = Tt0 * iab;

temp = (iab0(1,:) + j * iab0(2,:)) .* exp(-j*theta);
idq0 = [real(temp);imag(temp);iab0(3,:)];



