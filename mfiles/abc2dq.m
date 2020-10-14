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

function idq = abc2dq(i1,i2,i3,theta)

theta_i = 0;    % was for multi 3 phase

% abc -> alpha beta
T32 = 2/3 * [cos(theta_i)     cos(theta_i+2*pi/3)    cos(theta_i-2*pi/3)
    sin(theta_i)     sin(theta_i+2*pi/3)    sin(theta_i-2*pi/3)];
iab = T32 * [i1;i2;i3];

% alpha beta --> dq
temp = (iab(1,:) + j * iab(2,:)) .* exp(-j*theta);
idq = [real(temp);imag(temp)];


