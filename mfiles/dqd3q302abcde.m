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

function i12345 = dqd3q302abcde(idqd3q30,theta)

theta_i = theta;    %

% dq -> alpha beta
T5i = [cos(theta_i)         sin(theta_i)            cos(3*theta_i)              sin(3*theta_i)             1/sqrt(2)
    cos(theta_i-2*pi/5)     sin(theta_i-2*pi/5)     cos(3*(theta_i-2*pi/5))     sin(3*(theta_i-2*pi/5))    1/sqrt(2)
    cos(theta_i-4*pi/5)     sin(theta_i-4*pi/5)     cos(3*(theta_i-4*pi/5))     sin(3*(theta_i-4*pi/5))    1/sqrt(2)
    cos(theta_i+4*pi/5)     sin(theta_i+4*pi/5)     cos(3*(theta_i+4*pi/5))     sin(3*(theta_i+4*pi/5))    1/sqrt(2)
    cos(theta_i+2*pi/5)     sin(theta_i+2*pi/5)     cos(3*(theta_i+2*pi/5))     sin(3*(theta_i+2*pi/5))    1/sqrt(2)
    ];

% alpha beta -> abc
i12345 = T5i * idqd3q30;

