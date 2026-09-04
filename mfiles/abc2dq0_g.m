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

function idq0 = abc2dq0_g(n3ph, ia, ib, ic,theta,delta)

theta_i = theta;    %

for i=1:n3ph
    i123(3*i-2) = ia(i);
    i123(3*i-1) = ib(i);
    i123(3*i) = ic(i);
end

% abcdef -> dqdq00
% delta = 5*0.5*pi/3; % 0 D3P, pi/6 A6P, pi/3 S6P
alpha = 2*pi/3;

if(n3ph == 2)
    
    T3n = (2/3)*[cos(theta_i) cos(theta_i-alpha) cos(theta_i-2*alpha) cos(theta_i-delta) cos(theta_i-alpha-delta) cos(theta_i-2*alpha-delta);
            -sin(theta_i) -sin(theta_i-alpha) -sin(theta_i-2*alpha) -sin(theta_i-delta) -sin(theta_i-alpha-delta) -sin(theta_i-2*alpha-delta);
            cos(theta_i) cos(theta_i-alpha) cos(theta_i-2*alpha) -cos(theta_i-delta) -cos(theta_i-alpha-delta) -cos(theta_i-2*alpha-delta);
            -sin(theta_i) -sin(theta_i-alpha) -sin(theta_i-2*alpha) sin(theta_i-delta) sin(theta_i-alpha-delta) sin(theta_i-2*alpha-delta);        
            1 1 1 0 0 0;
            0 0 0 1 1 1];
else
    T3n = (2/3)*[cos(theta_i) cos(theta_i-alpha) cos(theta_i-2*alpha) cos(theta_i-delta) cos(theta_i-alpha-delta) cos(theta_i-2*alpha-delta) cos(theta_i-2*delta) cos(theta_i-alpha-2*delta) cos(theta_i-2*alpha-2*delta);
            -sin(theta_i) -sin(theta_i-alpha) -sin(theta_i-2*alpha) -sin(theta_i-delta) -sin(theta_i-alpha-delta) -sin(theta_i-2*alpha-delta) -sin(theta_i-2*delta) -sin(theta_i-alpha-2*delta) -sin(theta_i-2*alpha-2*delta);
            cos(theta_i) cos(theta_i-alpha) cos(theta_i-2*alpha) -cos(theta_i-delta) -cos(theta_i-alpha-delta) -cos(theta_i-2*alpha-delta) cos(theta_i-2*delta) cos(theta_i-alpha-2*delta) cos(theta_i-2*alpha-2*delta);
            -sin(theta_i) -sin(theta_i-alpha) -sin(theta_i-2*alpha) sin(theta_i-delta) sin(theta_i-alpha-delta) sin(theta_i-2*alpha-delta) -sin(theta_i-2*delta) -sin(theta_i-alpha-2*delta) -sin(theta_i-2*alpha-2*delta); 
            cos(theta_i) cos(theta_i-alpha) cos(theta_i-2*alpha) cos(theta_i-delta) cos(theta_i-alpha-delta) cos(theta_i-2*alpha-delta) -cos(theta_i-2*delta) -cos(theta_i-alpha-2*delta) -cos(theta_i-2*alpha-2*delta);
            -sin(theta_i) -sin(theta_i-alpha) -sin(theta_i-2*alpha) -sin(theta_i-delta) -sin(theta_i-alpha-delta) -sin(theta_i-2*alpha-delta) sin(theta_i-2*delta) sin(theta_i-alpha-2*delta) sin(theta_i-2*alpha-2*delta);        
            1 1 1 0 0 0 0 0 0;
            0 0 0 1 1 1 0 0 0;
            0 0 0 0 0 0 1 1 1];
end

% alpha beta -> abc
idq0 = T3n * i123';



