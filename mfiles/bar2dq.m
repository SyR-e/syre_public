% Copyright 2022
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

% trasformazione abc -> dq
% fbar righe: k barra
% fbar colonne: tempo o posizione

function idq = bar2dq(fbar,theta,N)


% zero is on bar1
ph = (0:N-1) * 360/N;
cos_ph = cosd(ph);
sin_ph = sind(ph);

TN2 = 2/N * [cos_ph ; sin_ph];

% 123 -> alpha beta
iab = TN2 * fbar;
% dq -> alpha beta
temp = (iab(1) + j * iab(2)) * exp(-j*theta);

idq = [real(temp) imag(temp)];
