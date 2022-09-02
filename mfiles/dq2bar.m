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

% trasformazione dq -> abc

function ibar = dq2bar(id,iq,theta,N)


% zero is on bar1
ph = (0:N-1) * 360/N;
cos_ph = cosd(ph);
sin_ph = sind(ph);

T2N = [cos_ph' sin_ph'];

% dq -> alpha beta
iab = (id + j * iq) .* exp(j*theta);
% alpha beta -> 123
ibar = T2N * [real(iab);imag(iab)];
