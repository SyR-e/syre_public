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

% rot_point - rotate the point (xi,yi) by the angle "angle"

function [xo,yo] = rot_point(xi,yi,angle)

xo = real((xi + 1j*yi).* exp(1j*angle));
yo = imag((xi + 1j*yi).* exp(1j*angle));
