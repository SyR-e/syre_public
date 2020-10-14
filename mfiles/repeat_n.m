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

function s_out = repeat_n(s_in,n)

s_out = s_in;

for j = 1:(n-1)
    if size(s_in,1) == 1
        s_out = [s_out s_in(2:end)];
    else
        s_out = [s_out s_in(:,2:end)];
    end
end
