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

%retta tangente ad una circonferenza di centro alpha e beta ed una
%retta passante per Po=(x0,y0) e coeff_angolare=45°.
%retta nella forma:
function [a b c]=retta_tg_ad_una_circonferenza(alpha,beta,x0,y0)
a=x0-alpha;
b=y0-beta;
c=-(x0*(x0-alpha)+y0*(y0-beta));