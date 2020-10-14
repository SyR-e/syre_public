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

%restituisce i parametri a b c di una retta passante per 2 pti:
function [a,b,c]=retta_per_2pti(x1,y1,x2,y2)
a=y1-y2;
b=x2-x1;
c=x1.*(y2-y1)-y1.*(x2-x1);