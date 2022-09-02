% Copyright 2021
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

function [map] = SDE_scaleL(map,kL)

map.fd    = map.fd*kL;
map.fq    = map.fq*kL;
map.Lbase = map.Lbase*kL;
map.Lmd   = map.Lbase*kL;
map.Rs    = map.Rs*kL;
map.T     = map.T*kL;
map.mPM   = map.mPM*kL;
map.mCu   = map.mCu*kL;

map.dataSet.StackLength = map.dataSet.StackLength*kL;

