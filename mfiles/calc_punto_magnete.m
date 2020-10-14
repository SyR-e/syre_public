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

% calc_punto_magnete.m
% calcola l'intersezione tra circonferenza centrata in x0,0 di raggio r1 e
% la retta che delimita il magnete, con fase gammam
% input: raggi r ed r1
% output: x,y del punto

function [x,y] = calc_punto_magnete(r1,gammam,x0)

x = x0 - r1 * cos(gammam);
y = r1 * sin(gammam);