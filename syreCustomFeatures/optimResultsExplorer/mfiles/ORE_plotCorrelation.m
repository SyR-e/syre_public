% Copyright 2026
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

function ORE_plotCorrelation(optRes)

hfig = figure();
screenPos=get(groot,'ScreenSize')/get(groot,'ScreenPixelsPerInch')*2.54; % cm
figSetting(screenPos(3)*0.9,screenPos(4)*0.9);

front = [optRes.varParetoFront optRes.objParetoFront];
designIdx = 1:1:length(optRes.varNames);
objIdx    = designIdx(end)+[1:1:length(optRes.objNames)];
designBounds = optRes.varBounds;
feasibilityLimits = optRes.objLimits;
feasibilityLimits = [-inf(size(feasibilityLimits)) feasibilityLimits];
varNames = horzcat(optRes.varNames,optRes.objNames);

customCorrPlot(hfig,front,designIdx,objIdx,designBounds,feasibilityLimits,varNames)