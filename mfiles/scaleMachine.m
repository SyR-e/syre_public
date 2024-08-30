% Copyright 2021
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

function [dataSet] = scaleMachine(dataSet,preValue)

ScaleFactor = dataSet.StatorOuterRadius/preValue;

% Main Data
dataSet.ShaftRadius     = dataSet.ShaftRadius*ScaleFactor;
dataSet.AirGapRadius    = dataSet.AirGapRadius*ScaleFactor;
dataSet.AirGapThickness = dataSet.AirGapThickness*ScaleFactor;

% Geometry
dataSet.ToothLength    = dataSet.ToothLength*ScaleFactor;
dataSet.ToothWidth     = dataSet.ToothWidth*ScaleFactor;
dataSet.ToothTangDepth = dataSet.ToothTangDepth*ScaleFactor;
dataSet.FilletCorner   = dataSet.FilletCorner*ScaleFactor;
dataSet.RadShiftInner  = dataSet.RadShiftInner*ScaleFactor;

% Options
dataSet.TanRibEdit      = dataSet.TanRibEdit*ScaleFactor;
dataSet.RotorFilletTan1 = dataSet.RotorFilletTan1*ScaleFactor;
dataSet.RotorFilletTan2 = dataSet.RotorFilletTan2*ScaleFactor;
dataSet.RadRibEdit      = dataSet.RadRibEdit*ScaleFactor;
dataSet.RotorFilletIn   = dataSet.RotorFilletIn*ScaleFactor;
dataSet.RotorFilletOut  = dataSet.RotorFilletOut*ScaleFactor;

dataSet.SleeveThickness = dataSet.SleeveThickness*ScaleFactor;

dataSet.Mesh      = dataSet.Mesh*ScaleFactor;
% dataSet.Mesh_MOOA = dataSet.Mesh_MOOA*ScaleFactor;
dataSet.MinMechTol = dataSet.MinMechTol*ScaleFactor ;
% dataSet.OverSpeed = dataSet.OverSpeed*ScaleFactor;
% dataSet.RadRibCheck = 1;
dataSet.OverSpeed = dataSet.OverSpeed/ScaleFactor;

% Materials
dataSet.PMdim = dataSet.PMdim*ScaleFactor;


end