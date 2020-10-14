% Copyright 2019
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

function SimulateOCT(varargin)

% script for simulating an existing machine (OCT or Matlab without GUI)
% Uses matlabpool (parfor)

%   Key INPUTs: CurrLoPP: current to be simulated
%               GammaPP: current phase angle
%               tempPP: PM temperature --> PM remanence 
%               NumOfRotPosPP: # simulated positions
%               AngularSpanPP: angular span of simulation
%               NumGrid: number of points in [0 Imax] for the single machine post-processing
%               EvalType: type of simulation (see description)
%=========================================================================

% EvalType = singt (single point)
% simulates single or multiple (id,iq) conditions
% example inputs:
% single condition: CurrLoPP = 1, GammaPP = 45
% multiple points:  CurrLoPP = [1 1.5 2], gamma = [45 45 45]

% EvalType = singm 
% calculates the flux map(id,iq) of an existing machine
% regular grid of (id,iq) combinations, d,q flux linkages versus id,iq
% example of inputs:
% CurrLoPP = 1, NumGrid = 10

% EvalType = flxdn
% same as singt + calculates Bgap (B at the airgap)

% EvalType = force
% same as singt + calculates magnetic pull (rotor to stator force)

% EvalType = izero
% same as singt, with 10% zero sequence current superimposed to id, iq

% EvalType = idemag
% evaluaters the demagnetization current limit (steady state)

% EvalType = ichval
% evaluates the characteristic current iteratively (PM only)

% input dialog box, script run out of GUI
[filemot, pathname, fltidx]=uigetfile(' *.fem', 'Pick a motor');
load(strrep(filemot,'.fem','.mat'));

dataIn = dataSet;

if isoctave()
    temp = inputdlg({'current load [p.u.]';'gamma [deg]';'Br [T]'; ...
        'number of rotor positions';'angular span (elt. deg.)'; ...
        'points in [0 Imax]';'Evaluation type'},'INPUT',...
        1,{dataSet.CurrLoPP;dataSet.GammaPP;dataSet.tempPP;dataSet.NumOfRotPosPP;dataSet.AngularSpanPP;dataSet.NumGrid;dataSet.EvalType});
else
    temp = inputdlg({'current load [p.u.]';'gamma [deg]';'PM temperature [°C]'; ...
        'number of rotor positions';'angular span (elt. deg.)'; ...
        'points in [0 Imax]';'Evaluation type'},'INPUT', ...
        1,{num2str(dataSet.CurrLoPP);num2str(dataSet.GammaPP);num2str(dataSet.tempPP);num2str(dataSet.NumOfRotPosPP);num2str(dataSet.AngularSpanPP);num2str(dataSet.NumGrid);dataSet.EvalType});
end

dataIn.CurrLoPP = eval(cell2mat(temp(1)));         % current to be simulated (p.u.)
dataIn.GammaPP = eval(cell2mat(temp(2)));          % current phase angle
dataIn.tempPP = eval(cell2mat(temp(3)));           % PM temperature
dataIn.NumOfRotPosPP = eval(cell2mat(temp(4)));    % # simulated positions
dataIn.AngularSpanPP = eval(cell2mat(temp(5)));    % angular span of simulation
dataIn.NumGrid = eval(cell2mat(temp(6)));          % number of points in [0 Imax] for the single machine post-processing
dataIn.EvalType = (cell2mat(temp(7)));

% translate tempPP into BrPP (PM remanence)
mat=material_properties_layer(dataSet.FluxBarrierMaterial);
if isfield(mat,'temp')
    dataIn.BrPP=interp1(mat.temp.temp,mat.temp.Br,dataIn.tempPP);
else
    warning('This PM material does not have temperature data!!!')
    if isfield(mat,'LayerMag')
        dataIn.BrPP = mat.LayerMag.Br;
    else
        dataIn.BrPP = 0;
    end
end

clc;

dataIn.currentpathname = pathname;
dataIn.currentfilename = strrep(filemot,'.fem','.mat');
[dataIn.RatedCurrent,~] = calc_io(geo,per);
dataIn.SimulatedCurrent = dataIn.CurrLoPP*dataIn.RatedCurrent;

switch dataIn.EvalType
    case {'singt','flxdn','izero','force'}
        eval_operatingPoint(dataIn);
    case 'singm'
        eval_fluxMap(dataIn);
    case 'idemag'
        eval_idemag(dataIn);
    case 'ichval'
        eval_ich(dataIn);
end




