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

function eval_fluxMapMN(dataIn,varargin)

%% singm mode simulation in MagNet (MN)
% regular grid of (id,iq) combinations, builds the magnetic model
% (d,q flux linkages over id,iq)

% Open matlabpool manually prior to execution



pathname=dataIn.currentpathname;
%  filemot = strrep(dataIn.currentfilename,'.mat','.fem');
filemot= dataIn.currentfilename;
load([dataIn.currentpathname dataIn.currentfilename]);

CurrLoPP = dataIn.CurrLoPP;
GammaPP  = dataIn.GammaPP;
BrPP = dataIn.BrPP;
NumOfRotPosPP = dataIn.NumOfRotPosPP;
AngularSpanPP = dataIn.AngularSpanPP;
NumGrid = dataIn.NumGrid;
per.EvalSpeed = dataIn.EvalSpeed;


clc;
syreRoot = fileparts(which('GUI_Syre.mlapp'));
current_path = syreRoot;

overload_temp =  CurrLoPP;   % current to be simulated
gamma_temp = GammaPP;        % current phase angle
Br = BrPP;                   % remanence of all barriers magnets

if (strcmp(dataIn.EvalType,'singt')||strcmp(dataIn.EvalType,'singtIron'))
    eval_type='singt';
elseif (strcmp(dataIn.EvalType,'singm')||strcmp(dataIn.EvalType,'singmIron'))
    eval_type='singm';
else
    error('Evaluation type not supported in MagNet!')
end

per.overload=CurrLoPP;
per.BrPP=BrPP;

per.nsim_singt = NumOfRotPosPP;       % # simulated positions
per.delta_sim_singt = AngularSpanPP;  % angular span of simulation

% Magnet coercivity (QUESTO NON SERVE IN MAGNET - CONTANO SOLO LA TEMPERATURA DEL MAGNETE E LA DEFINIZIONE DEL MATERIALE)
% Hc = 1/(4e-7*pi)*Br;    % PM coercivity
% mat.LayerMag.Hc = Hc;


% flux map over a grid of id,iq combinations
n_grid = NumGrid;     % number of points in [0 Imax] for the single machine post-processing
iAmp = overload_temp*calc_io(geo,per);
switch length(iAmp)
    case 1  % square domain
        if strcmp(geo.RotType,'SPM')
            idvect = linspace(-iAmp,0,n_grid);
            iqvect = linspace(0,iAmp,n_grid);
        else
            idvect = linspace(0,iAmp,n_grid);
            iqvect = linspace(0,iAmp,n_grid);
        end
    case 2  % rectangular domain
        if strcmp(geo.RotType,'SPM')
            idvect = linspace(-iAmp(1),0,n_grid);
            iqvect = linspace(0,iAmp(2),n_grid);
        else
            idvect = linspace(0,iAmp(1),n_grid);
            iqvect = linspace(0,iAmp(2),n_grid);
        end
    case 4 % CurrLoPP = [IdMin IdMax IqMin IqMax]
        idvect=linspace(iAmp(1),iAmp(2),n_grid);
        iqvect=linspace(iAmp(3),iAmp(4),n_grid);
end

[F_map,OUT] = eval_FdFq_tables_in_MN(geo,per,mat,dataIn,idvect,iqvect,eval_type,pathname, filemot);

% builds a new folder for each id, iq simulation
Idstr=num2str(max(abs(idvect)),3); Idstr = strrep(Idstr,'.','A');
Iqstr=num2str(max(abs(iqvect)),3); Iqstr = strrep(Iqstr,'.','A');

if isempty(strfind(Idstr, 'A'))
    Idstr = [Idstr 'A'];
end
if isempty(strfind(Iqstr, 'A'))
    Iqstr = [Iqstr 'A'];
end

resFolder = [filemot(1:end-4) '_results\FEA results\'];
if ~exist([pathname resFolder],'dir')
    mkdir([pathname resFolder]);
end

NewDir=[pathname resFolder filemot(1:end-4) '_F_map_' Idstr 'x' Iqstr '_MN'];
mkdir(NewDir);
NewDir=[NewDir '\'];
if isoctave()            %OCT
    file_name1= strcat(NewDir,'F_map','.mat');
    save('-v7', file_name1,'F_map','OUT');
    clear file_name1
else
    save([NewDir,'F_map','.mat'],'F_map','OUT');
end

% interp and then plots the magnetic curves
filename=filemot(1:end-4);
%         plot_singm;
plot_singm(F_map,NewDir,filename);
save([NewDir,'fdfq_idiq_n256.mat'],'dataSet','geo','per','mat','-append'); 


