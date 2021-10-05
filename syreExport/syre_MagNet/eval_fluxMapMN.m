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
filemot= dataIn.currentfilename;
load([dataIn.currentpathname dataIn.currentfilename]);

RatedCurrent = dataIn.RatedCurrent;
CurrLoPP = dataIn.CurrLoPP;
SimulatedCurrent = dataIn.SimulatedCurrent;
GammaPP  = dataIn.GammaPP;
BrPP = dataIn.BrPP;
NumOfRotPosPP = dataIn.NumOfRotPosPP;
AngularSpanPP = dataIn.AngularSpanPP;
NumGrid = dataIn.NumGrid;
MapQuadrants = dataIn.MapQuadrants;
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
% n_grid = NumGrid;     % number of points in [0 Imax] for the single machine post-processing
% iAmp = overload_temp*calc_io(geo,per);
% switch length(iAmp)
%     case 1  % square domain
%         if strcmp(geo.RotType,'SPM')
%             idvect = linspace(-iAmp,0,n_grid);
%             iqvect = linspace(0,iAmp,n_grid);
%         else
%             idvect = linspace(0,iAmp,n_grid);
%             iqvect = linspace(0,iAmp,n_grid);
%         end
%     case 2  % rectangular domain
%         if strcmp(geo.RotType,'SPM')
%             idvect = linspace(-iAmp(1),0,n_grid);
%             iqvect = linspace(0,iAmp(2),n_grid);
%         else
%             idvect = linspace(0,iAmp(1),n_grid);
%             iqvect = linspace(0,iAmp(2),n_grid);
%         end
%     case 4 % CurrLoPP = [IdMin IdMax IqMin IqMax]
%         idvect=linspace(iAmp(1),iAmp(2),n_grid);
%         iqvect=linspace(iAmp(3),iAmp(4),n_grid);
% end

switch MapQuadrants
    case 1
        if (strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype'))
            idvect = linspace(-SimulatedCurrent,0,NumGrid);
            iqvect = linspace(0,SimulatedCurrent,NumGrid);
        else
            idvect = linspace(0,SimulatedCurrent,NumGrid);
            iqvect = linspace(0,SimulatedCurrent,NumGrid);
        end
    case 2
        if (strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype'))
            idvect = linspace(-SimulatedCurrent,SimulatedCurrent,NumGrid+NumGrid-1);
            iqvect = linspace(0,SimulatedCurrent,NumGrid);
        else
            idvect = linspace(0,SimulatedCurrent,NumGrid);
            iqvect = linspace(-SimulatedCurrent,SimulatedCurrent,NumGrid+NumGrid-1);
        end
    case 4
        if (strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype'))
            idvect = linspace(-SimulatedCurrent,SimulatedCurrent,NumGrid+NumGrid-1);
            iqvect = linspace(-SimulatedCurrent,SimulatedCurrent,NumGrid+NumGrid-1);
        else
            idvect = linspace(-SimulatedCurrent,SimulatedCurrent,NumGrid+NumGrid-1);
            iqvect = linspace(-SimulatedCurrent,SimulatedCurrent,NumGrid+NumGrid-1);
        end
end

% ----------------------------------------------------------------
%% temp
[Id,Iq] = meshgrid(idvect,iqvect);

I = Id + 1i * Iq;
iAmp = abs(I);
gamma = angle(I) * 180/pi;

mat = mat;  %#ok parfor compatibility (do not comment)
geo = geo;  %#ok parfor compatibility (do not comment)

% create vector for parpool efficiency
[nR,nC] = size(iAmp);
iVect = iAmp(:);
gVect = gamma(:);

for ii=1:length(iVect)
    disp(['Evaluation of position I:',num2str(iVect(ii)),' gamma:',num2str(gVect(ii))]);
    perTmp = per;
    perTmp.gamma    = gVect(ii);
    perTmp.overload = iVect(ii)/per.i0;
    [~,~,OUT{ii},~] = MNfitness([],geo,perTmp,mat,eval_type,pathname,filemot);
end

Id = Id(:);
Iq = Iq(:);
Fd   = zeros(size(Id));
Fq   = zeros(size(Id));
T    = zeros(size(Id));
dT   = zeros(size(Id));
dTpp = zeros(size(Id));
if isfield(OUT{1},'Pfes_h')
    Pfes_h = zeros(size(Id));
    Pfes_c = zeros(size(Id));
    Pfer_h = zeros(size(Id));
    Pfer_c = zeros(size(Id));
    Ppm    = zeros(size(Id));
    velDim = OUT{1,1}.velDim;
end

for ii=1:length(Id)
    Fd(ii)   = OUT{ii}.fd;
    Fq(ii)   = OUT{ii}.fq;
    T(ii)    = OUT{ii}.T;
    dT(ii)   = OUT{ii}.dT;
    dTpp(ii) = OUT{ii}.dTpp;
    SOL{ii}  = OUT{ii}.SOL;
    
    if isfield(OUT{ii},'Pfes_h')
        Pfes_h(ii) = OUT{ii}.Pfes_h;
        Pfes_c(ii) = OUT{ii}.Pfes_c;
        Pfer_h(ii) = OUT{ii}.Pfer_h;
        Pfer_c(ii) = OUT{ii}.Pfer_c;
        Ppm(ii)    = OUT{ii}.Ppm;
    end
end

Id   = reshape(Id,[nR,nC]);
Iq   = reshape(Iq,[nR,nC]);
Fd   = reshape(Fd,[nR,nC]);
Fq   = reshape(Fq,[nR,nC]);
T    = reshape(T,[nR,nC]);
dT   = reshape(dT,[nR,nC]);
dTpp = reshape(dTpp,[nR,nC]);
SOL  = reshape(SOL,[nR,nC]);

if isfield(OUT{1},'Pfes_h')
    Pfes_h = reshape(Pfes_h,[nR,nC]);
    Pfes_c = reshape(Pfes_c,[nR,nC]);
    Pfer_h = reshape(Pfer_h,[nR,nC]);
    Pfer_c = reshape(Pfer_c,[nR,nC]);
    Ppm    = reshape(Ppm,[nR,nC]);
end

F_map.Id   = Id;
F_map.Iq   = Iq;
F_map.Fd   = Fd;
F_map.Fq   = Fq;
F_map.T    = T;
F_map.dT   = dT;
F_map.dTpp = dTpp;

if exist('Pfes_h')
    F_map.Pfes_h = Pfes_h;
    F_map.Pfes_c = Pfes_c;
    F_map.Pfer_h = Pfer_h;
    F_map.Pfer_c = Pfer_c;
    F_map.Ppm    = Ppm;
    F_map.Pfe    = Pfes_h+Pfes_c+Pfer_h+Pfer_c;
    F_map.velDim = velDim;
end

% results folder
Idstr=num2str(max(abs(idvect)),3); Idstr = strrep(Idstr,'.','A');
Iqstr=num2str(max(abs(iqvect)),3); Iqstr = strrep(Iqstr,'.','A');

if ~contains(Idstr, 'A')
    Idstr = [Idstr 'A'];
end
if ~contains(Iqstr, 'A')
    Iqstr = [Iqstr 'A'];
end

% NewDir=[pathname,[filemot(1:end-4) '_F_map_' Idstr 'x' Iqstr]];

resFolder = [filemot(1:end-4) '_results\FEA results\'];
if ~exist([pathname resFolder],'dir')
    mkdir([pathname resFolder]);
end

% NewDir = [pathname resFolder filemot(1:end-4) '_F_map_' Idstr 'x' Iqstr '_' int2str(per.tempPP) 'deg'];
NewDir = [pathname resFolder 'F_map_' Idstr 'x' Iqstr '_' int2str(per.tempPP) 'deg_' int2str(MapQuadrants) 'Q_MN'];
if exist('Pfes_h')
    nStr = int2str(per.EvalSpeed);
    nStr = strrep(nStr,'.','rpm');
    if ~strcmpi(nStr,'rpm')
        nStr = [nStr 'rpm'];
    end
    NewDir = [NewDir '_' nStr '_ironLoss'];
end
mkdir(NewDir);
NewDir=[NewDir '\'];
if isoctave()            %OCT
    file_name1= strcat(NewDir,'F_map','.mat');
    save('-v7', file_name1,'F_map','SOL');
    clear file_name1
else
    save([NewDir,'F_map','.mat'],'F_map','SOL','dataSet','geo','per','mat');
end

% interp and then plots the magnetic curves
[pathname,filename,ext] = fileparts(filemot);
plot_singm(F_map,NewDir,filename);

% add motor information to flux map files
dataSet.RatedCurrent     = RatedCurrent;
dataSet.CurrLoPP         = CurrLoPP;
dataSet.SimulatedCurrent = SimulatedCurrent;
dataSet.BrPP             = BrPP;
dataSet.NumOfRotPosPP    = NumOfRotPosPP;
dataSet.AngularSpanPP    = AngularSpanPP;
dataSet.NumGrid          = NumGrid;
dataSet.EvalSpeed        = per.EvalSpeed;

save([NewDir,'F_map','.mat'],'dataSet','geo','per','mat','-append');
save([NewDir,'fdfq_idiq_n256.mat'],'dataSet','geo','per','mat','-append'); 
