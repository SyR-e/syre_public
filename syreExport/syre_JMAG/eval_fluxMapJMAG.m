 %  Copyright 2024
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

function [NewDir] = eval_fluxMapJMAG(dataIn)

% calculates the flux map(id,iq) of an existing machine
% regular grid of (id,iq) combinations, d,q flux linkages versus id,iq
% example of inputs:
% CurrLoPP = 1, NumGrid = 10

% Uses matlabpool (parfor)

% Key INPUTs: CurrLoPP: current to be simulated
%             BrPP: remanence of all barriers magnets
%             NumOfRotPosPP: # simulated positions
%             AngularSpanPP: angular span of simulation
%             NumGrid: number of points in [0 Imax] for the single machine post-processing
%=========================================================================

if isfield(dataIn,'flagHWCMap')
    dataIntmp = dataIn;
end
pathname=dataIn.currentpathname;
filemot = strrep(dataIn.currentfilename,'.mat','.jmag');
load([dataIn.currentpathname dataIn.currentfilename]);

if isfield(dataIn,'flagHWCMap')
    dataIn = dataIntmp;
end

RatedCurrent = dataIn.RatedCurrent;
CurrLoPP = dataIn.CurrLoPP;
%SimulatedCurrent = dataIn.SimulatedCurrent;
SimulatedCurrent = RatedCurrent*CurrLoPP;
% GammaPP  = dataIn.GammaPP;
BrPP = dataIn.BrPP;
NumOfRotPosPP = dataIn.NumOfRotPosPP;
AngularSpanPP = dataIn.AngularSpanPP;
NumGrid = dataIn.NumGrid;
tempPP = dataIn.tempPP;

per.flag3phaseSet = dataIn.Active3PhaseSets;

per.EvalSpeed = dataIn.EvalSpeed;
per = rmfield(per,'custom_act');

MapQuadrants = dataIn.MapQuadrants;

clc;

if ~isfield(geo,'axisType')
    if strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype')
        geo.axisType = 'PM';
    else
        geo.axisType = 'SR';
    end
end

if ~strcmp(geo.axisType,dataIn.axisType)
    %geo.axisType = dataIn.axisType;
    if strcmp(dataIn.axisType,'PM')
        geo.th0 = geo.th0 - 90;
    else
        geo.th0 = geo.th0 + 90;
    end
end


eval_type = dataIn.EvalType;

per.overload        = CurrLoPP;
per.i0              = RatedCurrent;
per.BrPP            = BrPP;
per.nsim_singt      = NumOfRotPosPP;  % # simulated positions
per.delta_sim_singt = AngularSpanPP;  % angular span of simulation
per.tempPP          = tempPP;
%per.offset          = 0;

% iAmp = dataIn.SimulatedCurrent;

switch MapQuadrants
    case 1
        %if (strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype')) || strcmp(dataSet.axisType,'PM')
        if strcmp(dataIn.axisType,'PM')
            idvect = linspace(-SimulatedCurrent,0,NumGrid);
            iqvect = linspace(0,SimulatedCurrent,NumGrid);
        else
            idvect = linspace(0,SimulatedCurrent,NumGrid);
            iqvect = linspace(0,SimulatedCurrent,NumGrid);
        end
    case 2
        if strcmp(dataIn.axisType,'PM')
            idvect = linspace(-SimulatedCurrent,SimulatedCurrent,NumGrid+NumGrid-1);
            iqvect = linspace(0,SimulatedCurrent,NumGrid);
        else
            idvect = linspace(0,SimulatedCurrent,NumGrid);
            iqvect = linspace(-SimulatedCurrent,SimulatedCurrent,NumGrid+NumGrid-1);
        end
    case 4
        idvect = linspace(-SimulatedCurrent,SimulatedCurrent,NumGrid+NumGrid-1);
        iqvect = linspace(-SimulatedCurrent,SimulatedCurrent,NumGrid+NumGrid-1);
end

if isfield(dataIn,'flagHWCMap')
    CurrStep = dataIn.SimulatedCurrent/NumGrid;
    NumGridPM = ceil(dataIn.SimulatedCurrent_HWC/CurrStep);
    if strcmp(dataIn.axisType,'PM')
        idvect = linspace(-abs(dataIn.SimulatedCurrent_HWC),SimulatedCurrent,NumGrid+NumGridPM-1);
%         idvect = linspace(-SimulatedCurrent,dataIn.SimulatedCurrent_HWC,NumGrid+NumGridPM-1);
        iqvect = linspace(0,SimulatedCurrent,NumGrid);
    else
        idvect = linspace(0,SimulatedCurrent,NumGrid);
        iqvect = linspace(-SimulatedCurrent,abs(dataIn.SimulatedCurrent_HWC),NumGrid+NumGridPM-1);
    end
end


if strcmp(dataIn.TypeOfRotor,'IM')
    % remove the axis from the identification
    iStep = SimulatedCurrent/(NumGrid-1);
    idvect(idvect==0) = iStep/10;
    iqvect(iqvect==0) = iStep/10;
end


[Id,Iq] = meshgrid(idvect,iqvect);

I = Id + 1i * Iq;
iAmp = abs(I);
gamma = angle(I) * 180/pi;
% [i0,~]=calc_io(geo,per);

filemot = strrep(filemot,'.mat','.fem');

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
    [~,~,OUT{ii},~] = JMAGfitness([],geo,perTmp,mat,pathname,filemot);
end

Id = Id(:);
Iq = Iq(:);
Fd   = zeros(size(Id));
Fq   = zeros(size(Id));
T    = zeros(size(Id));
dT   = zeros(size(Id));
dTpp = zeros(size(Id));
% We   = zeros(size(Id));
% Wc   = zeros(size(Id));
SOL  = cell(size(Id));
if isfield(OUT{1},'Pfes_h')
    Pfes_h = zeros(size(Id));
    Pfes_c = zeros(size(Id));
    Pfer_h = zeros(size(Id));
    Pfer_c = zeros(size(Id));
    Ppm    = zeros(size(Id));
    velDim = OUT{1,1}.velDim;
end
if isfield(OUT{1},'IM')
    Ir      = zeros(size(Id));
    Fdr     = zeros(size(Id));
    Fqr     = zeros(size(Id));
    kr      = zeros(size(Id));
    Ibar    = cell(size(Id));
    Vbar    = cell(size(Id));
    Fbar    = cell(size(Id));
    FbarTot = cell(size(Id));
end

for ii=1:length(Id)
    Fd(ii)   = OUT{ii}.fd;
    Fq(ii)   = OUT{ii}.fq;
    T(ii)    = OUT{ii}.T;
    dT(ii)   = OUT{ii}.dT;
    dTpp(ii) = OUT{ii}.dTpp;
%     We(ii)   = OUT{ii}.We;
%     Wc(ii)   = OUT{ii}.Wc;
    SOL{ii}  = OUT{ii}.SOL;
    
    if isfield(OUT{ii},'Pfes_h')
        Pfes_h(ii) = OUT{ii}.Pfes_h;
        Pfes_c(ii) = OUT{ii}.Pfes_c;
        Pfer_h(ii) = OUT{ii}.Pfer_h;
        Pfer_c(ii) = OUT{ii}.Pfer_c;
        Ppm(ii)    = OUT{ii}.Ppm;
    end

    if isfield(OUT{ii},'IM')
        Ir(ii)      = OUT{ii}.IM.ir;
        Fdr(ii)     = OUT{ii}.IM.fdr;
        Fqr(ii)     = OUT{ii}.IM.fqr;
        kr(ii)      = OUT{ii}.IM.kr;
        Ibar{ii}    = OUT{ii}.IM.Ibar;
        Vbar{ii}    = OUT{ii}.IM.Vbar;
        Fbar{ii}    = OUT{ii}.IM.Fbar;
        FbarTot{ii} = OUT{ii}.IM.FbarTot;
    end
end

Id   = reshape(Id,[nR,nC]);
Iq   = reshape(Iq,[nR,nC]);
Fd   = reshape(Fd,[nR,nC]);
Fq   = reshape(Fq,[nR,nC]);
T    = reshape(T,[nR,nC]);
dT   = reshape(dT,[nR,nC]);
dTpp = reshape(dTpp,[nR,nC]);
% We   = reshape(We,[nR,nC]);
% Wc   = reshape(Wc,[nR,nC]);
SOL  = reshape(SOL,[nR,nC]);

if isfield(OUT{1},'Pfes_h')
    Pfes_h = reshape(Pfes_h,[nR,nC]);
    Pfes_c = reshape(Pfes_c,[nR,nC]);
    Pfer_h = reshape(Pfer_h,[nR,nC]);
    Pfer_c = reshape(Pfer_c,[nR,nC]);
    Ppm    = reshape(Ppm,[nR,nC]);
end

if isfield(OUT{1},'IM')
    Ir      = reshape(Ir,[nR,nC]);
    Fdr     = reshape(Fdr,[nR,nC]);
    Fqr     = reshape(Fqr,[nR,nC]);
    kr      = reshape(kr,[nR,nC]);
    Ibar    = reshape(Ibar,[nR,nC]);
    Vbar    = reshape(Vbar,[nR,nC]);
    Fbar    = reshape(Fbar,[nR,nC]);
    FbarTot = reshape(FbarTot,[nR,nC]);
end


F_map.Id   = Id;
F_map.Iq   = Iq;
F_map.Fd   = Fd;
F_map.Fq   = Fq;
F_map.T    = T;
F_map.dT   = dT;
F_map.dTpp = dTpp;
% F_map.We   = We;
% F_map.Wc   = Wc;

if exist('Pfes_h','var')
    F_map.Pfes_h = Pfes_h;
    F_map.Pfes_c = Pfes_c;
    F_map.Pfer_h = Pfer_h;
    F_map.Pfer_c = Pfer_c;
    F_map.Ppm    = Ppm;
    F_map.Pfe    = Pfes_h+Pfes_c+Pfer_h+Pfer_c;
    F_map.velDim = velDim;
end

if exist('Ir','var')
    F_map.IM.Ir  = Ir;
    F_map.IM.Fdr = Fdr;
    F_map.IM.Fqr = Fqr;
    F_map.IM.kr  = kr;

    F_map.bar.I    = Ibar;
    F_map.bar.V    = Vbar;
    F_map.bar.F    = Fbar;
    F_map.bar.Ftot = FbarTot;

    [F_map] = elab_Fmap_IM(F_map,geo,per);

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
NewDir = [pathname resFolder 'F_map_' Idstr 'x' Iqstr '_' int2str(per.tempPP) 'deg_' int2str(MapQuadrants) 'Q'];
if sum(per.flag3phaseSet)~=geo.win.n3phase
    NewDir = [NewDir '_' mat2str(per.flag3phaseSet)];
end

if strcmp(dataIn.EvalType,'singmIron')
    nStr = int2str(per.EvalSpeed);
    nStr = strrep(nStr,'.','rpm');
    if ~strcmpi(nStr,'rpm')
        nStr = [nStr 'rpm'];
    end
    NewDir = [NewDir '_' nStr '_ironLoss'];
end

NewDir = [NewDir '_JMAG'];

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
plot_singm(F_map,NewDir);

% add motor information to flux map files
dataSet.RatedCurrent     = RatedCurrent;
dataSet.CurrLoPP         = CurrLoPP;
dataSet.SimulatedCurrent = SimulatedCurrent;
dataSet.BrPP             = BrPP;
dataSet.NumOfRotPosPP    = NumOfRotPosPP;
dataSet.AngularSpanPP    = AngularSpanPP;
dataSet.NumGrid          = NumGrid;
dataSet.EvalSpeed        = per.EvalSpeed;
dataSet.axisType         = dataIn.axisType;

save([NewDir,'F_map','.mat'],'dataSet','geo','per','mat','-append');
save([NewDir,'fdfq_idiq_n256.mat'],'dataSet','geo','per','mat','-append'); 

if nargout()==0
    clear NewDir
end

