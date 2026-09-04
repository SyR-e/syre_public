% Copyright 2026
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

function [dataSet,model] = eval_vonMisesStress_COMSOL(dataSet,model)

%  ========================================================================
%  Inizializzazione
%  ========================================================================

import com.comsol.model.*
import com.comsol.model.util.*

disp('COMSOL structural analysis...')

info = dataSet.infoComsol;
selDom = info.tags.selection.domain;
selBnd = info.tags.selection.boundary;

compTag = info.tags.component;
geomTag = info.tags.geometry;
meshTag = info.tags.mesh;

comp = model.component(compTag);

rpm = dataSet.EvalSpeed;                           % Velocita meccanica [rpm]
meshK = dataSet.Mesh;                               % Livello mesh COMSOL [1-9]   ------------------------------- LIMITARE min(max())
nDig = 16;                                         % Cifre conversione COMSOL    ------------------------------- Prendere info salvata da infoComsol
rpmC = [num2str(rpm,nDig) '[rpm]'];                % Velocita per COMSOL

nMaxStress = 1;                                    % Numero massimi nodali mostrati
defScale = 50;                                     % Scala grafica deformazione
stressDec = 0;                                     % Decimali etichette stress
nMaxStress = max(1,round(nMaxStress));

structDomTag  = 'domStructR';                      % Domini strutturali rotore
structSideTag = 'bndStructSideR';                  % Lati strutturali radiali
physTag       = 'solid';                           % Solid Mechanics
symTag        = 'sym1';                            % Simmetria
rotTag        = 'rotf1';                           % Rotating Frame
stdTag        = 'stdStruct';                       % Studio strutturale
statTag       = 'statStruct';                      % Step stazionario
dsetTag       = 'dsetStruct';                      % Dataset soluzione
pgTag         = 'pgStress';                        % Plot von Mises
surfTag       = 'surfStress';                      % Superficie von Mises


%  ========================================================================
%  Selezioni strutturali del rotore
%  ========================================================================

structIn = {selDom.feRot};

if isfield(selDom,'cuRot')
    structIn{end+1} = selDom.cuRot;
end

if strcmp(info.machineType,'EESM') && isfield(selDom,'fillRot')
    structIn{end+1} = selDom.fillRot;
end

comp.selection.create(structDomTag,'Union');
comp.selection(structDomTag).label('Rotor structural domains');
comp.selection(structDomTag).set('entitydim',2);
comp.selection(structDomTag).set('input',structIn);

dStruct = double(comp.selection(structDomTag).entities());
dStruct = unique(dStruct(:)');

bStruct = double(mphgetadj(model,geomTag,'boundary','domain',dStruct));
bSideR = double(comp.selection(selBnd.sideRot).entities());
bStruct = unique(bStruct(:)');
bSideR = unique(bSideR(:)');
bSym = intersect(bStruct,bSideR);

comp.selection.create(structSideTag,'Explicit').geom(1);
comp.selection(structSideTag).label('Rotor structural sector sides');
comp.selection(structSideTag).set(bSym);

model.nodeGroup(info.tags.group.domain).add('selection',structDomTag);
model.nodeGroup(info.tags.group.boundary).add('selection',structSideTag);


%  ========================================================================
%  Fisica meccanica
%  ========================================================================

comp.physics.create(physTag,'SolidMechanics',geomTag);
comp.physics(physTag).label('Rotor structural mechanics');
comp.physics(physTag).selection().named(structDomTag);

comp.physics(physTag).create(symTag,'SymmetrySolid',1);
comp.physics(physTag).feature(symTag).label('Sector symmetry');
comp.physics(physTag).feature(symTag).selection().named(structSideTag);

comp.physics(physTag).create(rotTag,'RotatingFrame',2);
comp.physics(physTag).feature(rotTag).label('Centrifugal load');
comp.physics(physTag).feature(rotTag).selection().named(structDomTag);
comp.physics(physTag).feature(rotTag).set('RotationalFrequency','RevolutionsPerTime');
comp.physics(physTag).feature(rotTag).set('rpt',rpmC);


%  ========================================================================
%  Mesh e studio
%  ========================================================================

comp.mesh(meshTag).automatic(true);
comp.mesh(meshTag).autoMeshSize(meshK);
comp.mesh(meshTag).run;

model.study.create(stdTag);
model.study(stdTag).label(['Structural analysis - ' num2str(rpm) ' rpm']);
model.study(stdTag).create(statTag,'Stationary');
model.study(stdTag).setGenPlots(false);
model.study(stdTag).run;

solTags = model.study(stdTag).getSolverSequences('SolverSequence');
if iscell(solTags)
    solTag = char(solTags{end});
else
    solTag = char(solTags(end));
end

model.result.dataset.create(dsetTag,'Solution');
model.result.dataset(dsetTag).label('Structural solution');
model.result.dataset(dsetTag).set('solution',solTag);
model.result.dataset(dsetTag).set('comp',compTag);
model.result.dataset(dsetTag).set('frametype','material');
model.result.dataset(dsetTag).selection().named(structDomTag);


%  ========================================================================
%  Individuazione dei punti di massimo stress
%  ========================================================================

pd = mpheval(model,{'solid.misesGp','solid.disp','u','v','X','Y'}, ...
    'dataset',dsetTag, ...
    'selection',structDomTag, ...
    'edim','domain', ...
    'pattern','lagrange', ...
    'refine',1, ...
    'smooth','internal', ...
    'complexout','off', ...
    'unit',{'MPa','mm','mm','mm','mm','mm'});

sN  = real(pd.d1(:));                             % Von Mises [MPa]
uN  = real(pd.d2(:));                             % Spostamento totale [mm]
uxN = real(pd.d3(:));                             % Spostamento x [mm]
uyN = real(pd.d4(:));                             % Spostamento y [mm]
xN  = real(pd.d5(:));                             % Coordinata materiale x [mm]
yN  = real(pd.d6(:));                             % Coordinata materiale y [mm]

valid = isfinite(sN) & isfinite(uN) & isfinite(uxN) & ...
        isfinite(uyN) & isfinite(xN) & isfinite(yN);

sN  = sN(valid);
uN  = uN(valid);
uxN = uxN(valid);
uyN = uyN(valid);
xN  = xN(valid);
yN  = yN(valid);

if isempty(sN)
    error('eval_vonMisesStress_COMSOL:EmptyStressData', ...
        'Nessun valore di stress valido trovato nei domini strutturali.');
end

[~,ordStress] = sort(sN,'descend');

nStressMax = min(nMaxStress,numel(ordStress));
idxStressMax = zeros(nStressMax,1);
tolNode = max(1e-9,1e-8*info.geometry.rotorOuterRadius);

nFound = 0;

for ii = 1:numel(ordStress)

    idx = ordStress(ii);

    if nFound == 0
        isNewNode = true;
    else
        dist = hypot(xN(idx)-xN(idxStressMax(1:nFound)), ...
                     yN(idx)-yN(idxStressMax(1:nFound)));

        isNewNode = all(dist > tolNode);
    end

    if isNewNode
        nFound = nFound+1;
        idxStressMax(nFound) = idx;
    end

    if nFound == nStressMax
        break
    end

end

idxStressMax = idxStressMax(1:nFound);
nStressMax = nFound;

stressMax = sN(idxStressMax);
xStressMax = xN(idxStressMax);
yStressMax = yN(idxStressMax);
uxStressMax = uxN(idxStressMax);
uyStressMax = uyN(idxStressMax);

[maxDisp,iMaxDisp] = max(uN);

maxStress = stressMax(1);
maxStressPos = [xStressMax(1) yStressMax(1)];
maxDispPos = [xN(iMaxDisp) yN(iMaxDisp)];


%  ========================================================================
%  Plot dei risultati
%  ========================================================================

model.result.create(pgTag,'PlotGroup2D');
model.result(pgTag).label(['Von Mises stress - deformation x' num2str(defScale)]);
model.result(pgTag).set('data',dsetTag);
model.result(pgTag).set('frametype','material');

model.result(pgTag).create(surfTag,'Surface');
model.result(pgTag).feature(surfTag).label('Von Mises stress');
model.result(pgTag).feature(surfTag).set('expr',{'solid.misesGp'});
model.result(pgTag).feature(surfTag).set('unit','MPa');
model.result(pgTag).feature(surfTag).set('coloring','colortable');
model.result(pgTag).feature(surfTag).set('colortable','Prism');
model.result(pgTag).feature(surfTag).set('colorlegend','on');
model.result(pgTag).feature(surfTag).set('colorscalemode','linear');
model.result(pgTag).feature(surfTag).set('resolution','custom');
model.result(pgTag).feature(surfTag).set('refine',2);
model.result(pgTag).feature(surfTag).set('smooth','internal');

model.result(pgTag).feature(surfTag).create('sel1','Selection');
model.result(pgTag).feature(surfTag).feature('sel1').selection().named(structDomTag);

model.result(pgTag).feature(surfTag).create('def1','Deform');
model.result(pgTag).feature(surfTag).feature('def1').set('expr',{'u','v'});
model.result(pgTag).feature(surfTag).feature('def1').set('unit',{'mm','mm'});
model.result(pgTag).feature(surfTag).feature('def1').set('scaleactive','on');
model.result(pgTag).feature(surfTag).feature('def1').set('scale',defScale);

maxStressTags = cell(1,nStressMax);

for ii = 1:nStressMax

    maxStressTags{ii} = ['maxStress' num2str(ii)];

    xPlot = xStressMax(ii) + defScale*uxStressMax(ii);
    yPlot = yStressMax(ii) + defScale*uyStressMax(ii);

    stressText = sprintf('%.*f MPa',stressDec,stressMax(ii));

    model.result(pgTag).create(maxStressTags{ii},'AnnotationData');
    model.result(pgTag).feature(maxStressTags{ii}).label(['Maximum stress ' num2str(ii)]);
    model.result(pgTag).feature(maxStressTags{ii}).set('pos',[xPlot yPlot]);
    model.result(pgTag).feature(maxStressTags{ii}).set('text',stressText);
    model.result(pgTag).feature(maxStressTags{ii}).set('color','black');
    model.result(pgTag).feature(maxStressTags{ii}).set('backgroundcolor','none');
    model.result(pgTag).feature(maxStressTags{ii}).set('showpoint',true);
    model.result(pgTag).feature(maxStressTags{ii}).set('pointradius',5);
    model.result(pgTag).feature(maxStressTags{ii}).set('showframe',false);
    model.result(pgTag).feature(maxStressTags{ii}).set('titletype','none');

end

model.result(pgTag).run;

fig = figure;
ax = axes('Parent',fig);

mphplot(model,pgTag,'parent',ax,'rangenum',1,'run','off');

% Limiti della scala di stress
% stressLim = [min(sN) max(sN)];

% if stressLim(1) == stressLim(2)
%     stressLim = stressLim + [-1 1]*max(abs(stressLim(1))*1e-6,eps);
% end

% Colormap COMSOL
% colormap(ax,mphcolortable(model,'Prism'));
% clim(ax,stressLim);

% Colorbar MATLAB
cb = colorbar(ax,'eastoutside');
% cb.Label.String = 'von Mises stress (MPa)';
% cb.Label.Interpreter = 'none';
cb.Label.String = '(MPa)';

% ax.Position = [0.08 0.10 0.76 0.84];


%  ========================================================================
%  Aggiornamento infoComsol e salvataggio
%  ========================================================================

dataSet.infoComsol.tags.selection.domain.structRot = structDomTag;
dataSet.infoComsol.tags.selection.boundary.structSideRot = structSideTag;
dataSet.infoComsol.tags.physics.structural = physTag;
dataSet.infoComsol.tags.study.structural = stdTag;
dataSet.infoComsol.tags.solution.structural = solTag;
dataSet.infoComsol.tags.dataset.structural = dsetTag;
dataSet.infoComsol.tags.result.vonMises = pgTag;
dataSet.infoComsol.tags.result.maxStressMarkers = maxStressTags;

dataSet.infoComsol.structural.EvalSpeed = rpm;
dataSet.infoComsol.structural.meshK = meshK;
dataSet.infoComsol.structural.deformationScale = defScale;

dataSet.infoComsol.structural.nMaxStressRequested = nMaxStress;
dataSet.infoComsol.structural.nMaxStressFound = nStressMax;
dataSet.infoComsol.structural.maxStressMPa = stressMax(:)';
dataSet.infoComsol.structural.maxStressPositionMM = [xStressMax yStressMax];
dataSet.infoComsol.structural.maxStressDisplacementMM = [uxStressMax uyStressMax];
dataSet.infoComsol.structural.maxStressPlotPositionMM = [xStressMax+defScale*uxStressMax, yStressMax+defScale*uyStressMax];

dataSet.infoComsol.structural.maxElementVonMisesMPa = maxStress;
dataSet.infoComsol.structural.maxElementPositionMM = maxStressPos;
dataSet.infoComsol.structural.maxDisplacementMM = maxDisp;
dataSet.infoComsol.structural.maxDisplacementPositionMM = maxDispPos;

[saveDir,fileName,~] = fileparts(info.files.model);
speedName = strrep(num2str(rpm,'%g'),'.','p');
structName = [fileName '_structural_' speedName 'rpm'];
structFile = fullfile(saveDir,[structName '.mph']);

model.modelPath(saveDir);
model.label([structName '.mph']);

dataSet.infoComsol.files.structural = structFile;

mphsave(model,structFile);

disp('Done!')

for ii = 1:nStressMax

    fprintf('\nMaximum stress point %d: %.1f MPa at [%.1f, %.1f] mm\n', ii,stressMax(ii),xStressMax(ii),yStressMax(ii));

end

fprintf('Maximum displacement:   %.4f mm at [%.1f, %.1f] mm\n\n', maxDisp,maxDispPos(1),maxDispPos(2));

fprintf('Structural model saved: %s\n',structFile);


end
