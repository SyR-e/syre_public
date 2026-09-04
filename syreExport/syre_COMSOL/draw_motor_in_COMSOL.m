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

function [model,geo] = draw_motor_in_COMSOL(geo,mat,pathIn,nameIn)

%  ========================================================================
%  Inizializzazione del modello
%  ========================================================================

import com.comsol.model.*
import com.comsol.model.util.*

[~,fileName,~] = fileparts(nameIn);

comsolDir = fullfile(pathIn,[fileName '_Comsol\']);

if exist(comsolDir,'dir') ~= 7
    mkdir(comsolDir);
end

modelFile     = fullfile(comsolDir,[fileName '.mph']);
rotorDxfFile  = fullfile(comsolDir,[fileName '_rot.dxf' ]);
statorDxfFile = fullfile(comsolDir,[fileName '_stat.dxf']);

% Tag principali del modello COMSOL
modelTag     = 'MotorModel';
% modelTag     = ['MotorModel' num2str(geo.Acoilf)];
componentTag = 'comp1';
geometryTag  = 'geom1';
meshTag      = 'mesh1';

% Tag delle feature geometriche
rotorImportTag  = 'impRot';
splitRotorTag   = 'splRot';
if strcmp(geo.RotType,'EESM')
    unionCoilFillTag = 'uniCuFillR';
end
statorImportTag = 'impStat';



%  ========================================================================
%  Creazione del modello Comsol
%  ========================================================================

model = ModelUtil.create(modelTag);

model.modelPath(pathIn);
model.label([fileName '.mph']);

model.component.create(componentTag,true);

model.component(componentTag).geom.create(geometryTag,2);
model.component(componentTag).geom(geometryTag).lengthUnit('mm');

% Il nodo mesh viene creato ma non viene costruito.
% Le funzioni fisiche successive definiranno la mesh appropriata.
model.component(componentTag).mesh.create(meshTag);


%  ========================================================================
%  Costruzione delle due metà di traferro ed esportazione .dxf
%  ========================================================================

sectorAngle = pi/geo.p*geo.ps;

r_ro = geo.r;
r_ag = geo.r + geo.g/2;
r_si = geo.r + geo.g;

materialCodes;

idxMatAirRot = max(geo.rotor(:,end)) + 1;

geo.rotor_n = [r_ro, 0, r_ag, 0, NaN, NaN, 0, codMatAirRot, idxMatAirRot;
               0,    0, r_ag, 0, r_ag*cos(sectorAngle), r_ag*sin(sectorAngle), 1, codMatAirRot, idxMatAirRot;
               r_ro*cos(sectorAngle), r_ro*sin(sectorAngle), r_ag*cos(sectorAngle), r_ag*sin(sectorAngle), NaN, NaN, 0, codMatAirRot, idxMatAirRot];


idxMatAirSta = max(geo.stator(:,end)) + 1;

geo.stator_n = [r_ag, 0, r_si, 0, NaN, NaN, 0, codMatAirSta, idxMatAirSta;
                0, 0, r_ag, 0, r_ag*cos(sectorAngle), r_ag*sin(sectorAngle), 1, codMatAirSta, idxMatAirSta;
                r_si*cos(sectorAngle), r_si*sin(sectorAngle), r_ag*cos(sectorAngle), r_ag*sin(sectorAngle), NaN, NaN, 0, codMatAirSta, idxMatAirSta];


% Geometrie complete da esportare
comsolRotor =  [geo.rotor;
                geo.rotor_n];

comsolStator = [geo.stator;
                geo.stator_n];


syreToDxf(NaN, comsolRotor,  comsolDir, [fileName '_rot.dxf' ]);
syreToDxf(comsolStator, NaN, comsolDir, [fileName '_stat.dxf']);


%  ========================================================================
%  Importazione .dxf su modello Comsol
%  ========================================================================

geom = model.component(componentTag).geom(geometryTag);
geom.create(rotorImportTag,'Import');
geom.feature(rotorImportTag).label('Rotor geometry');
geom.feature(rotorImportTag).set('type', 'dxf');
geom.feature(rotorImportTag).set('filename', rotorDxfFile);
geom.feature(rotorImportTag).set('selresult','on');
geom.feature(rotorImportTag).set('selindividual','on');
geom.feature(rotorImportTag).importData();

geom.create(splitRotorTag,'Split');
geom.feature(splitRotorTag).selection('input').set({rotorImportTag});
geom.run(splitRotorTag);

if strcmp(geo.RotType,'EESM')

    matRot = geo.BLKLABELS.rotore.xy(:,1:3);

    pCuR = matRot(matRot(:,3)==codMatCuRot,1:2);
    pFillR = matRot(matRot(:,3)==codMatAirRot,1:2);

    pCuR = unique(pCuR(all(isfinite(pCuR),2),:),'rows');
    pFillR = unique(pFillR(all(isfinite(pFillR),2),:),'rows');

    rObjSel = max(1e-7,1e-6*geo.r);
    objSel = cell(1,size(pCuR,1)+size(pFillR,1));
    kk = 0;

    for ii = 1:size(pCuR,1)

        kk = kk+1;
        tag = ['selCuObj' num2str(ii)];

        geom.create(tag,'DiskSelection');
        geom.feature(tag).set('entitydim',-1);
        geom.feature(tag).set('posx',pCuR(ii,1));
        geom.feature(tag).set('posy',pCuR(ii,2));
        geom.feature(tag).set('r',rObjSel);
        geom.feature(tag).set('condition','intersects');
        geom.feature(tag).set('selshow','off');

        objSel{kk} = tag;

    end

    for ii = 1:size(pFillR,1)

        kk = kk+1;
        tag = ['selFillObj' num2str(ii)];

        geom.create(tag,'DiskSelection');
        geom.feature(tag).set('entitydim',-1);
        geom.feature(tag).set('posx',pFillR(ii,1));
        geom.feature(tag).set('posy',pFillR(ii,2));
        geom.feature(tag).set('r',rObjSel);
        geom.feature(tag).set('condition','intersects');
        geom.feature(tag).set('selshow','off');

        objSel{kk} = tag;

    end

    geom.create('selCuFillObj','UnionSelection');
    geom.feature('selCuFillObj').set('entitydim',-1);
    geom.feature('selCuFillObj').set('input',objSel);
    geom.feature('selCuFillObj').set('selshow','off');

    geom.create(unionCoilFillTag,'Union');
    geom.feature(unionCoilFillTag).label('Rotor winding and fill');
    geom.feature(unionCoilFillTag).selection('input').named('selCuFillObj');
    geom.feature(unionCoilFillTag).set('intbnd','on');
    geom.run(unionCoilFillTag);

end

geom.create(statorImportTag,'Import');
geom.feature(statorImportTag).label('Stator geometry');
geom.feature(statorImportTag).set('type', 'dxf');
geom.feature(statorImportTag).set('filename', statorDxfFile);
geom.feature(statorImportTag).importData();

geom.feature('fin').set('action','assembly');
geom.feature('fin').set('pairtype','contact');
geom.run('fin');


%  ========================================================================
%  MEMORIZZAZIONE DELLE INFORMAZIONI COMSOL IN geo
%  ========================================================================

geo.infoComsol = struct();
geo.infoComsol.machineType = geo.RotType;
geo.infoComsol.modelName = fileName;

geo.infoComsol.tags = struct();
geo.infoComsol.tags.model = modelTag;
geo.infoComsol.tags.component = componentTag;
geo.infoComsol.tags.geometry = geometryTag;
geo.infoComsol.tags.mesh = meshTag;


% Contenitori che verranno riempiti nei blocchi successivi
geo.infoComsol.tags.selection = struct();
geo.infoComsol.tags.selection.domain = struct();
geo.infoComsol.tags.selection.boundary = struct();
geo.infoComsol.tags.material = struct();


% ------------------------------------------------------------------------
% File e cartelle
% ------------------------------------------------------------------------

geo.infoComsol.files = struct();
geo.infoComsol.files.inputDirectory = pathIn;
geo.infoComsol.files.comsolDirectory = comsolDir;
geo.infoComsol.files.model = modelFile;
geo.infoComsol.files.rotorDxf = rotorDxfFile;
geo.infoComsol.files.statorDxf = statorDxfFile;


% ------------------------------------------------------------------------
% Informazioni geometriche
% ------------------------------------------------------------------------

geo.infoComsol.geometry = struct();
geo.infoComsol.geometry.dimension = 2; % Modello 2D (to  be verified)
geo.infoComsol.geometry.lengthUnit = 'mm';
geo.infoComsol.geometry.sectorAngleDeg = sectorAngle*180/pi;
geo.infoComsol.geometry.rotorOuterRadius = r_ro;
geo.infoComsol.geometry.airgapMiddleRadius = r_ag;
geo.infoComsol.geometry.statorInnerRadius = r_si;

geo.infoComsol.geometry.numberOfDomains = double(geom.getNDomains());
geo.infoComsol.geometry.numberOfBoundaries = double(geom.getNBoundaries());


%  ========================================================================
%  Selezioni dei domini
%  ========================================================================

comp = model.component(componentTag);

matR = geo.BLKLABELS.rotore.xy(:,1:3);
matS = geo.BLKLABELS.statore.xy(:,1:3);

matDom = [matR; matS];

% matDom = matDom(all(isfinite(matDom(:,1:2)),2),:); In caso di punti non validi

comp.selection.create('findDom','Disk');
findDom = comp.selection('findDom');
findDom.set('entitydim',2);

dFillR = [];
dFillS = [];
dCuS = [];
dFeS = [];
dFeR = [];
dPM = [];
dShaft = [];
dCuR = [];
dSleeve = [];

for ii = 1:size(matDom,1)

    findDom.set('posx',matDom(ii,1));
    findDom.set('posy',matDom(ii,2));

    id = double(findDom.entities());
    id = id(:)';

    switch matDom(ii,3)

        case codMatAirRot
            dFillR = [dFillR id];

        case codMatAirSta
            dFillS = [dFillS id];

        case codMatCu
            dCuS = [dCuS id];

        case codMatFeSta
            dFeS = [dFeS id];

        case codMatFeRot
            dFeR = [dFeR id];

        case codMatBar
            dPM = [dPM id];

        case codMatShaft
            dShaft = [dShaft id];

        case codMatCuRot
            dCuR = [dCuR id];

        case codMatSleeve
            dSleeve = [dSleeve id];

    end

end

thAg = sectorAngle/2;
rAgR = (r_ro+r_ag)/2;
rAgS = (r_ag+r_si)/2;

findDom.set('posx',rAgR*cos(thAg));
findDom.set('posy',rAgR*sin(thAg));
dAgR = double(findDom.entities());
dAgR = dAgR(:)';

findDom.set('posx',rAgS*cos(thAg));
findDom.set('posy',rAgS*sin(thAg));
dAgS = double(findDom.entities());
dAgS = dAgS(:)';

dFillR = unique(dFillR);
dFillS = unique(dFillS);
dCuS = unique(dCuS);
dFeS = unique(dFeS);
dFeR = unique(dFeR);
dPM = unique(dPM);
dShaft = unique(dShaft);
dCuR = unique(dCuR);
dSleeve = unique(dSleeve);
dAgR = unique(dAgR);
dAgS = unique(dAgS);

dAg = unique([dAgR dAgS]);
dRot = unique([dFillR dAgR dFeR dPM dShaft dCuR dSleeve]);
dSta = unique([dFillS dAgS dCuS dFeS]);

selDom = {
    'airGap',   'domAg',     'Air gap',             dAg;
    'airGapRot','domAgR',    'Rotor air gap',       dAgR;
    'airGapSta','domAgS',    'Stator air gap',      dAgS;

    'rotor',     'domRot',     'Rotor domains',            dRot;
    'feRot',     'domFeR',     'Rotor iron',               dFeR;
    'cuRot',     'domCuR',     'Rotor conductors',         dCuR;
    'PMRot',     'domPMR',     'Rotor permanent magnets',  dPM;
    'fillRot',   'domFillR',   'Rotor fill',               dFillR;
    'sleeveRot', 'domSleeveR', 'Rotor sleeve',             dSleeve;
    'shaft',     'domShaft',   'Shaft',                    dShaft;
    
    'stator',   'domSta',    'Stator domains',      dSta;
    'feSta',    'domFeS',    'Stator iron',         dFeS;
    'cuSta',    'domCuS',    'Stator conductors',   dCuS;
    'fillSta',   'domFillS', 'Stator fill',         dFillS;
    };

for ii = 1:size(selDom,1)
    if ~isempty(selDom{ii,4})

        fld = selDom{ii,1};
        tag = selDom{ii,2};
        comp.selection.create(tag,'Explicit').geom(2);
        comp.selection(tag).label(selDom{ii,3});
        comp.selection(tag).set(selDom{ii,4});

        geo.infoComsol.tags.selection.domain.(fld) = tag;
    end

end

comp.selection.remove('findDom');


%  ========================================================================
%  Creazione e assegnazione dei materiali
%  ========================================================================

sel = geo.infoComsol.tags.selection.domain;

% ------------------------------------------------------------------------
% Definizione proprietà materiali
% ------------------------------------------------------------------------

nDig     = 16;                                      % Cifre conversione COMSOL
dimDom   = 2;                                       % Dimensione domini
matClass = 'Common';                                % Classe materiale COMSOL
matNS    = 'nonSolid';                              % Tipo materiale non solido
famAir   = 'air';                                   % Famiglia COMSOL aria
famCu    = 'copper';                                % Famiglia COMSOL rame
famFill  = 'plastic';                               % Famiglia COMSOL filler

pgDef = 'def';                                      % Gruppo proprietà base
pgEnu = 'Enu';                                      % Gruppo proprietà elastiche
pgBH  = 'BHCurve';                                  % Gruppo curva BH
pgRes = 'linzRes';                                  % Gruppo resistività lineare
pgPM  = 'RemanentFluxDensity';                      % Gruppo magneti permanenti
   
typeEnu  = 'Enu';                                   % Tipo gruppo elastico
descEnu  = 'Young''s_modulus_and_Poisson''s_ratio'; % Descrizione gruppo elastico
typeBH   = 'B-H Curve';                             % Tipo gruppo BH
typeInt  = 'Interpolation';                         % Tipo interpolazione BH
typePM   = 'Remanent flux density';                 % Tipo gruppo PM
typeRes  = 'Linearized resistivity';                % Tipo gruppo resistivita
selUnion = 'Union';                                 % Tipo selection aggregata

funBH  = 'BH';                                      % Nome funzione BH
bhExt  = 'linear';                                  % Estrapolazione BH
bhUnit = 'T';                                       % Unità B
hUnit  = 'A/m';                                     % Unità H
bhB    = [funBH '(normHin)'];                       % Relazione B(H)
bhH    = [funBH '_inv(normBin)'];                   % Relazione H(B)
bhW    = [funBH '_prim(normHin)'];                  % Energia magnetica
inTemp = 'temperature';                             % Input temperatura
inH    = 'magneticfield';                           % Input campo magnetico
inB    = 'magneticfluxdensity';                     % Input induzione magnetica

uRho   = '[kg/m^3]';                                % u.d.m. densità
uSig   = '[S/m]';                                   % u.d.m. conducibilità
uE     = '[Pa]';                                    % u.d.m. modulo elastico
uCp    = '[J/(kg*K)]';                              % u.d.m. calore specifico
uK     = '[W/(m*K)]';                               % u.d.m. conducibilità termica
uRhoE  = '[ohm*m]';                                 % u.d.m. resistività
uAlpha = '[1/K]';                                   % u.d.m. coefficiente termico
uTemp  = '[K]';                                     % u.d.m. temperatura
uBr    = '[T]';                                     % u.d.m. induzione residua


% Air gap
tagAir  = 'matAir';                                 % Tag COMSOL
nameAir = [char(mat.SlotAir.MatName) ' - air gap']; % Nome materiale
rhoAir  = mat.SlotAir.kgm3;                         % Densità [kg/m^3]
sigAir  = mat.SlotAir.sigma;                        % Conducibilità [S/m]
muAir   = mat.LayerAir.mu;                          % Permeabilità relativa [-]
epsAir  = 1;                                        % Permittività [-], esterno
rhoAirC = [num2str(rhoAir,nDig) uRho];              % numero + u.d.m. in stringa per Comsol
sigAirC = [num2str(sigAir,nDig) uSig];
muAirC  = num2str(muAir,nDig);
epsAirC = num2str(epsAir,nDig);

% ROTORE ------------------------------------------------------------------
% Ferro di rotore
tagFeR  = 'matFeR';                                 % Tag COMSOL
nameFeR = [char(mat.Rotor.MatName) ' - rotor'];     % Nome materiale
rhoFeR  = mat.Rotor.kgm3;                           % Densità [kg/m^3]
EFeR    = mat.Rotor.E*1e9;                          % Modulo di Young [Pa]  MM -> 185e9
nuFeR   = 0.28;                                     % Poisson [-], esterno
sigFeR  = 0;                                        % Conducibilità [S/m], laminato
epsFeR  = 1;                                        % Permittività [-], esterno
rhoFeRC = [num2str(rhoFeR,nDig) uRho];              % numero + u.d.m. in stringa per Comsol
EFeRC   = [num2str(EFeR,nDig) uE];
nuFeRC  = num2str(nuFeR,nDig);
sigFeRC = [num2str(sigFeR,nDig) uSig];
epsFeRC = num2str(epsFeR,nDig);
bhFeR = num2cell([mat.Rotor.BH(:,2) mat.Rotor.BH(:,1)]);
for ii = 1:numel(bhFeR)
    bhFeR{ii} = num2str(bhFeR{ii},nDig);
end
%bhFeR   = vertcat(bhFeR{:});

% Conduttori di rotore
tagCuR    = 'matCuR';                               % Tag COMSOL
nameCuR   = [char(mat.BarCond.MatName) ' - rotor']; % Nome materiale
rhoCuR    = mat.BarCond.kgm3*geo.win.kcuf;          % Densità equivalente [kg/m^3]  MM -> 5480
sigCuR    = mat.BarCond.sigma;                      % Conducibilità [S/m]
rho0CuR   = 1/sigCuR;                               % Resistività a Tref [ohm*m]
alphaCuR  = mat.BarCond.alpha;                      % Coefficiente termico [1/K]
TrefCuR   = 293.15;                                 % Temperatura riferimento [K], esterno
ECuR      = mat.BarCond.E*10^6;                     % Modulo equivalente [Pa], esterno -> se non esiste 5.34e9
nuCuR     = mat.BarCond.nu;                         % Poisson [-], esterno -> se non esiste 0.34345
muCuR     = 1;                                      % Permeabilità [-], esterno
epsCuR    = 1;                                      % Permittività [-], esterno
cpCuR     = 385;                                    % Calore specifico [J/(kg*K)], esterno
kCuR      = 400;                                    % Conducibilità termica [W/(m*K)], esterno
emisCuR   = 0.5;                                    % Emissività [-], esterno
rhoCuRC   = [num2str(rhoCuR,nDig) uRho];            % numero + u.d.m. in stringa per Comsol
sigCuRC   = [num2str(sigCuR,nDig) uSig];
rho0CuRC  = [num2str(rho0CuR,nDig) uRhoE];
alphaCuRC = [num2str(alphaCuR,nDig) uAlpha];
TrefCuRC  = [num2str(TrefCuR,nDig) uTemp];
ECuRC     = [num2str(ECuR,nDig) uE];
nuCuRC    = num2str(nuCuR,nDig);
muCuRC    = num2str(muCuR,nDig);
epsCuRC   = num2str(epsCuR,nDig);
cpCuRC    = [num2str(cpCuR,nDig) uCp];
kCuRC     = [num2str(kCuR,nDig) uK];
emisCuRC  = num2str(emisCuR,nDig);

% Magneti di rotore
tagPMR  = 'matPMR';                                 % Tag COMSOL
namePMR = [char(mat.LayerMag.MatName) ' - rotor'];  % Nome materiale
rhoPMR  = mat.LayerMag.kgm3;                        % Densità [kg/m^3]
sigPMR  = mat.LayerMag.sigmaPM;                     % Conducibilità [S/m]
muPMR   = mat.LayerMag.mu;                          % Permeabilità recoil [-]
BrPMR   = mat.LayerMag.Br;                          % Induzione residua [T]
HcPMR   = mat.LayerMag.Hc;                          % Campo coercitivo [A/m], futuro
EPMR    = [];                                       % Modulo di Young, assente in SyR-e
nuPMR   = [];                                       % Poisson, assente in SyR-e
epsPMR  = 1;                                        % Permittività [-], esterno
cpPMR   = 400;                                      % Calore specifico [J/(kg*K)], esterno
kPMR    = 7;                                        % Conducibilità termica [W/(m*K)], esterno
rhoPMRC = [num2str(rhoPMR,nDig) uRho];              % numero + u.d.m. in stringa per Comsol
sigPMRC = [num2str(sigPMR,nDig) uSig];
epsPMRC = num2str(epsPMR,nDig);
cpPMRC  = [num2str(cpPMR,nDig) uCp];
kPMRC   = [num2str(kPMR,nDig) uK];
zPMRC   = num2str(0,nDig);
muPMRC  = num2str(muPMR,nDig);
muPMRT  = {muPMRC,zPMRC,zPMRC,zPMRC,muPMRC,zPMRC,zPMRC,zPMRC,muPMRC};
BrPMRC  = {[num2str(BrPMR,nDig) uBr]};
EPMRC   = '';
nuPMRC  = '';
if ~isempty(EPMR) && ~isempty(nuPMR)
    EPMRC  = [num2str(EPMR,nDig) uE];
    nuPMRC = num2str(nuPMR,nDig);
end

% Filler di rotore
tagFillR = 'matFillR';                          % Tag COMSOL
if strcmp(geo.RotType,'EESM')
    nameFillR = 'Rotor fill';                   % Nome materiale
    rhoFillR  = 1000;                           % Densità [kg/m^3], esterno
    EFillR    = 5e9;                            % Modulo di Young [Pa], esterno
    nuFillR   = 0.34;                           % Poisson [-], esterno
    sigFillR  = 0;                              % Conducibilità [S/m], esterno
    muFillR   = 1;                              % Permeabilità [-], esterno
    epsFillR  = 1;                              % Permittività [-], esterno
    fillRNS   = false;                          % Filler solido
else
    nameFillR = [char(mat.LayerAir.MatName) ' - rotor fill']; % Nome materiale
    rhoFillR  = mat.LayerAir.kgm3;                            % Densità [kg/m^3]
    EFillR    = [];                                           % Non richiesto per aria
    nuFillR   = [];                                           % Non richiesto per aria
    sigFillR  = mat.SlotAir.sigma;                            % Conducibilità [S/m]
    muFillR   = mat.LayerAir.mu;                              % Permeabilità relativa [-]
    epsFillR  = 1;                                            % Permittività [-], esterno
    fillRNS   = true;                                         % Dominio non solido
end             
rhoFillRC = [num2str(rhoFillR,nDig) uRho];                    % numero + u.d.m. in stringa per Comsol
sigFillRC = [num2str(sigFillR,nDig) uSig];
muFillRC = num2str(muFillR,nDig);
epsFillRC = num2str(epsFillR,nDig);
if ~isempty(EFillR)
    EFillRC = [num2str(EFillR,nDig) uE];
    nuFillRC = num2str(nuFillR,nDig);
end

% Sleeve di rotore
tagSleeveR  = 'matSleeveR';                                 % Tag COMSOL
nameSleeveR = [char(mat.Sleeve.MatName) ' - rotor sleeve']; % Nome materiale
rhoSleeveR  = mat.Sleeve.kgm3;                              % Densità [kg/m^3]
ESleeveR    = mat.Sleeve.E*1e9;                             % Modulo di Young [Pa]
nuSleeveR   = [];                                           % Poisson, assente in SyR-e
sigSleeveR  = 0;                                            % Conducibilità [S/m], esterno
muSleeveR   = 1;                                            % Permeabilità [-], esterno
epsSleeveR  = 1;                                            % Permittività [-], esterno
rhoSleeveRC = [num2str(rhoSleeveR,nDig) uRho];              % numero + u.d.m. in stringa per Comsol
ESleeveRC   = [num2str(ESleeveR,nDig) uE];
sigSleeveRC = [num2str(sigSleeveR,nDig) uSig];
muSleeveRC  = num2str(muSleeveR,nDig);
epsSleeveRC = num2str(epsSleeveR,nDig);
if ~isempty(nuSleeveR)
    nuSleeveRC = num2str(nuSleeveR,nDig);
end

% Shaft
tagShaft  = 'matShaft';                          % Tag COMSOL
nameShaft = char(mat.Shaft.MatName);             % Nome materiale
rhoShaft  = mat.Shaft.kgm3;                      % Densità [kg/m^3]
EShaft    = [];                                  % Modulo di Young, assente in SyR-e
nuShaft   = [];                                  % Poisson, assente in SyR-e
sigShaft  = 0;                                   % Conducibilità [S/m], esterno
epsShaft  = 1;                                   % Permittività [-], esterno
shaftNS   = strcmpi(nameShaft,'ShaftAir');       % Albero non solido
rhoShaftC = [num2str(rhoShaft,nDig) uRho];       % numero + u.d.m. in stringa per Comsol
sigShaftC = [num2str(sigShaft,nDig) uSig];
epsShaftC = num2str(epsShaft,nDig);
if shaftNS
    muShaft  = mat.LayerAir.mu;                  % Permeabilità relativa [-]
    muShaftC = num2str(muShaft,nDig);
else
        bhShaft = num2cell([mat.Shaft.BH(:,2) mat.Shaft.BH(:,1)]);
    for ii = 1:numel(bhShaft)
        bhShaft{ii} = num2str(bhShaft{ii},nDig);
    end
    %bhShaft = vertcat(bhShaft{:});
end
if ~isempty(EShaft) && ~isempty(nuShaft)
    EShaftC  = [num2str(EShaft,nDig) uE];
    nuShaftC = num2str(nuShaft,nDig);
end

%STATORE ------------------------------------------------------------------
% Ferro di statore
tagFeS  = 'matFeS';                               % Tag COMSOL
nameFeS = [char(mat.Stator.MatName) ' - stator']; % Nome materiale
rhoFeS  = mat.Stator.kgm3;                        % Densità [kg/m^3]
EFeS    = mat.Stator.E*1e9;                       % Modulo di Young [Pa]
nuFeS   = 0.28;                                   % Poisson [-], esterno
sigFeS  = 0;                                      % Conducibilità [S/m], laminato
epsFeS  = 1;                                      % Permittività [-], esterno
rhoFeSC = [num2str(rhoFeS,nDig) uRho];            % numero + u.d.m. in stringa per Comsol
EFeSC   = [num2str(EFeS,nDig) uE];
nuFeSC  = num2str(nuFeS,nDig);
sigFeSC = [num2str(sigFeS,nDig) uSig];
epsFeSC = num2str(epsFeS,nDig);
bhFeS = num2cell([mat.Stator.BH(:,2) mat.Stator.BH(:,1)]);
for ii = 1:numel(bhFeS)
    bhFeS{ii} = num2str(bhFeS{ii},nDig);
end

% Conduttori di statore
tagCuS    = 'matCuS';                                 % Tag COMSOL
nameCuS   = [char(mat.SlotCond.MatName) ' - stator']; % Nome materiale
rhoCuS    = mat.SlotCond.kgm3*geo.win.kcu;            % Densità equivalente [kg/m^3]
sigCuS    = mat.SlotCond.sigma;                       % Conducibilità [S/m]
rho0CuS   = 1/sigCuS;                                 % Resistività a Tref [ohm*m]
alphaCuS  = mat.SlotCond.alpha;                       % Coefficiente termico [1/K]
TrefCuS   = 293.15;                                   % Temperatura riferimento [K], esterno
ECuS      = 126e9;                                    % Modulo di Young [Pa], esterno
nuCuS     = 0.34;                                     % Poisson [-], esterno
muCuS     = 1;                                        % Permeabilità [-], esterno
epsCuS    = 1;                                        % Permittività [-], esterno
cpCuS     = 385;                                      % Calore specifico [J/(kg*K)], esterno
kCuS      = 400;                                      % Conducibilità termica [W/(m*K)], esterno
emisCuS   = 0.5;                                      % Emissività [-], esterno
rhoCuSC   = [num2str(rhoCuS,nDig) uRho];              % numero + u.d.m. in stringa per Comsol
sigCuSC   = [num2str(sigCuS,nDig) uSig];
rho0CuSC  = [num2str(rho0CuS,nDig) uRhoE];
alphaCuSC = [num2str(alphaCuS,nDig) uAlpha];
TrefCuSC  = [num2str(TrefCuS,nDig) uTemp];
ECuSC     = [num2str(ECuS,nDig) uE];
nuCuSC    = num2str(nuCuS,nDig);
muCuSC    = num2str(muCuS,nDig);
epsCuSC   = num2str(epsCuS,nDig);
cpCuSC    = [num2str(cpCuS,nDig) uCp];
kCuSC     = [num2str(kCuS,nDig) uK];
emisCuSC  = num2str(emisCuS,nDig);

% Filler di statore
tagFillS  = 'matFillS';                                   % Tag COMSOL
nameFillS = [char(mat.SlotAir.MatName) ' - stator fill']; % Nome materiale
rhoFillS  = mat.SlotAir.kgm3;                             % Densità [kg/m^3]
sigFillS  = mat.SlotAir.sigma;                            % Conducibilità [S/m]
muFillS   = mat.LayerAir.mu;                              % Permeabilità relativa [-]
epsFillS  = 1;                                            % Permittività [-], esterno
rhoFillSC = [num2str(rhoFillS,nDig) uRho];
sigFillSC = [num2str(sigFillS,nDig) uSig];
muFillSC  = num2str(muFillS,nDig);
epsFillSC = num2str(epsFillS,nDig);

% Selection aggregata dei domini d'aria
tagAirDom  = 'domAir';                          % Tag selection
nameAirDom = 'Air domains';                     % Nome selection
airIn = {sel.airGap};
if isfield(sel,'fillSta')
    airIn{end+1} = sel.fillSta;
end
if isfield(sel,'fillRot') && fillRNS
    airIn{end+1} = sel.fillRot;
end
if isfield(sel,'shaft') && shaftNS
    airIn{end+1} = sel.shaft;
end
comp.selection.create(tagAirDom,selUnion);
comp.selection(tagAirDom).label(nameAirDom);
comp.selection(tagAirDom).set('entitydim',dimDom);
comp.selection(tagAirDom).set('input',airIn);
geo.infoComsol.tags.selection.domain.air = tagAirDom;

% ------------------------------------------------------------------------
% Creazione materiali in COMSOL e assegnazione ai domini
% ------------------------------------------------------------------------

% Air gap
comp.material.create(tagAir,matClass);
comp.material(tagAir).label(nameAir);
comp.material(tagAir).set('family',famAir);
comp.material(tagAir).materialType(matNS);
comp.material(tagAir).propertyGroup(pgDef).set('density',rhoAirC);
comp.material(tagAir).propertyGroup(pgDef).set('electricconductivity',sigAirC);
comp.material(tagAir).propertyGroup(pgDef).set('relpermeability',muAirC);
comp.material(tagAir).propertyGroup(pgDef).set('relpermittivity',epsAirC);
comp.material(tagAir).selection().named(sel.airGap);
geo.infoComsol.tags.material.airGap = tagAir;

% Ferro di rotore
if isfield(sel,'feRot')
    comp.material.create(tagFeR,matClass);
    comp.material(tagFeR).label(nameFeR);
    comp.material(tagFeR).propertyGroup(pgDef).set('density',rhoFeRC);
    comp.material(tagFeR).propertyGroup(pgDef).set('electricconductivity',sigFeRC);
    comp.material(tagFeR).propertyGroup(pgDef).set('relpermittivity',epsFeRC);
    comp.material(tagFeR).propertyGroup().create(pgEnu,typeEnu,descEnu);
    comp.material(tagFeR).propertyGroup(pgEnu).set('E',EFeRC);
    comp.material(tagFeR).propertyGroup(pgEnu).set('nu',nuFeRC);
    comp.material(tagFeR).propertyGroup().create(pgBH,typeBH);
    comp.material(tagFeR).propertyGroup(pgBH).func().create(funBH,typeInt);
    comp.material(tagFeR).propertyGroup(pgBH).func(funBH).set('table',bhFeR);
    comp.material(tagFeR).propertyGroup(pgBH).func(funBH).set('extrap',bhExt);
    comp.material(tagFeR).propertyGroup(pgBH).func(funBH).set('fununit',bhUnit);
    comp.material(tagFeR).propertyGroup(pgBH).func(funBH).set('argunit',hUnit);
    comp.material(tagFeR).propertyGroup(pgBH).func(funBH).set('defineinv',true);
    comp.material(tagFeR).propertyGroup(pgBH).func(funBH).set('defineprimfun',true);
    comp.material(tagFeR).propertyGroup(pgBH).set('normB',bhB);
    comp.material(tagFeR).propertyGroup(pgBH).set('normH',bhH);
    comp.material(tagFeR).propertyGroup(pgBH).set('Wpm',bhW);
    comp.material(tagFeR).propertyGroup(pgBH).addInput(inH);
    comp.material(tagFeR).propertyGroup(pgBH).addInput(inB);
    comp.material(tagFeR).selection().named(sel.feRot);
    geo.infoComsol.tags.material.feRot = tagFeR;
end

% Conduttori di rotore
if isfield(sel,'cuRot')
    comp.material.create(tagCuR,matClass);
    comp.material(tagCuR).label(nameCuR);
    comp.material(tagCuR).set('family',famCu);
    comp.material(tagCuR).propertyGroup(pgDef).set('density',rhoCuRC);
    comp.material(tagCuR).propertyGroup(pgDef).set('electricconductivity',sigCuRC);
    comp.material(tagCuR).propertyGroup(pgDef).set('relpermeability',muCuRC);
    comp.material(tagCuR).propertyGroup(pgDef).set('relpermittivity',epsCuRC);
    comp.material(tagCuR).propertyGroup(pgDef).set('heatcapacity',cpCuRC);
    comp.material(tagCuR).propertyGroup(pgDef).set('thermalconductivity',kCuRC);
    comp.material(tagCuR).propertyGroup(pgDef).set('emissivity',emisCuRC);
    comp.material(tagCuR).propertyGroup().create(pgEnu,typeEnu,descEnu);
    comp.material(tagCuR).propertyGroup(pgEnu).set('E',ECuRC);
    comp.material(tagCuR).propertyGroup(pgEnu).set('nu',nuCuRC);
    comp.material(tagCuR).propertyGroup().create(pgRes,typeRes);
    comp.material(tagCuR).propertyGroup(pgRes).set('rho0',rho0CuRC);
    comp.material(tagCuR).propertyGroup(pgRes).set('alpha',alphaCuRC);
    comp.material(tagCuR).propertyGroup(pgRes).set('Tref',TrefCuRC);
    comp.material(tagCuR).propertyGroup(pgRes).addInput(inTemp);
    comp.material(tagCuR).selection().named(sel.cuRot);
    geo.infoComsol.tags.material.cuRot = tagCuR;
end

% Magneti di rotore
if isfield(sel,'PMRot')
    comp.material.create(tagPMR,matClass);
    comp.material(tagPMR).label(namePMR);
    comp.material(tagPMR).propertyGroup(pgDef).set('density',rhoPMRC);
    comp.material(tagPMR).propertyGroup(pgDef).set('electricconductivity',sigPMRC);
    comp.material(tagPMR).propertyGroup(pgDef).set('relpermittivity',epsPMRC);
    comp.material(tagPMR).propertyGroup(pgDef).set('heatcapacity',cpPMRC);
    comp.material(tagPMR).propertyGroup(pgDef).set('thermalconductivity',kPMRC);
    comp.material(tagPMR).propertyGroup().create(pgPM,typePM);
    comp.material(tagPMR).propertyGroup(pgPM).set('murec',muPMRT);
    comp.material(tagPMR).propertyGroup(pgPM).set('normBr',BrPMRC);
    if ~isempty(EPMR) && ~isempty(nuPMR)
        comp.material(tagPMR).propertyGroup().create(pgEnu,typeEnu,descEnu);
        comp.material(tagPMR).propertyGroup(pgEnu).set('E',EPMRC);
        comp.material(tagPMR).propertyGroup(pgEnu).set('nu',nuPMRC);
    end
    comp.material(tagPMR).selection().named(sel.PMRot);
    geo.infoComsol.tags.material.PMRot = tagPMR;
end

% Filler di rotore
if isfield(sel,'fillRot')
    comp.material.create(tagFillR,matClass);
    comp.material(tagFillR).label(nameFillR);
    comp.material(tagFillR).propertyGroup(pgDef).set('density',rhoFillRC);
    comp.material(tagFillR).propertyGroup(pgDef).set('electricconductivity',sigFillRC);
    comp.material(tagFillR).propertyGroup(pgDef).set('relpermeability',muFillRC);
    comp.material(tagFillR).propertyGroup(pgDef).set('relpermittivity',epsFillRC);
    if fillRNS
        comp.material(tagFillR).set('family',famAir);
        comp.material(tagFillR).materialType(matNS);
    else
        comp.material(tagFillR).set('family',famFill);
        comp.material(tagFillR).propertyGroup().create(pgEnu,typeEnu,descEnu);
        comp.material(tagFillR).propertyGroup(pgEnu).set('E',EFillRC);
        comp.material(tagFillR).propertyGroup(pgEnu).set('nu',nuFillRC);
    end
    comp.material(tagFillR).selection().named(sel.fillRot);
    geo.infoComsol.tags.material.fillRot = tagFillR;
end

% Sleeve di rotore
if isfield(sel,'sleeveRot')
    comp.material.create(tagSleeveR,matClass);
    comp.material(tagSleeveR).label(nameSleeveR);
    comp.material(tagSleeveR).propertyGroup(pgDef).set('density',rhoSleeveRC);
    comp.material(tagSleeveR).propertyGroup(pgDef).set('electricconductivity',sigSleeveRC);
    comp.material(tagSleeveR).propertyGroup(pgDef).set('relpermeability',muSleeveRC);
    comp.material(tagSleeveR).propertyGroup(pgDef).set('relpermittivity',epsSleeveRC);
    if ~isempty(nuSleeveR)
        comp.material(tagSleeveR).propertyGroup().create(pgEnu,typeEnu,descEnu);
        comp.material(tagSleeveR).propertyGroup(pgEnu).set('E',ESleeveRC);
        comp.material(tagSleeveR).propertyGroup(pgEnu).set('nu',nuSleeveRC);
    end
    comp.material(tagSleeveR).selection().named(sel.sleeveRot);
    geo.infoComsol.tags.material.sleeveRot = tagSleeveR;
end

% Shaft
if isfield(sel,'shaft')
    comp.material.create(tagShaft,matClass);
    comp.material(tagShaft).label(nameShaft);
    comp.material(tagShaft).propertyGroup(pgDef).set('density',rhoShaftC);
    comp.material(tagShaft).propertyGroup(pgDef).set('electricconductivity',sigShaftC);
    comp.material(tagShaft).propertyGroup(pgDef).set('relpermittivity',epsShaftC);
    if shaftNS
        comp.material(tagShaft).set('family',famAir);
        comp.material(tagShaft).materialType(matNS);
        comp.material(tagShaft).propertyGroup(pgDef).set('relpermeability',muShaftC);
    else
        comp.material(tagShaft).propertyGroup().create(pgBH,typeBH);
        comp.material(tagShaft).propertyGroup(pgBH).func().create(funBH,typeInt);
        comp.material(tagShaft).propertyGroup(pgBH).func(funBH).set('table',bhShaft);
        comp.material(tagShaft).propertyGroup(pgBH).func(funBH).set('extrap',bhExt);
        comp.material(tagShaft).propertyGroup(pgBH).func(funBH).set('fununit',bhUnit);
        comp.material(tagShaft).propertyGroup(pgBH).func(funBH).set('argunit',hUnit);
        comp.material(tagShaft).propertyGroup(pgBH).func(funBH).set('defineinv',true);
        comp.material(tagShaft).propertyGroup(pgBH).func(funBH).set('defineprimfun',true);
        comp.material(tagShaft).propertyGroup(pgBH).set('normB',bhB);
        comp.material(tagShaft).propertyGroup(pgBH).set('normH',bhH);
        comp.material(tagShaft).propertyGroup(pgBH).set('Wpm',bhW);
        comp.material(tagShaft).propertyGroup(pgBH).addInput(inH);
        comp.material(tagShaft).propertyGroup(pgBH).addInput(inB);
        if ~isempty(EShaft) && ~isempty(nuShaft)
            comp.material(tagShaft).propertyGroup().create(pgEnu,typeEnu,descEnu);
            comp.material(tagShaft).propertyGroup(pgEnu).set('E',EShaftC);
            comp.material(tagShaft).propertyGroup(pgEnu).set('nu',nuShaftC);
        end
    end
    comp.material(tagShaft).selection().named(sel.shaft);
    geo.infoComsol.tags.material.shaft = tagShaft;
end

% Ferro di statore
if isfield(sel,'feSta')
    comp.material.create(tagFeS,matClass);
    comp.material(tagFeS).label(nameFeS);
    comp.material(tagFeS).propertyGroup(pgDef).set('density',rhoFeSC);
    comp.material(tagFeS).propertyGroup(pgDef).set('electricconductivity',sigFeSC);
    comp.material(tagFeS).propertyGroup(pgDef).set('relpermittivity',epsFeSC);
    comp.material(tagFeS).propertyGroup().create(pgEnu,typeEnu,descEnu);
    comp.material(tagFeS).propertyGroup(pgEnu).set('E',EFeSC);
    comp.material(tagFeS).propertyGroup(pgEnu).set('nu',nuFeSC);
    comp.material(tagFeS).propertyGroup().create(pgBH,typeBH);
    comp.material(tagFeS).propertyGroup(pgBH).func().create(funBH,typeInt);
    comp.material(tagFeS).propertyGroup(pgBH).func(funBH).set('table',bhFeS);
    comp.material(tagFeS).propertyGroup(pgBH).func(funBH).set('extrap',bhExt);
    comp.material(tagFeS).propertyGroup(pgBH).func(funBH).set('fununit',bhUnit);
    comp.material(tagFeS).propertyGroup(pgBH).func(funBH).set('argunit',hUnit);
    comp.material(tagFeS).propertyGroup(pgBH).func(funBH).set('defineinv',true);
    comp.material(tagFeS).propertyGroup(pgBH).func(funBH).set('defineprimfun',true);
    comp.material(tagFeS).propertyGroup(pgBH).set('normB',bhB);
    comp.material(tagFeS).propertyGroup(pgBH).set('normH',bhH);
    comp.material(tagFeS).propertyGroup(pgBH).set('Wpm',bhW);
    comp.material(tagFeS).propertyGroup(pgBH).addInput(inH);
    comp.material(tagFeS).propertyGroup(pgBH).addInput(inB);
    comp.material(tagFeS).selection().named(sel.feSta);
    geo.infoComsol.tags.material.feSta = tagFeS;
end

% Conduttori di statore
if isfield(sel,'cuSta')
    comp.material.create(tagCuS,matClass);
    comp.material(tagCuS).label(nameCuS);
    comp.material(tagCuS).set('family',famCu);
    comp.material(tagCuS).propertyGroup(pgDef).set('density',rhoCuSC);
    comp.material(tagCuS).propertyGroup(pgDef).set('electricconductivity',sigCuSC);
    comp.material(tagCuS).propertyGroup(pgDef).set('relpermeability',muCuSC);
    comp.material(tagCuS).propertyGroup(pgDef).set('relpermittivity',epsCuSC);
    comp.material(tagCuS).propertyGroup(pgDef).set('heatcapacity',cpCuSC);
    comp.material(tagCuS).propertyGroup(pgDef).set('thermalconductivity',kCuSC);
    comp.material(tagCuS).propertyGroup(pgDef).set('emissivity',emisCuSC);
    comp.material(tagCuS).propertyGroup().create(pgEnu,typeEnu,descEnu);
    comp.material(tagCuS).propertyGroup(pgEnu).set('E',ECuSC);
    comp.material(tagCuS).propertyGroup(pgEnu).set('nu',nuCuSC);
    comp.material(tagCuS).propertyGroup().create(pgRes,typeRes);
    comp.material(tagCuS).propertyGroup(pgRes).set('rho0',rho0CuSC);
    comp.material(tagCuS).propertyGroup(pgRes).set('alpha',alphaCuSC);
    comp.material(tagCuS).propertyGroup(pgRes).set('Tref',TrefCuSC);
    comp.material(tagCuS).propertyGroup(pgRes).addInput(inTemp);
    comp.material(tagCuS).selection().named(sel.cuSta);
    geo.infoComsol.tags.material.cuSta = tagCuS;
end

% Filler di statore
if isfield(sel,'fillSta')
    comp.material.create(tagFillS,matClass);
    comp.material(tagFillS).label(nameFillS);
    comp.material(tagFillS).set('family',famAir);
    comp.material(tagFillS).materialType(matNS);
    comp.material(tagFillS).propertyGroup(pgDef).set('density',rhoFillSC);
    comp.material(tagFillS).propertyGroup(pgDef).set('electricconductivity',sigFillSC);
    comp.material(tagFillS).propertyGroup(pgDef).set('relpermeability',muFillSC);
    comp.material(tagFillS).propertyGroup(pgDef).set('relpermittivity',epsFillSC);
    comp.material(tagFillS).selection().named(sel.fillSta);
    geo.infoComsol.tags.material.fillSta = tagFillS;
end


%  ========================================================================
%  Selezioni dei bordi
%  ========================================================================

tolBnd = max(1e-6,1e-5*r_si);

bAdjR = double(mphgetadj(model,geometryTag,'boundary','domain',dRot));
bAdjS = double(mphgetadj(model,geometryTag,'boundary','domain',dSta));

bAdjR = unique(bAdjR(:)');
bAdjS = unique(bAdjS(:)');

bR1 = [];
bR2 = [];
bS1 = [];
bS2 = [];


% Lati radiali del rotore
if geo.ps < 2*geo.p

    for ii = 1:numel(bAdjR)

        id = bAdjR(ii);

        dAdj = double(mphgetadj(model,geometryTag,'domain','boundary',id));
        dAdj = unique(dAdj(:)');

        % I lati di settore sono bordi esterni, non interfacce tra domini
        if numel(dAdj) ~= 1
            continue
        end

        pnt = double(mphgetadj(model,geometryTag,'point','boundary',id));
        pnt = unique(pnt(:)');

        xy = zeros(numel(pnt),2);

        for jj = 1:numel(pnt)
            coord = double(mphgetcoords(model,geometryTag,'point',pnt(jj)));
            coord = coord(:);
            xy(jj,:) = coord(1:2)';
        end

        dist1 = abs(xy(:,2));
        proj1 = xy(:,1);

        dist2 = abs(xy(:,1)*sin(sectorAngle)-xy(:,2)*cos(sectorAngle));
        proj2 = xy(:,1)*cos(sectorAngle)+xy(:,2)*sin(sectorAngle);

        if all(dist1 <= tolBnd) && all(proj1 >= -tolBnd)
            bR1 = [bR1 id];
        elseif all(dist2 <= tolBnd) && all(proj2 >= -tolBnd)
            bR2 = [bR2 id];
        end

    end


    % Lati radiali dello statore
    for ii = 1:numel(bAdjS)

        id = bAdjS(ii);

        dAdj = double(mphgetadj(model,geometryTag,'domain','boundary',id));
        dAdj = unique(dAdj(:)');

        if numel(dAdj) ~= 1
            continue
        end

        pnt = double(mphgetadj(model,geometryTag,'point','boundary',id));
        pnt = unique(pnt(:)');

        xy = zeros(numel(pnt),2);

        for jj = 1:numel(pnt)
            coord = double(mphgetcoords(model,geometryTag,'point',pnt(jj)));
            coord = coord(:);
            xy(jj,:) = coord(1:2)';
        end

        dist1 = abs(xy(:,2));
        proj1 = xy(:,1);

        dist2 = abs(xy(:,1)*sin(sectorAngle)-xy(:,2)*cos(sectorAngle));
        proj2 = xy(:,1)*cos(sectorAngle)+xy(:,2)*sin(sectorAngle);

        if all(dist1 <= tolBnd) && all(proj1 >= -tolBnd)
            bS1 = [bS1 id];
        elseif all(dist2 <= tolBnd) && all(proj2 >= -tolBnd)
            bS2 = [bS2 id];
        end

    end

end

bR1 = unique(bR1);
bR2 = unique(bR2);
bS1 = unique(bS1);
bS2 = unique(bS2);

bRSides = unique([bR1 bR2]);
bSSides = unique([bS1 bS2]);
bSides = unique([bRSides bSSides]);


% Interfaccia tra le due metà del traferro
bAdjAgR = double(mphgetadj(model,geometryTag,'boundary','domain',dAgR));
bAdjAgS = double(mphgetadj(model,geometryTag,'boundary','domain',dAgS));

bAdjAgR = unique(bAdjAgR(:)');
bAdjAgS = unique(bAdjAgS(:)');

bAgR = [];
bAgS = [];

for ii = 1:numel(bAdjAgR)

    id = bAdjAgR(ii);

    pnt = double(mphgetadj(model,geometryTag,'point','boundary',id));
    pnt = unique(pnt(:)');

    xy = zeros(numel(pnt),2);

    for jj = 1:numel(pnt)
        coord = double(mphgetcoords(model,geometryTag,'point',pnt(jj)));
        coord = coord(:);
        xy(jj,:) = coord(1:2)';
    end

    radius = hypot(xy(:,1),xy(:,2));

    if all(abs(radius-r_ag) <= tolBnd)
        bAgR = [bAgR id];
    end

end

for ii = 1:numel(bAdjAgS)

    id = bAdjAgS(ii);

    pnt = double(mphgetadj(model,geometryTag,'point','boundary',id));
    pnt = unique(pnt(:)');

    xy = zeros(numel(pnt),2);

    for jj = 1:numel(pnt)
        coord = double(mphgetcoords(model,geometryTag,'point',pnt(jj)));
        coord = coord(:);
        xy(jj,:) = coord(1:2)';
    end

    radius = hypot(xy(:,1),xy(:,2));

    if all(abs(radius-r_ag) <= tolBnd)
        bAgS = [bAgS id];
    end

end

bAgR = unique(bAgR);
bAgS = unique(bAgS);
bAg = unique([bAgR bAgS]);


% Bore rotorico e statorico
bBoreR = setdiff(bAdjAgR,unique([bAgR bRSides]));
bBoreS = setdiff(bAdjAgS,unique([bAgS bSSides]));


% Bordo esterno dello statore
bndSta = geo.BLKLABELS.statore.boundary(:,1:3);
pOutS = bndSta(all(isfinite(bndSta(:,1:2)),2) & bndSta(:,3)==0,1:2);

comp.selection.create('findBnd','Disk');
findBnd = comp.selection('findBnd');
findBnd.set('entitydim',1);
findBnd.set('r',max(1e-6,geo.g/100));

bOutS = [];

for ii = 1:size(pOutS,1)
    findBnd.set('posx',pOutS(ii,1));
    findBnd.set('posy',pOutS(ii,2));
    id = double(findBnd.entities());
    bOutS = [bOutS id(:)'];
end

bOutS = setdiff(unique(bOutS),unique([bSSides bBoreS bAgS]));


% Interfaccia dello shaft
bShaft = [];

if ~isempty(dShaft)
    bShaft = double(mphgetadj(model,geometryTag,'boundary','domain',dShaft));
    bShaft = setdiff(unique(bShaft(:)'),bRSides);
end


% Interfaccia interna dello sleeve
bSleeveIn = [];

if ~isempty(dSleeve)
    bSleeveIn = double(mphgetadj(model,geometryTag,'boundary','domain',dSleeve));
    bSleeveIn = setdiff(unique(bSleeveIn(:)'),unique([bRSides bBoreR]));
end


% Interfaccia tra polo rotorico e bobina rotorica
bPoleCoilFe = [];
bPoleCoilCu = [];

if ~isempty(dFeR) && ~isempty(dCuR)

    bAdjFeR = double(mphgetadj(model,geometryTag,'boundary','domain',dFeR));
    bAdjCuR = double(mphgetadj(model,geometryTag,'boundary','domain',dCuR));

    bAdjFeR = unique(bAdjFeR(:)');
    bAdjCuR = unique(bAdjCuR(:)');

end


% Contatto tra bobina rotorica e filler rotorico
bCuFillR = [];

if ~isempty(dCuR) && ~isempty(dFillR)

    if ~exist('bAdjCuR','var')
        bAdjCuR = double(mphgetadj(model,geometryTag,'boundary','domain',dCuR));
        bAdjCuR = unique(bAdjCuR(:)');
    end

    bAdjFillR = double(mphgetadj(model,geometryTag,'boundary','domain',dFillR));
    bAdjFillR = unique(bAdjFillR(:)');

    bCuFillR = intersect(bAdjCuR,bAdjFillR);

end

% Lato bobina: esclude filler e lati radiali
if ~isempty(bAdjCuR)
    bPoleCoilCu = setdiff(bAdjCuR,unique([bCuFillR bRSides]));
end

% Lato ferro: esclude bordi che non possono appartenere al contatto
if ~isempty(bAdjFeR)
    bPoleCoilFe = setdiff(bAdjFeR,unique([bRSides bBoreR bShaft bSleeveIn bAgR]));
end


% Creazione delle selections
selBnd = {
    'sideRot1',    'bndR1',        'Rotor sector side 1',             bR1;
    'sideRot2',    'bndR2',        'Rotor sector side 2',             bR2;
    'sideRot',     'bndRSides',    'Rotor sector sides',              bRSides;
    'boreRot',     'bndBoreR',     'Rotor bore',                      bBoreR;
    'shaft',       'bndShaft',     'Shaft interface',                 bShaft;
    'sleeveInRot', 'bndSleeveIn',  'Rotor sleeve inner interface',    bSleeveIn;
    'poleCoilRotFe', 'bndPoleCoilFe', 'Rotor pole-winding interface - iron',    bPoleCoilFe;
    'poleCoilRotCu', 'bndPoleCoilCu', 'Rotor pole-winding interface - winding', bPoleCoilCu;
    'coilFillRot', 'bndCuFillR',   'Rotor winding-fill interface',    bCuFillR;

    'sideSta1',    'bndS1',        'Stator sector side 1',            bS1;
    'sideSta2',    'bndS2',        'Stator sector side 2',            bS2;
    'sideSta',     'bndSSides',    'Stator sector sides',             bSSides;
    'boreSta',     'bndBoreS',     'Stator bore',                     bBoreS;
    'outerSta',    'bndOutS',      'Stator outer boundary',           bOutS;

    'sideAll',     'bndSides',     'All sector sides',                bSides;
    'airGapRot',   'bndAgR',       'Rotor air-gap interface',         bAgR;
    'airGapSta',   'bndAgS',       'Stator air-gap interface',        bAgS;
    'airGap',      'bndAg',        'Air-gap interfaces',              bAg;
    };

for ii = 1:size(selBnd,1)

    if ~isempty(selBnd{ii,4})

        fld = selBnd{ii,1};
        tag = selBnd{ii,2};

        comp.selection.create(tag,'Explicit').geom(1);
        comp.selection(tag).label(selBnd{ii,3});
        comp.selection(tag).set(selBnd{ii,4});

        geo.infoComsol.tags.selection.boundary.(fld) = tag;

    end

end

comp.selection.remove('findBnd');


%  ========================================================================
%  Raggruppamento delle selections
%  ========================================================================

grpDomTag = 'grpDom';
grpBndTag = 'grpBnd';

model.nodeGroup.create(grpDomTag,'Definitions',componentTag);
model.nodeGroup(grpDomTag).set('type','selection');
model.nodeGroup(grpDomTag).label('Domain selections');

fldDom = fieldnames(geo.infoComsol.tags.selection.domain);

for ii = 1:numel(fldDom)
    tag = geo.infoComsol.tags.selection.domain.(fldDom{ii});
    model.nodeGroup(grpDomTag).add('selection',tag);
end

model.nodeGroup.create(grpBndTag,'Definitions',componentTag);
model.nodeGroup(grpBndTag).set('type','selection');
model.nodeGroup(grpBndTag).label('Boundary selections');

fldBnd = fieldnames(geo.infoComsol.tags.selection.boundary);

for ii = 1:numel(fldBnd)
    tag = geo.infoComsol.tags.selection.boundary.(fldBnd{ii});
    model.nodeGroup(grpBndTag).add('selection',tag);
end


%  ========================================================================
%  Controllo finale dei materiali
%  ========================================================================

matDom = [dAg dFeR dCuR dPM dFillR dSleeve dShaft dFeS dCuS dFillS];
matDom = matDom(:)';

allDom = 1:double(geom.getNDomains());

matOverlap = numel(matDom) ~= numel(unique(matDom));
missDom = setdiff(allDom,unique(matDom));


%  ========================================================================
%  Aggiornamento infoComsol
%  ========================================================================

geo.infoComsol.tags.group.domain = grpDomTag;
geo.infoComsol.tags.group.boundary = grpBndTag;

geo.infoComsol.geometry.numberOfDomains = double(geom.getNDomains());
geo.infoComsol.geometry.numberOfBoundaries = double(geom.getNBoundaries());

geo.infoComsol.check.materialOverlap = matOverlap;
geo.infoComsol.check.missingMaterialDomains = missDom;
geo.infoComsol.check.materialsComplete = ~matOverlap && isempty(missDom);
geo.infoComsol.check.numberOfDomainSelections = numel(fldDom);
geo.infoComsol.check.numberOfBoundarySelections = numel(fldBnd);

geo.infoComsol.files.model = modelFile;

if matOverlap
    error('draw_motor_in_COMSOL:MaterialOverlap','Uno o più domini risultano assegnati a materiali differenti.');
end

if ~isempty(missDom)
    error('draw_motor_in_COMSOL:MissingMaterial','I domini %s non hanno un materiale assegnato.',mat2str(missDom));
end

mphsave(model,modelFile);
%mphlaunch(model)

end