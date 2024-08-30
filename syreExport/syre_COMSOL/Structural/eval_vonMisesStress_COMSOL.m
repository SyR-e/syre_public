% Copyright 2024
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

function [newDir] = eval_vonMisesStress_COMSOL(dataSet)

pathname = dataSet.currentpathname;
filename = dataSet.currentfilename;
load([pathname filename])
[~, file_name, ~] = fileparts(filename);

csv_dir = 'C:\Users\S296193\Desktop\syre_developers_20240610\syreCustomFeatures\COMSOL_functions\Structural\csv';

Stress_csv = fullfile(csv_dir, [file_name, '_Stress-csv']);

% materialCodes;

evalSpeed = dataSet.EvalSpeed;
if evalSpeed==0
    error('Please set a speed higher than zero!!!')
end

clc

resFolder = [filename(1:end-4) '_results\FEA results\'];
if ~exist([pathname resFolder],'dir')
    mkdir([pathname resFolder]);
end

newDir = ['structural_' int2str(evalSpeed) 'rpm\'];

newDir = [pathname resFolder newDir];
mkdir(newDir);

import com.comsol.model.*
import com.comsol.model.util.*
model = mphopen([pathname filename(1:end-4), '.mph']);

% Setting component2
model.component().create('comp2', true);
model.component('comp2').geom().create('geom2', 2);
model.component('comp2').mesh().create('mesh2');
model.component('comp2').geom('geom2').lengthUnit('mm');

% Setting geometry
magn_dxf = fullfile('C:\Users\S296193\Desktop\syre_developers_20240610\motorExamples\syreDefaultMotor_dxf\', [file_name, '_magn.dxf']);
rot_mecc_dxf = fullfile('C:\Users\S296193\Desktop\syre_developers_20240610\motorExamples\syreDefaultMotor_dxf\', [file_name, '_rot_mecc.dxf']);
flux_barr_dxf = fullfile('C:\Users\S296193\Desktop\syre_developers_20240610\motorExamples\syreDefaultMotor_dxf\', [file_name, '_flux_barr.dxf']);

% magn_dxf = fullfile(dxf_dir, [file_name, '_magn.dxf']);
% rot_mecc_dxf = fullfile(dxf_dir, [file_name, '_rot_mecc.dxf']);
% flux_barr_dxf = fullfile(dxf_dir, [file_name, '_flux_barr.dxf']);

model.component('comp2').geom('geom2').create('imp1', 'Import');
model.component('comp2').geom('geom2').feature('imp1').set('filename', rot_mecc_dxf);
model.component('comp2').geom('geom2').run('imp1');
model.component('comp2').geom('geom2').create('imp2', 'Import');
model.component('comp2').geom('geom2').feature('imp2').set('filename', flux_barr_dxf);
model.component('comp2').geom('geom2').run('imp2');
model.component('comp2').geom('geom2').create('dif1', 'Difference');
model.component('comp2').geom('geom2').feature('dif1').selection('input').set('imp1');
model.component('comp2').geom('geom2').feature('dif1').selection('input2').set('imp2');
model.component('comp2').geom('geom2').run('dif1');
model.component('comp2').geom('geom2').create('imp3', 'Import');
model.component('comp2').geom('geom2').feature('imp3').set('filename', magn_dxf);
model.component('comp2').geom('geom2').run('imp3');
model.component('comp2').geom('geom2').feature('fin').set('action', 'assembly');
model.component('comp2').geom('geom2').run('fin');

% Setting selections
model.component('comp2').selection().create('disk3', 'Disk');
disk3 = model.component('comp2').selection('disk3');
disk3 = model.component('comp2').selection('disk3').set('entitydim', 2);

model.component('comp2').selection().create('disk4', 'Disk');
disk4 = model.component('comp2').selection('disk4');
disk4 = model.component('comp2').selection('disk4').set('entitydim', 1);
disk4 = model.component('comp2').selection('disk4').set('r', 0.1);

% Defining study
model.study().create('std2');
model.study('std2').create('stat', 'Stationary');
model.study('std2').feature('stat').setSolveFor('/physics/rmm', true);
model.study('std2').feature('stat').setSolveFor('/physics/cir', true);
model.study('std2').feature('stat').setEntry('activate', 'cir', false);
model.study('std2').feature('stat').setEntry('activate', 'rmm', false);
model.component('comp2').physics().create('solid', 'SolidMechanics', 'geom2');
model.component('comp2').physics('solid').prop('d').set('d', 0.001);
% model.study('std1').feature('stat').setSolveFor('/physics/solid', true);
% model.study('std1').feature('time').setSolveFor('/physics/solid', true);
% model.study('std1').feature('emloss').setSolveFor('/physics/solid', true);
model.study('std2').feature('stat').setSolveFor('/physics/solid', true);
% model.study('std1').feature('stat').setEntry('activate', 'solid', false);
% model.study('std1').feature('time').setEntry('activate', 'solid', false);

% definizione materiali

%BH Curve Ferro
BH_curve = [
    0, 0;
    20, 0.064356436;
    22.5, 0.0817476280729167;
    25, 0.0995255775833333;
    27.5, 0.11807704196875;
    30, 0.137788778666667;
    32.5, 0.159008869369792;
    35, 0.181930692791667;
    37.5, 0.206708951901042;
    40, 0.233498349666667;
    42.5260416666667, 0.262646967765625;
    45.2083333333333, 0.295276402708333;
    48.203125, 0.332701629713542;
    51.6666666666667, 0.376237624;
    55.7291666666667, 0.426657900398437;
    60.4166666666667, 0.4825701321875;
    65.7291666666667, 0.542040532257812;
    71.6666666666667, 0.6031353135;
    78.3854166666667, 0.664475041114583;
    86.6666666666667, 0.726897689541667;
    97.4479166666667, 0.79179558553125;
    111.666666666667, 0.860561055833333;
    129.921875, 0.933593749796875;
    151.458333333333, 1.00732260716667;
    175.182291666667, 1.07718389028646;
    200, 1.1386138615;
    225.260416666667, 1.18827351500781;
    252.083333333333, 1.2277227724375;
    282.03125, 1.25974628727344;
    316.666666666667, 1.287128713;
    357.03125, 1.31216481035677;
    402.083333333333, 1.33518976910417;
    450.260416666667, 1.35604888625781;
    500, 1.37458745883333;
    550.260416666667, 1.39075391916667;
    602.083333333333, 1.404909240875;
    657.03125, 1.41751753289583;
    716.666666666667, 1.42904290416667;
    783.854166666667, 1.4400139231875;
    866.666666666667, 1.45121699670833;
    974.479166666667, 1.46350299104167;
    1116.66666666667, 1.4777227725;
    1299.21875, 1.49437912571615;
    1514.58333333333, 1.51258250860417;
    1751.82291666667, 1.53109529739844;
    2000, 1.54867986833333;
    2252.60416666667, 1.56447246315885;
    2520.83333333333, 1.5791047856875;
    2820.3125, 1.5935824052474;
    3166.66666666667, 1.60891089116667;
    3570.3125, 1.6257606229974;
    4020.83333333333, 1.6434612211875;
    4502.60416666667, 1.66100711640885;
    5000, 1.67739273933333;
    5502.60416666667, 1.69190903466667;
    6020.83333333333, 1.70503300325;
    6570.3125, 1.71753815995833;
    7166.66666666667, 1.73019801966667;
    7838.54166666667, 1.7437474215;
    8666.66666666667, 1.75876650158333;
    9744.79166666667, 1.77579672029167;
    11166.6666666667, 1.795379538;
    12992.1875, 1.81770833340625;
    15145.8333333333, 1.8415841585;
    17518.2291666667, 1.86545998359375;
    20000, 1.887788779;
    22552.0833333333, 1.9073844886224;
    25416.6666666667, 1.92450495072917;
    28906.25, 1.93976897717969;
    33333.3333333333, 1.95379537983333;
    38843.7447916667, 1.96701632833879;
    44916.625, 1.97911742350197;
    50864.4427083333, 1.9895976239181;
    55999.6666666667, 1.99795588818242;
    60093.015625, 2.00423728388627;
    64748.2083333333, 2.01067131460551;
    72027.2135416667, 2.02003359291211;
    83992, 2.03509973137807;
    102079.71875, 2.05782942911622;
    125228.25, 2.08691873140287;
    151750.65625, 2.12024777005519;
    179960, 2.15569667689034
];

Bf_max = max(BH_curve(:,2));

BH_curve_cell = arrayfun(@(x, y) {num2str(x), num2str(y)}, BH_curve(:, 1), BH_curve(:, 2), 'UniformOutput', false);    % Converte la matrice BH_curve in una cella di stringhe con due colonne
BH_curve_table = vertcat(BH_curve_cell{:});                                                                            % Converte la cella di stringhe in una matrice di stringhe

% N52 PM
model.component('comp2').material().create('mat5', 'Common');
model.component('comp2').material('mat5').propertyGroup().create('RemanentFluxDensity', 'Remanent flux density');
model.component('comp2').material('mat5').label('N52 (Sintered NdFeB)');
model.component('comp2').material('mat5').set('family', 'chrome');
model.component('comp2').material('mat5').propertyGroup('def').set('electricconductivity', {'1/1.4[uohm*m]', '0', '0', '0', '1/1.4[uohm*m]', '0', '0', '0', '1/1.4[uohm*m]'});
model.component('comp2').material('mat5').propertyGroup('def').set('relpermittivity', {'1', '0', '0', '0', '1', '0', '0', '0', '1'});
model.component('comp2').material('mat5').propertyGroup('RemanentFluxDensity').set('murec', {'1.05', '0', '0', '0', '1.05', '0', '0', '0', '1.05'});
model.component('comp2').material('mat5').propertyGroup('RemanentFluxDensity').set('normBr', '1.44[T]');
model.component('comp2').material('mat5').propertyGroup().create('Enu', 'Youngs_modulus_and_Poissons_ratio');
model.component('comp2').material('mat5').propertyGroup('Enu').set('E', '175e9');
model.component('comp2').material('mat5').propertyGroup('Enu').set('nu', '0.24');
model.component('comp2').material('mat5').propertyGroup('def').set('density', '7550');

% NGO 35PN270
model.component('comp2').material().create('mat6', 'Common');
model.component('comp2').material('mat6').propertyGroup().create('BHCurve', 'B-H Curve');
model.component('comp2').material('mat6').propertyGroup('BHCurve').func().create('BH', 'Interpolation');
model.component('comp2').material('mat6').label('Silicon Steel NGO 35PN270');
model.component('comp2').material('mat6').propertyGroup('def').set('electricconductivity', {'1.851852[MS/m]', '0', '0', '0', '1.851852[MS/m]', '0', '0', '0', '1.851852[MS/m]'});
model.component('comp2').material('mat6').propertyGroup('def').set('relpermittivity', {'1[1]', '0', '0', '0', '1[1]', '0', '0', '0', '1[1]'});
model.component('comp2').material('mat6').propertyGroup('BHCurve').label('B-H Curve');
model.component('comp2').material('mat6').propertyGroup('BHCurve').func('BH').label('Interpolation 1');
model.component('comp2').material('mat6').propertyGroup('BHCurve').func('BH').set('table', BH_curve_table);
model.component('comp2').material('mat6').propertyGroup('BHCurve').func('BH').set('extrap', 'linear');
model.component('comp2').material('mat6').propertyGroup('BHCurve').func('BH').set('fununit', {'T'});
model.component('comp2').material('mat6').propertyGroup('BHCurve').func('BH').set('argunit', {'A/m'});
model.component('comp2').material('mat6').propertyGroup('BHCurve').func('BH').set('defineinv', true);
model.component('comp2').material('mat6').propertyGroup('BHCurve').func('BH').set('defineprimfun', true);
model.component('comp2').material('mat6').propertyGroup('BHCurve').set('normB', 'BH(normHin)');
model.component('comp2').material('mat6').propertyGroup('BHCurve').set('normH', 'BH_inv(normBin)');
model.component('comp2').material('mat6').propertyGroup('BHCurve').set('Wpm', 'BH_prim(normHin)');
model.component('comp2').material('mat6').propertyGroup('BHCurve').descr('normHin', 'Magnetic field norm');
model.component('comp2').material('mat6').propertyGroup('BHCurve').descr('normBin', 'Magnetic flux density norm');
model.component('comp2').material('mat6').propertyGroup('BHCurve').addInput('magneticfield');
model.component('comp2').material('mat6').propertyGroup('BHCurve').addInput('magneticfluxdensity');
model.component('comp2').material('mat6').propertyGroup().create('Enu', 'Youngs_modulus_and_Poissons_ratio');
model.component('comp2').material('mat6').propertyGroup('Enu').set('E', '200e9');
model.component('comp2').material('mat6').propertyGroup('Enu').set('nu', '0.3');
model.component('comp2').material('mat6').propertyGroup('def').set('density', '7650');

% Assegnazione materiali
selNumber_fer = [];
fer = [];

for kk = 1:size(geo.BLKLABELS.rotore.xy, 1)
    x = geo.BLKLABELS.rotore.xy(kk, 1);
    y = geo.BLKLABELS.rotore.xy(kk, 2);
    disk3.set('posx', x);
    disk3.set('posx', y);
    selNumber_fer = disk3.entities();
    if geo.BLKLABELS.rotore.xy(kk, 3) == 5
       fer = horzcat(fer, selNumber_fer);
    end
end

model.component('comp2').material('mat6').selection().set(fer);

% Assegnazione boundary conditions
% setting rotating frame
model.component('comp2').physics('solid').create('rotf1', 'RotatingFrame', 2);
model.component('comp2').physics('solid').feature('rotf1').set('SpinSoftening', false);
model.component('comp2').physics('solid').feature('rotf1').set('Ovm', dataSet.EvalSpeed/30*pi);

% setting sector symmetry
model.component('comp2').physics('solid').create('sym1', 'SymmetrySolid', 1);

simm_bound = [];

for kk = 1:size(geo.BLKLABELS.rotore.boundary, 1)
    x = geo.BLKLABELS.rotore.boundary(kk, 1);
    y = geo.BLKLABELS.rotore.boundary(kk, 2);
    disk4.set('posx', x);
    disk4.set('posy', y);
    selNumber = disk4.entities();
    simm_bound = horzcat(simm_bound, selNumber);
end

model.component('comp2').physics('solid').feature('sym1').selection().set(simm_bound);

selNumber_shaft = [];

raggio = geo.Ar;
alpha_mecc = pi/geo.p/2;
x_mecc_m = raggio*cos(alpha_mecc);
y_mecc_m = raggio*sin(alpha_mecc);
disk4.set('posx', x_mecc_m);
disk4.set('posy', y_mecc_m);
selNumber_shaft = disk4.entities();

model.component('comp2').physics('solid').create('fix1', 'Fixed', 1);
model.component('comp2').physics('solid').feature('fix1').selection().set(selNumber_shaft);

% mphrun(model)

model.study('std2').run();      % run singolo studio (senza barra di progresso)

%% ====== Post-processing ====== %% 

% nel caso in cui si apra un file già risolto, bisogna cambiare i numeri
% dei plotgroup (pg), pg1 diventa pg5 e così via e verificare che i dataset
% siano quelli corretti

width    = 12;
height   = 10;
textSize = 12;

% stress von Mises deform
model.result().create('pg1', 'PlotGroup2D');
model.result('pg1').set('data', 'dset2');
model.result('pg1').set('defaultPlotID', 'stress');
model.result('pg1').label('Stress (solid)');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').set('edges', false);
model.result('pg1').create('surf1', 'Surface');
model.result('pg1').feature('surf1').set('expr', {'solid.misesGp'});
model.result('pg1').feature('surf1').set('threshold', 'manual');
model.result('pg1').feature('surf1').set('thresholdvalue', 0.2);
model.result('pg1').feature('surf1').set('colortable', 'Rainbow');
model.result('pg1').feature('surf1').set('colortabletrans', 'none');
model.result('pg1').feature('surf1').set('colorscalemode', 'linear');
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg1').feature('surf1').set('colortable', 'Prism');
model.result('pg1').feature('surf1').create('def1', 'Deform');
model.result('pg1').feature('surf1').feature('def1').set('scaleactive', true);
model.result('pg1').feature('surf1').feature('def1').set('scale', 100);
model.result('pg1').feature('surf1').feature('def1').set('expr', {'u', 'v'});
model.result('pg1').feature('surf1').feature('def1').set('descr', 'Displacement field');
model.result('pg1').feature('surf1').set('unit', 'MPa');
model.result('pg1').feature('surf1').set('colortable', 'RainbowLightClassic');

% Points exceeding the yield limit
model.result().export().create('plot1', 'pg1', 'surf1', 'Plot');
model.result().export('plot1').set('filename', Stress_csv);
model.result().export('plot1').run();

Stress = readmatrix(fullfile(Stress_csv));

[out] = eval_maxStress_COMSOL(Stress,mat);

figure()
% set(gcf,'color','w')
mphplot(model, 'pg1', 'rangenum', 1);
% xlabel('[mm]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
% ylabel('[mm]', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');
% set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
% set(gca, 'GridLineStyle', ':');
set(gcf,'defaultTextInterpreter','Latex');
set(gcf,'defaultLegendInterpreter','Latex');
set(gcf,'defaultAxesTickLabelInterpreter','Latex');

set(gcf,'defaultAxesLineWidth',1);
set(gcf,'defaultLineLineWidth',1.5);

set(gcf,'defaultAxesGridLineStyle',':');
set(gcf,'defaultAxesXGrid','on');
set(gcf,'defaultAxesYGrid','on');
set(gcf,'defaultAxesZGrid','on');
set(gcf,'defaultAxesXColor',0*[1 1 1]);
set(gcf,'defaultAxesYColor',0*[1 1 1]);
set(gcf,'defaultAxesZColor',0*[1 1 1]);
%set(gcf,'defaultAxesGridAlpha',1);
%set(gcf,'defaultAxesLayer','top');
set(gcf,'defaultAxesBox','on');
set(gcf,'defaultAxesNextPlot','add');

% Define text size
textSize = 12;  % Or set it to your desired value

set(gcf,'defaultAxesFontSize',textSize);
set(gcf,'defaultTextFontSize',textSize);
set(gcf,'defaultLegendFontSize',textSize);
set(gcf,'defaultAxesFontSizeMode','manual');

set(gcf,'defaultTextFontSizeMode','manual');
set(gcf,'defaultLegendFontSizeMode','manual');
set(gcf,'defaultAxesLabelFontSizeMultiplier',1);
set(gcf,'defaultAxesTitleFontSizeMultiplier',1);

set(gcf,'defaultAxesFontName','Times');
set(gcf,'defaultTextFontName','Times');

screenPos=get(groot,'ScreenSize')/get(groot,'ScreenPixelsPerInch')*2.54; % cm
figPos(1)=screenPos(3)/2-width/2;
figPos(2)=screenPos(4)/2-height/2;
figPos(3)=width;
figPos(4)=height;

colors{1} = [0.0 0.0 1.0];
colors{2} = [1.0 0.0 0.0];
colors{3} = [0.0 0.8 0.0];
colors{4} = [1.0 0.5 0.0];
colors{5} = [0.0 0.8 0.8];
colors{6} = [0.8 0.0 0.8];
colors{7} = [1.0 0.8 0.0];

set(gcf,'defaultAxesColorOrder',[colors{1};colors{2};colors{3};colors{4};colors{5};colors{6};colors{7}]);

% set(gcf,'Renderer','painters');
set(gcf,'Units','centimeters');
set(gcf,'Position',figPos);
set(gcf,'Color',[1 1 1]);

% Set the colormap to 'turbo'
colormap('turbo');

title('von Mises Stress [MPa] - deformation scale = 100', 'FontWeight', 'normal', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');

% von Mises Stress
model.result().create("pg3", "PlotGroup2D");
model.result("pg3").set("data", "dset2");
model.result("pg3").create("surf1", "Surface");
model.result('pg3').set('edges', false);
model.result("pg3").feature("surf1").set("expr", "solid.misesGp");
model.result("pg3").feature("surf1").set("descr", "von Mises stress");
model.result("pg3").feature("surf1").set("unit", "MPa");
model.result("pg3").feature("surf1").set("colortable", "RainbowLightClassic");

figure()
% set(gcf,'color','w')
mphplot(model, 'pg3', 'rangenum', 1);
% xlabel('[mm]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
% ylabel('[mm]', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');
% set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
% set(gca, 'GridLineStyle', ':');
set(gcf,'defaultTextInterpreter','Latex');
set(gcf,'defaultLegendInterpreter','Latex');
set(gcf,'defaultAxesTickLabelInterpreter','Latex');

set(gcf,'defaultAxesLineWidth',1);
set(gcf,'defaultLineLineWidth',1.5);

set(gcf,'defaultAxesGridLineStyle',':');
set(gcf,'defaultAxesXGrid','on');
set(gcf,'defaultAxesYGrid','on');
set(gcf,'defaultAxesZGrid','on');
set(gcf,'defaultAxesXColor',0*[1 1 1]);
set(gcf,'defaultAxesYColor',0*[1 1 1]);
set(gcf,'defaultAxesZColor',0*[1 1 1]);
%set(gcf,'defaultAxesGridAlpha',1);
%set(gcf,'defaultAxesLayer','top');
set(gcf,'defaultAxesBox','on');
set(gcf,'defaultAxesNextPlot','add');

% Define text size
textSize = 12;  % Or set it to your desired value

set(gcf,'defaultAxesFontSize',textSize);
set(gcf,'defaultTextFontSize',textSize);
set(gcf,'defaultLegendFontSize',textSize);
set(gcf,'defaultAxesFontSizeMode','manual');

set(gcf,'defaultTextFontSizeMode','manual');
set(gcf,'defaultLegendFontSizeMode','manual');
set(gcf,'defaultAxesLabelFontSizeMultiplier',1);
set(gcf,'defaultAxesTitleFontSizeMultiplier',1);

set(gcf,'defaultAxesFontName','Times');
set(gcf,'defaultTextFontName','Times');

screenPos=get(groot,'ScreenSize')/get(groot,'ScreenPixelsPerInch')*2.54; % cm
figPos(1)=screenPos(3)/2-width/2;
figPos(2)=screenPos(4)/2-height/2;
figPos(3)=width;
figPos(4)=height;

colors{1} = [0.0 0.0 1.0];
colors{2} = [1.0 0.0 0.0];
colors{3} = [0.0 0.8 0.0];
colors{4} = [1.0 0.5 0.0];
colors{5} = [0.0 0.8 0.8];
colors{6} = [0.8 0.0 0.8];
colors{7} = [1.0 0.8 0.0];

set(gcf,'defaultAxesColorOrder',[colors{1};colors{2};colors{3};colors{4};colors{5};colors{6};colors{7}]);

% set(gcf,'Renderer','painters');
set(gcf,'Units','centimeters');
set(gcf,'Position',figPos);
set(gcf,'Color',[1 1 1]);

% Set the colormap to 'turbo'
colormap('turbo');

title('von Mises Stress [MPa]', 'FontWeight', 'normal', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');

% displacement
model.result().create('pg2', 'PlotGroup2D');
model.result('pg2').set('data', 'dset2');
model.result('pg2').label('Displacement');
model.result('pg2').set('frametype', 'spatial');
model.result('pg2').set('edges', false);
model.result('pg2').create('surf1', 'Surface');
model.result('pg2').feature('surf1').create('def1', 'Deform');
model.result('pg2').feature('surf1').feature('def1').set('scaleactive', true);
model.result('pg2').feature('surf1').feature('def1').set('scale', 1);
model.result('pg2').feature('surf1').set('colortable', 'RainbowLightClassic');

figure()
% set(gcf,'color','w')
mphplot(model, 'pg2', 'rangenum', 1);
% xlabel('[mm]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
% ylabel('[mm]', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');
% set(gca, 'TickLabelInterpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
% set(gca, 'GridLineStyle', ':');
set(gcf,'defaultTextInterpreter','Latex');
set(gcf,'defaultLegendInterpreter','Latex');
set(gcf,'defaultAxesTickLabelInterpreter','Latex');

set(gcf,'defaultAxesLineWidth',1);
set(gcf,'defaultLineLineWidth',1.5);

set(gcf,'defaultAxesGridLineStyle',':');
set(gcf,'defaultAxesXGrid','on');
set(gcf,'defaultAxesYGrid','on');
set(gcf,'defaultAxesZGrid','on');
set(gcf,'defaultAxesXColor',0*[1 1 1]);
set(gcf,'defaultAxesYColor',0*[1 1 1]);
set(gcf,'defaultAxesZColor',0*[1 1 1]);
%set(gcf,'defaultAxesGridAlpha',1);
%set(gcf,'defaultAxesLayer','top');
set(gcf,'defaultAxesBox','on');
set(gcf,'defaultAxesNextPlot','add');

% Define text size
textSize = 12;  % Or set it to your desired value

set(gcf,'defaultAxesFontSize',textSize);
set(gcf,'defaultTextFontSize',textSize);
set(gcf,'defaultLegendFontSize',textSize);
set(gcf,'defaultAxesFontSizeMode','manual');

set(gcf,'defaultTextFontSizeMode','manual');
set(gcf,'defaultLegendFontSizeMode','manual');
set(gcf,'defaultAxesLabelFontSizeMultiplier',1);
set(gcf,'defaultAxesTitleFontSizeMultiplier',1);

set(gcf,'defaultAxesFontName','Times');
set(gcf,'defaultTextFontName','Times');

screenPos=get(groot,'ScreenSize')/get(groot,'ScreenPixelsPerInch')*2.54; % cm
figPos(1)=screenPos(3)/2-width/2;
figPos(2)=screenPos(4)/2-height/2;
figPos(3)=width;
figPos(4)=height;

colors{1} = [0.0 0.0 1.0];
colors{2} = [1.0 0.0 0.0];
colors{3} = [0.0 0.8 0.0];
colors{4} = [1.0 0.5 0.0];
colors{5} = [0.0 0.8 0.8];
colors{6} = [0.8 0.0 0.8];
colors{7} = [1.0 0.8 0.0];

set(gcf,'defaultAxesColorOrder',[colors{1};colors{2};colors{3};colors{4};colors{5};colors{6};colors{7}]);

% set(gcf,'Renderer','painters');
set(gcf,'Units','centimeters');
set(gcf,'Position',figPos);
set(gcf,'Color',[1 1 1]);

% Set the colormap to 'turbo'
colormap('turbo');

title('Displacement [mm]', 'FontWeight', 'normal', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');

% %Point exceeding yield limit
% if any(out.stressrot(:) ~= 0)
% figure ()
% figSetting()
% xlabel('[mm]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
% ylabel('[mm]', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');
% title('Von Mises Stress Limit Exceeded', 'FontWeight', 'normal', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');
% GUI_Plot_Machine(gca, geo.rotor(1:end-3, :))
% plot(out.stressrot(:,1), out.stressrot(:,2), 'ro', 'MarkerSize', 1, 'LineWidth', 2);
% end

mphsave(model, [pathname filename(1:end-4), '_str_solved.mph'])

