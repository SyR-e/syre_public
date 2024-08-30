% Copyright 2024
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

function [geo,mat] = draw_motor_in_COMSOL(geo,mat,pathIn,nameIn)

% Connessione a COMSOL e Apertura modello default
import com.comsol.model.*
import com.comsol.model.util.*
model = mphopen('Test_auto.mph');

if ~isfile([pathIn nameIn(1:end-4),'_dxf'])
    button='Yes';
else
    button = questdlg('DXF existing. Replace it?','SELECT','Yes','No','Yes');
end

% if strcmp(button,'Yes')
%     syreToDxf(geo.stator, geo.rotor,pathname, filename);
% end

% Estrazione del nome del file senza estensione
[~, file_name, ~] = fileparts(nameIn);

%Caricamento dei file
load([pathIn nameIn(1:end-4),'.mat']);

% Caricamento dei file DXF
dxf_dir = 'C:\Users\S296193\Desktop\syre_developers_20240610\motorExamples\syreDefaultMotor_dxf\';
rot_dxf = fullfile(dxf_dir, [file_name, '_rot.dxf']);
stat_dxf = fullfile(dxf_dir, [file_name, '_stat.dxf']);

% ============== Import e Costruzione della Geometria ============== %

geom = model.component('comp1').geom('geom1');
model.component('comp1').geom('geom1').create('imp1', 'Import');
model.component('comp1').geom('geom1').feature('imp1').set('filename', rot_dxf); 
model.component('comp1').geom('geom1').feature('imp1').importData();
model.component('comp1').geom('geom1').create('imp2', 'Import');
model.component('comp1').geom('geom1').feature('imp2').set('filename', stat_dxf);
model.component('comp1').geom('geom1').feature('imp2').importData();
model.component('comp1').geom('geom1').feature('fin').set('action', 'assembly');
model.component('comp1').geom('geom1').run('fin');
model.component('comp1').geom('geom1').lengthUnit('mm');

% Implementazione lunghezza assiale
model.component('comp1').physics('rmm').prop('d').set('d', 'l [mm]');

% Implementazione Air Gap
%rotore
xre2_n = geo.r + geo.g/2;
yre2_n = 0;
xre3_n = (geo.r + geo.g/2)*cos(pi/geo.p*geo.ps);
yre3_n = (geo.r + geo.g/2)*sin(pi/geo.p*geo.ps);

xra2_n = geo.r;
yra2_n = 0;
xra3_n = geo.r*cos(pi/geo.p*geo.ps);
yra3_n = geo.r*sin(pi/geo.p*geo.ps);

index_el_r = max(geo.rotor(:,9)) + 1;
air_r = 1;

%statore
xse1_n = geo.r + geo.g;
yse1_n = 0;
xse2_n = (geo.r + geo.g)*cos(pi/geo.p*geo.ps);
yse2_n = (geo.r + geo.g)*sin(pi/geo.p*geo.ps);

xsi1_n = geo.r + geo.g/2;
ysi1_n = 0;
xsi2_n = (geo.r + geo.g/2)*cos(pi/geo.p*geo.ps);
ysi2_n = (geo.r + geo.g/2)*sin(pi/geo.p*geo.ps);

index_el_s = max(geo.stator(:,9)) + 1;
air_s = 2;

%assegnazione coordinate
geo_stator = [];
geo_rotor = [];
geo.stator_n = [xsi1_n, ysi1_n, xse1_n, yse1_n, NaN, NaN, 0, air_s, index_el_s; 
                0, 0, xse1_n, yse1_n, xse2_n, yse2_n, 1, air_s, index_el_s;
                xse2_n, yse2_n, xsi2_n, ysi2_n, NaN, NaN, 0, air_s, index_el_s
                ];
geo.rotor_n = [xra2_n, yra2_n, xre2_n, yre2_n, NaN, NaN, 0, air_r, index_el_r;
               0, 0, xre2_n, yre2_n, xre3_n, yre3_n, 1, air_r, index_el_r;
               xre3_n, yre3_n, xra3_n, yra3_n, NaN, NaN, 0, air_r, index_el_r;];

geo_stator = [geo.stator; geo.stator_n];
geo_rotor = [geo.rotor; geo.rotor_n];

% Ndomains = max(geo.stator(:,end))+(geo.Qs)+(geo.Qs-1)+max(geo.rotor(:,end));
Ndomains_n = max(geo_stator(:,end))+(geo.Qs)+(geo.Qs-1)+max(geo_rotor(:,end));

% Correzione matrici boundaries
stat_boundary = [];
rot_boundary = [];

%rotore
x_bnd_r_1 = (xra2_n + xre2_n)/2;
y_bnd_r_1 = (yra2_n + yre2_n)/2;
x_bnd_r_2 = (xre3_n + xra3_n)/2;
y_bnd_r_2 = (yre3_n + yra3_n)/2;

%statore
x_bnd_s_1 = (xsi1_n + xse1_n)/2;
y_bnd_s_1 = (ysi1_n + yse1_n)/2;
x_bnd_s_2 = (xse2_n + xsi2_n)/2;
y_bnd_s_2 = (yse2_n + ysi2_n)/2;

%assegnazione coordinate
gm.s.boundary = [x_bnd_s_1, y_bnd_s_1, 10;
                 x_bnd_s_2, y_bnd_s_2, 10
                 ];
gm.r.boundary = [x_bnd_r_1, y_bnd_r_1, 10;
                 x_bnd_r_2, y_bnd_r_2, 10
                 ];

stat_boundary = [geo.BLKLABELS.statore.boundary; gm.s.boundary];
rot_boundary = [geo.BLKLABELS.rotore.boundary; gm.r.boundary];

% Correzione matrici materiali (BLKLABELS)
stat_mat = [];
rot_mat = [];

%definizione punto medio
xm_arc = geo.r*cos(pi/geo.p); 
ym_arc = geo.r*sin(pi/geo.p);

%rotore
xm_rot = xm_arc + geo.g/4;
ym_rot = ym_arc + geo.g/4;

%statore
xm_stat = xm_arc + (3*geo.g/4);
ym_stat = ym_arc + (3*geo.g/4);

stat_mat = [geo.BLKLABELS.statore.xy; xm_stat, ym_stat, 2, 1.6667, 1];
rot_mat = [geo.BLKLABELS.rotore.xy; xm_rot, ym_rot, 1, 1.6667, 1, 0, 0, 0];

% ============== Assegnazione dei Materiali ============== %

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

% Air
model.component('comp1').material().create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup('def').func().create('eta', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func().create('Cp', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func().create('rho', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func().create('k', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func().create('cs', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func().create('an1', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func().create('an2', 'Analytic');
model.component('comp1').material('mat1').propertyGroup().create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat1').propertyGroup().create('NonlinearModel', 'Nonlinear model');
model.component('comp1').material('mat1').propertyGroup().create('idealGas', 'Ideal gas');
model.component('comp1').material('mat1').propertyGroup('idealGas').func().create('Cp', 'Piecewise');
model.component('comp1').material('mat1').label('Air');
model.component('comp1').material('mat1').set('family', 'air');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('pieces', {'200.0', '1600.0', '-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('fununit', 'Pa*s');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('pieces', {'200.0', '1600.0', '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('fununit', 'J/(kg*K)');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('expr', 'pA*0.02897/R_const[K*mol/J]/T');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('args', {'pA', 'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('fununit', 'kg/m^3');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('argunit', {'Pa', 'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('plotargs', {'pA', '101325', '101325'; 'T', '273.15', '293.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('pieces', {'200.0', '1600.0', '-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('fununit', 'W/(m*K)');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('expr', 'sqrt(1.4*R_const[K*mol/J]/0.02897*T)');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('fununit', 'm/s');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('argunit', {'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('plotargs', {'T', '273.15', '373.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('funcname', 'alpha_p');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('expr', '-1/rho(pA,T)*d(rho(pA,T),T)');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('args', {'pA', 'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('fununit', '1/K');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('argunit', {'Pa', 'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('plotargs', {'pA', '101325', '101325'; 'T', '273.15', '373.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('funcname', 'muB');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('expr', '0.6*eta(T)');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('fununit', 'Pa*s');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('argunit', {'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('plotargs', {'T', '200', '1600'});
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', '');
model.component('comp1').material('mat1').propertyGroup('def').set('molarmass', '');
model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', '');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', {'alpha_p(pA,T)', '0', '0', '0', 'alpha_p(pA,T)', '0', '0', '0', 'alpha_p(pA,T)'});
model.component('comp1').material('mat1').propertyGroup('def').set('molarmass', '0.02897[kg/mol]');
model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', 'muB(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1', '0', '0', '0', '1', '0', '0', '0', '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'1', '0', '0', '0', '1', '0', '0', '0', '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('ratioofspecificheat', '1.4');
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'0[S/m]', '0', '0', '0', '0[S/m]', '0', '0', '0', '0[S/m]'});
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', 'Cp(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rho(pA,T)');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'k(T)', '0', '0', '0', 'k(T)', '0', '0', '0', 'k(T)'});
model.component('comp1').material('mat1').propertyGroup('def').set('soundspeed', 'cs(T)');
model.component('comp1').material('mat1').propertyGroup('def').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('def').addInput('pressure');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'1', '0', '0', '0', '1', '0', '0', '0', '1'});
model.component('comp1').material('mat1').propertyGroup('NonlinearModel').set('BA', '(def.gamma+1)/2');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').label('Piecewise 2');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('pieces', {'200.0', '1600.0', '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('fununit', 'J/(kg*K)');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('Rs', 'R_const/Mn');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('heatcapacity', 'Cp(T)');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('ratioofspecificheat', '1.4');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('molarmass', '0.02897');
model.component('comp1').material('mat1').propertyGroup('idealGas').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('idealGas').addInput('pressure');
model.component('comp1').material('mat1').materialType('nonSolid');
model.component('comp1').material('mat1').set('family', 'air');

% NGO 35PN270
model.component('comp1').material().create('mat2', 'Common');
model.component('comp1').material('mat2').propertyGroup().create('BHCurve', 'B-H Curve');
model.component('comp1').material('mat2').propertyGroup('BHCurve').func().create('BH', 'Interpolation');
model.component('comp1').material('mat2').label('Silicon Steel NGO 35PN270');
model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', {'0'});
model.component('comp1').material('mat2').propertyGroup('def').set('relpermittivity', {'1[1]', '0', '0', '0', '1[1]', '0', '0', '0', '1[1]'});
model.component('comp1').material('mat2').propertyGroup('BHCurve').label('B-H Curve');
model.component('comp1').material('mat2').propertyGroup('BHCurve').func('BH').label('Interpolation 1');
model.component('comp1').material('mat2').propertyGroup('BHCurve').func('BH').set('table', BH_curve_table);
model.component('comp1').material('mat2').propertyGroup('BHCurve').func('BH').set('extrap', 'linear');
model.component('comp1').material('mat2').propertyGroup('BHCurve').func('BH').set('fununit', 'T');
model.component('comp1').material('mat2').propertyGroup('BHCurve').func('BH').set('argunit', 'A/m');
model.component('comp1').material('mat2').propertyGroup('BHCurve').func('BH').set('defineinv', true);
model.component('comp1').material('mat2').propertyGroup('BHCurve').func('BH').set('defineprimfun', true);
model.component('comp1').material('mat2').propertyGroup('BHCurve').set('normB', 'BH(normHin)');
model.component('comp1').material('mat2').propertyGroup('BHCurve').set('normH', 'BH_inv(normBin)');
model.component('comp1').material('mat2').propertyGroup('BHCurve').set('Wpm', 'BH_prim(normHin)');
model.component('comp1').material('mat2').propertyGroup('BHCurve').descr('normHin', 'Magnetic field norm');
model.component('comp1').material('mat2').propertyGroup('BHCurve').descr('normBin', 'Magnetic flux density norm');
model.component('comp1').material('mat2').propertyGroup('BHCurve').addInput('magneticfield');
model.component('comp1').material('mat2').propertyGroup('BHCurve').addInput('magneticfluxdensity');
model.component('comp1').material('mat2').set('family', 'plastic');

% Copper
model.component('comp1').material().create('mat3', 'Common');
model.component('comp1').material('mat3').propertyGroup().create('Enu', 'Young''s modulus and Poisson''s ratio');
model.component('comp1').material('mat3').propertyGroup().create('linzRes', 'Linearized resistivity');
model.component('comp1').material('mat3').label('Copper');
model.component('comp1').material('mat3').set('family', 'copper');
model.component('comp1').material('mat3').propertyGroup('def').set('relpermeability', {'1', '0', '0', '0', '1', '0', '0', '0', '1'});
model.component('comp1').material('mat3').propertyGroup('def').set('electricconductivity', {'5.998e7[S/m]', '0', '0', '0', '5.998e7[S/m]', '0', '0', '0', '5.998e7[S/m]'});
model.component('comp1').material('mat3').propertyGroup('def').set('heatcapacity', '385[J/(kg*K)]');
model.component('comp1').material('mat3').propertyGroup('def').set('relpermittivity', {'1', '0', '0', '0', '1', '0', '0', '0', '1'});
model.component('comp1').material('mat3').propertyGroup('def').set('emissivity', '0.5');
model.component('comp1').material('mat3').propertyGroup('def').set('density', '8940[kg/m^3]');
model.component('comp1').material('mat3').propertyGroup('def').set('thermalconductivity', {'400[W/(m*K)]', '0', '0', '0', '400[W/(m*K)]', '0', '0', '0', '400[W/(m*K)]'});
model.component('comp1').material('mat3').propertyGroup('Enu').set('E', '126e9[Pa]');
model.component('comp1').material('mat3').propertyGroup('Enu').set('nu', '0.34');
model.component('comp1').material('mat3').propertyGroup('linzRes').set('rho0', '1.667e-8[ohm*m]');
model.component('comp1').material('mat3').propertyGroup('linzRes').set('alpha', '3.862e-3[1/K]');
model.component('comp1').material('mat3').propertyGroup('linzRes').set('Tref', '293.15[K]');
model.component('comp1').material('mat3').propertyGroup('linzRes').addInput('temperature');
model.component('comp1').material('mat3').set('family', 'copper');

% N52 PM
model.component('comp1').material().create('mat4', 'Common');
model.component('comp1').material('mat4').propertyGroup().create('RemanentFluxDensity', 'Remanent flux density');
model.component('comp1').material('mat4').label('N52 (Sintered NdFeB)');
model.component('comp1').material('mat4').set('family', 'chrome');
model.component('comp1').material('mat4').propertyGroup('def').set('electricconductivity', {'1/1.4[uohm*m]', '0', '0', '0', '1/1.4[uohm*m]', '0', '0', '0', '1/1.4[uohm*m]'});
model.component('comp1').material('mat4').propertyGroup('def').set('relpermittivity', {'1', '0', '0', '0', '1', '0', '0', '0', '1'});
model.component('comp1').material('mat4').propertyGroup('RemanentFluxDensity').set('murec', {'1.05', '0', '0', '0', '1.05', '0', '0', '0', '1.05'});
model.component('comp1').material('mat4').propertyGroup('RemanentFluxDensity').set('normBr', '1.44[T]');
model.component('comp1').material('mat4').set('family', 'chrome');

% Definizione selection inspector
model.component('comp1').selection().create('disk1', 'Disk');
model.component('comp1').selection().create('disk2', 'Disk');

% Assegnazione default Aria 
sel = 1:Ndomains_n;
model.component('comp1').material('mat1').selection().set(sel);
model.component('comp1').physics('rmm').feature('al1').label("Ampere's Law - Air");
model.component('comp1').physics('rmm').feature('al1').create('loss1', 'LossCalculation', 2);

% Vettori per assegnazione materiali
tmp = [];
tmp_rot = rot_mat(:, 1:3);
tmp_stat = stat_mat(:, 1:3);
tmp = [tmp_rot; tmp_stat]; 
disk1 = model.component('comp1').selection('disk1');

% Vettori per settare i domini
A = []; 
B = [];
C = [];

% Assegnazione materiali
for kk = 1:Ndomains_n
    x = tmp(kk,1);
    y = tmp(kk,2);
    disk1.set('posx', x);
    disk1.set('posy', y);
    selNumber = disk1.entities();
    switch tmp(kk, 3)
        case 3
            A = horzcat(A, selNumber);
        case {4, 5}
            B = horzcat(B, selNumber);
        case 6
            C = horzcat(C, selNumber);
    end
end

model.component('comp1').material('mat3').selection().set(A);
model.component('comp1').material('mat2').selection().set(B);
model.component('comp1').material('mat4').selection().set(C);

% ============== Fitting Steinmetz ============== %

% Calcolo tempo di simulazione

p = geo.p;                                  % paia poli macchina
w = per.EvalSpeed*pi/30;                    % velocità di rotazione [rad/s]
freq = w*p/2/pi;                            % frequenza di alimentazione [Hz]

ff = linspace(50, freq, 51);
Bf = linspace(0, Bf_max, 51);
kh = mat.Rotor.kh;
ke = mat.Rotor.ke;
alpha = mat.Rotor.alpha;
beta = mat.Rotor.beta;
pfe_f = kh .*(ff .^alpha) .*(Bf .^beta) +  ke .*(ff .^2) .*(Bf .^2);

[fitresult, gof] = createFit(ff, Bf, pfe_f);
KH = fitresult.KH;
ALPHA = fitresult.ALPHA;
BETA = fitresult.BETA;

% ============== Assegnazione Condizioni a Contorno ============== %

% Definizione boundary condition PM
tmp_rot_righe = size(rot_mat, 1);
AM = [];
BC_pm = [];

for kk = 1:tmp_rot_righe
    x = tmp(kk,1);
    y = tmp(kk,2);
    disk1.set('posx', x);
    disk1.set('posy', y);
    selNumber = disk1.entities();
    if rot_mat(kk, 3)==6
        BC_pm = horzcat(BC_pm, selNumber);
        xm = rot_mat(kk, 6);
        ym = rot_mat(kk, 7);
        zm = rot_mat(kk, 8);
        AM = [xm, ym, zm];
        alnumber_m = ['al_m' num2str(kk+1)];
        model.component('comp1').physics('rmm').create(alnumber_m, 'AmperesLaw', 2);
        model.component('comp1').physics('rmm').feature(alnumber_m).selection().set(selNumber);
        model.component('comp1').physics('rmm').feature(alnumber_m).set('ConstitutiveRelationBH', 'RemanentFluxDensity');
        model.component('comp1').physics('rmm').feature(alnumber_m).set('e_crel_BH_RemanentFluxDensity', AM);
        model.component('comp1').physics('rmm').feature(alnumber_m).create('loss1', 'LossCalculation', 2);
    end
end

% Memorizzazione domini barriere di flusso rotore
Bar = [];

for kk = 1:tmp_rot_righe
    x = tmp(kk,1);
    y = tmp(kk,2);
    disk1.set('posx', x);
    disk1.set('posy', y);
    selNumber = disk1.entities();
    if rot_mat(kk, 3)==1
        Bar = horzcat(Bar, selNumber);
    end
end

% Definizione boundary condition su ferro di statore
BC_fe_s = [];

for kk = 1:size(tmp_stat, 1)
    x = tmp_stat(kk,1);
    y = tmp_stat(kk,2);
    disk1.set('posx', x);
    disk1.set('posy', y);
    selNumber = disk1.entities();
    if tmp_stat(kk, 3)==4
       BC_fe_s = horzcat(BC_fe_s, selNumber);
    end
end

model.component('comp1').physics('rmm').create('al_fs', 'AmperesLaw', 2);
model.component('comp1').physics('rmm').feature('al_fs').selection().set(BC_fe_s);
model.component('comp1').physics('rmm').feature('al_fs').set('ConstitutiveRelationBH', 'BHCurve');
model.component('comp1').physics('rmm').feature('al_fs').create('loss_f', 'LossCalculation', 2);
model.component('comp1').physics('rmm').feature('al_fs').feature('loss_f').set('LossModel', 'Steinmetz');
model.component('comp1').physics('rmm').feature('al_fs').label("Ampere's Law - Stator");
model.component('comp1').physics('rmm').feature('al_fs').feature('loss_f').set('kh_steinmetz', KH);
model.component('comp1').physics('rmm').feature('al_fs').feature('loss_f').set('alpha', ALPHA);
model.component('comp1').physics('rmm').feature('al_fs').feature('loss_f').set('beta_steinmetz', BETA);

% Definizione boundary condition su ferro di rotore
BC_fe_r = [];

for kk = 1:size(tmp_rot, 1)
    x = tmp_rot(kk,1);
    y = tmp_rot(kk,2);
    disk1.set('posx', x);
    disk1.set('posy', y);
    selNumber = disk1.entities();
    if tmp_rot(kk,3)==5
       BC_fe_r = horzcat(BC_fe_r, selNumber);
    end
end

model.component('comp1').physics('rmm').create('al_fr', 'AmperesLaw', 2);
model.component('comp1').physics('rmm').feature('al_fr').selection().set(BC_fe_r);
model.component('comp1').physics('rmm').feature('al_fr').set('ConstitutiveRelationBH', 'BHCurve');
model.component('comp1').physics('rmm').feature('al_fr').create('loss_f', 'LossCalculation', 2);
model.component('comp1').physics('rmm').feature('al_fr').feature('loss_f').set('LossModel', 'Steinmetz');
model.component('comp1').physics('rmm').feature('al_fr').label("Ampere's Law - Rotor");
model.component('comp1').physics('rmm').feature('al_fs').feature('loss_f').set('kh_steinmetz', KH);
model.component('comp1').physics('rmm').feature('al_fs').feature('loss_f').set('alpha', ALPHA);
model.component('comp1').physics('rmm').feature('al_fs').feature('loss_f').set('beta_steinmetz', BETA);

% Definizione coils
avv = geo.win.avv;
[num_righe_avv, num_colonne_avv] = size(avv);
coil1 = [];
coil2 = [];
coil3 = [];

num_colonne_A = size(A, 2);    % Calcolo numero di colonne di A
n = numel(A)/3;                % Calcolo numero domini per fase

% Estrazione delle tre parti di A
parte1 = A(:, 1:n);
parte2 = A(:, n+1:2*n);
parte3 = A(:, 2*n+1:end);

% Reshape delle parti in matrici 2x3
matrice1 = reshape(parte1, 2, []);
matrice2 = reshape(parte2, 2, []);
matrice3 = reshape(parte3, 2, []);

% Creazione della matrice ordinata correttamente Ac
Ac = horzcat(matrice1, matrice2, matrice3); 

model.component('comp1').physics('rmm').create('coil1', 'Coil', 2);
model.component('comp1').physics('rmm').feature('coil1').label('Phase 1');
model.component('comp1').physics('rmm').feature('coil1').set('ConductorModel', 'Multi');
model.component('comp1').physics('rmm').feature('coil1').set('coilGroup', true);
model.component('comp1').physics('rmm').feature('coil1').set('CoilExcitation', 'CircuitCurrent');
model.component('comp1').physics('rmm').feature('coil1').set('N', {num2str(geo.win.Nbob*2*geo.p)});
model.component('comp1').physics('rmm').feature('coil1').set('AreaFrom', 'FillingFactor');
model.component('comp1').physics('rmm').feature('coil1').set('FillingFactor', {num2str(geo.win.kcu)});
model.component('comp1').physics('rmm').create('coil2', 'Coil', 2);
model.component('comp1').physics('rmm').feature('coil2').label('Phase 2');
model.component('comp1').physics('rmm').feature('coil2').set('ConductorModel', 'Multi');
model.component('comp1').physics('rmm').feature('coil2').set('coilGroup', true);
model.component('comp1').physics('rmm').feature('coil2').set('CoilExcitation', 'CircuitCurrent');
model.component('comp1').physics('rmm').feature('coil2').set('N', {num2str(geo.win.Nbob*2*geo.p)});
model.component('comp1').physics('rmm').feature('coil2').set('AreaFrom', 'FillingFactor');
model.component('comp1').physics('rmm').feature('coil2').set('FillingFactor', {num2str(geo.win.kcu)});
model.component('comp1').physics('rmm').create('coil3', 'Coil', 2);
model.component('comp1').physics('rmm').feature('coil3').label('Phase 3');
model.component('comp1').physics('rmm').feature('coil3').set('ConductorModel', 'Multi');
model.component('comp1').physics('rmm').feature('coil3').set('coilGroup', true);
model.component('comp1').physics('rmm').feature('coil3').set('CoilExcitation', 'CircuitCurrent');
model.component('comp1').physics('rmm').feature('coil3').set('N', {num2str(geo.win.Nbob*2*geo.p)});
model.component('comp1').physics('rmm').feature('coil3').set('AreaFrom', 'FillingFactor');
model.component('comp1').physics('rmm').feature('coil3').set('FillingFactor', {num2str(geo.win.kcu)});

% Ciclo per assegnare i domini alle coil in base ai valori di avv
for i = 1:num_righe_avv
    for c = 1:num_colonne_avv
        switch avv(i, c)
            case 1
                coil1 = horzcat(coil1, Ac(i, c));
            case -3
                coil3 = horzcat(coil3, Ac(i, c));
            case 2
                coil2 = horzcat(coil2, Ac(i, c));
        end
    end
end

% Assegna i domini alle bobine
model.component('comp1').physics('rmm').feature('coil1').selection().set(coil1);
model.component('comp1').physics('rmm').feature('coil2').selection().set(coil2);
model.component('comp1').physics('rmm').feature('coil3').selection().set(coil3);

% Assegna i domini alla bobina con corrente inversa 
model.component('comp1').physics('rmm').feature('coil3').create('rcd1', 'ReverseCoilGroupDomain', 2);
model.component('comp1').physics('rmm').feature('coil3').feature('rcd1').selection().set(coil3);
model.component('comp1').physics('rmm').feature('coil3').feature('rcd1').label('Reverse Current Phase 3');

%Assegna calcolo perdite alle bobine
model.component('comp1').physics('rmm').feature('coil1').create('loss1', 'LossCalculation', 2);
model.component('comp1').physics('rmm').feature('coil2').create('loss1', 'LossCalculation', 2);
model.component('comp1').physics('rmm').feature('coil3').create('loss1', 'LossCalculation', 2);

% Vettori assegnazione periodic condition sui bordi
tmp = [];
tmp_rot = rot_boundary(:, 1:3);
tmp_stat = stat_boundary(:, 1:3);
tmp_righe_s = size(tmp_stat, 1);
tmp_righe_r = size(tmp_rot, 1);
AP_s = [];
AP_r = [];  

disk2 = model.component('comp1').selection('disk2').set('entitydim', 1);
disk2 = model.component('comp1').selection('disk2').set('r', 0.1);

% Definizione periodic condition statore (continuità/antiperiodicità)
for kk = 1:tmp_righe_s
    x = tmp_stat(kk,1);
    y = tmp_stat(kk,2);
    disk2.set('posx', x);
    disk2.set('posy', y);
    selNumber = disk2.entities();       
    if tmp_stat(kk,3)==10
        AP_s = horzcat(AP_s, selNumber);
    end
end

model.component('comp1').physics('rmm').create('pc1', 'PeriodicCondition', 1);
model.component('comp1').physics('rmm').feature('pc1').selection().set(AP_s);
if mod(geo.ps, 2)==0
model.component('comp1').physics('rmm').feature('pc1').set('PeriodicType', 'Continuity');
else
model.component('comp1').physics('rmm').feature('pc1').set('PeriodicType', 'AntiPeriodicity');
end
model.component('comp1').physics('rmm').feature('pc1').label("Periodic Condition - Stator");

% Definizione periodic condition rotore (continuità/antiperiodicità)
for kk = 1:tmp_righe_r
    x = tmp_rot(kk,1);
    y = tmp_rot(kk,2);
    disk2.set('posx', x);
    disk2.set('posy', y);
    selNumber = disk2.entities();       
    if tmp_rot(kk,3)==10
        AP_r = horzcat(AP_r, selNumber);
    end
end

model.component('comp1').physics('rmm').create('pc2', 'PeriodicCondition', 1);
model.component('comp1').physics('rmm').feature('pc2').selection().set(AP_r);
if mod(geo.ps, 2)==0
model.component('comp1').physics('rmm').feature('pc2').set('PeriodicType', 'Continuity');
else
model.component('comp1').physics('rmm').feature('pc2').set('PeriodicType', 'AntiPeriodicity');
end
model.component('comp1').physics('rmm').feature('pc2').label("Periodic Condition - Rotor");

% Definizione di Sector Symmetry su Air Gap
model.component('comp1').physics('rmm').create('ssc1', 'SectorSymmetry', 1);
model.component('comp1').physics('rmm').feature('ssc1').set('pairs', 'ap1');
if mod(geo.ps, 2)==0
model.component('comp1').physics('rmm').feature('ssc1').set('nsector', geo.p);
model.component('comp1').physics('rmm').feature('ssc1').set('PeriodicType', 'Continuity');
else
model.component('comp1').physics('rmm').feature('ssc1').set('nsector', geo.p*2);
model.component('comp1').physics('rmm').feature('ssc1').set('PeriodicType', 'AntiPeriodicity');
end
model.component('comp1').physics('rmm').feature('ssc1').set('constraintOptions', 'weakConstraints');

% Definizione Arkkio Torque Calculation
model.component('comp1').physics('rmm').create('ark1', 'ArkkioTorqueCalculation', 2);

% ============== Definizione Moving Mesh ============== %
T = [];
tmp = rot_mat;

for kk = 1:tmp_rot_righe        
    x = tmp(kk,1);
    y = tmp(kk,2);
    disk1.set('posx', x);
    disk1.set('posy', y);   
    selNumber = disk1.entities();
    T = horzcat(T, selNumber);
end

model.component('comp1').common().create('rot1', 'RotatingDomain');
model.component('comp1').common('rot1').selection().set(T);
model.component('comp1').common('rot1').set('rotationType', 'rotationalVelocity');
model.component('comp1').common('rot1').set('rotationalVelocityExpression', 'constantAngularVelocity');
model.component('comp1').common('rot1').set('angularVelocity', w);

% ============== Definizione Circuito ============== %

% Calcolo correnti dq
iAmp = per.overload*per.i0;         % ampiezza corrente (tiene conto di eventuale sovraccarico)
gamma = dataSet.GammaPP;                         % angolo del vettore I rispetto all'asse d [deg]
theta_i = (geo.th0 + gamma)*pi/180;             % Angolo corrente [rad]

id = iAmp*cos(theta_i);       
iq = iAmp*sin(theta_i);       
Imod = abs(id + 1i*iq);            % Modulo corrente [A]
Iarg = angle(Imod(end));           % Angolo corrente [rad]

model.component('comp1').physics().create('cir', 'Circuit', 'geom1');
model.component('comp1').physics('cir').create('I1', 'CurrentSourceCircuit', -1);
model.component('comp1').physics('cir').create('I2', 'CurrentSourceCircuit', -1);
model.component('comp1').physics('cir').create('I3', 'CurrentSourceCircuit', -1);
model.component('comp1').physics('cir').create('termI1', 'ModelTerminalIV', -1);
model.component('comp1').physics('cir').create('termI2', 'ModelTerminalIV', -1);
model.component('comp1').physics('cir').create('termI3', 'ModelTerminalIV', -1);
model.component('comp1').physics('cir').feature('I2').setIndex('Connections', 1, 0, 0);
model.component('comp1').physics('cir').feature('I2').setIndex('Connections', 3, 1, 0);
model.component('comp1').physics('cir').feature('I3').setIndex('Connections', 1, 0, 0);
model.component('comp1').physics('cir').feature('I3').setIndex('Connections', 4, 1, 0);
model.component('comp1').physics('cir').feature('termI1').set('Connections', 2);
model.component('comp1').physics('cir').feature('termI2').set('Connections', 3);
model.component('comp1').physics('cir').feature('termI3').set('Connections', 4);
model.component('comp1').physics('cir').create('R1', 'Resistor', -1);
model.component('comp1').physics('cir').feature('R1').set('R', '10000 [Ω]');
model.component('comp1').physics('cir').feature('R1').setIndex('Connections', 1, 0, 0);
model.component('comp1').physics('cir').feature('R1').setIndex('Connections', 0, 1, 0);
model.component('comp1').physics('cir').feature('I1').set('sourceType', 'SineSource');
model.component('comp1').physics('cir').feature('I2').set('sourceType', 'SineSource');
model.component('comp1').physics('cir').feature('I3').set('sourceType', 'SineSource');
model.component('comp1').physics('cir').feature('I1').set('value', [num2str(Imod) ' [A]']);
model.component('comp1').physics('cir').feature('I1').set('freq', [num2str(freq) ' [Hz]']);
model.component('comp1').physics('cir').feature('I1').set('phase', [num2str(theta_i)]);
model.component('comp1').physics('cir').feature('I2').set('value', [num2str(Imod) ' [A]']);
model.component('comp1').physics('cir').feature('I2').set('freq', [num2str(freq) ' [Hz]']);
model.component('comp1').physics('cir').feature('I2').set('phase', [num2str(theta_i) ' - 120*pi/180']);
model.component('comp1').physics('cir').feature('I3').set('value', [num2str(Imod) ' [A]']);
model.component('comp1').physics('cir').feature('I3').set('freq', [num2str(freq) ' [Hz]']);
model.component('comp1').physics('cir').feature('I3').set('phase', [num2str(theta_i) ' + 120*pi/180']);

% ============== Costruzione Mesh ============== %

model.component('comp1').mesh('mesh1').autoMeshSize(6);    % mesh size (1-10) dalla più fitta alla meno fitta
model.component('comp1').mesh('mesh1').run();

% ============== Salvataggio Modello Inizzializzato ============== %

geo.BC_fe_s = BC_fe_s;
geo.BC_fe_r = BC_fe_r;
geo.BC_pm = BC_pm;
geo.Bar = Bar;

mphsave(model, [pathIn nameIn(1:end-4), '.mph']);
