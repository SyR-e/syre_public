clear

% Connessione a COMSOL e apertura modello default
import com.comsol.model.*
import com.comsol.model.util.*
model = mphopen('Test_auto.mph');

% Caricamento file .mat

file_path = 'C:\Users\S296193\Dropbox (Politecnico Di Torino Studenti)\Tesi Magistrale\syre_developers_20240410\motorExamples\syreDefaultMotor.mat';
load(file_path);
[~, file_name, ~] = fileparts(file_path);
rot_dxf = fullfile('C:\Users\S296193\Dropbox (Politecnico Di Torino Studenti)\Tesi Magistrale\syre_developers_20240410\COMSOL\Tests\dxf\', [file_name, '_rot.dxf']);
stat_dxf = fullfile('C:\Users\S296193\Dropbox (Politecnico Di Torino Studenti)\Tesi Magistrale\syre_developers_20240410\COMSOL\Tests\dxf\', [file_name, '_stat.dxf']);

% torque_csv = fullfile('C:\Users\S296193\Desktop\', [file_name, '_Torque.csv']);
% flux_csv = fullfile('C:\Users\S296193\Desktop\', [file_name, '_Flux.csv']);
% current_csv = fullfile('C:\Users\S296193\Desktop\', [file_name, '_Current.csv']);

% Definizione del percorso del file di input
pathname = 'C:\Users\S296193\Dropbox (Politecnico Di Torino Studenti)\Tesi Magistrale\syre_developers_20240410\motorExamples\';
filename = 'syreDefaultMotor.mat';

% Caricamento del file
load(fullfile(pathname, filename));

% Estrazione del nome del file senza estensione
[~, file_name, ~] = fileparts(filename);

% Caricamento dei file DXF
dxf_dir = 'C:\Users\S296193\Dropbox (Politecnico Di Torino Studenti)\Tesi Magistrale\syre_developers_20240410\COMSOL\Tests\dxf\';
rot_dxf = fullfile(dxf_dir, [file_name, '_rot.dxf']);
stat_dxf = fullfile(dxf_dir, [file_name, '_stat.dxf']);


% Definizione della directory di output per i file CSV
csv_dir = 'C:\Users\S296193\Desktop\';

% Creazione dei percorsi completi per i file CSV
torque_csv = fullfile(csv_dir, [file_name, '_Torque.csv']);
flux_csv = fullfile(csv_dir, [file_name, '_Flux.csv']);
current_csv = fullfile(csv_dir, [file_name, '_Current.csv']);
Pfe_r_csv = fullfile(csv_dir, [file_name, '_Pfe_r.csv']);
Pfe_s_csv = fullfile(csv_dir, [file_name, '_Pfe_s.csv']);
PM_losses_csv = fullfile(csv_dir, [file_name, '_PM_losses.csv']);
PjrBar_csv = fullfile(csv_dir, [file_name, '_PjrBar.csv']);
Stress_csv = fullfile(csv_dir, [file_name, '_Stress-csv']);

% Verifica geometria del rotore

% if strcmp(geo.RotType,'SPM') || strcmp(geo.RotType,'Vtype')
%     geo.axisType = 'PM';
%     geo.th0 = geo.th0 - 90;
% else
%     geo.axisType = 'SR';
%     geo.th0 = geo.th0 + 90;
% end
    
% Import e costruzione della geometria

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

% Create the selection inspector
model.component('comp1').selection().create('disk1', 'Disk');
model.component('comp1').selection().create('disk2', 'Disk');

% Material Definition

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
% plot(BH_curve(:,1), BH_curve(:,2), 'r', 'LineWidth', 2);

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

% Implementazione lunghezza assiale
model.component('comp1').physics('rmm').prop('d').set('d', 'l');

% Impostazione discretizzazione vettore potenziale magnetico (lineare = 1, quadratica = 2, cubica = 3)
model.component('comp1').physics('rmm').prop('ShapeProperty').set('order_magneticvectorpotential', 2);

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

% Assegnazione default Aria 
sel = 1:Ndomains_n;
% mphsave(model, 'rotto.mph');
model.component('comp1').material('mat1').selection().set(sel);
model.component('comp1').physics('rmm').feature('al1').label("Ampere's Law - Air");
model.component('comp1').physics('rmm').feature('al1').create('loss1', 'LossCalculation', 2);

% Vettori per assegnazione materiali
tmp = [];
tmp_rot = rot_mat(:, 1:3);
tmp_stat = stat_mat(:, 1:3);
tmp = [tmp_rot; tmp_stat]; 
disk1 = model.component('comp1').selection('disk1');
% disk1 = model.component('comp1').selection('disk1').set('entitydim', 1);
% disk1 = model.component('comp1').selection('disk1').set('r', 0);

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
theta_e = 2;                                % intervallo angolare elettrico [deg]
w = per.EvalSpeed*pi/30;                    % velocità di rotazione [rad/s]
theta = theta_e*pi/180/p;                   % angolo meccanico [rad]
ts = theta / w;                             % vettore di time steps [s]
t = ts*360/theta_e;                         % time step [s]
freq = w*p/2/pi;                            % frequenza di alimentazione [Hz]
theta_m_tot = w*t*180/pi;                   % angolo meccanico di simulazione [deg]
theta_e_tot = 2*pi*freq*t*180/pi;           % angolo elettrico di simulazione [deg]

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

% evenuale implementazione 'conducting magnet'
% model.component("comp1").physics("rmm").create("cmag1", "ConductingMagnet", 2);
% model.component("comp1").physics("rmm").feature("cmag1").selection().set(5);
% model.component("comp1").physics("rmm").feature("cmag1").feature("north1").selection().set(11);
% model.component("comp1").physics("rmm").feature("cmag1").feature("south1").selection().set(7);
% model.component("comp1").physics("rmm").create("cmag2", "ConductingMagnet", 2);
% model.component("comp1").physics("rmm").feature("cmag2").selection().set(7);
% model.component("comp1").physics("rmm").feature("cmag2").feature("north1").selection().set(19);
% model.component("comp1").physics("rmm").feature("cmag2").feature("south1").selection().set(15);

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

% % Definizione force calculation (rotore)
% T = [];
% 
% for kk = 1:tmp_rot_righe        
%     x = tmp(kk,1);
%     y = tmp(kk,2);
%     disk1.set('posx', x);
%     disk1.set('posy', y);
%     selNumber = disk1.entities();
%     T = horzcat(T, selNumber);
% end
% 
% model.component('comp1').physics('rmm').create('fcal1', 'ForceCalculation', 2);
% model.component('comp1').physics('rmm').feature('fcal1').selection().set(T);

% Vettori assegnazione periodic condition sui bordi
tmp = [];
tmp_rot = rot_boundary(:, 1:3);
tmp_stat = stat_boundary(:, 1:3);
tmp_righe_s = size(tmp_stat, 1);
tmp_righe_r = size(tmp_rot, 1);
% tmp_stat_a = gm.s.boundary;
% tmp_rot_a = gm.r.boundary; 
% tmp_righe_as = size(tmp_stat_a, 1);
% tmp_righe_ar = size(tmp_rot_a, 1);
AP_s = [];
AP_r = [];  
% AP_as = [];
% AP_ar = [];

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

% Definizione periodic condition aria statore (continuità/antiperiodicità)
% for kk = 1:tmp_righe_as
%     x = tmp_stat_a(kk,1);
%     y = tmp_stat_a(kk,2);
%     disk2.set('posx', x);
%     disk2.set('posy', y);
%     selNumber = disk2.entities();       
%     if tmp_stat_a(kk,3)==10
%         AP_as = horzcat(AP_as, selNumber);
%     end
% end
% 
% model.component('comp1').physics('rmm').create('pc2', 'PeriodicCondition', 1);
% model.component('comp1').physics('rmm').feature('pc2').selection().set(AP_as);
% if mod(geo.ps, 2)==0
% model.component('comp1').physics('rmm').feature('pc2').set('PeriodicType', 'Continuity');
% else
% model.component('comp1').physics('rmm').feature('pc2').set('PeriodicType', 'AntiPeriodicity');
% end
% model.component('comp1').physics('rmm').feature('pc2').label("Periodic Condition - Stator Air");

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

% Definizione periodic condition aria rotore (continuità/antiperiodicità)
% for kk = 1:tmp_righe_ar
%     x = tmp_rot_a(kk,1);
%     y = tmp_rot_a(kk,2);
%     disk2.set('posx', x);
%     disk2.set('posy', y);
%     selNumber = disk2.entities();       
%     if tmp_rot_a(kk,3)==10
%         AP_ar = horzcat(AP_ar, selNumber);
%     end
% end
% 
% model.component('comp1').physics('rmm').create('pc4', 'PeriodicCondition', 1);
% model.component('comp1').physics('rmm').feature('pc4').selection().set(AP_ar);
% if mod(geo.ps, 2)==0
% model.component('comp1').physics('rmm').feature('pc4').set('PeriodicType', 'Continuity');
% else
% model.component('comp1').physics('rmm').feature('pc4').set('PeriodicType', 'AntiPeriodicity');
% end
% model.component('comp1').physics('rmm').feature('pc4').label("Periodic Condition - Rotor Air");

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

% Definizione moving mesh
T = [];
tmp = rot_mat;

% % Calcolo tempo di simulazione
% 
% p = geo.p;                                  % paia poli macchina
% theta_e = 2;                                % intervallo angolare elettrico [deg]
% w = per.EvalSpeed*pi/30;                    % velocità di rotazione [rad/s]
% theta = theta_e*pi/180/p;                   % angolo meccanico [rad]
% ts = theta / w;                             % vettore di time steps [s]
% t = ts*360/theta_e;                         % time step [s]
% freq = w*p/2/pi;                            % frequenza di alimentazione [Hz]
% theta_m_tot = w*t*180/pi;                   % angolo meccanico di simulazione [deg]
% theta_e_tot = 2*pi*freq*t*180/pi;           % angolo elettrico di simulazione [deg]

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

% Definizione Circuito

% Calcolo correnti dq
iAmp = 802.1/6; %per.overload*per.i0/(2*geo.p/geo.ps);         % ampiezza corrente (tiene conto di eventuale sovraccarico)
gamma = 148; %per.GammaPP;                      % angolo del vettore I rispetto all'asse d [deg] (default = 55°, corrisponde all'asse_d considerando l'offset th0)
theta_i = (gamma+geo.th0)*pi/180;             % Angolo corrente [rad] 
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


% import correnti dq
% iAmp = per.overload*per.i0;        % ampiezza corrente (tiene conto di eventuale sovraccarico)
% gamma = dataSet.GammaPP;           % angolo del vettore I rispetto all'asse d [deg]
% th0 = geo.th0;                     % offset [deg]
% 
% id = iAmp*cos(gamma*180/pi);    
% iq = iAmp*sin(gamma*180/pi);    
% Imod = abs(id + 1i* iq);           % Modulo corrente [A]
% theta_i = (th0 + gamma)*pi/180;    % Angolo corrente [rad]
% 
% theta = teta_e*pi/180;          
% pos_0 = 0;
% pos_Mov = t * per.EvalSpeed * pi/30 + pos_0;      % rotor position [rad mec]
% 
% 
% % dq2abc
% i123 = dq2abc(id,iq,theta);
% theta_i = th0 * pi/180 + pos_Mov * geo.p;


% Costruzione mesh
model.component('comp1').mesh('mesh1').autoMeshSize(6);    % mesh size (1-10) dalla più fitta alla meno fitta
model.component('comp1').mesh('mesh1').run();

% Definizione studio
model.param().set('t', '0 [s]');                                          % impostazione tempo 0 s per la simulazione 
model.param().set('l', [num2str(geo.l) ' [mm]']);                         % impostazione lunghezza assiale l [mm]
model.study('std1').setGenIntermediatePlots(true);
% Time dependent
model.study('std1').create('time', 'Transient');
model.study('std1').feature('time').set('tlist', ['range(0,' num2str(ts) ',' num2str(t) ')']);
model.study('std1').feature('time').set('plot', true);
% TimetoFrequency Losses
model.study('std1').create('emloss', 'TimeToFrequencyLosses');
model.study('std1').feature('emloss').set('endtmethod', 'userdef');
% w_loss = 3*(w/2/pi);
model.study('std1').feature('emloss').set('lossstarttime', '0');
model.study('std1').feature('emloss').set('lossendtime', num2str(t));
% Studio Circuito
model.study('std1').feature('stat').setSolveFor('/physics/cir', true);
model.study('std1').feature('time').setSolveFor('/physics/cir', true);
model.study('std1').feature('emloss').setSolveFor('/physics/cir', true);

% Solver Settings
model.sol().create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature('fc1').set('jtech', 'once');
model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 100);
model.sol('sol1').feature('s1').feature('fc1').set('ntolfact', 1);
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('s1').feature('fc1').set('jtech', 'once');
model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 100);
model.sol('sol1').feature('s1').feature('fc1').set('ntolfact', 1);
model.sol('sol1').feature('s1').feature().remove('fcDef');
model.sol('sol1').create('su1', 'StoreSolution');
model.sol('sol1').create('st2', 'StudyStep');
model.sol('sol1').feature('st2').set('study', 'std1');
model.sol('sol1').feature('st2').set('studystep', 'time');
model.sol('sol1').create('v2', 'Variables');
model.sol('sol1').feature('v2').set('initmethod', 'sol');
model.sol('sol1').feature('v2').set('initsol', 'sol1');
model.sol('sol1').feature('v2').set('initsoluse', 'sol2');
model.sol('sol1').feature('v2').set('notsolmethod', 'sol');
model.sol('sol1').feature('v2').set('notsol', 'sol1');
model.sol('sol1').feature('v2').set('notsoluse', 'sol2');
model.sol('sol1').feature('v2').set('control', 'time');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').set('tlist', ['range(0,' num2str(ts) ',' num2str(t) ')']);
model.sol('sol1').feature('t1').set('plot', 'on');
model.sol('sol1').feature('t1').set('plotgroup', 'Default');
model.sol('sol1').feature('t1').set('plotfreq', 'tout');
model.sol('sol1').feature('t1').set('probesel', 'all');
model.sol('sol1').feature('t1').set('probes', {});
model.sol('sol1').feature('t1').set('probefreq', 'tsteps');
model.sol('sol1').feature('t1').set('tout', 'tstepsclosest');
model.sol('sol1').feature('t1').set('atolglobalvaluemethod', 'factor');
model.sol('sol1').feature('t1').set('reacf', true);
model.sol('sol1').feature('t1').set('storeudot', true);
model.sol('sol1').feature('t1').set('tstepsbdf', 'strict');
model.sol('sol1').feature('t1').set('endtimeinterpolation', true);
model.sol('sol1').feature('t1').set('estrat', 'exclude');
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature('fc1').set('jtech', 'once');
model.sol('sol1').feature('t1').feature('fc1').set('maxiter', 25);
model.sol('sol1').feature('t1').feature('fc1').set('ntolfact', 0.2);
model.sol('sol1').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('t1').feature('fc1').set('jtech', 'once');
model.sol('sol1').feature('t1').feature('fc1').set('maxiter', 25);
model.sol('sol1').feature('t1').feature('fc1').set('ntolfact', 0.2);
model.sol('sol1').feature('t1').feature().remove('fcDef');
model.sol('sol1').create('su2', 'StoreSolution');
model.sol('sol1').create('st3', 'StudyStep');
model.sol('sol1').feature('st3').set('study', 'std1');
model.sol('sol1').feature('st3').set('studystep', 'emloss');
model.sol('sol1').create('v3', 'Variables');
model.sol('sol1').feature('v3').set('initmethod', 'sol');
model.sol('sol1').feature('v3').set('initsol', 'sol1');
model.sol('sol1').feature('v3').set('initsoluse', 'sol3');
model.sol('sol1').feature('v3').set('notsolmethod', 'sol');
model.sol('sol1').feature('v3').set('notsol', 'sol1');
model.sol('sol1').feature('v3').set('control', 'emloss');
model.sol('sol1').create('fft1', 'FFT');
model.sol('sol1').feature('fft1').set('control', 'emloss');
model.study('std1').feature('emloss').set('fftstarttime', '0');
model.study('std1').feature('emloss').set('fftendtime', num2str(t));
model.study('std1').feature('emloss').set('fftmaxfreq', ['6/(' num2str(t) '-0)']);
model.sol('sol1').create('su3', 'StoreSolution');
model.sol('sol1').create('cms1', 'CombineSolution');
model.sol('sol1').feature('cms1').set('soloper', 'gensum');
model.sol('sol1').feature('cms1').setIndex('gensumexpressionactive', 'on', 7);
model.sol('sol1').feature('cms1').setIndex('gensumexpression', 'comp1.rmm.Qloss', 7);
model.sol('sol1').feature('cms1').set('control', 'emloss');
model.sol('sol1').feature('v3').set('notsolnum', 'auto');
model.sol('sol1').feature('v3').set('notsolvertype', 'solnum');
model.sol('sol1').attach('std1');
model.sol('sol1').feature('t1').set('maxorder', 2);
model.sol('sol1').feature('t1').feature('aDef').set('cachepattern', false);
model.sol('sol1').feature('t1').feature('fc1').set('dtech', 'auto');
model.sol('sol1').feature('t1').feature('fc1').set('rstep', 4);

% Impostazione bobine nel circuit interface (NB: non spostare i comandi da qui!!)
model.component('comp1').physics('cir').feature('termI1').set('V_src', 'root.comp1.rmm.VCoil_1');
model.component('comp1').physics('cir').feature('termI2').set('V_src', 'root.comp1.rmm.VCoil_2');
model.component('comp1').physics('cir').feature('termI3').set('V_src', 'root.comp1.rmm.VCoil_3');

% mphrun(model);                  % run intero modello (con barra di progresso)

% model.study('std1').run();      % run singolo studio (senza barra di progresso)

%%----------------Studio Meccanico------------------%%
model.component().create('comp2', true);
model.component('comp2').geom().create('geom2', 2);
model.component('comp2').mesh().create('mesh2');
model.component('comp2').geom('geom2').lengthUnit('mm');

load(fullfile(pathname, filename));

magn_dxf = fullfile('C:\Users\S296193\Dropbox (Politecnico Di Torino Studenti)\Tesi Magistrale\syre_developers_20240410\COMSOL\Tests\dxf\', [file_name, '_magn.dxf']);
rot_mecc_dxf = fullfile('C:\Users\S296193\Dropbox (Politecnico Di Torino Studenti)\Tesi Magistrale\syre_developers_20240410\COMSOL\Tests\dxf\', [file_name, '_rot_mecc.dxf']);
flux_barr_dxf = fullfile('C:\Users\S296193\Dropbox (Politecnico Di Torino Studenti)\Tesi Magistrale\syre_developers_20240410\COMSOL\Tests\dxf\', [file_name, '_flux_barr.dxf']);

magn_dxf = fullfile(dxf_dir, [file_name, '_magn.dxf']);
rot_mecc_dxf = fullfile(dxf_dir, [file_name, '_rot_mecc.dxf']);
flux_barr_dxf = fullfile(dxf_dir, [file_name, '_flux_barr.dxf']);

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

model.component('comp2').selection().create('disk3', 'Disk');
disk3 = model.component('comp2').selection('disk3');
disk3 = model.component('comp2').selection('disk3').set('entitydim', 2);

model.component('comp2').selection().create('disk4', 'Disk');
disk4 = model.component('comp2').selection('disk4');
disk4 = model.component('comp2').selection('disk4').set('entitydim', 1);
disk4 = model.component('comp2').selection('disk4').set('r', 0.1);

model.study().create('std2');
model.study('std2').create('stat', 'Stationary');
model.study('std2').feature('stat').setSolveFor('/physics/rmm', true);
model.study('std2').feature('stat').setSolveFor('/physics/cir', true);
model.study('std2').feature('stat').setEntry('activate', 'cir', false);
model.study('std2').feature('stat').setEntry('activate', 'rmm', false);
model.component('comp2').physics().create('solid', 'SolidMechanics', 'geom2');
model.component('comp2').physics('solid').prop('ShapeProperty').set('order_displacement', '2s');
model.component('comp2').physics('solid').prop('d').set('d', 0.35*1e-3);
model.study('std1').feature('stat').setSolveFor('/physics/solid', true);
model.study('std1').feature('time').setSolveFor('/physics/solid', true);
model.study('std1').feature('emloss').setSolveFor('/physics/solid', true);
model.study('std2').feature('stat').setSolveFor('/physics/solid', true);
model.study('std1').feature('stat').setEntry('activate', 'solid', false);
model.study('std1').feature('time').setEntry('activate', 'solid', false);

% definizione materiali

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
    disk3.set('posy', y);
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
model.component('comp2').physics('solid').feature('rotf1').set('Ovm', 18100*pi/30);

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

%Mesh settings
model.component('comp2').mesh('mesh2').autoMeshSize(1);

% model.study('std2').run();      % run singolo studio (senza barra di progresso)

mphrun(model);                  % run intero modello (con barra di progresso)

%%
%%-----------Post-processing meccanico--------------%%

width    = 12;
height   = 10;
textSize = 12;

% stress von Mises deform
model.result('pg5').set('data', 'dset6');
model.result('pg5').set('defaultPlotID', 'stress');
model.result('pg5').label('Stress (solid)');
model.result('pg5').set('frametype', 'spatial');
model.result('pg5').set('edges', false);
model.result('pg5').feature('surf1').set('expr', {'solid.misesGp'});
model.result('pg5').feature('surf1').set('threshold', 'manual');
model.result('pg5').feature('surf1').set('thresholdvalue', 0.2);
model.result('pg5').feature('surf1').set('colortable', 'Rainbow');
model.result('pg5').feature('surf1').set('colortabletrans', 'none');
model.result('pg5').feature('surf1').set('colorscalemode', 'linear');
model.result('pg5').feature('surf1').set('resolution', 'normal');
model.result('pg5').feature('surf1').set('colortable', 'Prism');
model.result('pg5').feature('surf1').feature('def').set('scaleactive', true);
model.result('pg5').feature('surf1').feature('def').set('scale', 100);
model.result('pg5').feature('surf1').feature('def').set('expr', {'u', 'v'});
model.result('pg5').feature('surf1').feature('def').set('descr', 'Displacement field');
model.result('pg5').feature('surf1').set('unit', 'MPa');
model.result('pg5').feature('surf1').set('colortable', 'RainbowLightClassic');

figure()
% set(gcf,'color','w')
mphplot(model, 'pg5', 'rangenum', 1);
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
model.result().create("pg7", "PlotGroup2D");
model.result("pg7").set("data", "dset6");
model.result("pg7").create("surf1", "Surface");
model.result('pg7').set('edges', false);
model.result("pg7").feature("surf1").set("expr", "solid.misesGp");
model.result("pg7").feature("surf1").set("descr", "von Mises stress");
model.result("pg7").feature("surf1").set("unit", "MPa");
model.result("pg7").feature("surf1").set("colortable", "RainbowLightClassic");

figure()
% set(gcf,'color','w')
mphplot(model, 'pg7', 'rangenum', 1);
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
model.result().create('pg6', 'PlotGroup2D');
model.result('pg6').set('data', 'dset6');
model.result('pg6').label('Displacement');
model.result('pg6').set('frametype', 'spatial');
model.result('pg6').set('edges', false);
model.result('pg6').create('surf1', 'Surface');
model.result('pg6').feature('surf1').create('def1', 'Deform');
model.result('pg6').feature('surf1').feature('def1').set('scaleactive', true);
model.result('pg6').feature('surf1').feature('def1').set('scale', 1);
model.result('pg6').feature('surf1').set('colortable', 'RainbowLightClassic');

figure()
% set(gcf,'color','w')
mphplot(model, 'pg6', 'rangenum', 1);
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


% Points exceeding the yield limit
model.result().export().create('plot1', 'pg5', 'surf1', 'Plot');
model.result().export('plot1').set('filename', Stress_csv);
model.result().export('plot1').run();

Stress = readmatrix(fullfile(Stress_csv));

stress_rot =[];
x_stress = [];
y_stress = [];

for kk = 1:size(Stress, 1)
    if Stress(kk, 3) > mat.Rotor.sigma_max
       x_stress = Stress(kk, 1); 
       y_stress = Stress(kk, 2);    
    end
    stress_rot = [stress_rot; x_stress, y_stress];
end

if any(stress_rot(:) ~= 0)
figure ()
figSetting()
xlabel('[mm]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
ylabel('[mm]', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');
title('Von Mises Stress Limit Exceeded', 'FontWeight', 'normal', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');
GUI_Plot_Machine(gca, geo.rotor(1:end-3, :))
plot(stress_rot(:,1), stress_rot(:,2), 'ro', 'MarkerSize', 1, 'LineWidth', 2);
end

%% Post-processing

% gamma_p = [0:360, 361]; % Vettore da 0 a 360 gradi con step di 1 grado
% id_p = iAmp*cos(gamma_p*pi/180);    
% iq_p = iAmp*sin(gamma_p*pi/180); 
% 
% figure (10); 
% hold on;
% plot(id_p, iq_p, 'k-'); 
% quiver(0, 0, real(Imod(end)*cos((gamma-55)*pi/180)), real(Imod(end)*sin((gamma-55)*pi/180)), 'r', 'LineWidth', 2);
% xlabel('id');
% ylabel('iq');
% title('I nel piano dq'); 
% grid on;
% hold off;

% Estrazione flussi concatenati
% model.result().numerical().create('gev1', 'EvalGlobal');
% model.result().numerical('gev1').set('data', 'dset3');
% model.result().numerical('gev1').set('expr', {'rmm.PhiCoil_1', 'rmm.PhiCoil_2', 'rmm.PhiCoil_3'});
% model.result().numerical('gev1').set('descr', {'Lambda 1', 'Lambda 2', 'Lambda 3'});
% model.result().table().create('tbl1', 'Table');
% model.result().table('tbl1').comments('Global Evaluation 1');
% model.result().numerical('gev1').set('table', 'tbl1');
% model.result().numerical('gev1').setResult();
% model.result().table('tbl1').set('storetable', 'onfile');
% model.result().table('tbl1').set('filename', 'C:\\Users\\S296193\\Dropbox (Politecnico Di Torino Studenti)\\Tesi Magistrale\\SyR-e\\SyR-e\\concatenated_flux.csv');

% concatenated_flux = readmatrix('C:\\Users\\S296193\\Dropbox (Politecnico Di Torino Studenti)\\Tesi Magistrale\\SyR-e\\SyR-e\\concatenated_flux.csv');
% tempo = [concatenated_flux(:, 1)]';
% lambda_1 = [concatenated_flux(:, 2)]'; 
% lambda_3 = [concatenated_flux(:, 3)]'; 
% lambda_2 = [concatenated_flux(:, 4)]'; 
% theta_l = tempo*2*pi*freq/p;    
% 
% lambda_dq = abc2dq(lambda_1, lambda_3, lambda_2, theta_l);
% lambda_d = lambda_dq(1, :);
% lambda_q = lambda_dq(2, :);
% 
% Lambda = abs(lambda_d + 1i*lambda_q);
% Arg_L = angle(Lambda(end));
% 
% theta_pe = [0:theta_e_tot, theta_e_tot + 1];                         % Vettore da 0 a theta_e_tot [deg] con step di 1 grado
% ld_p = Lambda(end)*cos(theta_pe*pi/180);    
% lq_p = Lambda(end)*sin(theta_pe*pi/180);

% Calcolo delle componenti x e y della freccia
% arrow_length = real(Lambda(end));               % Lunghezza della freccia
% arrow_angle = Arg_L;                            % Angolo della freccia in radianti
% arrow_x = arrow_length * cos(arrow_angle);     
% arrow_y = arrow_length * sin(arrow_angle);     
% 
% figure (2); 
% hold on;
% plot(ld_p, lq_p, 'k-'); 
% quiver(0, 0, arrow_x', arrow_y, 'r', 'LineWidth', 2);
% xlabel('\lambda d');
% ylabel('\lambda q');
% title('\lambda nel piano dq'); 
% grid on;
% hold off;
% 
% figure (4);
% hold on;
% plot(theta_pe, ld_p, 'r');
% hold on;
% plot(theta_pe, lq_p, 'b');
% xlabel('\theta e');
% ylabel('\lambda d (r) \lambda q (b)');
% title('\lambda dq in funzione di \theta e'); 
% grid on;
% hold off;

% % data = mpheval(model,'rmm.normB','dataset','dset1');
% % mphplot(data, 'rangenum', 1, 'colortable', 'Prism', 'mesh', 'off')
% % hold on
% % mphgeom(model, 'geom1', 'facemode', 'off')
% 
% % data = mpheval(model,'rmm.normB','dataset','dset1');
% Munit = mphevalpointmatrix(model,'rmm.normB','dataset','dset1');   % matrice di B durante tutta la simulazione
% model.result('pg1').feature('surf1').stepFirst(0);                 % set a 0s della geometria per la figura 
% model.result('pg1').feature('surf1').set('colortable', 'Prism');   % colortable
% 
% Perdite nel ferro
% figure(1)
% mphplot(model, 'pg1', 'rangenum', 1);                              % plot del plotgroup1

% Induzione magnetica
% model.result().create('pg2', 'PlotGroup2D');
% model.result('pg2').set('data', 'dset3');
% model.result('pg2').stepFirst(0);
% model.result('pg2').label('Magnetic Flux Density');
% model.result('pg2').create('surf2', 'Surface');
% model.result('pg2').feature('surf2').set('data', 'dset3');
% model.result('pg2').feature('surf2').stepFirst(0);
% model.result('pg2').feature('surf2').set('colortable', 'Prism');
% model.result('pg2').create('str2', 'Streamline');
% model.result('pg2').feature('str2').set('data', 'dset3');
% model.result('pg2').feature('str2').stepFirst(0);
% model.result('pg2').feature('str2').set('posmethod', 'uniform');
% model.result('pg2').feature('str2').set('udist', 0.015);

% figure(2)
% mphplot(model, 'pg2', 'rangenum', 1);                              % plot del plotgroup2

% model.result('pg3').stepFirst(0);
% 
% figure(3)
% mphplot(model, 'pg3', 'rangenum', 1);                              % plot del plotgroup3 (Induzione magnetica in rotazione)

% Flussi concatenati singoli
% model.result().create('pg4', 'PlotGroup1D');
% model.result('pg4').create('ptgr1', 'PointGraph');
% model.result('pg4').feature('ptgr1').set('markerpos', 'datapoints');
% model.result('pg4').feature('ptgr1').set('linewidth', 'preference');
% model.result('pg4').feature('ptgr1').selection().set(99);
% model.result('pg4').set('data', 'dset3');
% model.result('pg4').feature('ptgr1').set('expr', 'rmm.PhiCoil_1');
% model.result('pg4').feature('ptgr1').set('descr', 'Coil concatenated flux');
% model.result('pg4').feature().duplicate('ptgr2', 'ptgr1');
% model.result('pg4').feature('ptgr2').set('expr', 'rmm.PhiCoil_2');
% model.result('pg4').feature().duplicate('ptgr3', 'ptgr2');
% model.result('pg4').feature('ptgr3').set('expr', 'rmm.PhiCoil_3');

% figure(4)
% mphplot(model, 'pg4', 'rangenum', 1);                              % plot del plotgroup4 (perdite nel ferro in funzione di f)

% % Flussi concatenati 
% model.result().create('pg5', 'PlotGroup1D');
% model.result('pg5').set('data', 'dset3');
% model.result('pg5').create('glob1', 'Global');
% model.result('pg5').feature('glob1').set('markerpos', 'datapoints');
% model.result('pg5').feature('glob1').set('linewidth', 'preference');
% model.result('pg5').feature('glob1').set('expr', {'rmm.PhiCoil_1', 'rmm.PhiCoil_2', 'rmm.PhiCoil_3'});
% model.result('pg5').feature('glob1').set('descr', {'Coil concatenated flux', 'Coil concatenated flux', 'Coil concatenated flux'});
% 
% figure(5)
% mphplot(model, 'pg5', 'rangenum', 1);                              % plot del plotgroup5

% % Correnti
% model.result().create('pg6', 'PlotGroup1D');
% model.result('pg6').set('data', 'dset3');
% model.result('pg6').create('glob1', 'Global');
% model.result('pg6').feature('glob1').set('markerpos', 'datapoints');
% model.result('pg6').feature('glob1').set('linewidth', 'preference');
% model.result('pg6').feature('glob1').set('expr', {'rmm.ICoil_1', 'rmm.ICoil_2', 'rmm.ICoil_3'});
% model.result('pg6').feature('glob1').set('descr', {'Coil current', 'Coil current', 'Coil current'});
% 
% figure(6)
% mphplot(model, 'pg6', 'rangenum', 1);                              % plot del plotgroup6

% Coppia Motrice (Arkkio)
% model.result().create('pg7', 'PlotGroup1D');
% model.result('pg7').create('glob1', 'Global');
% model.result('pg7').feature('glob1').set('markerpos', 'datapoints');
% model.result('pg7').feature('glob1').set('linewidth', 'preference');
% model.result('pg7').set('data', 'dset3');
% model.result('pg7').feature('glob1').set('expr', {'rmm.Tark_1'});
% model.result('pg7').feature('glob1').set('descr', {'Axial torque'});
% model.result('pg7').feature('glob1').set('unit', {'N*m'});
% 
% figure(7)
% mphplot(model, 'pg7', 'rangenum', 1);                              % plot del plotgroup7 

%% Post-processing Elettromagnetico (Creazione Dataset)

% Coppia Motrice (Arkkio)
model.result().numerical().create('gev1', 'EvalGlobal');
model.result().table().create('tbl1', 'Table');
model.result().export().create('tbl1', 'Table');

% Flussi concatenati dq
model.result().numerical().create('gev2', 'EvalGlobal');
model.result().table().create('tbl2', 'Table');
model.result().export().create('tbl2', 'Table');

% Correnti abc
model.result().numerical().create('gev3', 'EvalGlobal');
model.result().table().create('tbl3', 'Table');
model.result().export().create('tbl3', 'Table');

% Pfe rotore
model.result().numerical().create('int1', 'IntSurface');
model.result().table().create('tbl4', 'Table');
model.result().export().create('tbl4', 'Table');

% Pfe statore
model.result().numerical().create('int2', 'IntSurface');
model.result().table().create('tbl5', 'Table');
model.result().export().create('tbl5', 'Table');

% PM losses
model.result().numerical().create('int3', 'IntSurface');
model.result().table().create('tbl6', 'Table');
model.result().export().create('tbl6', 'Table');

% PM losses
model.result().numerical().create('int4', 'IntSurface');
model.result().table().create('tbl7', 'Table');
model.result().export().create('tbl7', 'Table');

pause(0.1)

%% Post-processing Elettromagnetico (Salvataggio Risultati)

% Coppia Motrice (Arkkio)
model.result().numerical('gev1').set('data', 'dset3');
model.result().numerical('gev1').set('expr', {'rmm.Tark_1'});
model.result().numerical('gev1').set('descr', {'Axial Torque'});
model.result().table('tbl1').comments('Global Evaluation 1');
model.result().numerical('gev1').set('table', 'tbl1');
model.result().numerical('gev1').setResult();
model.result().export('tbl1').set('filename', torque_csv);
model.result().export('tbl1').set('table', 'tbl1');
model.result().export('tbl1').run();

% Flussi concatenati dq
model.result().numerical('gev2').set('data', 'dset3');
model.result().numerical('gev2').set('expr', {'rmm.PhiCoil_1', 'rmm.PhiCoil_2', 'rmm.PhiCoil_3'});
model.result().numerical('gev2').set('descr', {'Lambda 1', 'Lambda 2', 'Lambda 3'});
model.result().table('tbl2').comments('Global Evaluation 2');
model.result().numerical('gev2').set('table', 'tbl2');
model.result().numerical('gev2').setResult();
model.result().export('tbl2').set('filename', flux_csv);
model.result().export('tbl2').set('table', 'tbl2');
model.result().export('tbl2').run();

% Correnti abc
model.result().numerical('gev3').set('data', 'dset3');
model.result().numerical('gev3').set('expr', {'rmm.ICoil_1', 'rmm.ICoil_2', 'rmm.ICoil_3'});
model.result().numerical('gev3').set('descr', {'Coil current', 'Coil current', 'Coil current'});
model.result().table('tbl3').comments('Global Evaluation 3');
model.result().numerical('gev3').set('table', 'tbl3');
model.result().numerical('gev3').setResult();
model.result().export('tbl3').set('filename', current_csv);
model.result().export('tbl3').set('table', 'tbl3');
model.result().export('tbl3').run();

% Pfe rotore
model.result().numerical('int1').set('intvolume', true);
model.result().numerical('int1').set('expr', 'rmm.Qh');
model.result().numerical('int1').set('descr', 'Volumetric loss density, electromagnetic');
model.result().numerical('int1').setIndex('expr', 'rmm.Qh*l', 0);
model.result().numerical('int1').selection().set(BC_fe_r);
model.result().table('tbl4').comments('Surface Integration 1');
model.result().numerical('int1').set('table', 'tbl4');
model.result().numerical('int1').setResult();
model.result().export('tbl4').set('filename', Pfe_r_csv);
model.result().export('tbl4').set('table', 'tbl4');
model.result().export('tbl4').run();

% Pfe statore
model.result().numerical('int2').set('intvolume', true);
model.result().numerical('int2').set('expr', 'rmm.Qh');
model.result().numerical('int2').set('descr', 'Volumetric loss density, electromagnetic');
model.result().numerical('int2').setIndex('expr', 'rmm.Qh*l', 0);
model.result().numerical('int2').selection().set(BC_fe_s);
model.result().table('tbl5').comments('Surface Integration 1');
model.result().numerical('int2').set('table', 'tbl5');
model.result().numerical('int2').setResult();
model.result().export('tbl5').set('filename', Pfe_s_csv);
model.result().export('tbl5').set('table', 'tbl5');
model.result().export('tbl5').run();

% PM losses
model.result().numerical('int3').set('intvolume', true);
model.result().numerical('int3').set('expr', 'rmm.Qh');
model.result().numerical('int3').set('descr', 'Volumetric loss density, electromagnetic');
model.result().numerical('int3').setIndex('expr', 'rmm.Qh*l', 0);
model.result().numerical('int3').selection().set(BC_pm);
model.result().table('tbl6').comments('Surface Integration 1');
model.result().numerical('int3').set('table', 'tbl6');
model.result().numerical('int3').setResult();
model.result().export('tbl6').set('filename', PM_losses_csv);
model.result().export('tbl6').set('table', 'tbl6');
model.result().export('tbl6').run();

% Pj barriere di flusso rotore
model.result().numerical('int4').set('intvolume', true);
model.result().numerical('int4').set('expr', 'rmm.Qh');
model.result().numerical('int4').set('descr', 'Volumetric loss density, electromagnetic');
model.result().numerical('int4').setIndex('expr', 'rmm.Qh*l', 0);
model.result().numerical('int4').selection().set(Bar);
model.result().table('tbl7').comments('Surface Integration 1');
model.result().numerical('int4').set('table', 'tbl7');
model.result().numerical('int4').setResult();
model.result().export('tbl7').set('filename', PjrBar_csv);
model.result().export('tbl7').set('table', 'tbl7');
model.result().export('tbl7').run();

% Saving Updated Model

mphsave(model, 'Test_auto_model_updated_dm_multifisico.mph');

% pause(0.1);

%% Post-processing Elettromagnetico (Estrazione Grafici)

% Coppia Motrice (Arkkio)
axial_torque = readmatrix(fullfile(torque_csv));
tempo = [axial_torque(:, 1)]';
Tor = [axial_torque(:, 2)]'; 
theta_elt = tempo*2*pi*freq+geo.th0*pi/180-pi/2;
theta_elt_deg = theta_elt*180/pi;

% Flussi concatenati 
concatenated_flux = readmatrix(flux_csv);
lambda_1 = [concatenated_flux(:, 2)]'; 
lambda_2 = [concatenated_flux(:, 3)]'; 
lambda_3 = [concatenated_flux(:, 4)]';

lambda_dq = abc2dq(lambda_1, lambda_2, lambda_3, theta_elt);
lambda_d = lambda_dq(1, :);
lambda_q = lambda_dq(2, :);

lambda_mod = abs(lambda_d + 1i*lambda_q);         % modulo flusso [Vs]
lambda_arg = angle(lambda_mod(end));              % angolo flusso [rad]

% Correnti
currents = readmatrix(current_csv);
tempo_i = [currents(:, 1)]';
i_1 = [currents(:, 2)]'*2*geo.p/geo.ps; 
i_2 = [currents(:, 3)]'*2*geo.p/geo.ps; 
i_3 = [currents(:, 4)]'*2*geo.p/geo.ps; 
theta_i_deg = tempo_i*2*pi*freq/pi*180+geo.th0-90;

i_dq = abc2dq(i_1, i_2, i_3, theta_elt);
i_d = i_dq(1, :);
i_q = i_dq(2, :);

iqm = mean(i_q);
idm = mean(i_d);

lqm = mean(lambda_q);
ldm = mean(lambda_d);

% IPF (internal power factor)
IPF = sin(atan2(i_q,i_d)-atan2(lambda_q,lambda_d));
% IPF = sin(Iarg - lambda_arg);

% Plot Torque - IPF
figure();

set(gcf,'color','w');
subplot(2,1,1);
hold on;
plot((theta_elt_deg - theta_elt_deg(1)), Tor, 'r', 'LineWidth', 1.5);
mean_Tor = mean(Tor);
grid on;
set(gca, 'GridLineStyle', ':');
xlabel('$\theta$ [elt deg]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
ylabel('[Nm]', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');
set(gca, 'XLim', [0 360], 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XTick', 0:60:360, 'TickLabelInterpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
title(['Mean Torque = ', num2str(mean_Tor), ' Nm'], 'FontWeight', 'normal', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');
box on;

subplot(2,1,2);
hold on;
plot((theta_elt_deg - theta_elt_deg(1)), IPF, 'r', 'LineWidth', 1.5);
mean_IPF = mean(IPF);
grid on;
set(gca, 'GridLineStyle', ':');
xlabel('$\theta$ [elt deg]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
ylabel('IPF', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');
set(gca, 'XLim', [0 360], 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XTick', 0:60:360, 'TickLabelInterpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
title(['Mean IPF = ', num2str(mean_IPF)], 'FontWeight', 'normal', 'FontName', 'Times', 'FontSize', 12, 'Interpreter', 'latex');
box on;

% Plot concatenated flux dq
figure();

set(gcf,'color','w');
subplot(2,1,1);                    % Plot lambda_d
hold on;
plot((theta_elt_deg - theta_elt_deg(1)), lambda_d, 'r', 'LineWidth', 1.5);
mean_lambda_d = mean(lambda_d);
grid on;
set(gca, 'GridLineStyle', ':');
xlabel('$\theta$ [elt deg]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
ylabel('$\lambda_d$ [Vs]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XLim', [0 360], 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XTick', 0:60:360, 'TickLabelInterpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);      
title(['Mean $\lambda_d$ = ', num2str(mean_lambda_d), ' Vs'], 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
box on;

subplot(2,1,2);                    % Plot lambda_q
hold on;
plot((theta_elt_deg - theta_elt_deg(1)), lambda_q, 'r', 'LineWidth', 1.5);
mean_lambda_q = mean(lambda_q);
grid on;
set(gca, 'GridLineStyle', ':');
xlabel('$\theta$ [elt deg]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
ylabel('$\lambda_q$ [Vs]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XLim', [0 360], 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XTick', 0:60:360, 'TickLabelInterpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
title(['Mean $\lambda_q$ = ', num2str(mean_lambda_q), 'Vs'], 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
box on;

% Plot concatenated flux abc - currents abc
figure();

set(gcf,'color','w');

subplot(2,1,1);                                    % Plot lambda abc
plot((theta_elt_deg - theta_elt_deg(1)), lambda_1, 'b', 'LineWidth', 1.5);
hold on;
plot((theta_elt_deg - theta_elt_deg(1)), lambda_2, 'r', 'LineWidth', 1.5);
hold on;
plot((theta_elt_deg - theta_elt_deg(1)), lambda_3, 'g', 'LineWidth', 1.5);
grid on;
set(gca, 'GridLineStyle', ':');
xlabel('$\theta$ [elt deg]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
ylabel('$\lambda_{abc}$ [Vs]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XLim', [0 360], 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XTick', 0:60:360, 'TickLabelInterpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
title('Phase flux linkages', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontName', 'Times', 'FontSize', 12);
box on;

subplot(2,1,2);                                    % Plot i abc
plot((theta_i_deg - theta_i_deg(1)), i_1, 'b', 'LineWidth', 1.5);
hold on;
plot((theta_i_deg - theta_i_deg(1)), i_2, 'r', 'LineWidth', 1.5);
hold on;
plot((theta_i_deg - theta_i_deg(1)), i_3, 'g', 'LineWidth', 1.5);
grid on;
set(gca, 'GridLineStyle', ':');
xlabel('$\theta$ [elt deg]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
ylabel('$i_{abc}$ [A]', 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XLim', [0 360], 'FontName', 'Times', 'FontSize', 12);
set(gca, 'XTick', 0:60:360, 'TickLabelInterpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
title('Phase currents', 'Interpreter', 'latex', 'FontWeight', 'normal', 'FontName', 'Times', 'FontSize', 12);
box on;
legend({'$i_a$', '$i_b$', '$i_c$'}, 'Interpreter', 'latex', 'FontName', 'Times', 'FontSize', 12);
