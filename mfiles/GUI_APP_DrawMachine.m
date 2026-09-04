% Copyright 2019
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

function app = GUI_APP_DrawMachine(app)

% flag_plot = 'Y';
h = app.AxisGeometry;
dataSet = app.dataSet;

cla(h);


if app.dataSet.custom
    app = customDraw(app);
    dataSet = app.dataSet;
end 

[~, ~, geo,per,mat] = data0(dataSet);
[geo,gamma,mat] = interpretRQ(geo.RQ,geo,mat);

geo.x0 = geo.r/cos(pi/2/geo.p);
% geo.x0 = (geo.r-geo.hs)/cos(pi/2/geo.p);

fem.res = 0;
fem.res_traf = 0;

% nodes
[rotor,~,geo] = ROTmatr(geo,fem,mat);
[geo,stator,~] = STATmatr(geo,fem);


if strcmp(geo.RotType, 'Hybrid')
    warning('Model under development ')
    GUI_Plot_Hybrid(h,rotor,stator);
else
    GUI_Plot_Machine(h,rotor);
    GUI_Plot_Machine(h,stator);
end


% Axis limits (to center the figure)
set(h,'dataAspectRatio',[1 1 1]);
xMax = geo.R*1.05;
if geo.ps>=geo.p
    xMin = -xMax;
else
    xMin = geo.R*cos(pi/geo.p*geo.ps)*1.05;
end

if xMin>0
    xMin = -geo.R*0.05;
end
set(h,'XLim',[xMin xMax]);

% end winding inductance
per.Lend = calc_Lend(geo);
dataSet.Lend = per.Lend;

geo.lend = calc_endTurnLength(geo);

% Rated current computation (thermal model)
if  isnan(dataSet.SimIth0)
    per = calc_i0(geo,per,mat);
    per.tempcuest = temp_est_simpleMod(geo,per);
    
    dataSet.AdmiJouleLosses     = per.Loss;
    dataSet.ThermalLoadKj       = per.kj;
    dataSet.CurrentDensity      = per.J;
    dataSet.RatedCurrent        = per.i0;
    dataSet.Rs                  = per.Rs;
    dataSet.SimulatedCurrent    = per.overload*per.i0;
    dataSet.EstimatedCopperTemp = per.tempcuest;
end
geo.lend = calc_endTurnLength(geo);
per.i0 = dataSet.RatedCurrent;
per.Lend = calc_Lend(geo);
dataSet.EndWindingsLength = geo.lend;
dataSet.Lend = per.Lend;

[~,geo] = calc_endTurnFieldLength(geo);
dataSet.RotorEndWindingsLength = geo.lendf;
per = calc_if(geo,per,mat);
dataSet.RotorCurrentDensity = per.Jf;
dataSet.RatedFieldCurrent = per.if0;
dataSet.Rf = per.Rf;
dataSet.RotorCurrentDensity = per.Jf;
per.if = per.if0;

% Mass and Inertia computation
geo.pShape = dataSet.pShape;
geo.mCu = calcMassCu(geo,mat);
geo.mPM = calcMassPM(geo,mat);
geo.mAl = calcMassAl(geo,mat);
[geo.mFeS,geo.mFeR] = calcMassFe(geo,mat);
geo.J = calcRotorInertia(geo,mat);


dataSet.MassWinding = geo.mCu;
dataSet.MassMagnet = geo.mPM;
dataSet.MassStatorIron = geo.mFeS;
dataSet.MassRotorIron  = geo.mFeR;
dataSet.MassRotorBar   = geo.mAl;
dataSet.RotorInertia   = geo.J;

% Refresh display
dataSet.SlotWidth = geo.st;
dataSet.SlotConductorParallel = geo.win.pCond;
dataSet.DepthOfBarrier = round(geo.dx,2);
dataSet.HCpu = round(geo.hc_pu,2);
dataSet.betaPMshape = round(geo.betaPMshape,2);
dataSet.YokeLength = geo.ly;
dataSet.RotorSlotFilletTop    = geo.IM.filletTop;
dataSet.RotorSlotFilletBottom = geo.IM.filletBot;

% dataSet.CurrentDensity = per.i0*geo.win.Nbob*2/(sqrt(2)*geo.Aslot*geo.win.kcu);
%dataSet.CurrentDensity = per.i0/(sqrt(2)*geo.Aslot*geo.win.kcu*geo.win.pCond);

dataSet.pontRangEdit = geo.pontRang;
% dataSet.pontRoffsetEdit = geo.pontRoffset;
dataSet.RadRibSplit = geo.radial_ribs_split;
dataSet.RotorFilletIn = geo.RotorFillet1;
dataSet.RotorFilletOut = geo.RotorFillet2;
dataSet.RadShiftInner = geo.dxIB; 
dataSet.NarrowFactor = geo.kOB;

dataSet.CentralShrink = geo.hcShrink;

% dalpha = geo.dalpha;            % barriers ends (deg)
% hc = geo.hc;                    % barriers hieghts (mm)

% set(app.EstimatedCoppTemp,'String',num2str(dataSet.EstimatedCopperTemp));
% set(app.CalculatedRatedCurrent,'String',num2str(dataSet.RatedCurrent));
% set(app.CurrentPP,'String',num2str(dataSet.RatedCurrent));
% set(app.Rsedit,'String',num2str(dataSet.Rs));
if (~strcmp(geo.RotType, 'SPM') && ~strcmp(geo.RotType,'IM') && ~strcmp(geo.RotType,'EESM') && ~strcmp(geo.RotType,'Hybrid'))
    dataSet.ALPHAdeg = round(100*geo.dalpha)/100;
    dataSet.HCmm = round(100*geo.hc)/100;
%     set(app.AlphadegreeEdit,'String',mat2str(dataSet.ALPHAdeg));
%     set(app.hcmmEdit,'String',mat2str(dataSet.HCmm));
elseif strcmp(geo.RotType,'SPM')
    dataSet.ALPHAdeg = round(100*geo.dalpha*2)/100;
    dataSet.HCmm     = round(100*geo.hc)/100;
elseif strcmp(geo.RotType, 'Spoke-type')
    dataSet.TanRibEdit = geo.pontT;
elseif(strcmp(geo.RotType, 'EESM'))
    dataSet.ALPHAdeg = round(geo.dalpha,2);
    dataSet.ALPHApu  = round(geo.dalpha_pu,3);
    dataSet.YokeWidth         = geo.lyr;
    dataSet.PoleBodyHeight    = geo.hpb;
    dataSet.PoleHeadHeight    = geo.hph;
    dataSet.PoleAnglepu       = geo.dalpha_pu;
    dataSet.PoleWidth         = geo.wp;
    dataSet.CoilWidth         = geo.wb;
    dataSet.CoilHeight        = geo.hb;
    dataSet.PoleRotHeadAngle  = geo.thHead_deg;
    dataSet.PoleRotHeadFillet = geo.r_fillet;
    dataSet.PoleRotYokeAngle  = geo.thYoke_deg;
    dataSet.PoleRotYokeFillet = geo.r_bfillet;

%==============================================================
elseif strcmp(geo.RotType, 'Hybrid')
    dataSet.ALPHAdeg = round(100*geo.dalpha)/100;
    dataSet.YokeWidth         = geo.lyr;
    dataSet.PoleBodyHeight    = geo.hpb;
    dataSet.PoleHeadHeight    = geo.hph;
    dataSet.PoleAnglepu       = geo.dalpha_pu;
    dataSet.PoleWidth         = geo.wp;
    dataSet.CoilWidth         = geo.wb;
    dataSet.CoilHeight        = geo.hb;
    dataSet.PoleRotHeadAngle  = geo.thHead_deg;
    dataSet.PoleRotHeadFillet = geo.r_fillet; 

    % Parametri specifici del Magnete per test_VFM
    dataSet.MagWidth  = geo.w_pm;
    % dataSet.MagHeight = geo.h_pm;
    % dataSet.MagOffset = geo.o_pm;  
    dataSet.MagnetOffset = geo.o_pm;
    dataSet.MagnetLength = geo.h_pm;
    dataSet.RadRibEdit = geo.w_bridge_pm;

    

    if isfield(geo, 'AreaPM')
        dataSet.AreaPM = geo.AreaPM;
    end
%=================================================================
end

dataSet.RadRibEdit = round(geo.pontR*100)/100;
% set(app.RadRibEdit,'String',mat2str(tmp));
dataSet.FilletCorner = geo.SFR;
dataSet.RotorFilletTan1  = geo.RotorFilletTan1;
dataSet.RotorFilletTan2  = geo.RotorFilletTan2;
% set(app.FillCorSlotEdit,'String',mat2str(round(100*geo.SFR)/100));
dataSet.ShaftRadius = floor(geo.Ar*100)/100;

dataSet.PMdim   = geo.PMdim;
dataSet.PMclear = geo.PMclear;
dataSet.PMNc    = geo.PMNc;
dataSet.PMNa    = geo.PMNa;
if isfield(geo,'AreaC')
    dataSet.PMdimPU = [geo.AreaC;geo.AreaE]./[geo.AreaCMax;geo.AreaEMax];
    dataSet.PMdimPU(isnan(dataSet.PMdimPU))=0;
else
    dataSet.PMdimPU = nan(2,geo.nlay);
end

%AGGIORNAMENTO PARAMETRI IN GUI

if strcmp(geo.RotType, 'Hybrid')
    dataSet.ALPHAdeg = round(100*geo.dalpha)/100;
    dataSet.YokeWidth         =  round(100*geo.lyr)/100;
    dataSet.PoleBodyHeight    =  round(100*geo.hpb)/100;
    dataSet.PoleHeadHeight    =  round(100*geo.hph)/100;
    dataSet.PoleAnglepu       =  round(100*geo.dalpha_pu)/100;
    dataSet.PoleWidth         =  round(100*geo.wp)/100;
    dataSet.CoilWidth         =  round(100*geo.wb)/100;
    dataSet.CoilHeight        =  round(100*geo.hb)/100;
    dataSet.PoleRotHeadAngle  =  round(100*geo.thHead_deg)/100;
    dataSet.PoleRotHeadFillet =  round(100*geo.r_fillet)/100; 

    % Parametri specifici del Magnete per test_VFM
    dataSet.MagWidth  =  round(100*geo.w_pm)/100;
    % dataSet.MagHeight = geo.h_pm;
    % dataSet.MagOffset = geo.o_pm;  
    dataSet.MagnetOffset =  round(100*geo.o_pm)/100;
    dataSet.MagnetLength =  round(100*geo.h_pm)/100;
    dataSet.RadRibEdit =  round(100*geo.w_bridge_pm)/100;

end

app.dataSet = dataSet;



