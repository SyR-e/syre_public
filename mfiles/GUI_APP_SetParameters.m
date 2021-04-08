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

function GUI_APP_SetParameters(app)

dataSet = app.dataSet;

set(app.currentMotFileName,'Value',dataSet.currentfilename);

% Main data panel
set(app.PolePairsEdit,'Value',num2str(dataSet.NumOfPolePairs));
set(app.NumOfSlotsEdit,'Value',num2str(dataSet.NumOfSlots));
set(app.GapThiEdit,'Value',num2str(dataSet.AirGapThickness));
set(app.StatorOuterRadEdit,'Value',num2str(dataSet.StatorOuterRadius));
set(app.AirGapRadiusEdit,'Value',num2str(dataSet.AirGapRadius));
set(app.ShaftRadEdit,'Value',num2str(dataSet.ShaftRadius));
set(app.StackLenghtEdit,'Value',num2str(dataSet.StackLength));
set(app.TypeOfRotorList,'Value',dataSet.TypeOfRotor);
% Preliminary design
set(app.bRangeEdit,'Value',mat2str(dataSet.bRange));
set(app.xRangeEdit,'Value',mat2str(dataSet.xRange));

set(app.BfeEdit,'Value',num2str(dataSet.Bfe));
set(app.ktEdit,'Value',num2str(dataSet.kt));
set(app.FEAfixPopUp,'Value',int2str(dataSet.FEAfixN));
set(app.CurrLoXBEdit,'Value',num2str(dataSet.CurrLoPP));

% Stator panel
set(app.ToothLengEdit,'Value',num2str(dataSet.ToothLength));
set(app.ToothWidthEdit,'Value',num2str(dataSet.ToothWidth));
if ~dataSet.ParallelSlotCheck
    set(app.SlotshapeDropDown,'Value','Trapezoidal');
else
    set(app.SlotshapeDropDown,'Value','Rectangular');
end
set(app.StatorSlotOpeEdit,'Value',num2str(dataSet.StatorSlotOpen));
set(app.ToothTanDepEdit,'Value',num2str(dataSet.ToothTangDepth));
set(app.ToothTangAngleEdit,'Value',num2str(dataSet.ToothTangAngle));
set(app.FillCorSlotEdit,'Value',num2str(dataSet.FilletCorner));
%% Rotor panel
if strcmp(dataSet.TypeOfRotor,'SPM')
    set(app.NumberOfLayersEdit,'Enable','off');
    set(app.AlphapuEdit,'Enable','off','Value','NaN');
    set(app.AlphadegreeEdit,'Enable','on','Editable','on','Value',mat2str(dataSet.AngleSpanOfPM));
    set(app.hcpuEdit,'Enable','off','Value','NaN');
    set(app.hcmmEdit,'Enable','on','Editable','on','Value',num2str(dataSet.ThicknessOfPM));
    set(app.DxEdit,'Enable','on','Value',num2str(dataSet.DepthOfBarrier));
    set(app.BetaEdit,'Enable','on','Value',mat2str(dataSet.betaPMshape));
    set(app.ThetaFBSEdit,'Enable','off');
    set(app.TanRibEdit,'Enable','off','Value',mat2str(dataSet.TanRibEdit));
    set(app.RadRibEdit,'Enable','off','Value',mat2str(dataSet.RadRibEdit));
    set(app.SplitRibsEditField,'Enable','off','Value','0');
    set(app.RadRibCheck,'Enable','off','Value',dataSet.RadRibCheck);
    set(app.RadialRibsAngleEditField,'Enable','off','Value',mat2str(dataSet.pontRangEdit));
    set(app.CentralBarriersShrinkEdit,'Enable','off','Value',mat2str(dataSet.CentralShrink));
    set(app.RadialRibsOffsetEditField,'Enable','off','Value',mat2str(dataSet.pontRoffsetEdit));
    set(app.RotorFilletRadInEditField,'Enable','off','Value',mat2str(dataSet.RotorFilletIn));
    set(app.RotorFilletRadOutEditField,'Enable','off','Value',mat2str(dataSet.RotorFilletOut));
elseif strcmp(dataSet.TypeOfRotor,'Vtype')
    set(app.NumberOfLayersEdit,'Enable','on','Value',num2str(dataSet.NumOfLayers));
    set(app.AlphapuEdit,'Enable','on','Value',mat2str(dataSet.ALPHApu));
    set(app.AlphadegreeEdit,'Enable','on','Editable','off','Value',mat2str(dataSet.ALPHAdeg));
    set(app.hcpuEdit,'Enable','on','Value',mat2str(dataSet.HCpu));
    set(app.hcmmEdit,'Enable','on','Editable','off','Value',mat2str(dataSet.HCmm));
    set(app.DxEdit,'Enable','on','Value',mat2str(dataSet.DepthOfBarrier));
    set(app.BetaEdit,'Enable','on','Value',mat2str(dataSet.betaPMshape));
    set(app.CentralBarriersShrinkEdit,'Enable','off','Value',mat2str(dataSet.CentralShrink));
    set(app.NarrowFactorEdit,'Enable','off','Value',mat2str(dataSet.NarrowFactor));
    set(app.RadShiftInnerEdit,'Enable','off','Value',mat2str(dataSet.RadShiftInner));
    set(app.ThetaFBSEdit,'Enable','off');
    set(app.TanRibEdit,'Enable','on','Value',mat2str(dataSet.TanRibEdit));
    set(app.RadRibEdit,'Enable','on','Value',mat2str(dataSet.RadRibEdit));
    set(app.SplitRibsEditField,'Enable','off','Value','0');
    set(app.RadRibCheck,'Enable','on','Value',dataSet.RadRibCheck);
    if dataSet.RadRibCheck
        set(app.RadRibEdit,'Enable','on')
    else
        set(app.RadRibEdit,'Enable','off')
    end
    set(app.RadialRibsAngleEditField,'Enable','off','Value',mat2str(dataSet.pontRangEdit));
    set(app.RadialRibsOffsetEditField,'Enable','off','Value',mat2str(dataSet.pontRoffsetEdit));
    set(app.RotorFilletRadInEditField,'Enable','on','Value',mat2str(dataSet.RotorFilletIn));
    set(app.RotorFilletRadOutEditField,'Enable','on','Value',mat2str(dataSet.RotorFilletOut));
else
    set(app.NumberOfLayersEdit,'Enable','on','Value',num2str(dataSet.NumOfLayers));
    set(app.AlphapuEdit,'Enable','on','Value',mat2str(dataSet.ALPHApu));
    set(app.AlphadegreeEdit,'Enable','on','Editable','off','Value',mat2str(dataSet.ALPHAdeg));
    set(app.hcpuEdit,'Enable','on','Value',mat2str(dataSet.HCpu));
    set(app.hcmmEdit,'Enable','on','Editable','off','Value',mat2str(dataSet.HCmm));
    if strcmp(dataSet.TypeOfRotor,'ISeg')
        set(app.DxEdit,'Enable','off','Value',mat2str(dataSet.DepthOfBarrier));
    else
        set(app.DxEdit,'Enable','on','Value',mat2str(dataSet.DepthOfBarrier));
    end
    set(app.BetaEdit,'Enable','off','Value',mat2str(dataSet.betaPMshape));
    if (strcmp(dataSet.TypeOfRotor,'Circular')||strcmp(dataSet.TypeOfRotor,'Seg'))
        set(app.ThetaFBSEdit,'Enable','on','Value',num2str(dataSet.thetaFBS))
    else
        set(app.ThetaFBSEdit,'Enable','off','Value',num2str(dataSet.thetaFBS))
    end
    set(app.TanRibEdit,'Enable','on','Value',mat2str(dataSet.TanRibEdit));
    set(app.RadRibEdit,'Enable','on','Value',mat2str(dataSet.RadRibEdit));
    
    if strcmp(dataSet.TypeOfRotor,'Seg')
        set(app.CentralBarriersShrinkEdit,'Enable','on','Value',mat2str(dataSet.CentralShrink));
        set(app.SplitRibsEditField,'Enable','on','Value',mat2str(dataSet.RadRibSplit));
        set(app.NarrowFactorEdit,'Enable','on','Value',mat2str(dataSet.NarrowFactor));
        set(app.RadShiftInnerEdit,'Enable','on','Value',mat2str(dataSet.RadShiftInner));
        set(app.RotorFilletRadInEditField,'Enable','on','Value',mat2str(dataSet.RotorFilletIn));
        set(app.RotorFilletRadOutEditField,'Enable','on','Value',mat2str(dataSet.RotorFilletOut));
    else
        set(app.SplitRibsEditField,'Enable','off','Value',mat2str(dataSet.RadRibSplit));
        set(app.CentralBarriersShrinkEdit,'Enable','off','Value',mat2str(dataSet.CentralShrink));
        set(app.NarrowFactorEdit,'Enable','off','Value',mat2str(dataSet.NarrowFactor));
        set(app.RadShiftInnerEdit,'Enable','off','Value',mat2str(dataSet.RadShiftInner));
        set(app.RotorFilletRadInEditField,'Enable','off','Value',mat2str(dataSet.RotorFilletIn));
        set(app.RotorFilletRadOutEditField,'Enable','off','Value',mat2str(dataSet.RotorFilletOut));
    end
    
    
    if any(dataSet.RadRibSplit)
        set(app.RadialRibsAngleEditField,'Enable','on','Value',mat2str(dataSet.pontRangEdit));
        set(app.RadialRibsOffsetEditField,'Enable','on','Value',mat2str(dataSet.pontRoffsetEdit));
    else
        set(app.RadialRibsAngleEditField,'Enable','off','Value',mat2str(dataSet.pontRangEdit));
        set(app.RadialRibsOffsetEditField,'Enable','off','Value',mat2str(dataSet.pontRoffsetEdit));
    end
    
    if strcmp(dataSet.TypeOfRotor,'Seg')
        set(app.RotorFilletTangEditField,'Enable','on','Value',mat2str(dataSet.RotorFillet));
    else
        set(app.RotorFilletTangEditField,'Enable','off','Value','0');
    end
    
    set(app.RadRibCheck,'Enable','on','Value',dataSet.RadRibCheck);
    if dataSet.RadRibCheck
        set(app.RadRibEdit,'Enable','on')
    else
        set(app.RadRibEdit,'Enable','off')
    end
end

% Options panel
set(app.ThermalLoadKj,'Value',num2str(dataSet.ThermalLoadKj));
set(app.CopperTempEdit,'Value',num2str(dataSet.TargetCopperTemp));
set(app.HousingTempEdit,'Value',num2str(dataSet.HousingTemp));
set(app.EstimatedCoppTemp,'Enable','on','Value',num2str(dataSet.EstimatedCopperTemp),'Editable','off');
set(app.CalculatedRatedCurrent,'Enable','on','Value',num2str(dataSet.RatedCurrent),'Editable','off');
set(app.JouleLossesEdit,'Enable','on','Value',num2str(dataSet.AdmiJouleLosses),'Editable','off');
set(app.CurrentdensityEdit,'Enable','on','Value',num2str(dataSet.CurrentDensity),'Editable','off');
set(app.Rsedit,'Enable','on','Value',num2str(dataSet.Rs),'Editable','off');
set(app.OverSpeedEdit,'Value',num2str(dataSet.OverSpeed));
set(app.MeshEdit,'Value',num2str(dataSet.Mesh));
set(app.MeshMOOAEdit,'Value',num2str(dataSet.Mesh_MOOA));
set(app.MecTolerEdit,'Value',num2str(dataSet.MinMechTol));


% Windings panel
set(app.SlotFillFacEdit,'Value',num2str(dataSet.SlotFillFactor));
if ~dataSet.SlotLayerPosCheck
    set(app.SlotlayerposDropDown,'Value','Stacked');
else
    set(app.SlotlayerposDropDown,'Value','Side-by-side');
end
set(app.TurnsSeriesEdit,'Value',num2str(dataSet.TurnsInSeries));
set(app.SlotSimulEdit,'Value',num2str(dataSet.Qs));
set(app.PitchWindEdit,'Value',num2str(dataSet.PitchShortFac));
set(app.Num3PhaseCircuitEdit,'Value',num2str(dataSet.Num3PhaseCircuit));
if dataSet.NumOfSlots<1
    set(app.PitchWindEdit,'Enable','off'); % concentrated winding
else
    set(app.PitchWindEdit,'Enable','on'); % distributed winding
end

% Winding table
columnName = cell(1,floor(dataSet.Qs));
wCols      = cell(1,floor(dataSet.Qs));
for ii = 1 : floor(dataSet.Qs)
    columnName{ii} = ['Slot n° ',num2str(ii)];
    wCols{ii} = 50;
end
rowName{1} = 'Layer 1';
rowName{2} = 'Layer 2';
set(app.WinTable,'rowname',rowName,'columnname',columnName,'Data',dataSet.WinMatr(:,1:floor(dataSet.Qs)),'ColumnWidth',wCols);

% slot model
set(app.ConductorTypeEdit,'Value',dataSet.SlotConductorType);
set(app.ConductorInsulationEdit,'Value',num2str(dataSet.SlotConductorInsulation));
%set(app.ConductorShapeEdit,'Value',num2str(dataSet.SlotConductorShape));
set(app.ConductorRadiusEdit,'Value',num2str(dataSet.SlotConductorRadius));
set(app.ConductorWidthEdit,'Value',num2str(dataSet.SlotConductorWidth));
set(app.ConductorHeightEdit,'Value',num2str(dataSet.SlotConductorHeight));
set(app.ConductorNumberEdit,'Value',num2str(dataSet.SlotConductorNumber));
set(app.SlotModelFrequencyEdit,'Value',mat2str(dataSet.SlotConductorFrequency));
set(app.SlotModelTemperatureEdit,'Value',mat2str(round(dataSet.SlotConductorTemperature,1)));

% Material panel
tmp = material_properties_conductor('0');
set(app.StatorslotmaterialDropDown,'Items',tmp.MatList);
tmp = material_properties_layer('0');
set(app.FluxbarriermaterialDropDown,'Items',tmp.MatList);
tmp = material_properties_iron('0');
set(app.StatorcorematerialDropDown,'Items',tmp.MatList);
set(app.RotorcorematerialDropDown,'Items',tmp.MatList);
set(app.ShaftmaterialDropDown,'Items',['Air' tmp.MatList]);

set(app.StatorslotmaterialDropDown,'Value',dataSet.SlotMaterial);
set(app.StatorcorematerialDropDown,'Value',dataSet.StatorMaterial);
set(app.RotorcorematerialDropDown,'Value',dataSet.RotorMaterial);
set(app.FluxbarriermaterialDropDown,'Value',dataSet.FluxBarrierMaterial);
set(app.ShaftmaterialDropDown,'Value',dataSet.ShaftMaterial);

set(app.MassWindingEditField,'Enable','on','Editable','off','Value',num2str(dataSet.MassWinding));
set(app.MassStatorIronEditField,'Enable','on','Editable','off','Value',num2str(dataSet.MassStatorIron));
set(app.MassRotorIronEditField,'Enable','on','Editable','off','Value',num2str(dataSet.MassRotorIron));
set(app.MassFluxBarrierEditField,'Enable','on','Editable','off','Value',num2str(dataSet.MassMagnet));
set(app.MassTotalEditField,'Enable','on','Editable','off','Value',num2str(dataSet.MassWinding+dataSet.MassStatorIron+dataSet.MassRotorIron+dataSet.MassMagnet));
set(app.RotorInertiaEditField,'Enable','on','Editable','off','Value',num2str(dataSet.RotorInertia));

% PM design panel
nameCols = cell(1,size(dataSet.PMdim,2));
wCols    = cell(1,size(dataSet.PMdim,2));
for ii=1:length(nameCols)
    nameCols{ii} = ['Layer ' int2str(ii)];
    wCols{ii} = 75;
end

nameRows{1} = 'Central width [mm]';
nameRows{2} = 'External width [mm]';
set(app.PMdimTable,'ColumnName',nameCols,'RowName',nameRows,'Data',round(dataSet.PMdim,4),'ColumnWidth',wCols);

nameCols = cell(1,size(dataSet.PMdim,2));
wCols    = cell(1,size(dataSet.PMdim,2));
for ii=1:length(nameCols)
    nameCols{ii} = ['Layer ' int2str(ii)];
    wCols{ii} = 75;
end
nameRows{1} = 'Central PM gap [mm]';
nameRows{2} = 'External PM gap [mm]';
set(app.PMclearanceTable,'ColumnName',nameCols,'RowName',nameRows,'Data',round(dataSet.PMclear,4),'ColumnWidth',wCols);

set(app.PMtemperatureEdit,'Value',num2str(dataSet.PMtemp));
set(app.CarCurEdit,'Value',num2str(dataSet.CurrPM,2));
set(app.BrPMEdit,'Value',num2str(dataSet.Br));
if strcmp(dataSet.FluxBarrierMaterial,'Air')
    set(app.PMdimTable,'Enable','off','ColumnEditable',false(1))
    set(app.PMclearanceTable,'Enable','off','ColumnEditable',false(1))
    set(app.PMDesignPush,'Enable','off')
    set(app.PMtemperatureEdit,'Enable','off')
    set(app.BrPMEdit,'Enable','off')
    set(app.CarCurEdit,'Enable','off')
else
    set(app.PMtemperatureEdit,'Enable','on')
    set(app.BrPMEdit,'Enable','on')
    if (strcmp(dataSet.TypeOfRotor,'SPM')||strcmp(dataSet.TypeOfRotor,'Fluid'))
        set(app.PMdimTable,'Enable','off','ColumnEditable',false(1))
        set(app.PMclearanceTable,'Enable','off','ColumnEditable',false(1))
        set(app.PMDesignPush,'Enable','off')
        set(app.CarCurEdit,'Enable','off')
    else
        set(app.PMdimTable,'Enable','on','ColumnEditable',true(1))
        set(app.PMclearanceTable,'Enable','on','ColumnEditable',true(1))
        set(app.PMDesignPush,'Enable','on')
        set(app.CarCurEdit,'Enable','on')
    end
end


% Optimization panel
set(app.MaxGenEdit,'Value',num2str(dataSet.MaxGen));
set(app.XPopEdit,'Value',num2str(dataSet.XPop));
set(app.SimPosMOEdit,'Value',num2str(dataSet.SimPoMOOA));
set(app.RotorPosiMOEdit,'Value',num2str(dataSet.RotPoMOOA));
set(app.SimPosFinerEdit,'Value',num2str(dataSet.SimPoFine));
set(app.RotorPosFinerEdit,'Value',num2str(dataSet.RotPoFine));
% opt variables and bounds
if strcmp(dataSet.TypeOfRotor,'SPM')
    set(app.Dalpha1BouCheck,'Enable','off','Value',0);
    set(app.Alpha1BouEdit,'Enable','off');
    set(app.DalphaBouCheck,'Enable','off','Value',0);
    set(app.DeltaAlphaBouEdit,'Enable','off','Value',mat2str(dataSet.DeltaAlphaBou));
else
    set(app.Dalpha1BouCheck,'Enable','on','Value',dataSet.Dalpha1BouCheck);
    set(app.DalphaBouCheck,'Enable','on','Value',dataSet.DalphaBouCheck);
    if dataSet.Dalpha1BouCheck
        set(app.Alpha1BouEdit,'Enable','on');
    else
        set(app.Alpha1BouEdit,'Enable','off');
    end
    if dataSet.DalphaBouCheck
        set(app.DeltaAlphaBouEdit,'Enable','on');
    else
        set(app.DeltaAlphaBouEdit,'Enable','off');
    end
end

set(app.hcBouCheck,'Enable','on','Value',dataSet.hcBouCheck);
if dataSet.hcBouCheck
    set(app.hcBouEdit,'Enable','on');
else
    set(app.hcBouEdit,'Enable','off');
end

if (strcmp(dataSet.TypeOfRotor,'ISeg')||strcmp(dataSet.TypeOfRotor,'SPM'))
    set(app.DxBouCheck,'Enable','off','Value',0);
    set(app.DfeBouEdit,'Enable','off');
else
    set(app.DxBouCheck,'Enable','on','Value',dataSet.DxBouCheck);
    if dataSet.DxBouCheck
        set(app.DfeBouEdit,'Enable','on');
    else
        set(app.DfeBouEdit,'Enable','off');
    end
end

set(app.GapBouCheck,'Enable','on','Value',dataSet.GapBouCheck);
if dataSet.GapBouCheck
    set(app.GapBouEdit,'Enable','on');
else
    set(app.GapBouEdit,'Enable','off');
end

if strcmp(dataSet.FluxBarrierMaterial,'Air')
    set(app.BrBouCheck,'Enable','off','Value',0);
    set(app.BrBouEdit,'Enable','off');
else
    set(app.BrBouCheck,'Enable','on','Value',dataSet.BrBouCheck);
    if dataSet.BrBouCheck
        set(app.BrBouEdit,'Enable','on');
    else
        set(app.BrBouEdit,'Enable','off');
    end
end

set(app.AirgapRadiusBouCheck,'Enable','on','Value',dataSet.AirgapRadiusBouCheck);
if dataSet.AirgapRadiusBouCheck
    set(app.AirgapRadiusBouEdit,'Enable','on');
else
    set(app.AirgapRadiusBouEdit,'Enable','off');
end

set(app.ToothWidthBouCheck,'Enable','on','Value',dataSet.ToothWidthBouCheck);
if dataSet.ToothWidthBouCheck
    set(app.ToothWidthBouEdit,'Enable','on');
else
    set(app.ToothWidthBouEdit,'Enable','off');
end

set(app.ToothLengthBouCheck,'Enable','on','Value',dataSet.ToothLengthBouCheck);
if dataSet.ToothLengthBouCheck
    set(app.ToothLenBouEdit,'Enable','on');
else
    set(app.ToothLenBouEdit,'Enable','off');
end

set(app.StatorSlotOpenBouCheck,'Enable','on','Value',dataSet.StatorSlotOpenBouCheck);
if dataSet.StatorSlotOpenBouCheck
    set(app.StatorSlotOpenBouEdit,'Enable','on');
else
    set(app.StatorSlotOpenBouEdit,'Enable','off');
end

set(app.ToothTangDepthBouCheck,'Enable','on','Value',dataSet.ToothTangDepthBouCheck);
if dataSet.ToothTangDepthBouCheck
    set(app.ToothTangDepthBouEdit,'Enable','on');
else
    set(app.ToothTangDepthBouEdit,'Enable','off');
end

if (strcmp(dataSet.TypeOfRotor,'SPM')||strcmp(dataSet.TypeOfRotor,'Vtype'))
    set(app.BetaPMshapeBouCheck,'Enable','on','Value',dataSet.BetaPMshapeBouCheck);
    if dataSet.BetaPMshapeBouCheck
        set(app.BetaPMshapeBouEdit,'Enable','on');
    else
        set(app.BetaPMshapeBouEdit,'Enable','off');
    end
else
    set(app.BetaPMshapeBouCheck,'Enable','off','Value',0);
    set(app.BetaPMshapeBouEdit,'Enable','off');
end

if (strcmp(dataSet.TypeOfRotor,'Circular')||strcmp(dataSet.TypeOfRotor,'Seg'))
    set(app.ThetaFBSBouCheck,'Enable','on','Value',dataSet.ThetaFBSBouCheck);
    if dataSet.ThetaFBSBouCheck
        set(app.ThetaFBSBouEdit,'Enable','on');
    else
        set(app.ThetaFBSBouEdit,'Enable','off');
    end
else
    set(app.ThetaFBSBouCheck,'Enable','off','Value',0);
    set(app.ThetaFBSBouEdit,'Enable','off');
end

if (strcmp(dataSet.FluxBarrierMaterial,'Air')||strcmp(dataSet.TypeOfRotor,'Fluid')||strcmp(dataSet.TypeOfRotor,'SPM'))
    set(app.PMdimBouCheck,'Enable','off','Value',0);
    set(app.PMdimBouEdit,'Enable','off');
else
    set(app.PMdimBouCheck,'Enable','on','Value',dataSet.PMdimBouCheck);
    if dataSet.PMdimBouCheck
        set(app.PMdimBouEdit,'Enable','on');
    else
        set(app.PMdimBouEdit,'Enable','off');
    end
end

set(app.GammaBouCheck,'Enable','on','Value',dataSet.GammaBouCheck);
if dataSet.GammaBouCheck
    set(app.PhaseAngleCurrBouEdit,'Enable','on');
else
    set(app.PhaseAngleCurrBouEdit,'Enable','off');
end

app.Alpha1BouEdit.Value         = mat2str(dataSet.Alpha1Bou);
app.DeltaAlphaBouEdit.Value     = mat2str(dataSet.DeltaAlphaBou);
app.hcBouEdit.Value             = mat2str(dataSet.hcBou);
app.DfeBouEdit.Value            = mat2str(dataSet.DfeBou);
app.GapBouEdit.Value            = mat2str(dataSet.GapBou);
app.BrBouEdit.Value             = mat2str(dataSet.BrBou);
app.AirgapRadiusBouEdit.Value   = mat2str(dataSet.GapRadiusBou);
app.ToothWidthBouEdit.Value     = mat2str(dataSet.ToothWiBou);
app.ToothLenBouEdit.Value       = mat2str(dataSet.ToothLeBou);
app.StatorSlotOpenBouEdit.Value = mat2str(dataSet.StatorSlotOpenBou);
app.ToothTangDepthBouEdit.Value = mat2str(dataSet.ToothTangDepthBou);
app.BetaPMshapeBouEdit.Value    = mat2str(dataSet.BetaPMshapeBou);
app.ThetaFBSBouEdit.Value       = mat2str(dataSet.ThetaFBSBou);
app.PMdimBouEdit.Value          = mat2str(dataSet.PMdimBou,2);
app.PhaseAngleCurrBouEdit.Value = mat2str(dataSet.PhaseAngleCurrBou);

% opt objectives
set(app.TorqueOptCheck,'Value',dataSet.TorqueOptCheck);
if dataSet.TorqueOptCheck
    set(app.MinExpTorEdit,'Enable','on');
else
    set(app.MinExpTorEdit,'Enable','off');
end
set(app.TorRipOptCheck,'Value',dataSet.TorRipOptCheck);
if dataSet.TorRipOptCheck
    set(app.MaxExpeRippleTorEdit,'Enable','on');
else
    set(app.MaxExpeRippleTorEdit,'Enable','off');
end
set(app.MassCuOptCheck,'Value',dataSet.MassCuOptCheck);
if dataSet.MassCuOptCheck
    set(app.MassCuEdit,'Enable','on');
else
    set(app.MassCuEdit,'Enable','off');
end
set(app.MassPMOptCheck,'Value',dataSet.MassPMOptCheck);
if dataSet.MassPMOptCheck
    set(app.MassPMEdit,'Enable','on');
else
    set(app.MassPMEdit,'Enable','off');
end

app.MinExpTorEdit.Value        = num2str(dataSet.MinExpTorque);
app.MaxExpeRippleTorEdit.Value = num2str(dataSet.MaxRippleTorque);
app.MassCuEdit.Value           = num2str(dataSet.MaxCuMass);
app.MassPMEdit.Value           = num2str(dataSet.MaxPMMass);

set(app.CurrentOverLoadEdit,'Value',num2str(dataSet.CurrOverLoad));

% Post-processing panel
set(app.SpanEltPPEdit,'Value',num2str(dataSet.AngularSpanPP));
set(app.GammaPPEdit,'Value',mat2str(dataSet.GammaPP));
set(app.CurrLoPPEdit,'Value',mat2str(dataSet.CurrLoPP));
set(app.CurrentPP,'Enable','off','Value',mat2str(dataSet.CurrLoPP*dataSet.RatedCurrent));
set(app.NumOfRotorPosiPPEdit,'Value',num2str(dataSet.NumOfRotPosPP));
set(app.BrPPEdit,'Value',mat2str(dataSet.BrPP));
set(app.TempPPEdit,'Value',mat2str(dataSet.tempPP));
set(app.NGridPPEdit,'Value',num2str(dataSet.NumGrid));
set(app.EvaluatedSpeedEdit,'Value',mat2str(dataSet.EvalSpeed));

switch dataSet.EvalType
    case 'singt'
        set(app.EvalTypePopUp,'Value','Single Point');
        set(app.SpanEltPPEdit,'Enable','on');
        set(app.GammaPPEdit,'Enable','on');
        set(app.CurrLoPPEdit,'Enable','on');
        set(app.NumOfRotorPosiPPEdit,'Enable','on');
        set(app.BrPPEdit,'Enable','on');
        set(app.TempPPEdit,'Enable','on');
        set(app.NGridPPEdit,'Enable','off');
        set(app.EvaluatedSpeedEdit,'Enable','on');
        set(app.StartPProMagnetPush,'Enable','on');
    case 'singm'
        set(app.EvalTypePopUp,'Value','Flux Map');
        set(app.SpanEltPPEdit,'Enable','on');
        set(app.GammaPPEdit,'Enable','off');
        set(app.CurrLoPPEdit,'Enable','on');
        set(app.NumOfRotorPosiPPEdit,'Enable','on');
        set(app.BrPPEdit,'Enable','on');
        set(app.TempPPEdit,'Enable','on');
        set(app.NGridPPEdit,'Enable','on');
        set(app.EvaluatedSpeedEdit,'Enable','on');
        set(app.StartPProMagnetPush,'Enable','on');
    case 'demagArea'
        set(app.EvalTypePopUp,'Value','Demagnetization Analysis');
        set(app.SpanEltPPEdit,'Enable','off');
        set(app.GammaPPEdit,'Enable','off');
        set(app.CurrLoPPEdit,'Enable','on');
        set(app.NumOfRotorPosiPPEdit,'Enable','off');
        set(app.BrPPEdit,'Enable','on');
        set(app.TempPPEdit,'Enable','on');
        set(app.NGridPPEdit,'Enable','off');
        set(app.EvaluatedSpeedEdit,'Enable','off');
        set(app.StartPProMagnetPush,'Enable','off');
    case 'idemag'
        set(app.EvalTypePopUp,'Value','Demagnetization Curve');
        set(app.SpanEltPPEdit,'Enable','off');
        set(app.GammaPPEdit,'Enable','off');
        set(app.CurrLoPPEdit,'Enable','off');
        set(app.NumOfRotorPosiPPEdit,'Enable','off');
        set(app.BrPPEdit,'Enable','on');
        set(app.TempPPEdit,'Enable','on');
        set(app.NGridPPEdit,'Enable','off');
        set(app.EvaluatedSpeedEdit,'Enable','off');
        set(app.StartPProMagnetPush,'Enable','off');
    case 'ichval'
        set(app.EvalTypePopUp,'Value','Characteristic Current');
        set(app.SpanEltPPEdit,'Enable','on');
        set(app.GammaPPEdit,'Enable','off');
        set(app.CurrLoPPEdit,'Enable','off');
        set(app.NumOfRotorPosiPPEdit,'Enable','on');
        set(app.BrPPEdit,'Enable','on');
        set(app.TempPPEdit,'Enable','on');
        set(app.NGridPPEdit,'Enable','off');
        set(app.EvaluatedSpeedEdit,'Enable','off');
        set(app.StartPProMagnetPush,'Enable','off');
    case 'flxdn'
        set(app.EvalTypePopUp,'Value','Flux Density Analysis');
        set(app.SpanEltPPEdit,'Enable','on');
        set(app.GammaPPEdit,'Enable','on');
        set(app.CurrLoPPEdit,'Enable','on');
        set(app.NumOfRotorPosiPPEdit,'Enable','on');
        set(app.BrPPEdit,'Enable','on');
        set(app.TempPPEdit,'Enable','on');
        set(app.NGridPPEdit,'Enable','off');
        set(app.EvaluatedSpeedEdit,'Enable','off');
        set(app.StartPProMagnetPush,'Enable','off');
    case 'izero'
        set(app.EvalTypePopUp,'Value','Current Offset');
        set(app.SpanEltPPEdit,'Enable','on');
        set(app.GammaPPEdit,'Enable','on');
        set(app.CurrLoPPEdit,'Enable','on');
        set(app.NumOfRotorPosiPPEdit,'Enable','on');
        set(app.BrPPEdit,'Enable','on');
        set(app.TempPPEdit,'Enable','on');
        set(app.NGridPPEdit,'Enable','off');
        set(app.EvaluatedSpeedEdit,'Enable','off');
        set(app.StartPProMagnetPush,'Enable','off');
    case 'force'
        set(app.EvalTypePopUp,'Value','Airgap Force');
        set(app.SpanEltPPEdit,'Enable','on');
        set(app.GammaPPEdit,'Enable','on');
        set(app.CurrLoPPEdit,'Enable','on');
        set(app.NumOfRotorPosiPPEdit,'Enable','on');
        set(app.BrPPEdit,'Enable','on');
        set(app.TempPPEdit,'Enable','on');
        set(app.NGridPPEdit,'Enable','off');
        set(app.EvaluatedSpeedEdit,'Enable','off');
        set(app.StartPProMagnetPush,'Enable','off');
    case 'singtIron'
        set(app.EvalTypePopUp,'Value','Iron Loss - Single Point');
        set(app.SpanEltPPEdit,'Enable','off');
        set(app.GammaPPEdit,'Enable','on');
        set(app.CurrLoPPEdit,'Enable','on');
        set(app.NumOfRotorPosiPPEdit,'Enable','on');
        set(app.BrPPEdit,'Enable','on');
        set(app.TempPPEdit,'Enable','on');
        set(app.NGridPPEdit,'Enable','off');
        set(app.EvaluatedSpeedEdit,'Enable','on');
        set(app.StartPProMagnetPush,'Enable','on');
    case 'singmIron'
        set(app.EvalTypePopUp,'Value','Iron Loss - Flux Map');
        set(app.SpanEltPPEdit,'Enable','off');
        set(app.GammaPPEdit,'Enable','off');
        set(app.CurrLoPPEdit,'Enable','on');
        set(app.NumOfRotorPosiPPEdit,'Enable','on');
        set(app.BrPPEdit,'Enable','on');
        set(app.TempPPEdit,'Enable','on');
        set(app.NGridPPEdit,'Enable','on');
        set(app.EvaluatedSpeedEdit,'Enable','on');
        set(app.StartPProMagnetPush,'Enable','on');
    case 'structural'
        set(app.EvalTypePopUp,'Value','Structural Analysis');
        set(app.SpanEltPPEdit,'Enable','off');
        set(app.GammaPPEdit,'Enable','off');
        set(app.CurrLoPPEdit,'Enable','off');
        set(app.NumOfRotorPosiPPEdit,'Enable','off');
        set(app.BrPPEdit,'Enable','off');
        set(app.TempPPEdit,'Enable','off');
        set(app.NGridPPEdit,'Enable','off');
        set(app.EvaluatedSpeedEdit,'Enable','on');
        set(app.StartPProMagnetPush,'Enable','off');
end

%% Motor-CAD
%Build Thermal
set(app.InlettemperatureEditField,'Value',num2str(dataSet.InletTemperature));
set(app.HousingTypeDropDown,'Value',dataSet.HousingType);
set(app.FluidDropDown,'Value',dataSet.Fluid);
set(app.FluidFlowRateEditField,'Value',num2str(dataSet.FlowRate));
% Thermal Sim
set(app.ThermalevaluationtypeDropDown,'Value',dataSet.th_eval_type);
set(app.TransientperiodEditField,'Value',num2str(dataSet.TransientPeriod));
set(app.TimestepEditField,'Value',num2str(dataSet.TransientTimeStep));
set(app.MachinetemperatureEditField,'Value',num2str(dataSet.MachineTemperature));
set(app.EvaluatedSpeedMirrorEdit,'Value',num2str(dataSet.EvalSpeed));

%set(app.ThermalLoadKj_mirror,'Value',num2str(dataSet.ThermalLoadKj));
%set(app.CurrLoPPMirrorEdit,'Value',mat2str(dataSet.CurrLoPP));
set(app.CurrentPPMirror,'Value',mat2str(round(dataSet.SimulatedCurrent,2)));
set(app.GammaPPMirrorEdit,'Value',mat2str(dataSet.GammaPP));
set(app.TempPPMirrorEdit,'Value',num2str(dataSet.tempPP));


