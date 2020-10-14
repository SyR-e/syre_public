% Copyright 2020
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

function MMM_GUI_SetParameters(app)

motorModel = app.motorModel;
data = motorModel.data;
Tw   = motorModel.Tw;
motorModelUnScale = app.motorModelUnScale;
motorModelUnSkew  = app.motorModelUnSkew;



% main data
set(app.CurrentPathEditField,...
    'Value',motorModel.data.pathname,...
    'Enable','on',...
    'Editable','off');
% set(app.ModelSavedCheckBox,...
%     'Enable','off',...
%     'Value',motorModel.data.flagSave);

% motor data
set(app.MotornameEditField,'Value',data.motorName)
set(app.Numberof3phasesetsEditField,...
    'Value',int2str(data.n3phase),...
    'Enable','on',...
    'Editable','off');
set(app.MotortypeEditField,'Value',data.motorType);
set(app.AxistypeEditField,'Value',data.axisType);
set(app.RatedpowerEditField,'Value',num2str(data.P0));
set(app.RatedcurrentEditField,'Value',num2str(data.i0));
set(app.MaximumcurrentEditField,'Value',num2str(data.Imax));
set(app.DClinkvoltageEditField,'Value',num2str(data.Vdc));
set(app.RatedspeedEditField,'Value',num2str(data.n0));
set(app.MaximumspeedEditField,'Value',num2str(data.nmax));
set(app.PhaseresistanceEditField,'Value',num2str(data.Rs));
set(app.WindingtemperatureEditField,'Value',num2str(data.tempCu));
% set(app.PMtemperatureEditField,'Value',num2str(data.tempPM),...
%     'Enable','on',...
%     'Editable','off');

tmpCell = cell(length(data.tempVectPM)+1,1);
for ii=1:length(data.tempVectPM)
    tmpCell{ii} = int2str(data.tempVectPM(ii));
end
tmpCell{end} = 'Add';
tempPM = data.tempPM;
if isempty(tempPM)
    str = 'Add';
else
    str = int2str(tempPM);
end
set(app.PMtemperatureDropDown,'Items',tmpCell);
set(app.PMtemperatureDropDown,'Value',str);
set(app.ActivelengthEditField,'Value',num2str(data.l),...
    'Enable','on',...
    'Editable','off');
set(app.TurnsinseriesperphaseEditField,'Value',num2str(data.Ns),...
    'Enable','on',...
    'Editable','off');
set(app.InertiaEditField,'Value',num2str(data.J));

% Main
% model loaded
if ~isempty(motorModel.fdfq)
    set(app.dqModelCheckBox,...
        'Enable','on',...
        'Value',1);
    set(app.PlotDQmodelButton,'Enable','on')
    set(app.SaveDQmodelButton,'Enable','on')
    set(app.PrintDQmodelButton,'Enable','on')
    set(app.EvaluateMTPAButton,'Enable','on')
    set(app.EvalInductanceButton,'Enable','on')
    set(app.EvalInverseDQButton,'Enable','on')
    set(app.PlotTgammaButton,'Enable','on')
    set(app.MaxTwPush,'Enable','on')
    set(app.EvaluateShortCircuitTorqueButton,'Enable','on')
else
    set(app.dqModelCheckBox,...
        'Enable','off',...
        'Value',0);
    set(app.PlotDQmodelButton,'Enable','off')
    set(app.SaveDQmodelButton,'Enable','off')
    set(app.PrintDQmodelButton,'Enable','off')
    set(app.EvaluateMTPAButton,'Enable','off')
    set(app.EvalInductanceButton,'Enable','off')
    set(app.EvalInverseDQButton,'Enable','off')
    set(app.PlotTgammaButton,'Enable','off')
    set(app.MaxTwPush,'Enable','off')
    set(app.EvaluateShortCircuitTorqueButton,'Enable','off')
end
if ~isempty(motorModel.dqtMap)
    set(app.dqtMapModelCheckBox,...
        'Enable','on',...
        'Value',1);
    set(app.PlotDQTmodelButton,'Enable','on')
    set(app.SaveDQTmodelButton,'Enable','on')
    set(app.EvalInverseDQTButton,'Enable','on')
else
    set(app.dqtMapModelCheckBox,...
        'Enable','off',...
        'Value',0);
    set(app.PlotDQTmodelButton,'Enable','off')
    set(app.SaveDQTmodelButton,'Enable','off')
    set(app.EvalInverseDQTButton,'Enable','off')
end
if ~isempty(motorModel.ironLoss)
    set(app.IronLossModelCheckBox,...
        'Enable','on',...
        'Value',1);
    set(app.PlotIronLossModelButton,'Enable','on')
    set(app.SaveIronLossModelButton,'Enable','on')
else
    set(app.IronLossModelCheckBox,...
        'Enable','off',...
        'Value',0);
    set(app.PlotIronLossModelButton,'Enable','off')
    set(app.SaveIronLossModelButton,'Enable','off')
end
if ~isempty(motorModel.skinEffect)
    set(app.SkinEffectModelCheckBox,...
        'Enable','on',...
        'Value',1);
    set(app.PlotSkinEffectButton,'Enable','on')
    set(app.SaveSkinEffectButton,'Enable','on')
else
    set(app.SkinEffectModelCheckBox,...
        'Enable','off',...
        'Value',0);
    set(app.PlotSkinEffectButton,'Enable','off')
    set(app.SaveSkinEffectButton,'Enable','off')
end
% AOA (Admitted Operating Area) MTPA/MTPV
if ~isempty(motorModel.AOA)
    set(app.AOACheckBox,...
        'Enable','on',...
        'Value',1);
    set(app.PlotMTPAButton,'Enable','on')
    set(app.SaveMTPAButton,'Enable','on')
    set(app.PrintMTPAButton,'Enable','on')
    set(app.MTPAmethodDropDown,'Value',motorModel.AOA.method)
    set(app.OpLimPush,'Enable','on')
    set(app.RatingsEvalPush,'Enable','on')
else
    set(app.AOACheckBox,...
        'Enable','off',...
        'Value',0);
    set(app.PlotMTPAButton,'Enable','off')
    set(app.SaveMTPAButton,'Enable','off')
    set(app.PrintMTPAButton,'Enable','off')
    set(app.MTPAmethodDropDown,'Value','LUT')
    set(app.OpLimPush,'Enable','off')
    set(app.RatingsEvalPush,'Enable','off')
end
% inductance map
if ~isempty(motorModel.Inductance)
    set(app.InductanceMapsCheckBox,...
        'Enable','on',...
        'Value',1);
    set(app.PlotInductanceButton,'Enable','on')
    set(app.SaveInductanceButton,'Enable','on')
else
    set(app.InductanceMapsCheckBox,...
        'Enable','off',...
        'Value',0);
    set(app.PlotInductanceButton,'Enable','off')
    set(app.SaveInductanceButton,'Enable','off')
end
% Current angle curves
set(app.NumCurrLevelTgammaEditField,'Value',mat2str(data.nCurr));
% inverse models
if ~isempty(motorModel.idiq)
    set(app.InversedqCheckBox,...
        'Enable','on',...
        'Value',1);
    set(app.PlotInverseDQButton,'Enable','on')
    set(app.SaveInverseDQButton,'Enable','on')
else
    set(app.InversedqCheckBox,...
        'Enable','off',...
        'Value',0);
    set(app.PlotInverseDQButton,'Enable','off')
    set(app.SaveInverseDQButton,'Enable','off')
end
if ~isempty(motorModel.dqtMapF)
    set(app.InversedqtMapCheckBox,...
        'Enable','on',...
        'Value',1);
    set(app.PlotInverseDQTButton,'Enable','on')
    set(app.SaveInverseDQTButton,'Enable','on')
else
    set(app.InversedqtMapCheckBox,...
        'Enable','off',...
        'Value',0);
    set(app.PlotInverseDQTButton,'Enable','off')
    set(app.SaveInverseDQTButton,'Enable','off')
end

% dqtMap model
if ~isempty(motorModel.dqtMap)
    set(app.DQTmodelCheckBox,...
        'Enable','on',...
        'Value',1);
    set(app.DQTplotPush,'Enable','on')
    set(app.DQTsavePush,'Enable','on')
    set(app.DQTharmButton,'Enable','on')
    set(app.DQTsingtButton,'Enable','on')
else
    set(app.DQTmodelCheckBox,...
        'Enable','off',...
        'Value',0);
    set(app.DQTplotPush,'Enable','off')
    set(app.DQTsavePush,'Enable','off')
    set(app.DQTharmButton,'Enable','off')
    set(app.DQTsingtButton,'Enable','off')
end
set(app.DQTharmEditField,'Value',mat2str(motorModel.dqtElab.harmonic))

set(app.DQTgammaEditField,'Value',mat2str(motorModel.dqtElab.CurrAngle,3))
set(app.DQTcurrentPUEditField,'Value',mat2str(motorModel.dqtElab.CurrLoad))
set(app.DQTcurrentEditField,...
    'Enable','on',...
    'Editable','off',...
    'Value',mat2str(motorModel.dqtElab.CurrAmpl))

% Scale & Skew
% Scale
set(app.NewTurnsEditField,'Value',int2str(motorModel.scale.Ns));
set(app.NewLengthEditField,'Value',num2str(motorModel.scale.l));
set(app.NewLldEditField,'Value',num2str(motorModel.scale.Lld));
set(app.NewLlqEditField,'Value',num2str(motorModel.scale.Llq));

% Skew
set(app.SkewAngleEditField,'Value',num2str(motorModel.skew.thSkw));
set(app.SkewSliceEditField,'Value',int2str(motorModel.skew.nSlice));
set(app.SkewPointsEditField,'Value',int2str(motorModel.skew.nPoints));
set(app.dqtMapskewingevaluationEditField,...
    'Value','0%',...
    'Enable','on',...
    'Editable','off');

if ~isempty(motorModelUnScale)
    set(app.NewTurnsEditField,'Enable','on');
    set(app.NewLengthEditField,'Enable','on');
    set(app.NewLldEditField,'Enable','on');
    set(app.NewLlqEditField,'Enable','on');
    set(app.ScaleModelPush,'Enable','on');
    set(app.UnscaleModelPush,'Enable','on');
    set(app.SkewAngleEditField,'Enable','off');
    set(app.SkewSliceEditField,'Enable','off');
    set(app.SkewPointsEditField,'Enable','off');
    set(app.SkewModelPush,'Enable','off');
    set(app.UnskewModelPush,'Enable','off');
elseif ~isempty(motorModelUnSkew)
    set(app.NewTurnsEditField,'Enable','off');
    set(app.NewLengthEditField,'Enable','off');
    set(app.NewLldEditField,'Enable','off');
    set(app.NewLlqEditField,'Enable','off');
    set(app.ScaleModelPush,'Enable','off');
    set(app.UnscaleModelPush,'Enable','off');
    set(app.SkewAngleEditField,'Enable','on');
    set(app.SkewSliceEditField,'Enable','on');
    set(app.SkewPointsEditField,'Enable','on');
    set(app.SkewModelPush,'Enable','on');
    set(app.UnskewModelPush,'Enable','on');
else
    set(app.NewTurnsEditField,'Enable','on');
    set(app.NewLengthEditField,'Enable','on');
    set(app.NewLldEditField,'Enable','on');
    set(app.NewLlqEditField,'Enable','on');
    set(app.ScaleModelPush,'Enable','off');
    set(app.UnscaleModelPush,'Enable','off');
    set(app.SkewAngleEditField,'Enable','on');
    set(app.SkewSliceEditField,'Enable','on');
    set(app.SkewPointsEditField,'Enable','on');
    set(app.SkewModelPush,'Enable','off');
    set(app.UnskewModelPush,'Enable','off');
end

if (~isempty(motorModelUnScale)||~isempty(motorModelUnSkew))
    set(app.SavePush,'Enable','off')
    set(app.SaveDQmodelButton,'Enable','off')
    set(app.SaveDQTmodelButton,'Enable','off')
    set(app.SaveIronLossModelButton,'Enable','off')
    set(app.SaveSkinEffectButton,'Enable','off')
    set(app.SaveMTPAButton,'Enable','off')
    set(app.SaveInductanceButton,'Enable','off')
    set(app.SaveInverseDQButton,'Enable','off')
    set(app.SaveInverseDQTButton,'Enable','off')
    set(app.PrintDQmodelButton,'Enable','off')
    set(app.PrintMTPAButton,'Enable','off')
    set(app.PlotDQTmodelButton,'Enable','off')
    set(app.PlotInverseDQTButton,'Enable','off')
    set(app.DQTevalPush,'Enable','off')
    set(app.DQTplotPush,'Enable','off')
    set(app.DQTsavePush,'Enable','off')
    set(app.DQTsingtButton,'Enable','off')
    set(app.OpLimPush,'Enable','off')
    set(app.RatingsEvalPush,'Enable','off')
    set(app.MaxTwPush,'Enable','off')
    set(app.PMtemperatureDropDown,'Enable','off');
else
    set(app.SavePush,'Enable','on')
    set(app.DQTevalPush,'Enable','on')
    set(app.PMtemperatureDropDown,'Enable','on');
end


% Torque-speed
set(app.OpLimLevelsEditField,'Value',mat2str(data.nCurr))

set(app.TwMinSpeedEditField,'Value',num2str(Tw.nmin))
set(app.TwMaxSpeedEditField,'Value',num2str(Tw.nmax))
set(app.TwNumSpeedEditField,'Value',int2str(Tw.nstep))
set(app.TwMinTorqueEditField,'Value',num2str(Tw.Tmin))
set(app.TwMaxTorqueEditField,'Value',num2str(Tw.Tmax))
set(app.TwNumTorqueEditField,'Value',int2str(Tw.Tstep))
set(app.TwTemperatureEditField,'Value',num2str(Tw.temperature))
set(app.TwMechLossEditField,'Value',mat2str(Tw.MechLoss))
set(app.TwIronLossDropDown,'Value',Tw.IronLossFlag)
set(app.TwIronLossFactorEditField,'Value',num2str(Tw.IronLossFactor))
set(app.TwPMlossDropDown,'Value',Tw.PMLossFlag)
set(app.TwPMlossFactorEditField,'Value',num2str(Tw.PMLossFactor))
set(app.TwSkinEffectDropDown,'Value',Tw.SkinEffectFlag)
set(app.TwSkinEffectMethodDropDown,'Value',Tw.SkinEffectMethod)
set(app.TwControlDropDown,'Value',Tw.Control)
set(app.TwAxis,...
    'XLim',[Tw.nmin Tw.nmax],...
    'YLim',[Tw.Tmin Tw.Tmax])

% syreDrive
set(app.ControltypeDropDown,'Value',motorModel.SyreDrive.Ctrl_type);
set(app.FluxmapsmodelDropDown,'Value',motorModel.SyreDrive.FMapsModel);
set(app.ONthreasholdEditField,'Value',num2str(motorModel.SyreDrive.Converter.V0));
set(app.IncrementalresistanceEditField,'Value',num2str(motorModel.SyreDrive.Converter.Rd));
set(app.DeadtimeEditField,'Value',num2str(motorModel.SyreDrive.Converter.dT));


