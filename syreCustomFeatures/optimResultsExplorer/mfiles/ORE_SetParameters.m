% Copyright 2026
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

function ORE_SetParameters(app)

optRes = app.optRes;
OREsetup = app.OREsetup;

if isempty(optRes)
    set(app.DatanotloadedLabel,'BackgroundColor',[1 0 0],'Text','Data not loaded');
    set(app.xaxisDropDown,'Items',{'none'},'Value',OREsetup.xaxis)
    set(app.yaxisDropDown,'Items',{'none'},'Value',OREsetup.yaxis)
    set(app.zaxisDropDown,'Items',{'none'},'Value',OREsetup.zaxis)
    set(app.caxisDropDown,'Items',{'none'},'Value',OREsetup.caxis)
    set(app.RemoveunfeasibleButton,'Value',OREsetup.removeUnfeasible)
    set(app.Plot2DButton,'Enable','off')
    set(app.Plot3DButton,'Enable','off')
    set(app.CorrelationplotButton,'Enable','off')
    set(app.DrawmotorButton,'Enable','off')
    set(app.VariablesListBox,'Items',{});
    set(app.ObjectivesListBox,'Items',{});
else
    set(app.DatanotloadedLabel,'BackgroundColor',[0 0.8 0],'Text','Data loaded');
    vartot = horzcat('none','motorID',optRes.varNames,optRes.objNames);
    set(app.xaxisDropDown,'Items',vartot,'Value',OREsetup.xaxis)
    set(app.yaxisDropDown,'Items',vartot,'Value',OREsetup.yaxis)
    set(app.zaxisDropDown,'Items',vartot,'Value',OREsetup.zaxis)
    set(app.caxisDropDown,'Items',vartot,'Value',OREsetup.caxis)
    set(app.RemoveunfeasibleButton,'Value',OREsetup.removeUnfeasible)
    set(app.Plot2DButton,'Enable','on')
    set(app.Plot3DButton,'Enable','on')
    set(app.CorrelationplotButton,'Enable','on')
    set(app.DrawmotorButton,'Enable','on')
    set(app.VariablesListBox,'Items',optRes.varNames);
    set(app.ObjectivesListBox,'Items',optRes.objNames);

    ORE_plot(optRes,OREsetup,app.UIAxes);
end