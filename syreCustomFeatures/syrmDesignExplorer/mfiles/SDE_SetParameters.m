% Copyright 2021
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

function SDE_SetParameters(app)

map = app.map;

if ~isempty(map)
    cla(app.UIAxes);
    contour(app.UIAxes,map.xx,map.bb,map.T,'-r','LineWidth',1,'DisplayName','T [Nm]','ShowText','on');
    contour(app.UIAxes,map.xx,map.bb,map.PF,'-b','LineWidth',1,'DisplayName','cos\phi','ShowText','on');
    plot(app.UIAxes,map.xRaw,map.bRaw,'Color',[0 0.5 0],'LineStyle','none','Marker','o','MarkerFaceColor',[0 0.5 0],'DisplayName','FEAfix')
    plot(app.UIAxes,map.xx(isnan(map.T)),map.bb(isnan(map.T)),'rx','DisplayName','unfeasible','MarkerSize',8,'LineWidth',1.5)
    legend(app.UIAxes,'show','Location','northeast');
    plot(app.UIAxes,map.xSelect,map.bSelect,'ko','MarkerFaceColor','k','DisplayName','Selected Motor');
    
    set(app.UIAxes,'XLim',[min(map.xx,[],'all') max(map.xx,[],'all')],'YLim',[min(map.bb,[],'all') max(map.bb,[],'all')])

    set(app.AvailabledataListBox,'Items',map.dataAvailable)

    set(app.xEditField,'Value',num2str(map.xSelect,4),'Editable','on');
    set(app.bEditField,'Value',num2str(map.bSelect,4),'Editable','on');
    set(app.VdcEditField,'Value',int2str(map.Vdc),'Editable','on');
    set(app.NsEditField,'Value',int2str(map.geo.win.Ns),'Editable','on');
    set(app.StacklengthEditField,'Value',int2str(map.geo.l),'Editable','on');
    set(app.PlaneNameEditField,'Value',app.filename,'Editable','off');

    set(app.DrawButton,'Enable','on');
    set(app.SaveButton,'Enable','on');
    set(app.ExporttoWorkSpaceButton,'Enable','on');
    set(app.EvaluatefeasibilityareaButton,'Enable','on');

    if ~isempty(map.dataSelect)
        set(app.Plot2DButton,'Enable','on')
        set(app.Plot3DButton,'Enable','on')
        set(app.AvailabledataListBox,'Value',map.dataSelect)
    else
        set(app.Plot2DButton,'Enable','off')
        set(app.Plot3DButton,'Enable','off')
        set(app.AvailabledataListBox,'Value',{})
    end
else
    cla(app.UIAxes);
    
    set(app.PlaneNameEditField,'Value','Plane not loaded','Editable','off');

    set(app.xEditField,'Value',num2str(0,4),'Editable','off');
    set(app.bEditField,'Value',num2str(0,4),'Editable','off');
    set(app.VdcEditField,'Value',int2str(0),'Editable','off');
    set(app.NsEditField,'Value',int2str(0),'Editable','off');
    set(app.StacklengthEditField,'Value',int2str(0),'Editable','off');

    set(app.DrawButton,'Enable','off');
    set(app.SaveButton,'Enable','off');
    set(app.ExporttoWorkSpaceButton,'Enable','off');
    set(app.EvaluatefeasibilityareaButton,'Enable','off');
    set(app.Plot2DButton,'Enable','off')
    set(app.Plot3DButton,'Enable','off')
end