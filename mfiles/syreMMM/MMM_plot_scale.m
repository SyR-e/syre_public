% Copyright 2020
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

function MMM_plot_scale(hax,motorModel,ls,name)

Id = motorModel.FluxMap_dq.Id;
Iq = motorModel.FluxMap_dq.Iq;
Fd = motorModel.FluxMap_dq.Fd;
Fq = motorModel.FluxMap_dq.Fq;

% flusso d
idPlot = unique(Id(:));
iqPlot = min(abs(Iq(:)),[],'all')*ones(size(idPlot));
fdPlot = interp2(Id,Iq,Fd,idPlot,iqPlot);
plot(hax,idPlot,fdPlot,ls,'Linewidth',1.5,'DisplayName',[name ' - D'])
% flusso q
iqPlot = unique(Iq(:));
idPlot = min(abs(Id(:)))*ones(size(iqPlot));
fqPlot = interp2(Id,Iq,Fq,idPlot,iqPlot);
plot(hax,iqPlot,fqPlot,ls,'Linewidth',1.5,'DisplayName',[name ' - Q'])

legend(hax,'off')
legend(hax,'show','Location','best')

ylim(hax,'auto')
xlim(hax,'auto')



