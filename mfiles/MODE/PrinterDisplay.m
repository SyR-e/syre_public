% Copyright 2014
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
function [OUT,Dat]=PrinterDisplay(OUT,Dat)

disp('------------------------------------------------')
disp(['Generation: ' num2str(Dat.CounterGEN)]);
disp(['FEs: ' num2str(Dat.CounterFES)]);
disp(['Pareto Front Size: ' mat2str(size(OUT.PFront,1))]);
disp(['Evaluation Time: ' int2str(Dat.EvalTime(4)) ' h ' num2str(Dat.EvalTime(5)) ' min ' num2str(round(Dat.EvalTime(6))) ' sec']);
disp(['Actual time             : ' datestr(now())]);
disp(['End of evolution process: ' datestr(Dat.EndTime)]);
disp('------------------------------------------------')

if mod(Dat.CounterGEN,1)==0
    if Dat.NOBJ==3
        stem3(OUT.PFront(:,1),OUT.PFront(:,2),OUT.PFront(:,3),'*r');
        grid on; hold on; drawnow
    elseif Dat.NOBJ==2
        figSetting
        plot(OUT.PFront(:,1),OUT.PFront(:,2),'*r');
        grid on; hold on; drawnow
    elseif Dat.NOBJ==1
        figSetting
        plot(Dat.CounterGEN,log(min(OUT.PFront(:,1))),'*r');
        grid on; hold on; drawnow
        %     else
        %         figSetting
        %         plot(OUT.PFront(:,1),OUT.PFront(:,2),'*r');
        %         grid on; hold on; drawnow
        %         title('Pareto front on the first two objs only')
    end
end