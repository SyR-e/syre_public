% Copyright 2022
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

function [Ns,windingData] = calc_hairpinTurns(setup,inputCheck)

if nargin==0
    prompt = {'Number of slots','Number of poles','Number of 3phase set','Number of conductors per slot','Number of parallel path'};
    answer = {'54','6','1','[2 4 6 8 10]','[1 2 3 4 5 6]'};
    answer = inputdlg(prompt,'Input',1,answer);
    Q     = eval(answer{1});
    P     = eval(answer{2});
    n3ph  = eval(answer{3});
    nCond = eval(answer{4});
    nPar  = eval(answer{5});
    flagPlot = 1;
elseif isfield(setup,'currentfilename')
    dataSet = setup;
    Q     = 6*dataSet.Num3PhaseCircuit*dataSet.NumOfPolePairs*dataSet.NumOfSlots;
    n3ph  = dataSet.Num3PhaseCircuit;
    P     = 2*dataSet.NumOfPolePairs;
    nCond = [2 4 6 8 10];
    nPar  = [1 2 3 4 5 6];
    flagPlot = 0;
else
    Q     = setup.Q;
    n3ph  = setup.n3ph;
    P     = setup.P;
    nCond = setup.nCond;
    nPar  = setup.nPar;
    flagPlot = 0;
end

if nargin()==2
    if inputCheck
        prompt = {'Number of slots','Number of poles','Number of 3phase set','Number of conductors per slot','Number of parallel path'};
        answer = {int2str(Q),int2str(P),int2str(n3ph),'[2 4 6 8 10]','[1 2 3 4 5 6]'};
        answer = inputdlg(prompt,'Input',1,answer);
        Q     = eval(answer{1});
        P     = eval(answer{2});
        n3ph  = eval(answer{3});
        nCond = eval(answer{4});
        nPar  = eval(answer{5});
        flagPlot = 1;
    end
end

p = P/2;
q = Q/(P*3*n3ph);

[nCond,nPar] = meshgrid(nCond,nPar);

Ns = nCond*Q./(2*nPar*3*n3ph);

feasibility = nan(size(nCond));
feasibility(nPar==1)              = 1;
feasibility(nPar==2)              = 1;
feasibility(nPar==nCond/2)        = 1;
feasibility(rem(nPar,p)==0)       = 1;
feasibility(rem(nPar,q)==0)       = 1;
feasibility(rem(nCond/2,nPar)==0) = 1;

Ns = Ns.*feasibility;

NsGrid = Ns;
Ns    = Ns(:);
nCond = nCond(:);
nPar  = nPar(:);

[Ns,index] = sort(Ns);
nCond = nCond(index);
nPar  = nPar(index);

Ns(rem(Ns,1)~=0) = NaN;

nCond = nCond(~isnan(Ns));
nPar  = nPar(~isnan(Ns));
Ns    = Ns(~isnan(Ns));

if flagPlot
    % figure
    % figSetting()
    % view(3)
    % xlabel('Number of conductor per slot');
    % ylabel('Number of parallel path');
    % zlabel('Number of turns in series per phase')
    % stem3(nCond,nPar,Ns,'-o','LineWidth',1.5)
    
    figure
    figSetting()
    xlabel('Number of conductors per slot');
    ylabel('Number of parallel paths');
    title('Feasible turns in series per phase $N_s$')
    set(gca,'Xlim',[min(nCond)-1 max(nCond)+1],'XTick',unique(nCond))
    set(gca,'Ylim',[0 max(nPar)+1],'YTick',unique(nPar))
    scatter(nCond, nPar, 1000, Ns, 'filled','Marker','s');
    colormap turbo
    dN = (max(Ns)-min(Ns))/5;
    clim([min(Ns)-dN max(Ns)+dN])
    for ii = 1:length(Ns)
        text(nCond(ii), nPar(ii), int2str(Ns(ii)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle');
    end

    clc
    disp(['Possible configuration for ' int2str(Q) ' slot, ' int2str(n3ph*3) ' phase stator:'])
    for ii=1:length(Ns)
        disp(['Ns = ' int2str(Ns(ii)) ' --> ' int2str(nCond(ii)) ' conductor/slot, ' int2str(nPar(ii)) ' parallel path'])
    end
end

windingData.Ns    = Ns;
windingData.nCond = nCond;
windingData.nPar  = nPar;
%windingData.setup = setup;

if nargout()==0
    clear Ns windingData;
end



