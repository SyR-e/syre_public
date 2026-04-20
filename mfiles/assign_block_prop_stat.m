% Copyright 2014
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

function [geo,FemmProblem] = assign_block_prop_stat(BLKLABELS,geo,mat,fem,group,FemmProblem)
BLKLABELSstat=BLKLABELS.statore;

if isempty(FemmProblem)
    flag_xfemm = 0;
else
    flag_xfemm = 1;
end

n3phase=geo.win.n3phase; %AS create circuits in FEMM

if(isnan(n3phase))
    phase_name{1}=strcat('fase',num2str(1));
    phase_name{2}=strcat('fase',num2str(2));
    phase_name{3}=strcat('fase',num2str(3));
    phase_name{4}=strcat('fase',num2str(4));
    phase_name{5}=strcat('fase',num2str(5));
    phase_name_neg{1}=strcat('fase',num2str(1),'n');
    phase_name_neg{2}=strcat('fase',num2str(2),'n');
    phase_name_neg{3}=strcat('fase',num2str(3),'n');
    phase_name_neg{4}=strcat('fase',num2str(4),'n');
    phase_name_neg{5}=strcat('fase',num2str(5),'n');
    
    if ~flag_xfemm
        mi_addcircprop(phase_name{1}, 0, 1);
        mi_addcircprop(phase_name{2}, 0, 1);
        mi_addcircprop(phase_name{3}, 0, 1);
        mi_addcircprop(phase_name{4}, 0, 1);
        mi_addcircprop(phase_name{5}, 0, 1);
        mi_addcircprop(phase_name_neg{1}, 0, 1);
        mi_addcircprop(phase_name_neg{2}, 0, 1);
        mi_addcircprop(phase_name_neg{3}, 0, 1);
        mi_addcircprop(phase_name_neg{4}, 0, 1);
        mi_addcircprop(phase_name_neg{5}, 0, 1);
    else
        FemmProblem = addcircuit_mfemm(FemmProblem,phase_name{1},'CircType',1);
        FemmProblem = addcircuit_mfemm(FemmProblem,phase_name{2},'CircType',1);
        FemmProblem = addcircuit_mfemm(FemmProblem,phase_name{3},'CircType',1);
        FemmProblem = addcircuit_mfemm(FemmProblem,phase_name{4},'CircType',1);
        FemmProblem = addcircuit_mfemm(FemmProblem,phase_name{5},'CircType',1);
        FemmProblem = addcircuit_mfemm(FemmProblem,phase_name_neg{1},'CircType',1);
        FemmProblem = addcircuit_mfemm(FemmProblem,phase_name_neg{2},'CircType',1);
        FemmProblem = addcircuit_mfemm(FemmProblem,phase_name_neg{3},'CircType',1);
        FemmProblem = addcircuit_mfemm(FemmProblem,phase_name_neg{4},'CircType',1);
        FemmProblem = addcircuit_mfemm(FemmProblem,phase_name_neg{5},'CircType',1);
    end
else
    for ii=0:(n3phase-1)
        phase_name{3*ii+1}=strcat('fase',num2str(3*ii+1));
        phase_name{3*ii+2}=strcat('fase',num2str(3*ii+2));
        phase_name{3*ii+3}=strcat('fase',num2str(3*ii+3));
        phase_name_neg{3*ii+1}=strcat('fase',num2str(3*ii+1),'n');
        phase_name_neg{3*ii+2}=strcat('fase',num2str(3*ii+2),'n');
        phase_name_neg{3*ii+3}=strcat('fase',num2str(3*ii+3),'n');
        
        if ~flag_xfemm
            mi_addcircprop(phase_name{3*ii+1}, 0, 1);
            mi_addcircprop(phase_name{3*ii+2}, 0, 1);
            mi_addcircprop(phase_name{3*ii+3}, 0, 1);
            mi_addcircprop(phase_name_neg{3*ii+1}, 0, 1);
            mi_addcircprop(phase_name_neg{3*ii+2}, 0, 1);
            mi_addcircprop(phase_name_neg{3*ii+3}, 0, 1);
        else
            FemmProblem = addcircuit_mfemm(FemmProblem,phase_name{3*ii+1},'CircType',1);
            FemmProblem = addcircuit_mfemm(FemmProblem,phase_name{3*ii+2},'CircType',1);
            FemmProblem = addcircuit_mfemm(FemmProblem,phase_name{3*ii+3},'CircType',1);
            FemmProblem = addcircuit_mfemm(FemmProblem,phase_name_neg{3*ii+1},'CircType',1);
            FemmProblem = addcircuit_mfemm(FemmProblem,phase_name_neg{3*ii+2},'CircType',1);
            FemmProblem = addcircuit_mfemm(FemmProblem,phase_name_neg{3*ii+3},'CircType',1);
        end
    end 
end

% avv=geo.defaultavv;  %AS
% avv=geo.win.avv;
ss=1; %slot index
ll=1; %layer index

% avv = windingCheck(geo);
avv = geo.win.avv;


for kk=1:length(BLKLABELSstat.xy(:,1))
    if BLKLABELSstat.xy(kk,3)==2 %Air
        if ~flag_xfemm
            mi_addblocklabel(BLKLABELSstat.xy(kk,1),BLKLABELSstat.xy(kk,2));
            mi_selectlabel(BLKLABELSstat.xy(kk,1),BLKLABELSstat.xy(kk,2));
            mi_setblockprop('Air', 0, fem.res,'None', 0, 10, 0);
            % mi_setblockprop(BLKLABELS.materials{BLKLABELSstat.xy(kk,3)}, 0, fem.res,'None', 0, 10, 0);
            %mi_setblockprop('Air', 0, fem.res,'None', 0, group, 0);
            mi_clearselected;
        else
            FemmProblem = addblocklabel_mfemm(FemmProblem,...
                BLKLABELSstat.xy(kk,1),BLKLABELSstat.xy(kk,2),...
                'BlockType','Air',...
                'MaxArea',fem.res,...
                'InGroup',10);
        end
    elseif BLKLABELSstat.xy(kk,3)==3 %Copper
        if avv(ll,ss) > 0
            fase = ['fase' num2str(abs(avv(ll,ss)))];
            dir = 1;
        else
            %fase = ['fase' num2str(abs(avv(ll,ss))) 'n'];
            fase = ['fase' num2str(abs(avv(ll,ss)))];
            dir = -1;
        end
        
        ll=ll+1;
        if ll>size(avv,1)
            ll=1;
            ss=ss+1;
        end
        
        if ~flag_xfemm
            mi_addblocklabel(BLKLABELSstat.xy(kk,1),BLKLABELSstat.xy(kk,2));
            mi_selectlabel(BLKLABELSstat.xy(kk,1),BLKLABELSstat.xy(kk,2));
            mi_setblockprop(mat.SlotCond.MatName, 0, fem.res, fase, 0, group, dir);
            %mi_setblockprop('Air', 0, fem.res,'None', 0, group, 0);
            mi_clearselected;
        else
            FemmProblem = addblocklabel_mfemm(FemmProblem,...
                BLKLABELSstat.xy(kk,1),BLKLABELSstat.xy(kk,2),...
                'BlockType',mat.SlotCond.MatName,...
                'MaxArea',fem.res,...
                'InGroup',group,...
                'InCircuit',fase,...
                'Turns',dir);
        end
    elseif BLKLABELSstat.xy(kk,3)==4 %Iron stator
        if ~flag_xfemm
            mi_addblocklabel(BLKLABELSstat.xy(kk,1),BLKLABELSstat.xy(kk,2));
            mi_selectlabel(BLKLABELSstat.xy(kk,1),BLKLABELSstat.xy(kk,2));
            mi_setblockprop(mat.Stator.MatName, 0, fem.res,'None', 0, 12, 0);
            %mi_setblockprop('Air', 0, fem.res,'None', 0, group, 0);
            mi_clearselected;
        else
            FemmProblem = addblocklabel_mfemm(FemmProblem,...
                BLKLABELSstat.xy(kk,1),BLKLABELSstat.xy(kk,2),...
                'BlockType',mat.Stator.MatName,...
                'MaxArea',fem.res,...
                'InGroup',12);
        end
    end
end
