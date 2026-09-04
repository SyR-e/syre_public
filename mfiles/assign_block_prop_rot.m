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

function [geo,FemmProblem] = assign_block_prop_rot(BLKLABELS,geo,mat,fem,group,FemmProblem)

BLKLABELSrot=BLKLABELS.rotore;
Br = repmat(mat.LayerMag.Br,1,geo.ps);

if isempty(FemmProblem)
    flag_xfemm = 0;
else
    flag_xfemm = 1;
end

if ~flag_xfemm
    % pulisce le selezioni precedenti
    mi_clearselected
end
% Q = geo.ns*geo.p;                    % number of slots
Q = 6*geo.q*geo.p;                    % number of slots

% Assegna aria alle barriere di flux:
if ((max(geo.betaPMshape)~=0)&&strcmp('Circular',geo.RotType))
    tmp = [mat.LayerMag.Br mat.LayerMag.Br];
    Br = repmat(tmp,1,geo.ps);
elseif strcmp(geo.RotType,'SPM')||strcmp(geo.RotType,'Spoke-type')||strcmp(geo.RotType,'SPM-Halbach')|| strcmp(geo.RotType,'Hybrid')
    Br = mat.LayerMag.Br(1)*ones(size(BLKLABELSrot.xy(:,1)));
end

if strcmp(geo.RotType,'EESM') || strcmp(geo.RotType,'Hybrid')
    if ~flag_xfemm
        mi_addcircprop('field',0,1);
    else
        FemmProblem = addcircuit_mfemm(FemmProblem,'field','CircType',1);
    end
end

kk = 1; % tiene conto di quale magnete sto assegnando
bb = 1; % tiene conto di quale barra di rotore sto assegnando
for ii=1:length(BLKLABELSrot.xy(:,1))
    switch BLKLABELSrot.xy(ii,3)
        case 1  % Aria
            % if strcmp(geo.RotType,'EESM')||strcmp(geo.RotType,'Hybrid')
            %     groupAir = 23;
            % else
                groupAir = group;
            % end
            if ~flag_xfemm
                mi_addblocklabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                mi_selectlabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(ii,3)}, 0, fem.res,'None', 0, groupAir, 0);
                mi_setblockprop('Air', 0, fem.res,'None', 0, groupAir, 0);
                mi_clearselected;
            else
                FemmProblem = addblocklabel_mfemm(FemmProblem,...
                    BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2),...
                    'BlockType','Air',...
                    'MaxArea',fem.res,...
                    'InGroup',groupAir);
            end
        %Definiamo Ideal Label
        case 10  % ideal flux barrier
            if ~flag_xfemm
                mi_addblocklabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                mi_selectlabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                mi_setblockprop(BLKLABELS.materials{BLKLABELSrot.xy(ii,3)}, 0, fem.res,'None', 0, 500, 0);
                mi_setblockprop('IdealBarrier', 0, fem.res,'None', 0, 500, 0);
                mi_clearselected;
            else
                FemmProblem = addblocklabel_mfemm(FemmProblem,...
                    BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2),...
                    'BlockType','IdealBarrier',...
                    'MaxArea',fem.res,...
                    'InGroup',500);
            end
        case 6 % PM
            if isfield(mat.LayerMag,'BH')
                magdir=atan2(BLKLABELSrot.xy(ii,2),BLKLABELSrot.xy(ii,1))*180/pi;
                if ~flag_xfemm
                    mi_addblocklabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                    mi_selectlabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                    mi_setblockprop(mat.LayerMag.MatName, 0, fem.res,'None', magdir, 200, 0);
                    

                    if strcmp(geo.RotType,'Hybrid')
                        mi_setblockprop(mat.LayerMag.MatName, 0, fem.res,'None', magdir, 401, 0);
                    end

                    mi_clearselected;

                else
                    if strcmp(geo.RotType,'Hybrid')
                        %modifica alla stringa mag name per farla
                        %funzionare in simulated xdeg
                        % indMat = length(FemmProblem.Materials)+1;
                        % FemmProblem.Materials(indMat) = newmaterial_mfemm([mat.LayerMag.MatName '_' num2str(kk)]);
                        % FemmProblem.Materials(indMat).Mu_x  = mat.LayerMag.mu;
                        % FemmProblem.Materials(indMat).Mu_y  = mat.LayerMag.mu;
                        % FemmProblem.Materials(indMat).H_c   = Hc;
                        % FemmProblem.Materials(indMat).Sigma = mat.LayerMag.sigmaPM/1e6;
                        % FemmProblem = addblocklabel_mfemm(FemmProblem,...
                        %     BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2),...
                        %     'BlockType',[mat.LayerMag.MatName '_' num2str(kk)],...
                        %     'MaxArea',fem.res/geo.mesh_kpm,...
                        %     'InGroup',400+kk,...
                        %     'MagDir',magdir);
                        FemmProblem = addblocklabel_mfemm(FemmProblem,...
                            BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2),...
                            'BlockType',mat.LayerMag.MatName,...
                            'MaxArea',fem.res,...
                            'InGroup',400,...
                            'MagDir',magdir);
                    else
                    FemmProblem = addblocklabel_mfemm(FemmProblem,...
                        BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2),...
                        'BlockType',mat.LayerMag.MatName,...
                        'MaxArea',fem.res,...
                        'InGroup',200,...
                        'MagDir',magdir);
                    end
                end
            else

                Hc=1/(4e-7*pi*mat.LayerMag.mu)*Br(kk);
                magdir=atan2(BLKLABELSrot.xy(ii,7),BLKLABELSrot.xy(ii,6))*180/pi;

                if ~flag_xfemm
                    mi_addmaterial([mat.LayerMag.MatName '_' num2str(kk)], mat.LayerMag.mu, mat.LayerMag.mu, Hc, 0, mat.LayerMag.sigmaPM/1e6);
                    mi_addblocklabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                    mi_selectlabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                    if strcmp(geo.RotType,'Hybrid')
                        mi_setblockprop([mat.LayerMag.MatName '_' num2str(kk)], 0, fem.res/geo.mesh_kpm,'None', magdir, 400+kk, 0);
                    else
                        mi_setblockprop([mat.LayerMag.MatName '_' num2str(kk)], 0, fem.res/geo.mesh_kpm,'None', magdir, 200+kk, 0);
                    end
                    
                    mi_clearselected;
                else
                if strcmp(geo.RotType,'Hybrid')
                    indMat = length(FemmProblem.Materials)+1;
                    FemmProblem.Materials(indMat) = newmaterial_mfemm([mat.LayerMag.MatName '_' num2str(kk)]);
                    FemmProblem.Materials(indMat).Mu_x  = mat.LayerMag.mu;
                    FemmProblem.Materials(indMat).Mu_y  = mat.LayerMag.mu;
                    FemmProblem.Materials(indMat).H_c   = Hc;
                    FemmProblem.Materials(indMat).Sigma = mat.LayerMag.sigmaPM/1e6;
                    FemmProblem = addblocklabel_mfemm(FemmProblem,...
                        BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2),...
                        'BlockType',[mat.LayerMag.MatName '_' num2str(kk)],...
                        'MaxArea',fem.res/geo.mesh_kpm,...
                        'InGroup',400+kk,...
                        'MagDir',magdir);
                else
                    indMat = length(FemmProblem.Materials)+1;
                    FemmProblem.Materials(indMat) = newmaterial_mfemm([mat.LayerMag.MatName '_' num2str(kk)]);
                    FemmProblem.Materials(indMat).Mu_x  = mat.LayerMag.mu;
                    FemmProblem.Materials(indMat).Mu_y  = mat.LayerMag.mu;
                    FemmProblem.Materials(indMat).H_c   = Hc;
                    FemmProblem.Materials(indMat).Sigma = mat.LayerMag.sigmaPM/1e6;
                    FemmProblem = addblocklabel_mfemm(FemmProblem,...
                        BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2),...
                        'BlockType',[mat.LayerMag.MatName '_' num2str(kk)],...
                        'MaxArea',fem.res/geo.mesh_kpm,...
                        'InGroup',200+kk,...
                        'MagDir',magdir);
                end

                end
                kk=kk+1;
            end
        case 5 % Ferro rotore
            if ~flag_xfemm
                mi_addblocklabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                mi_selectlabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                mi_setblockprop(mat.Rotor.MatName, 0, fem.res,'None', 0, 22, 0);
                mi_clearselected;
            else
                FemmProblem = addblocklabel_mfemm(FemmProblem,...
                    BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2),...
                    'BlockType',mat.Rotor.MatName,...
                    'MaxArea',fem.res,...
                    'InGroup',22);
            end
        case 7 % shaft
            if ~flag_xfemm
                mi_addblocklabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                mi_selectlabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                if isequal(mat.Shaft.MatName,'ShaftAir')
                    mi_setblockprop('Air', 0, fem.res,'None', 0, group, 0);
                else
                    mi_setblockprop(mat.Shaft.MatName, 0, fem.res,'None', 0, group, 0);
                end
                mi_clearselected;
            else
                if isequal(mat.Shaft.MatName,'Air')
                    mat.Shaft.MatName = 'Air';
                end
                FemmProblem = addblocklabel_mfemm(FemmProblem,...
                    BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2),...
                    'BlockType',mat.Shaft.MatName,...
                    'MaxArea',fem.res,...
                    'InGroup',group);
            end
        case 8 % rotor bar
            if strcmp(geo.RotType,'IM')
                barName = ['bar' int2str(bb)];
                if ~flag_xfemm
                    mi_addcircprop(barName,0,1);
                    mi_addblocklabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                    mi_selectlabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                    mi_setblockprop(mat.BarCond.MatName,0,fem.res,barName,0,200+bb,1);
                    mi_clearselected;
                else
                    FemmProblem = addcircuit_mfemm(FemmProblem,...
                        barName,...
                        'CircType',1);
                    FemmProblem = addblocklabel_mfemm(FemmProblem,...
                        BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2),...
                        'BlockType',mat.BarCond.MatName,...
                        'MaxArea',fem.res,...
                        'InGroup',200+bb,...
                        'InCircuit',barName,...
                        'Turns',1);
                end
                bb = bb+1;
            else
                barName = 'field';
                if ~flag_xfemm
                    mi_addblocklabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                    mi_selectlabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                    mi_setblockprop(mat.BarCond.MatName,0,fem.res,barName,0,200+bb,BLKLABELSrot.xy(ii,end));
                    mi_clearselected;
                else
                    FemmProblem = addblocklabel_mfemm(FemmProblem,...
                        BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2),...
                        'BlockType',mat.BarCond.MatName,...
                        'MaxArea',fem.res,...
                        'InGroup',200+bb,...
                        'InCircuit',barName,...
                        'Turns',BLKLABELSrot.xy(ii,end));
                end
                bb = bb+1;
            end
        case 9 % sleeve
            if ~flag_xfemm
                mi_addblocklabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                mi_selectlabel(BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2));
                mi_setblockprop(mat.Sleeve.MatName,0,fem.res_traf,'None',0,199,0);
                mi_clearselected;
            else
                FemmProblem = addblocklabel_mfemm(FemmProblem,...
                    BLKLABELSrot.xy(ii,1),BLKLABELSrot.xy(ii,2),...
                        'BlockType',mat.Sleeve.MatName,...
                        'MaxArea',fem.res_traf,...
                        'InGroup',199);
            end
    end
end

