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

function BarCenter = defineBlockCenters(temp,fem,geo)

% coordinates of center points of all FEMM blocks,
% where block labels will be eventually placed

res=fem.res;
p = geo.p;

% Material codes: refer to BLKLABELS.materials positions
materialCodes;  % load material codes

xc   = temp.xc;
yc   = temp.yc(~isnan(xc));
xmag = temp.xmag(~isnan(xc(1,:)));
ymag = temp.ymag(~isnan(xc(1,:)));
zmag = temp.zmag(~isnan(xc(1,:)));
mirrorFlag  = temp.mirrorFlag(~isnan(xc));
mirrorFlagAir = temp.mirrorFlagAir;
xc   = xc(~isnan(xc));

if isfield(temp,'xair')
    xair = temp.xair;
    yair = temp.yair(~isnan(xair));
    mirrorFlagAir = mirrorFlagAir(~isnan(xair));
    xair = xair(~isnan(xair));
    
else
    xair = [];
    yair = [];
end

if ~isempty(xair)
    xc   = [xc(mirrorFlag==1)';   xc(mirrorFlag==1)';    xc(mirrorFlag==0)'];
    yc   = [yc(mirrorFlag==1)';   -yc(mirrorFlag==1)';   yc(mirrorFlag==0)'];
%     xair = [xair(mirrorFlag1==1)'; xair(mirrorFlag1==1)';  xair(mirrorFlag1==0)'];
%     yair = [yair(mirrorFlag1==1)'; -yair(mirrorFlag1==1)'; -yair(mirrorFlag1==0)'];
    xair = [xair(mirrorFlagAir==1)'; xair(mirrorFlagAir==1)'; xair(mirrorFlagAir==0)'];
    yair = [yair(mirrorFlagAir==1)'; -yair(mirrorFlagAir==1)'; yair(mirrorFlagAir==0)'];
    xmag = [xmag(mirrorFlag==1)'; xmag(mirrorFlag==1)';  xmag(mirrorFlag==0)'];
    ymag = [ymag(mirrorFlag==1)'; -ymag(mirrorFlag==1)'; -ymag(mirrorFlag==0)'];
    zmag = [zmag(mirrorFlag==1)'; zmag(mirrorFlag==1)';  zmag(mirrorFlag==0)'];
    if strcmp(geo.RotType, 'EESM') || strcmp(geo.RotType,'Hybrid')
        zmag = [1; -1]; % for EESM, it is the current direction in the rotor
        BarCenter = [
            xc,yc,codMatCuRot*ones(length(xc),1),res*ones(length(xc),1),1*ones(length(xc),1),xmag,ymag,zmag;
            xair,yair,codMatAirRot*ones(length(xair),1),res*ones(length(xair),1),1*ones(length(xair),1),zeros(length(xair),1),zeros(length(xair),1),zeros(length(xair),1)
            ];
    else
        BarCenter = [
            xc,yc,codMatBar*ones(length(xc),1),res*ones(length(xc),1),1*ones(length(xc),1),xmag,ymag,zmag;
            xair,yair,codMatAirRot*ones(length(xair),1),res*ones(length(xair),1),1*ones(length(xair),1),zeros(length(xair),1),zeros(length(xair),1),zeros(length(xair),1)
            ];
    end




        %===============================================================================
    if strcmp(geo.RotType, 'Hybrid') && isfield(temp, 'hasPM') && temp.hasPM

        % Calcola il centro del magnete - anche per EESM xc,yc sono
        % calcolati con il mean

        x_pm_center = mean(temp.xpm);
        y_pm_center = 0;

        % Direzione magnete: (1,0) = Radiale, lungo l'asse X
        MagDir_X = 1;
        MagDir_Y = 0;


        codMag = codMatBar; % Valore di fallback (spesso usato per PM nei test)
        GroupPM = 1;


        % Crea la riga per la nuova etichetta in teoria xc e xair etc..
        % sono uguali a 1 nel nostro design,però usiamo la stessa logica
        % delle bobine
        MagLabel = [x_pm_center, y_pm_center, codMatBar*ones(length(x_pm_center),1),res*ones(length(x_pm_center),1),GroupPM, MagDir_X, MagDir_Y, NaN];
        


        x_g_center = mean([temp.xpm(1) temp.xg(1)]);
        magnet_air_gap_Label= [x_g_center, y_pm_center, codMatAirRot*ones(length(x_g_center),1),res*ones(length(x_g_center),1),zeros(length(x_g_center),1),zeros(length(x_g_center),1),zeros(length(x_g_center),1),1];
        
        %Ideal Barrier
        x_ib_center = mean(temp.xib);
        y_ib_center = mean(temp.yib); 
        % disp(codMatIdealBarrier)
        ideal_barrier_Label_top = [x_ib_center, y_ib_center, codMatIdealBarrier*ones(length(x_ib_center),1), res*ones(length(x_ib_center),1), zeros(length(x_ib_center),1), zeros(length(x_ib_center),1), zeros(length(x_ib_center),1), 1];
        ideal_barrier_Label_bottom = [x_ib_center, -y_ib_center, codMatIdealBarrier*ones(length(x_ib_center),1), res*ones(length(x_ib_center),1), zeros(length(x_ib_center),1), zeros(length(x_ib_center),1), zeros(length(x_ib_center),1), 1];
        
        BarCenter = [BarCenter; MagLabel; magnet_air_gap_Label; ideal_barrier_Label_top; ideal_barrier_Label_bottom];
    end
    %==============================================================================

else
    xc   = [xc';xc'];
    yc   = [yc';-yc'];
    xmag = [xmag';xmag'];
    ymag = [ymag';-ymag'];
    zmag = [zmag';zmag'];
    BarCenter = [xc,yc,codMatBar*ones(length(xc),1),res*ones(length(xc),1),1*ones(length(xc),1),xmag,ymag,zmag];
end
