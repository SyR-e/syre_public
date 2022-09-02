% Copyright 2022
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

function [out] = MMM_eval_thermalLimit_LossMap_sing(mcad,loss,TCu_lim,Tpm_lim,speed)

invoke(mcad,'ShowThermalContext');
MagnetTemperature = nan;

if speed==0
    speed = speed+1;
end
invoke(mcad,'SetVariable','ShaftSpeed',speed);

Pfe_s = loss(:,2);
Pfe_s = Pfe_s(isfinite(Pfe_s));

%%% stator loss split between yoke and teeth
[~,wt]  = invoke(mcad,'GetVariable','Tooth_Width');
[~,sd]  = invoke(mcad,'GetVariable','Slot_Depth');
[~,rse] = invoke(mcad,'GetVariable','Stator_Lam_Dia');
[~,rsi] = invoke(mcad,'GetVariable','Stator_Bore');
wy = (rse/2-rsi/2-sd);
ratio_pfe = wt/wy;

flag   = 1;
flagPM = 0;
kk = 1;
jj = 1;

doFig = 0;
if doFig
    hfig(1)=figure();
    figSetting
    hax(1) = axes('OuterPosition',[0 0 1 1]);
    xlabel('$Loss$ [W]')
    ylabel('$T_{Cu}$ [$^\circ$C]')
    % yaxis([0 TCu_lim+20])

    hfig(2)=figure();
    figSetting
    hax(2) = axes('OuterPosition',[0 0 1 1]);
    xlabel('$Loss$ [W]')
    ylabel('$T_{PM}$ [$^\circ$C]')
    % yaxis([0 Tpm_lim+20])
end

while flag==1
    if kk==1
        index = floor(length(Pfe_s)/4);
    elseif kk==2
        index = floor(length(Pfe_s)/4*3);
    elseif kk>2 && kk<30
        %lossTest(kk-1) = loss(index,end);
        if flagPM
            if jj==2
                 lossObj = 2/3* lossTest(kk-1); 
            else
                [a,b,c] = retta_per_2pti(lossTest(kk-2),Tpm(jj-2),lossTest(kk-1),Tpm(jj-1));
                [m,q] = retta_abc2mq(a,b,c);
                lossObj = (Tpm_lim-q)/m;
            end
        else
            [a,b,c] = retta_per_2pti(lossTest(kk-2),WindingTemperature_Max(kk-2),lossTest(kk-1),WindingTemperature_Max(kk-1));
            [m,q] = retta_abc2mq(a,b,c);
            lossObj = (TCu_lim-q)/m;
        end
        if lossObj > 0
            [~,index] = min(abs(lossObj-loss(:,end)));
            lossTest(kk) = loss(index,end);
        end
    elseif kk>30
        disp('Iteration failed');
        break
    end

    flag = 0;
    %%% insert the loss in MCAD
    Pfe_init  = loss(index,2);
    Pfet_init = ratio_pfe*Pfe_init;
    Pfey_init = (1-ratio_pfe)*Pfe_init;
    Pfer_init = loss(index,3);
    Ppm_init  = loss(index,4);
    Pjs_init  = loss(index,5);
    invoke(mcad,'SetVariable','Stator_Iron_Loss_@Ref_Speed_[Tooth]', Pfet_init);
    invoke(mcad,'SetVariable','Stator_Iron_Loss_@Ref_Speed_[Back_Iron]', Pfey_init);
    invoke(mcad,'SetVariable','Magnet_Iron_Loss_@Ref_Speed', Ppm_init);
    invoke(mcad,'SetVariable','Armature_Copper_Loss_@Ref_Speed', Pjs_init);
    invoke(mcad,'SetVariable','Rotor_Iron_Loss_@Ref_Speed_[Embedded_Magnet_Pole]', Pfer_init/2);
    invoke(mcad,'SetVariable','Rotor_Iron_Loss_@Ref_Speed_[Back_Iron]', Pfer_init/2);

    invoke(mcad,'DoSteadyStateAnalysis');

    %read the temperature
    [~, WindingTemperature_Max(kk)] = invoke(mcad,'GetVariable','T_[Winding_Max]');
    %[~, MagnetTemperature(kk)]      = invoke(mcad,'GetVariable','Magnet_Temperature');

    %chech the temperature
    if (WindingTemperature_Max(kk)<TCu_lim-5 || WindingTemperature_Max(kk)>TCu_lim+5) && flagPM==0
        flag = 1;
    end

    if flag == 0
        %%% it can be optimized if MotorCAD will add the magnet temperature
        %%% activeX parameter
        invoke(mcad,'SetVariable','MagneticThermalCoupling',2);
        invoke(mcad,'DoMagneticCalculation');
        [~, MagnetTemperature] = invoke(mcad,'GetVariable','Magnet_Temperature');
        if MagnetTemperature>Tpm_lim && flagPM==0
            flag = 1;
            flagPM = 1;
            Tpm(jj) = MagnetTemperature;
            jj = jj + 1; 
        elseif  (MagnetTemperature<Tpm_lim-5 || MagnetTemperature>Tpm_lim+5) && flagPM==1
            flag = 1;
            %flagPM = 1;
            Tpm(jj) = MagnetTemperature;
            jj = jj + 1; 
        end   
    end

    if index==1
        disp(['No continous point feasible at ' num2str(speed) ' rpm']);
        break
    end

    lossTest(kk) = loss(index,end);

    %plot
    if doFig
        plot(hax(1),loss(index,end),WindingTemperature_Max(kk),'o','DisplayName', [num2str(kk) ' iteration'])
        hleg(1) = legend(hax(1),'show');
        set(hleg(1),'Location','northwest');
%         plot(hax(2),loss(index,end),MagnetTemperature,'o','DisplayName', [num2str(kk) ' iteration'])
%         hleg(2) = legend(hax(2),'show');
%         set(hleg(2),'Location','northwest');
    end

    kk = kk + 1;
end

if flag==0
    [~, WJ_out]  = invoke(mcad,'GetVariable','WJ_Fluid_Outlet_Temp_[Active]');
    [~, WJ_in]   = invoke(mcad,'GetVariable','WJ_Fluid_Inlet_Temperature');
    deltaWJ      = WJ_out - WJ_in;
    Loss_Tot     = loss(index,end);
else
    Loss_Tot                   = nan;
    deltaWJ                    = nan;
    WindingTemperature_Max(kk) = nan;
    MagnetTemperature          = nan;
end

%% Output
out.Ploss   = Loss_Tot;
out.deltaWJ = deltaWJ;
out.TempCu  = WindingTemperature_Max(end);
out.TempPM  = MagnetTemperature;
