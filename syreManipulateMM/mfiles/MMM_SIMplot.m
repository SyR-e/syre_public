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

function MMM_SIMplot(motorModel)
    Ctrl_type = motorModel.SyreDrive.Ctrl_type;
    
    Time            = motorModel.SyreDrive.simOut.time;
    Torque_ref      = motorModel.SyreDrive.simOut.torque(:,1);
    Torque          = motorModel.SyreDrive.simOut.torque(:,2);
    Torque_obs      = motorModel.SyreDrive.simOut.torque(:,3);
    Speed_ref       = motorModel.SyreDrive.simOut.speed(:,1);
    Speed           = motorModel.SyreDrive.simOut.speed(:,3);
    id_ref          = motorModel.SyreDrive.simOut.id(:,1);
    id              = motorModel.SyreDrive.simOut.iq(:,2);
    iq_ref          = motorModel.SyreDrive.simOut.iq(:,1);
    iq              = motorModel.SyreDrive.simOut.iq(:,2);
    lambdad_obs     = motorModel.SyreDrive.simOut.lambdad(:,1);
    lambdad         = motorModel.SyreDrive.simOut.lambdad(:,2);
    lambdaq_obs     = motorModel.SyreDrive.simOut.lambdaq(:,1);
    lambdaq         = motorModel.SyreDrive.simOut.lambdaq(:,2);
    
    % TORQUE    
    figure
    if strcmp(Ctrl_type,'Torque control')
        plot(Time,Torque_ref,'LineWidth',1);
        hold on
    end
    plot(Time,Torque,Time,Torque_obs,'LineWidth',1);
    xlabel('Time [s]')
    ylabel('Torque [Nm]')
    if strcmp(Ctrl_type,'Torque control')
        legend('Reference torque','Torque','Torque from flux observer','Location','best')
    else
        legend('Torque','Torque from flux observer','Location','best')
    end
    title('Torque')
    grid on
    
    % SPEED
    figure
    if strcmp(Ctrl_type,'Speed control')
        plot(Time,Speed_ref,'LineWidth',1)
        hold on
        legend('Reference speed')
    end
    plot(Time,Speed,'LineWidth',1);
    xlabel('Time [s]')
    ylabel('Speed [rpm]')
    if strcmp(Ctrl_type,'Speed control')
        legend('Reference speed','Speed','Location','best')
    else
        legend('Speed','Location','best')
    end
    title('Speed')
    grid on

    % dq CURRENTS
    figure
    subplot(2,1,1)
    plot(Time,id_ref,Time,id,'LineWidth',1)
    xlabel('Time [s]')
    ylabel('d Current [A]')
    legend('Reference d Current','d Current')
    grid on
    title('dq Currents')

    
    subplot(2,1,2)
    plot(Time,iq_ref,Time,iq,'LineWidth',1)
    xlabel('Time [s]')
    ylabel('q Current [A]')
    legend('Reference q Current','q Current')
    grid on
    
   % dq FLUXES
    figure
    subplot(2,1,1)
    plot(Time,lambdad_obs,Time,lambdad,'LineWidth',1)
    xlabel('Time [s]')
    ylabel('d Flux [A]')
    legend('Observed d Flux','d Flux')
    grid on
    title('dq Fluxes')

    
    subplot(2,1,2)
    plot(Time,lambdaq_obs,Time,lambdaq,'LineWidth',1)
    xlabel('Time [s]')
    ylabel('q Flux [A]')
    legend('Observed q Flux','q Flux')
    grid on
end