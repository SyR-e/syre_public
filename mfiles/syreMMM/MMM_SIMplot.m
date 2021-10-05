

function MMM_SIMplot(motorModel)
    Ctrl_type = motorModel.SyreDrive.Ctrl_type;
    
    Time                = motorModel.SyreDrive.simOut.Outputs.id.Time;
    Torque_ref          = motorModel.SyreDrive.simOut.Outputs.T_ext.Data;
    Torque              = motorModel.SyreDrive.simOut.Out_M.T_m.Data;
    Torque_obs          = motorModel.SyreDrive.simOut.Outputs.T_elt.Data;
    Torque_load         = motorModel.SyreDrive.simOut.Out_M.T_load.Data;
    Speed               = motorModel.SyreDrive.simOut.Out_M.n_m.Data;
    Speed_obs           = motorModel.SyreDrive.simOut.Outputs.wr.Data;
    Speed_ref           = motorModel.SyreDrive.simOut.Outputs.wr_ref.Data;
    pos_err             = motorModel.SyreDrive.simOut.Outputs.pos_err_real.Data;
    pos_err_signal      = motorModel.SyreDrive.simOut.Outputs.pos_err.Data;
    id_ref              = motorModel.SyreDrive.simOut.Outputs.id_ref.Data;
    id                  = motorModel.SyreDrive.simOut.Outputs.id.Data;
    iq_ref              = motorModel.SyreDrive.simOut.Outputs.iq_ref.Data;
    iq                  = motorModel.SyreDrive.simOut.Outputs.iq.Data;
    lambdad_obs         = motorModel.SyreDrive.simOut.Outputs.lambdaD.Data;
    lambdad_CM          = motorModel.SyreDrive.simOut.Outputs.lambdaD_CM.Data;
    lambdaq_obs         = motorModel.SyreDrive.simOut.Outputs.lambdaQ.Data;
    lambdaq_CM          = motorModel.SyreDrive.simOut.Outputs.lambdaQ_CM.Data;
    lambdad             = motorModel.SyreDrive.simOut.Out_M.lambda_dq.Data(:,1);
    lambdaq             = motorModel.SyreDrive.simOut.Out_M.lambda_dq.Data(:,2);
    f_w                 = motorModel.SyreDrive.simOut.Outputs.f_w.Data;
    iA                  = motorModel.SyreDrive.simOut.Outputs.iA.Data;
    iB                  = motorModel.SyreDrive.simOut.Outputs.iB.Data;
    iC                  = motorModel.SyreDrive.simOut.Outputs.iC.Data;
    theta_elt           = motorModel.SyreDrive.simOut.Outputs.theta_elt_meas.Data;
    
    %----------------Mechanical Quantities------------------------------
    
    figure
    figSetting
    subplot(4,1,1)
    if strcmp(Ctrl_type,'Torque control')
        plot(Time,Torque,Time,Torque_obs,Time,Torque_ref);
        xlabel('t (s)')
        ylabel('T (Nm)')
        legend('$T$','$\hat{T}$','$T^*$','Location','best')
        ylim([1.05*min(min(min(Torque),min(Torque_obs)),min(Torque_ref)) 1.05*max(max(max(Torque),max(Torque_obs)),max(Torque_ref))])
        xlim([0 Time(end)])
        title('Torque')
        grid on
    end
    if strcmp(Ctrl_type,'Speed control')
        plot(Time,Torque,Time,Torque_obs,Time,Torque_load);
        xlabel('t (s)')
        ylabel('T (Nm.)')
        ylim([-inf 1.05*max(max(max(Torque),max(Torque_obs)),max(Torque_load))])
        xlim([0 Time(end)])
        legend('$T$','$\hat{T}$','$T_L$','Location','best')
        title('Torque')
        grid on
    end
    if strcmp(Ctrl_type,'Current control')
        plot(Time,Torque,Time,Torque_obs);
        xlabel('t (s)')
        ylabel('T (Nm.)')
        ylim([-inf 1.05*max(max(max(Torque),max(Torque_obs)),max(Torque_load))])
        xlim([0 Time(end)])
        legend('$T$','$\hat{T}$','Location','best')
        title('Torque')
        grid on
    end
    
    % SPEED
    subplot(4,1,2)
    if strcmp(Ctrl_type,'Speed control')
        plot(Time,Speed_ref,Time,Speed,Time,Speed_obs)
        xlabel('t (s)')
        ylabel('$\omega$ (rpm)')
        legend('$\omega_r^*$','$\omega_r$','$\hat{\omega}_r$','Location','best')
        ylim([-inf 1.05*max(max(Speed_ref),max(Speed))])
        xlim([0 Time(end)])
        title('Speed')
        grid on
    end
    if strcmp(Ctrl_type,'Torque control')
        plot(Time,Speed,Time,Speed_obs)
        xlabel('t (s)')
        ylabel('$\omega$ (rpm)')
        ylim([1.05*min(min(Speed_ref),min(Speed)) 1.05*max(max(Speed_ref),max(Speed))])
        xlim([0 Time(end)])
        legend('$\omega_r$','$\hat{\omega}_r$','Location','best')
        title('Speed')
        grid on
    end
    if strcmp(Ctrl_type,'Current control')
        plot(Time,Speed,Time,Speed_obs)
        xlabel('t (s)')
        ylabel('$\omega$ (rpm)')
        ylim([1.05*min(min(Speed_ref),min(Speed)) 1.05*max(max(Speed_ref),max(Speed))])
        xlim([0 Time(end)])
        legend('$\omega_r$','$\hat{\omega}_r$','Location','best')
        title('Speed')
        grid on
    end
    
    % POSITION
    subplot(4,1,3)
    plot(Time, pos_err_signal, Time, pos_err)
    xlabel('Time (s)')
    ylabel('$\tilde{\theta}$ (deg)')
    ylim([-20 20]);
    xlim([0 Time(end)])
    legend('$\epsilon$','$\tilde{\theta}$','Location','best')
    grid on
    title('Position error')
    
    subplot(4,1,4);
    plot(Time, f_w)
    xlabel('t (s)')
    ylabel('$f_{\omega}$')
    ylim([-0.2 1.2])
    xlim([0 Time(end)])
    grid on
    title('Fusion coefficient')
    
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 8], 'PaperUnits', 'Inches');

    
    %-------------------------Electrical Quantities-----------------------
    
    figure
    figSetting
    subplot(4,1,1)
    plot(Time,id_ref,Time,id)
    xlabel('t (s)')
    ylabel('$i_d$ (A)')
    xlim([0 Time(end)])
    legend('$i_d^*$','$i_d$','Location','best')
    grid on
    title('d Currents')

    subplot(4,1,2)
    plot(Time,iq_ref,Time,iq)
    xlabel('t (s)')
    ylabel('$i_q$ (A)')
    legend('$i_q^*$','$i_q$','Location','best')
    xlim([0 Time(end)])
    grid on
    title('q Currents')

   % dq FLUXES
    subplot(4,1,3)
    plot(Time,lambdad,Time,lambdad_obs,Time,lambdad_CM)
    xlabel('t (s)')
    ylabel('$\lambda_d$ (Vs)')
    xlim([0 Time(end)])
    legend('$\lambda_d$','$\hat{\lambda}_d$','$\lambda^i_d$','Location','best')
    grid on
    title('d Fluxes')

    
    subplot(4,1,4)
    plot(Time,lambdaq,Time,lambdaq_obs,Time,lambdaq_CM)
    xlabel('t (s)')
    ylabel('$\lambda_q$ (Vs)')
    xlim([0 Time(end)])
    legend('$\lambda_q$','$\hat{\lambda}_q$','$\lambda^i_q$','Location','best')
    grid on
    title('q Fluxes')
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 8], 'PaperUnits', 'Inches');
    
    
    %----------------------------------Electrical Cycle---------------------
    
    t1 = Time(end)-60/Speed(end)*1/motorModel.data.p;
    
    figure
    figSetting
    subplot(4,1,1)
    plot(Time,iA,Time,iB,Time,iC)
    xlabel('t (s)')
    ylabel('$i_{dq}$ (A)')
    xlim([t1 Time(end)])
    grid on
    title('abc Currents')
    
    subplot(4,1,2)
    yyaxis left;
    plot(Time,id)
    xlim([t1 Time(end)]);
    ylim([min(id(Time>t1))-1 max(id(Time>t1))+2]);
    ylabel('$i_d$ (A)')
    yyaxis right;
    plot(Time,iq);
    ylabel('$i_q$ (A)')
    xlabel('t (s)')
    xlim([t1 Time(end)]);
    ylim([min(iq(Time>t1))-2 max(iq(Time>t1))+1]);
    grid on
    title('dq Currents')
    
    
    subplot(4,1,3)
    yyaxis left;
    plot(Time,lambdad)
    ylabel('$\lambda_d$ (Vs)')
    xlim([t1 Time(end)]);
    ylim([min(lambdad(Time>t1))-0.1 max(lambdad(Time>t1))+0.2]);
    yyaxis right;
    plot(Time,lambdaq);
    ylabel('$\lambda_q$ (Vs)')
    xlim([t1 Time(end)]);
    ylim([min(lambdaq(Time>t1))-0.2 max(lambdaq(Time>t1))+0.1]);
    xlabel('t (s)')
    grid on
    title('dq Fluxes')

    subplot(4,1,4)
    yyaxis left;
    plot(Time,theta_elt)
    xlabel('t (s)')
    ylabel('$\theta$ (deg)')
    xlim([t1 Time(end)])
    grid on
    yyaxis right;
    plot(Time,pos_err)
    xlabel('t (s)')
    ylabel('$\tilde{\theta}$ (deg)')
    xlim([t1 Time(end)])
    grid on
    title('Position (electrical)')
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 5, 8], 'PaperUnits', 'Inches');
    
end