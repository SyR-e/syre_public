%% VERSIONE 20 04 2018
%% end debug PostProcessing

Doc         = invoke(h.magnetHandler, 'getDocument');
View        = invoke(Doc, 'getCurrentView');
Solution    = invoke(Doc, 'getSolution');
calculator  = invoke(h.magnetHandler, 'getStackCalculator');

% read simulation time values
NTime       = invoke(Solution,'getFieldSolutionTimeInstants',1);  % numero istanti

Command     = 'TimeInstants = getDocument().getSolution().getFieldSolutionTimeInstants(1,time_instants)';
invoke(h.magnetHandler, 'processCommand', Command);

time = [];
for ij = 1:NTime
    invoke(h.magnetHandler, 'processCommand', ['Call setVariant(0, time_instants(' num2str(ij-1) '))']);
    time(ij) = invoke(h.magnetHandler, 'getVariant', 0);
end
%%

Nbodies = invoke(Doc,'getNumberOfBodies',1);
invoke(h.magnetHandler, 'processCommand', 'REDIM ArrayOfValues(0)');

% IPM(1) o SMPM (0)
% IPM = isfield(Mac,'tipo_strati');

%% 21 nov 2011

% if exist('n_mag_simulati')
%     n_mag = n_mag_simulati;
% else
%     if IPM
%         n_mag = size(Mac.magneti,1)/3;
%     else %if motor is an SPM model
%         n_mag = size(Mac.magneti,1)/4;
%     end
% end

Flux_1sim = [];
Torque_1sim = [];
pm_loss = zeros(n_mag_simulati,NTime);

Ncoils = invoke(Doc,'getNumberOfCoils');

%% MATTEO solo per S1 ed S4 settare Ncoils %% Solo per il caso S1 S4 Matteo RICORDATELO

% Ncoils=3;
Porz_stat=geo.Qs;
Nbarrier=0;
%%
rotor_barrier_loss=zeros(NTime,Nbarrier);

%%
for ii=1:NTime
    % flux, currents, resitance
    Flux_1sim_1ph = [];
    Torque_r = [];
    
    % problem ID
    Command = 'ReDim ProblemID(1)';
    invoke(h.magnetHandler, 'processCommand', Command);
    Command = 'ProblemID(0) = 1';
    invoke(h.magnetHandler, 'processCommand', Command);
    Command = ['ProblemID(1) = ' num2str(time(ii) * 1.001)];
    invoke(h.magnetHandler, 'processCommand', Command);
    
    for jj=1:Ncoils
        Command = ['Call getDocument.getSolution.getFluxLinkageThroughCoil(ProblemID, ' num2str(jj) ', magnitude, phase)'];
        invoke(h.magnetHandler, 'processCommand', Command);
        invoke(h.magnetHandler, 'processCommand', 'Call setVariant(0, magnitude)');
        x = invoke(h.magnetHandler, 'getVariant', 0);
        Flux_1sim_1ph = [Flux_1sim_1ph x];
    end
    
    % torque
    for jjj=1:Nbodies
        Command = ['Call getDocument.getSolution.getTorqueAboutOriginOnBody(ProblemID, ' num2str(jjj) ', torque_x, torque_y, torque_z)'];
        invoke(h.magnetHandler, 'processCommand', Command);
        invoke(h.magnetHandler, 'processCommand', 'Call setVariant(0, torque_z)');
        torque = invoke(h.magnetHandler, 'getVariant', 0);
        Torque_r = [Torque_r torque];
    end
    
    %     % co-energy
    %     Command = ['Coenergy = getDocument.getSolution.getCoenergy(ProblemID)'];
    %     invoke(h.magnetHandler, 'processCommand', Command);
    %     invoke(h.magnetHandler, 'processCommand', 'Call setVariant(0, Coenergy)');
    %     CoEn(ii) = invoke(h.magnetHandler, 'getVariant', 0);
    %
    
    %% load iron and pm losses (only for 360 deg simulation)
    if (ii == 1) && (xdeg == 360)
        
        % stator: hysteresis and classical
        %         h_loss_s = zeros(1,Mac.Qs);
        %         c_loss_s = zeros(1,Mac.Qs);
        
        for jjj = 1:geo.Qs
            %             for jjj = 1:Porz_stat
            Command = ['IronLoss = getDocument.getSolution.getIronLossInComponent (ProblemID,"statore_' num2str(jjj) '" , Losses)'];
            invoke(h.magnetHandler, 'processCommand', Command);
            invoke(h.magnetHandler, 'processCommand', 'Call setVariant(0, Losses)');
            test0 = invoke(h.magnetHandler, 'getVariant', 0);
            test0 = cell2mat(test0);
            h_loss_s(jjj) = test0(1);
            c_loss_s(jjj) = test0(2);
        end
        
        % rotor
        Command = ['IronLoss = getDocument.getSolution.getIronLossInComponent (ProblemID,"rotor" , Losses)'];
        invoke(h.magnetHandler, 'processCommand', Command);
        invoke(h.magnetHandler, 'processCommand', 'Call setVariant(0, Losses)');
        test0 = invoke(h.magnetHandler, 'getVariant', 0);
        test0 = cell2mat(test0);
        h_loss_r = test0(1);
        c_loss_r = test0(2);
    end
    
    %% PM loss
    for jjj = 1:n_mag_simulati
        Command = ['PowerLoss=getDocument().getSolution().getOhmicLossInConductor(ProblemID,"Magnet_Bar_' num2str(jjj) '")'];
        invoke(h.magnetHandler, 'processCommand', Command);
        invoke(h.magnetHandler, 'processCommand', 'Call setVariant(0, PowerLoss)');
        pm_loss(jjj,ii) = invoke(h.magnetHandler, 'getVariant', 0);
    end
    Flux_1sim = [Flux_1sim;Flux_1sim_1ph];
    Torque_1sim = [Torque_1sim;Torque_r];
    ii;
    %     keyboard
    %% Rotor conductor Losses:
    if (Nbarrier==0)
        jjj=1;
        rotor_barrier_loss(ii,jjj) = 0;
    else
        %         keyboard
        jk=1;
        for jjj=1:Nbarrier
            for jj=1:geo.ps
                Command = ['PowerLoss=getDocument().getSolution().getOhmicLossInConductor(ProblemID,"barrier#' num2str(jjj),'_',num2str(jj) '")'];
                invoke(h.magnetHandler, 'processCommand', Command);
                invoke(h.magnetHandler, 'processCommand', 'Call setVariant(0, PowerLoss)');
                rotor_barrier_loss(ii,jk) = invoke(h.magnetHandler, 'getVariant', 0);
                jk=jk+1;
            end
        end
    end
    
    
end
% keyboard
Torque_1sim = Torque_1sim * Q /geo.Qs;
Flux_1sim = Flux_1sim * Q /geo.Qs / N_parallel;
% CoEnergy = CoEn * Mac.Q /Mac.Qs;
pm_loss = pm_loss * Q /geo.Qs;
Ppm = mean(sum(pm_loss));
if (xdeg == 360)
    %     Pfes_h = mean(h_loss_s) * Mac.Q;
    %     Pfes_c = mean(c_loss_s) * Mac.Q;
    %     Pfer_h = h_loss_r * 2 * Mac.p;
    %     Pfer_c = c_loss_r * 2 * Mac.p;
    %%  Matteo modifica preliminare
    Period_mot=gcd(geo.Qs,geo.p);
    if (Porz_stat==1)
        Pfes_h = mean(h_loss_s) * Q /geo.Qs;
        Pfes_c = mean(c_loss_s) * Q /geo.Qs;
    else
        Pfes_h = mean(h_loss_s) * Q;
        Pfes_c = mean(c_loss_s) * Q;
        
    end
    %     Pfer_h = h_loss_r * 2 * Mac.p/Period_mot;     % usati nel caso delle macchine S1 ed S2 possono essere erronee come linee di prog
    %     Pfer_c = c_loss_r * 2 * Mac.p/Period_mot;
    
    Pfer_h = h_loss_r / geo.ps * 2 * geo.p;
    Pfer_c = c_loss_r / geo.ps * 2 * geo.p;
    
    PjrBar=sum(mean(rotor_barrier_loss))*2 * geo.p/geo.ps;
    
    %%
end

% --> Fluxabc, Fluxdq (Npos x 1)
RunS_FluxElab_MN_2
RunS_VoltageElab_MN

Torque_0 = -mean(Torque_1sim(:,1));
% SaveData - single simulation
tmp1 = id*ones(size(time,2),1);
tmp2 = iq*ones(size(time,2),1);

%per avere il grafico Torque coincidente come in FEMM devo:
%parto dalla seconda riga, metto la sesta riga nella prima riga,elimino
%l'ultima riga
Torque_1sim = Torque_1sim((2:(end)),:);
Torque_1sim = [Torque_1sim(NTime-1,:);Torque_1sim];
Torque_1sim = Torque_1sim(1:(end-1),:);

% Ft.values=[time(2:end)' Flux_1sim(2:end,:) Ea Eb Ec tmp1(2:end) tmp2(2:end) Fluxd Fluxq Torque_1sim((2:end),:) CoEnergy(2:end)'];
% Ft.values=[time(2:end)' Flux_1sim(2:end,:) Ea Eb Ec tmp1(2:end) tmp2(2:end) Fluxd Fluxq Torque_1sim((1:(end-1)),:)];
% if (xdeg == 360)
% %     save([CaseFileName(1:end-3) '.mat'],'Ft','Cas','Mac','Sim','pm_loss','Ppm','Pfes_h','Pfes_c','Pfer_h','Pfer_c','PjrBar');
%       save([CaseFileName(1:end-3) '.mat'],'pm_loss','Ppm','Pfes_h','Pfes_c','Pfer_h','Pfer_c','PjrBar');
% else
% %     save([CaseFileName(1:end-3) '.mat'],'Ft','Cas','Mac','Sim','pm_loss','Ppm');
%       save([CaseFileName(1:end-3) '.mat'],'pm_loss','Ppm');
% end

