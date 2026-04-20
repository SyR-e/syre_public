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

function [SOL] = simulate_xdeg(geo,per,mat,eval_type,pathname,filename)

% Input: nsim, id, iq;
% Output:
% - vector of the solutions SOL (nsim x 6);
% - Temp\sim_gamma_numerosimulazione.fem;
% - Temp\sim_gamma.mat (memorizza SOL)

% - Output structure: SOL

% NB: pathname with final slash
% eval type determines number of simulated positions

gamma      = per.gamma;
th0        = geo.th0;
p          = geo.p;
ps         = geo.ps;
n3phase    = geo.win.n3phase; %AS number of 3-phase circuits
if(isnan(n3phase))
    Nbob       = geo.win.Nbob;
    nphases = 5;
else
    Nbob       = geo.win.Nbob.*ones(1,n3phase);
end

% Nbob       = geo.win.Nbob*[0.5 1 1];
% warning('Test different Nbob')
l          = geo.l;
if strcmp(geo.RotType,'EESM')
    If = per.if;       % Field current
    Nf = geo.win.Nf;   % EESM number of turns per pole
end

if(~isnan(n3phase))
    if isfield(per,'flag3phaseSet')
        flag3phSet = per.flag3phaseSet;
    else
        flag3phSet = ones(1,n3phase);
    end
end

if isfield(geo,'slidingGap')
    flagSG=1;
else
    flagSG=0;
end

if geo.ps==2*geo.p
    flagSG=0;
end

custom_act = 0;

if isfield(per,'custom_act')
    if per.custom_act
        if(isnan(geo.win.n3phase))
            curr_pos = per.rotorPos;
            Ia         = per.custom_ia(1,:,curr_pos);
            Ib         = per.custom_ib(1,:,curr_pos);
            Ic         = per.custom_ic(1,:,curr_pos);
            Id         = per.custom_id(1,:,curr_pos);
            Ie         = per.custom_ie(1,:,curr_pos);
            rotorPos   = per.custom_rotorPos;
            custom_act = per.custom_act;
        else
            Ia         = per.custom_ia;
            Ib         = per.custom_ib;
            Ic         = per.custom_ic;
            %time       = per.custom_time;
            custom_act = per.custom_act;
        end
    end
end

if ~isfield(per,'offset')
    per.offset = 0;
end

if ~isfield(per,'iOffset')
    per.iOffset = 0;
end

flag_xfemm = geo.XFEMMsimulation;

switch eval_type
    case 'MO_OA' % optimization
        nsim = per.nsim_MOOA;
        xdeg = per.delta_sim_MOOA;
        sim_step=xdeg/(nsim+0.5);
        offset=sim_step*rand;   % during optimization, random position offset
        theta=offset:sim_step:xdeg+offset;
        thetaPark=th0(1)+[theta(1:nsim) theta(1)]; % disregard the last position
        iOffset = 0;
        iOffsetPU = 0;
    case 'MO_GA' % optimization
        nsim = per.nsim_MOOA;
        xdeg = per.delta_sim_MOOA;
        sim_step=xdeg/(nsim+0.5);
        offset=sim_step*rand;   % during optimization, random position offset
        theta=offset:sim_step:xdeg+offset;
        thetaPark=th0(1)+[theta(1:nsim) theta(1)]; % disregard the last position
        iOffset = 0;
        iOffsetPU = 0;
    case 'singt' % single working point (id,iq)
        xdeg=per.delta_sim_singt;
        nsim = round(per.nsim_singt*xdeg/per.delta_sim_singt);
        sim_step=xdeg/(nsim);               % during re-evaluation, regular position steps
        offset = per.offset;
        theta=offset:sim_step:xdeg+offset;
        if(custom_act)
            if(isnan(n3phase))
                thetaPark=0*th0(1)+0*offset+geo.p*rotorPos(curr_pos)*ones(1, nsim+1); % 
                wt = linspace(0,360,nsim);
            else
                thetaPark=th0(1)+[theta(1:nsim) theta(1)]; % disregard the last position
            end
        else
            thetaPark=th0(1)+[theta(1:nsim) theta(1)]; % disregard the last position
        end
        iOffset = per.iOffset;
        iOffsetPU = iOffset/per.i0;
    case 'singm' % flux map
        xdeg = per.delta_sim_singt;
        nsim = round(per.nsim_singt*xdeg/per.delta_sim_singt);
        sim_step=xdeg/(nsim);               % during re-evaluation, regular position steps
        offset=0;
        theta=offset:sim_step:xdeg+offset;
        thetaPark=th0(1)+[theta(1:nsim) theta(1)]; % disregard the last position
        iOffset = per.iOffset;
        iOffsetPU = iOffset/per.i0;
    case {'idemag','idemagmap','demagArea'} % demagnetization
        nsim = per.nsim_singt;
        nsim = 2;
        if nsim==2
            theta = 0.5*360/(6*geo.p*geo.q*geo.win.n3phase)*[0 1];
            thetaPark=th0(1)+theta;
        else
            xdeg = per.delta_sim_singt;
            nsim = round(per.nsim_singt*xdeg/per.delta_sim_singt);
            sim_step=xdeg/(nsim);
            theta=0:sim_step:xdeg;
            thetaPark=th0(1)+[theta(1:nsim) theta(1)]; % disregard the last position
        end
        
        iOffset = per.iOffset;
        iOffsetPU = iOffset/per.i0;
        tmp=geo.BLKLABELS.rotore.xy;
        tmp = tmp(tmp(:,3) == 6,:);
        BrDir = atan2(tmp(:,7),tmp(:,6));
        BrGro=200+[1:1:numel(BrDir)];
    case 'flxdn' % flux density analysis (airgap, tooth, stator yoke)
        xdeg = per.delta_sim_singt;
        nsim = round(per.nsim_singt*xdeg/per.delta_sim_singt);
        sim_step=xdeg/(nsim);               % during re-evaluation, regular position steps
        offset=0;
        theta=offset:sim_step:xdeg+offset;
        thetaPark=th0(1)+[theta(1:nsim) theta(1)]; % disregard the last position
        iOffset = per.iOffset;
        iOffsetPU = iOffset/per.i0;

        % initialize the additional variables for flux densities
        angIni  = 0;
        angFin  = 360/(2*p)*ps;
        angRes  = 3600/(2*p)*ps+2;
        angVect = linspace(angIni,angFin,angRes);
        angVect = angVect(2:end-1); % the first and the last point must be avoided because they are on the boundary

        radTooth = geo.r+geo.g+geo.lt/2;
        xTooth   = radTooth*cosd(angVect);
        yTooth   = radTooth*sind(angVect);

        geo.ly   = geo.R-geo.r-geo.g-geo.lt;
        radYoke  = geo.R-geo.ly/2;
        xYoke    = radYoke*cosd(angVect);
        yYoke    = radYoke*sind(angVect);

        Bg = zeros(angRes-2,nsim+1);
        Bt = zeros(angRes-2,nsim+1);
        By = zeros(angRes-2,nsim+1);
        Bg(:,1) = angVect';
        By(:,1) = angVect';
        Bt(:,1) = angVect';
    case 'izero'
        xdeg=per.delta_sim_singt;
        nsim = round(per.nsim_singt*xdeg/per.delta_sim_singt);
        sim_step=xdeg/(nsim);               % during re-evaluation, regular position steps
        offset=0;
        theta=offset:sim_step:xdeg+offset;
        thetaPark=th0(1)+[theta(1:nsim) theta(1)]; % disregard the last position
        iOffsetPU = 0.2;
        warning('Simulation mode to be removed')
    case 'force'
        xdeg = per.delta_sim_singt;
        nsim = round(per.nsim_singt*xdeg/per.delta_sim_singt);
        sim_step=xdeg/(nsim);               % during re-evaluation, regular position steps
        offset=0;
        theta=offset:sim_step:xdeg+offset;
        thetaPark=th0(1)+[theta(1:nsim) theta(1)]; % disregard the last position
        iOffset = per.iOffset;
        iOffsetPU = iOffset/per.i0;

        % initialize the additional variables for flux densities
        angIni  = 0;
        angFin  = 360/(2*p)*ps;
        angRes  = 1002;
        angRes  = 3600/(2*p)*ps+2;
        angRes  = 100*360/(6*geo.p*geo.q*geo.win.n3phase)*geo.Qs+2;
        angVect = linspace(angIni,angFin,angRes);
        %angVect = angVect(2:end-1); % the first and the last point must be avoided because they are on the boundary
        %angRef = cumsum(diff(angVect));
        angStp = angVect(2)-angVect(1);
        angRef = cumsum(diff(angVect))-angStp/2;

        rF = geo.r+geo.g/6;
        rF = geo.r+geo.g*5/6;
        xF = rF*cosd(angVect);
        yF = rF*sind(angVect);

        Fr = zeros(angRes-1,nsim+1);
        Ft = zeros(angRes-1,nsim+1);
        Fr(:,1) = angRef';
        Ft(:,1) = angRef';
        Fx(:,1) = angRef';
        Fy(:,1) = angRef';
    case {'singtIron','singmIron'} % simulation with iron loss
        xdeg = per.delta_sim_singt;
        nsim = round(per.nsim_singt*xdeg/per.delta_sim_singt);
        sim_step = xdeg/(nsim);               % during re-evaluation, regular position steps
        offset = per.offset;
        theta = offset:sim_step:xdeg+offset;
        thetaPark = th0(1)+[theta(1:nsim) theta(1)]; % disregard the last position
        iOffset = per.iOffset;
        iOffsetPU = iOffset/per.i0;
end

% evaluation of the phase current values for all positions to be simulated
% iAmp = per.overload*calc_io(geo,per);
iAmp = per.overload*per.i0;

iAmpCoil = iAmp.*Nbob;
if(isnan(n3phase))
    iOffCoil = iOffset.*Nbob;

    id1 = iAmpCoil.*cos(gamma*pi/180).*ones(1,1);
    iq1 = iAmpCoil.*sin(gamma*pi/180).*ones(1,1);
    id3 = iAmpCoil.*cos(gamma*pi/180).*zeros(1,1);
    iq3 = iAmpCoil.*sin(gamma*pi/180).*zeros(1,1);
    io = iOffCoil;
else
    iOffCoil = iOffset.*Nbob.*ones(1,n3phase);

    id = iAmpCoil.*cos(gamma*pi/180).*ones(1,n3phase);
    iq = iAmpCoil.*sin(gamma*pi/180).*ones(1,n3phase);
    io = iOffCoil.*ones(1,n3phase);
end

if(isnan(n3phase))
    i_tmp = zeros(nphases,nsim);   %matrix containing all phase current values for the simulated rotor position
else
    i_tmp = zeros(3*n3phase,nsim);   %matrix containing all phase current values for the simulated rotor position
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% ciclo for %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open and draw motor once, rotate and simulate nsim positions
if ~flag_xfemm
    openfemm(1);
    opendocument([pathname,filename]);
else
    FemmProblem = loadfemmfile(checkPathSyntax([pathname filename]));
    femmFileName = [pathname filename];
    FemmProblem.ProbInfo.SmartMesh = 0;
end

SOL.th = zeros(1,nsim);         % electrical angle in degree
SOL.T  = zeros(1,nsim);         % Torque

if(isnan(n3phase))

    SOL.id1 = zeros(1,nsim);         % d-axis current
    SOL.iq1 = zeros(1,nsim);         % q-axis current
    SOL.id3 = zeros(1,nsim);         % d3-axis current
    SOL.iq3 = zeros(1,nsim);         % q3-axis current
    SOL.i0 = zeros(1,nsim);         % 0-axis current (homopolar)
    SOL.fd1 = zeros(1,nsim);         % d-axis flux linkage
    SOL.fq1 = zeros(1,nsim);         % q-axis flux linkage
    SOL.fd3 = zeros(1,nsim);         % d3-axis flux linkage
    SOL.fq3 = zeros(1,nsim);         % q3-axis flux linkage
    SOL.f0 = zeros(1,nsim);         % 0-axis flux linkage (homopolar)

    SOL.ia = zeros(1,nsim);   % phase a current
    SOL.ib = zeros(1,nsim);   % phase b current
    SOL.ic = zeros(1,nsim);   % phase c current
    SOL.id = zeros(1,nsim);   % phase d current
    SOL.ie = zeros(1,nsim);   % phase e current
    SOL.fa = zeros(1,nsim);   % phase a flux linkage
    SOL.fb = zeros(1,nsim);   % phase b flux linkage
    SOL.fc = zeros(1,nsim);   % phase c flux linkage
    SOL.fd = zeros(1,nsim);   % phase d flux linkage
    SOL.fe = zeros(1,nsim);   % phase e flux linkage
else
    SOL.id = zeros(1,nsim);         % d-axis current
    SOL.iq = zeros(1,nsim);         % q-axis current
    SOL.i0 = zeros(1,nsim);         % 0-axis current (homopolar)
    SOL.fd = zeros(1,nsim);         % d-axis flux linkage
    SOL.fq = zeros(1,nsim);         % q-axis flux linkage
    SOL.f0 = zeros(1,nsim);         % 0-axis flux linkage (homopolar)
    SOL.ia = zeros(n3phase,nsim);   % phase a current
    SOL.ib = zeros(n3phase,nsim);   % phase b current
    SOL.ic = zeros(n3phase,nsim);   % phase c current
    SOL.fa = zeros(n3phase,nsim);   % phase a flux linkage
    SOL.fb = zeros(n3phase,nsim);   % phase b flux linkage
    SOL.fc = zeros(n3phase,nsim);   % phase c flux linkage
end
SOL.we = zeros(1,nsim);         % magnetic energy
SOL.wc = zeros(1,nsim);         % magnetic coenergy

if(isnan(n3phase))
    phase_name = cell(nphases,1);
    phase_name_neg = cell(nphases,1);
else    
    phase_name = cell(n3phase*3,1);
    phase_name_neg = cell(n3phase*3,1);
end

%% Custom Current
if custom_act
    sim_range = round(xdeg/360*length(Ia(1,:)));

    theta_custom = linspace(0,360,length(Ia));

    if(isnan(n3phase))
        Ia_i = zeros(1,nsim);
        Ib_i = zeros(1,nsim);
        Ic_i = zeros(1,nsim);
        Id_i = zeros(1,nsim);
        Ie_i = zeros(1,nsim);
    else
        Ia_i = zeros(n3phase,nsim);
        Ib_i = zeros(n3phase,nsim);
        Ic_i = zeros(n3phase,nsim);
    end

    if(isnan(n3phase))
        Ia_i(1,:) = interp1(theta_custom,Ia(1,:),wt(1:end));
        Ib_i(1,:) = interp1(theta_custom,Ib(1,:),wt(1:end));
        Ic_i(1,:) = interp1(theta_custom,Ic(1,:),wt(1:end));
        Id_i(1,:) = interp1(theta_custom,Id(1,:),wt(1:end));
        Ie_i(1,:) = interp1(theta_custom,Ie(1,:),wt(1:end));
    else
        for ik=1:(n3phase)
            %         Ia_i(ik,:) = interp1(1:sim_range,Ia(ik,1:sim_range),linspace(1,sim_range,nsim));
            %         Ib_i(ik,:) = interp1(1:sim_range,Ib(ik,1:sim_range),linspace(1,sim_range,nsim));
            %         Ic_i(ik,:) = interp1(1:sim_range,Ic(ik,1:sim_range),linspace(1,sim_range,nsim));
    
            Ia_i(ik,:) = interp1(theta_custom,Ia(ik,:),theta(1:end-1));
            Ib_i(ik,:) = interp1(theta_custom,Ib(ik,:),theta(1:end-1));
            Ic_i(ik,:) = interp1(theta_custom,Ic(ik,:),theta(1:end-1));
        end
    end

    if(isnan(n3phase))
        if(isscalar(rotorPos))

            figure
            figSetting
            title('Phase A - Current')

            plot(linspace(0,xdeg,nsim),Ia_i(1,:),'DisplayName',['Interpolated - ' num2str(1)] )
            plot(linspace(0,xdeg,sim_range),Ia(1,1:sim_range),'DisplayName',['Input - ' num2str(1)]')
            plot(linspace(0,xdeg,sim_range),Ib(1,1:sim_range),'DisplayName',['Input - ' num2str(2)]')
            plot(linspace(0,xdeg,sim_range),Ic(1,1:sim_range),'DisplayName',['Input - ' num2str(3)]')
            plot(linspace(0,xdeg,sim_range),Id(1,1:sim_range),'DisplayName',['Input - ' num2str(4)]')
            plot(linspace(0,xdeg,sim_range),Ie(1,1:sim_range),'DisplayName',['Input - ' num2str(5)]')
        end
    else

        figure
        figSetting
        title('Phase A - Current')

        for ik=1:(n3phase)
            plot(linspace(0,xdeg,nsim),Ia_i(ik,:),'DisplayName',['Interpolated - ' num2str(ik) ' 3phase set'] )
            plot(linspace(0,xdeg,sim_range),Ia(ik,1:sim_range),'DisplayName',['Input - ' num2str(ik) ' 3phase set']')
            %plot(linspace(0,xdeg,sim_range),Ib(1,1:sim_range),'DisplayName','B')
            %plot(linspace(0,xdeg,sim_range),Ic(1,1:sim_range),'DisplayName','C')
        end
    end

    if(isnan(n3phase))
        if(custom_act)   
            axis([0 xdeg -1.2*abs(max(max(Ia))) 1.2*abs(max(max(Ia)))])
        else
            axis([0 xdeg -1.2*max(max(Ia)) 1.2*max(max(Ia))])
        end
    else
        axis([0 xdeg -1.2*max(max(Ia)) 1.2*max(max(Ia))])
    end
    legend show

    if(isnan(n3phase))
        %%% Nbob rescaling
        Ia_i = Ia_i*Nbob;
        Ib_i = Ib_i*Nbob;
        Ic_i = Ic_i*Nbob;
        Id_i = Id_i*Nbob;
        Ie_i = Ie_i*Nbob;
    else
        %%% Nbob rescaling
        Ia_i = Ia_i*Nbob;
        Ib_i = Ib_i*Nbob;
        Ic_i = Ic_i*Nbob;
    end
end

%% Simulation
for jj = 1:nsim

    % flag_xfemm = 0;
    % openfemm();
    % opendocument(femmFileName);

    % assign the phase current values to the FEMM circuits
    if(isnan(n3phase))
       if custom_act
            i_tmp(1,jj) = Ia_i(1,jj);
            i_tmp(2,jj) = Ib_i(1,jj);
            i_tmp(3,jj) = Ic_i(1,jj);
            i_tmp(4,jj) = Id_i(1,jj);
            i_tmp(5,jj) = Ie_i(1,jj);

            idqd3q30 = abcde2dqd3q30([Ia_i(1,jj)/Nbob;Ib_i(1,jj)/Nbob;Ic_i(1,jj)/Nbob;Id_i(1,jj)/Nbob;Ie_i(1,jj)/Nbob],(thetaPark(jj)+(th0(1)-th0(1)))*pi/180);

            id1 = idqd3q30(1);
            iq1 = idqd3q30(2);
            id3 = idqd3q30(3);
            iq3 = idqd3q30(4);
            io = idqd3q30(5);
       else
            % i123 = dq2abc(id(ik+1),iq(ik+1),(thetaPark(jj)+(th0(ik+1)-th0(1)))*pi/180);      % each 3phase set has its own offset angle
            i12345 = dqd3q302abcde([id1(1);iq1(1);id3(1);iq3(1);io(1)],(thetaPark(jj)+(th0(1)-th0(1)))*pi/180);      % each 5phase set has its own offset angle
            % i_tmp((3*ik)+1,jj) = (i123(1)+iOffCoil(ik+1))*flag3phSet(ik+1);
            % i_tmp((3*ik)+2,jj) = (i123(2)+iOffCoil(ik+1))*flag3phSet(ik+1);
            % i_tmp((3*ik)+3,jj) = (i123(3)+iOffCoil(ik+1))*flag3phSet(ik+1);
            i_tmp(1,jj) = i12345(1);
            i_tmp(2,jj) = i12345(2);
            i_tmp(3,jj) = i12345(3);
            i_tmp(4,jj) = i12345(4);
            i_tmp(5,jj) = i12345(5);
        end

        phase_name{1}     = strcat('fase',num2str(1));
        phase_name{2}     = strcat('fase',num2str(2));
        phase_name{3}     = strcat('fase',num2str(3));
        phase_name{4}     = strcat('fase',num2str(4));
        phase_name{5}     = strcat('fase',num2str(5));
        phase_name_neg{1} = strcat('fase',num2str(1),'n');
        phase_name_neg{2} = strcat('fase',num2str(2),'n');
        phase_name_neg{3} = strcat('fase',num2str(3),'n');
        phase_name_neg{4} = strcat('fase',num2str(4),'n');
        phase_name_neg{5} = strcat('fase',num2str(5),'n');

        % change current value in FEMM
        if ~flag_xfemm
            mi_modifycircprop(phase_name{1}, 1,i_tmp(1,jj));
            mi_modifycircprop(phase_name{2}, 1,i_tmp(2,jj));
            mi_modifycircprop(phase_name{3}, 1,i_tmp(3,jj));
            mi_modifycircprop(phase_name{4}, 1,i_tmp(4,jj));
            mi_modifycircprop(phase_name{5}, 1,i_tmp(5,jj));
            mi_modifycircprop(phase_name_neg{1}, 1,-i_tmp(1,jj));
            mi_modifycircprop(phase_name_neg{2}, 1,-i_tmp(2,jj));
            mi_modifycircprop(phase_name_neg{3}, 1,-i_tmp(3,jj));
            mi_modifycircprop(phase_name_neg{4}, 1,-i_tmp(4,jj));
            mi_modifycircprop(phase_name_neg{5}, 1,-i_tmp(5,jj));
            if strcmp(geo.RotType,'EESM')
                mi_modifycircprop('field',1,If*Nf);
            end
        else
            FemmProblem = setcircuitcurrent(FemmProblem,phase_name{1},i_tmp(1,jj));
            FemmProblem = setcircuitcurrent(FemmProblem,phase_name{2},i_tmp(2,jj));
            FemmProblem = setcircuitcurrent(FemmProblem,phase_name{3},i_tmp(3,jj));
            FemmProblem = setcircuitcurrent(FemmProblem,phase_name{4},i_tmp(4,jj));
            FemmProblem = setcircuitcurrent(FemmProblem,phase_name{5},i_tmp(5,jj));
            FemmProblem = setcircuitcurrent(FemmProblem,phase_name_neg{1},-i_tmp(1,jj));
            FemmProblem = setcircuitcurrent(FemmProblem,phase_name_neg{2},-i_tmp(2,jj));
            FemmProblem = setcircuitcurrent(FemmProblem,phase_name_neg{3},-i_tmp(3,jj));
            FemmProblem = setcircuitcurrent(FemmProblem,phase_name_neg{4},-i_tmp(4,jj));
            FemmProblem = setcircuitcurrent(FemmProblem,phase_name_neg{5},-i_tmp(5,jj));
            if strcmp(geo.RotType,'EESM')
                FemmProblem = setcircuitcurrent(FemmProblem,'field',If*Nf);
            end
        end
    else
        for ik=0:(n3phase-1)
            if custom_act
                i_tmp((3*ik)+1,jj) = Ia_i(ik+1,jj);
                i_tmp((3*ik)+2,jj) = Ib_i(ik+1,jj);
                i_tmp((3*ik)+3,jj) = Ic_i(ik+1,jj);
            else
                % i123 = dq2abc(id(ik+1),iq(ik+1),(thetaPark(jj)+(th0(ik+1)-th0(1)))*pi/180);      % each 3phase set has its own offset angle
                i123 = dq02abc([id(ik+1);iq(ik+1);io(ik+1)],(thetaPark(jj)+(th0(ik+1)-th0(1)))*pi/180);      % each 3phase set has its own offset angle
                % i_tmp((3*ik)+1,jj) = (i123(1)+iOffCoil(ik+1))*flag3phSet(ik+1);
                % i_tmp((3*ik)+2,jj) = (i123(2)+iOffCoil(ik+1))*flag3phSet(ik+1);
                % i_tmp((3*ik)+3,jj) = (i123(3)+iOffCoil(ik+1))*flag3phSet(ik+1);
                i_tmp((3*ik)+1,jj) = i123(1)*flag3phSet(ik+1);
                i_tmp((3*ik)+2,jj) = i123(2)*flag3phSet(ik+1);
                i_tmp((3*ik)+3,jj) = i123(3)*flag3phSet(ik+1);
            end
    
    
            phase_name{3*ik+1}     = strcat('fase',num2str(3*ik+1));
            phase_name{3*ik+2}     = strcat('fase',num2str(3*ik+2));
            phase_name{3*ik+3}     = strcat('fase',num2str(3*ik+3));
            phase_name_neg{3*ik+1} = strcat('fase',num2str(3*ik+1),'n');
            phase_name_neg{3*ik+2} = strcat('fase',num2str(3*ik+2),'n');
            phase_name_neg{3*ik+3} = strcat('fase',num2str(3*ik+3),'n');
    
            % change current value in FEMM
            if ~flag_xfemm
                mi_modifycircprop(phase_name{3*ik+1}, 1,i_tmp((3*ik)+1,jj));
                mi_modifycircprop(phase_name{3*ik+2}, 1,i_tmp((3*ik)+2,jj));
                mi_modifycircprop(phase_name{3*ik+3}, 1,i_tmp((3*ik)+3,jj));
                mi_modifycircprop(phase_name_neg{3*ik+1}, 1,-i_tmp((3*ik)+1,jj));
                mi_modifycircprop(phase_name_neg{3*ik+2}, 1,-i_tmp((3*ik)+2,jj));
                mi_modifycircprop(phase_name_neg{3*ik+3}, 1,-i_tmp((3*ik)+3,jj));
                if strcmp(geo.RotType,'EESM')
                    mi_modifycircprop('field',1,If*Nf);
                end
            else
                FemmProblem = setcircuitcurrent(FemmProblem,phase_name{3*ik+1},i_tmp((3*ik)+1,jj));
                FemmProblem = setcircuitcurrent(FemmProblem,phase_name{3*ik+2},i_tmp((3*ik)+2,jj));
                FemmProblem = setcircuitcurrent(FemmProblem,phase_name{3*ik+3},i_tmp((3*ik)+3,jj));
                FemmProblem = setcircuitcurrent(FemmProblem,phase_name_neg{3*ik+1},-i_tmp((3*ik)+1,jj));
                FemmProblem = setcircuitcurrent(FemmProblem,phase_name_neg{3*ik+2},-i_tmp((3*ik)+2,jj));
                FemmProblem = setcircuitcurrent(FemmProblem,phase_name_neg{3*ik+3},-i_tmp((3*ik)+3,jj));
                if strcmp(geo.RotType,'EESM')
                    FemmProblem = setcircuitcurrent(FemmProblem,'field',If*Nf);
                end
            end
        end
    end

    if length(mat.LayerMag.Hc)==1
        numPM = sum(geo.BLKLABELS.rotore.xy(:,3)==6);
        Hc_vect = mat.LayerMag.Hc*ones(1,numPM);
    else
        %  Hc_vect=[Hc Hc];  // DUBBIO -- a cosa serve? 2018 07 26
    end
    if ~strcmp(mat.LayerMag.MatName,'Air')
        if ~flag_xfemm
            for ii = 1:length(Hc_vect)
                mi_modifymaterial([mat.LayerMag.MatName '_' num2str(ii)],3,Hc_vect(ii));
            end
        else
            for ii = 1:length(Hc_vect)
                FemmProblem = modifymaterial_mfemm(FemmProblem,[mat.LayerMag.MatName '_' num2str(ii)],'H_c',Hc_vect(ii));
            end
        end
    end

    if(isnan(n3phase))
        if(custom_act)
            theta_r = rotorPos(curr_pos) - (th0(1) + 0*per.offset)/p;
        else
            theta_r = theta(jj)/p;
        end
    else
        theta_r = theta(jj)/p;
    end
    if flagSG
        % sliding gap (since 2019)
        if ~flag_xfemm
            mi_modifyboundprop('AGap',10,theta_r);  % modify inner boundary angle
        else
            FemmProblem = modifyboundprop_mfemm(FemmProblem,'AGap','InnerAngle',theta_r);
        end
    else
        % before sliding gap was introduced - delete the airgap arc prior to moving the rotor
        %         mi_selectgroup(20), mi_deleteselectedarcsegments;
        % rotate the rotor
        tmp = geo.PMdim(:);
        if strcmp(geo.RotType,'Circular')
            nPM = geo.ps*numel(tmp(tmp~=0))*2*2;
        else
            nPM = geo.ps*numel(tmp(tmp~=0))*2;
        end


        tmp = geo.BLKLABELS.rotore.xy(:,3);
        nPM = sum(tmp==6);

        indexPM = 200+(1:1:nPM);

        mi_selectgroup(22), mi_selectgroup(2);
        if strcmp(geo.RotType,'EESM')
            for ii = 201:1:200+4*p
                mi_selectgroup(ii);
            end
        end
        for kk=1:length(indexPM)
            mi_selectgroup(indexPM(kk));
        end
        if jj > 1
            mi_moverotate(0,0,(thetaPark(jj) - thetaPark(jj-1))/p);
        else
            mi_moverotate(0,0,theta_r);
        end
        % redraw the airgap arc
        if (ps<2*p)
            draw_airgap_arc_with_mesh(geo,theta_r,geo.mesh_res)
        else
            %             draw_airgap_arc_with_mesh_fullMachine(geo,theta_r,geo.mesh_res)
        end
    end
    
    % mi_saveas(femmFileName);
    % closefemm();
    % flag_xfemm = 1;

    if ~flag_xfemm
        mi_analyze(1);
        mi_loadsolution;
    else
        if ~contains(eval_type,'singm')
            disp(['XFEMM: simulating point ' int2str(jj) ' of ' int2str(nsim)])
        end
        writefemmfile(checkPathSyntax(femmFileName), FemmProblem);
        femmFileName = fmesher(checkPathSyntax(femmFileName));
        ansFile = fsolver(femmFileName,0,1);
        myfpproc = fpproc();
        myfpproc.opendocument(ansFile);
    end


    % flag_xfemm = 0;
    % openfemm();
    % opendocument(femmFileName);
    % mi_loadsolution;

    % load phase flux linkages
    if ~flag_xfemm
        if(isnan(n3phase))
            temp_out = mo_getcircuitproperties(phase_name{1});
            temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{1});
            f(1) = temp_out(3) * 2 * p/ps; %ps number of poles in FEMM
            temp_out = mo_getcircuitproperties(phase_name{2});
            temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{2});
            f(2) = temp_out(3) * 2 * p/ps;
            temp_out = mo_getcircuitproperties(phase_name{3});
            temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{3});
            f(3) = temp_out(3) * 2 * p/ps;
            temp_out = mo_getcircuitproperties(phase_name{4});
            temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{4});
            f(4) = temp_out(3) * 2 * p/ps;
            temp_out = mo_getcircuitproperties(phase_name{5});
            temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{5});
            f(5) = temp_out(3) * 2 * p/ps;
        else
            for ii=0:(n3phase-1) %AS
                temp_out = mo_getcircuitproperties(phase_name{3*ii+1});
                temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{3*ii+1});
                f(3*ii+1) = temp_out(3) * 2 * p/ps; %ps number of poles in FEMM
                temp_out = mo_getcircuitproperties(phase_name{3*ii+2});
                temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{3*ii+2});
                f(3*ii+2) = temp_out(3) * 2 * p/ps;
                temp_out = mo_getcircuitproperties(phase_name{3*ii+3});
                temp_out = temp_out - mo_getcircuitproperties(phase_name_neg{3*ii+3});
                f(3*ii+3) = temp_out(3) * 2 * p/ps;
            end
        end
    else
        if(isnan(n3phase))
            temp_out = myfpproc.getcircuitprops(phase_name{1});
            temp_out = temp_out - myfpproc.getcircuitprops(phase_name_neg{1});
            f(1) = temp_out(3) * 2 * p/ps; %ps number of poles in FEMM
            temp_out = myfpproc.getcircuitprops(phase_name{2});
            temp_out = temp_out - myfpproc.getcircuitprops(phase_name_neg{2});
            f(2) = temp_out(3) * 2 * p/ps;
            temp_out = myfpproc.getcircuitprops(phase_name{3});
            temp_out = temp_out - myfpproc.getcircuitprops(phase_name_neg{3});
            f(3) = temp_out(3) * 2 * p/ps;
            temp_out = myfpproc.getcircuitprops(phase_name{4});
            temp_out = temp_out - myfpproc.getcircuitprops(phase_name_neg{4});
            f(4) = temp_out(3) * 2 * p/ps;
            temp_out = myfpproc.getcircuitprops(phase_name{5});
            temp_out = temp_out - myfpproc.getcircuitprops(phase_name_neg{5});
            f(5) = temp_out(3) * 2 * p/ps;
        else
            for ii=0:(n3phase-1) %AS
                temp_out = myfpproc.getcircuitprops(phase_name{3*ii+1});
                temp_out = temp_out - myfpproc.getcircuitprops(phase_name_neg{3*ii+1});
                f(3*ii+1) = temp_out(3) * 2 * p/ps; %ps number of poles in FEMM
                temp_out = myfpproc.getcircuitprops(phase_name{3*ii+2});
                temp_out = temp_out - myfpproc.getcircuitprops(phase_name_neg{3*ii+2});
                f(3*ii+2) = temp_out(3) * 2 * p/ps;
                temp_out = myfpproc.getcircuitprops(phase_name{3*ii+3});
                temp_out = temp_out - myfpproc.getcircuitprops(phase_name_neg{3*ii+3});
                f(3*ii+3) = temp_out(3) * 2 * p/ps;
            end
        end
    end

    if(isnan(n3phase))
        % fdq = abc2dq(f(3*ik+1),f(3*ik+2),f(3*ik+3),(thetaPark(jj)+(th0(ik+1)-th0(1)))*pi/180);
        fdqd3q30 = abcde2dqd3q30([f(1);f(2);f(3);f(4);f(5)],(thetaPark(jj)+(th0(1)-th0(1)))*pi/180);
        %         fdq = abc2dq(f(3*ik+1),f(3*ik+2),f(3*ik+3),(thetaPark(jj)-ik*60/n3phase)*pi/180);
        %         fdq = abc2dq(f(3*ik+1),f(3*ik+2),f(3*ik+3),thetaPark(jj)*pi/180,n3phase,ik);
        fd1_temp(1,jj) = fdqd3q30(1);
        fq1_temp(1,jj) = fdqd3q30(2);
        fd3_temp(1,jj) = fdqd3q30(3);
        fq3_temp(1,jj) = fdqd3q30(4);
        fo_temp(1,jj) = fdqd3q30(5);
        fa_temp(1,jj) = f(1);
        fb_temp(1,jj) = f(2);
        fc_temp(1,jj) = f(3);
        fd_temp(1,jj) = f(4);
        fe_temp(1,jj) = f(5);

        fd1 = mean(fd1_temp(:,jj));
        fq1 = mean(fq1_temp(:,jj));
        fd3 = mean(fd3_temp(:,jj));
        fq3 = mean(fq3_temp(:,jj));
        fo = mean(fo_temp(:,jj));
    else
        for ik=0:(n3phase-1) %AS
            % fdq = abc2dq(f(3*ik+1),f(3*ik+2),f(3*ik+3),(thetaPark(jj)+(th0(ik+1)-th0(1)))*pi/180);
            fdq0 = abc2dq0([f(3*ik+1);f(3*ik+2);f(3*ik+3)],(thetaPark(jj)+(th0(ik+1)-th0(1)))*pi/180);
            %         fdq = abc2dq(f(3*ik+1),f(3*ik+2),f(3*ik+3),(thetaPark(jj)-ik*60/n3phase)*pi/180);
            %         fdq = abc2dq(f(3*ik+1),f(3*ik+2),f(3*ik+3),thetaPark(jj)*pi/180,n3phase,ik);
            fd_temp(ik+1,jj) = fdq0(1);
            fq_temp(ik+1,jj) = fdq0(2);
            fo_temp(ik+1,jj) = fdq0(3);
            fa_temp(ik+1,jj) = f(3*ik+1);
            fb_temp(ik+1,jj) = f(3*ik+2);
            fc_temp(ik+1,jj) = f(3*ik+3);
        end

        fd = mean(fd_temp(:,jj));
        fq = mean(fq_temp(:,jj));
        fo = mean(fo_temp(:,jj));
    end

    % Torque computation. For old model, the rotor blocks are selected and
    % torque is computed as block integral. For new models (with
    % slidingGap), torque is directly computed from the airgap boundary
    if flagSG
        if ~flag_xfemm
            %mo_groupselectblock(2), mo_groupselectblock(22), mo_groupselectblock(200)
            if geo.ps==2*geo.p
                %             mo_groupselectblock(2), mo_groupselectblock(22), mo_groupselectblock(200)
                %             T=mo_blockintegral(22)*2*p/ps;
                %             mo_clearblock;
                for ii=1:length(geo.BLKLABELS.rotore.xy(:,1))
                    xB=geo.BLKLABELS.rotore.xy(ii,1);
                    yB=geo.BLKLABELS.rotore.xy(ii,2);
                    %         if ~flagSG
                    [xB,yB]=rot_point(xB,yB,theta_r*pi/180);
                    %         end
                    mo_selectblock(xB,yB);
                end
                T=mo_blockintegral(22)*2*p/ps;
                mo_clearblock;
            else
                T = mo_gapintegral('AGap',0);
            end
        else
            %mo_groupselectblock(2), mo_groupselectblock(22), mo_groupselectblock(200)
            if geo.ps==2*geo.p
                %             mo_groupselectblock(2), mo_groupselectblock(22), mo_groupselectblock(200)
                %             T=mo_blockintegral(22)*2*p/ps;
                %             mo_clearblock;
                for ii=1:length(geo.BLKLABELS.rotore.xy(:,1))
                    xB=geo.BLKLABELS.rotore.xy(ii,1);
                    yB=geo.BLKLABELS.rotore.xy(ii,2);
                    %         if ~flagSG
                    [xB,yB]=rot_point(xB,yB,theta_r*pi/180);
                    %         end
                    myfpproc.selectblock(xB,yB);
                end
                T=myfpproc.blockintegral(22)*2*p/ps;
                myfpproc.clearblock;
            else
                T = myfpproc.gapintegral('AGap',0);
            end
        end
    else
        if ~flag_xfemm
            for ii=1:length(geo.BLKLABELS.rotore.xy(:,1))
                xB=geo.BLKLABELS.rotore.xy(ii,1);
                yB=geo.BLKLABELS.rotore.xy(ii,2);
                %         if ~flagSG
                [xB,yB]=rot_point(xB,yB,theta_r*pi/180);
                %         end
                mo_selectblock(xB,yB);
            end
            T=mo_blockintegral(22)*2*p/ps;
            mo_clearblock;
        else
            for ii=1:length(geo.BLKLABELS.rotore.xy(:,1))
                xB=geo.BLKLABELS.rotore.xy(ii,1);
                yB=geo.BLKLABELS.rotore.xy(ii,2);
                %         if ~flagSG
                [xB,yB]=rot_point(xB,yB,theta_r*pi/180);
                %         end
                myfpproc.selectblock(xB,yB);
            end
            T=myfpproc.blockintegral(22)*2*p/ps;
            myfpproc.clearblock;
        end
    end
    
    if ~flag_xfemm
        mo_groupselectblock();
        we = mo_blockintegral(2)*2*p/ps;
        wc = mo_blockintegral(17)*2*p/ps;
        mo_clearblock();
    else
        %myfpproc.groupselectblock();
        % we = myfpproc.blockintegral(2)*2*p/ps;
        % wc = myfpproc.blockintegral(17)*2*p/ps;
        % myfpproc.clearblock();
        we = NaN;
        wc = NaN;
    end

    if(custom_act)
        if(isnan(n3phase))
            SOL.th(jj) = wt(jj);
        else
            SOL.th(jj) = thetaPark(jj);
        end
    else
        SOL.th(jj) = thetaPark(jj);
    end
    SOL.T(jj)  = T;
    SOL.we(jj) = we;
    SOL.wc(jj) = wc;

    if(isnan(n3phase))

        SOL.id1(jj) = mean(mean(id1)./Nbob); % Divide by Ns (simulation done with one turn per coil)
        SOL.iq1(jj) = mean(mean(iq1)./Nbob);
        SOL.id3(jj) = mean(mean(id3)./Nbob); % Divide by Ns (simulation done with one turn per coil)
        SOL.iq3(jj) = mean(mean(iq3)./Nbob);
        SOL.io(jj) = mean(mean(io)./Nbob);
        SOL.fd1(jj) = mean(fd1.*Nbob); % Times Ns
        SOL.fq1(jj) = mean(fq1.*Nbob);
        SOL.fd3(jj) = mean(fd3.*Nbob); % Times Ns
        SOL.fq3(jj) = mean(fq3.*Nbob);
        SOL.fo(jj) = mean(fo.*Nbob);

        SOL.ia(1,jj) = i_tmp(1,jj)/Nbob(1);
        SOL.ib(1,jj) = i_tmp(2,jj)/Nbob(1);
        SOL.ic(1,jj) = i_tmp(3,jj)/Nbob(1);
        SOL.id(1,jj) = i_tmp(4,jj)/Nbob(1);
        SOL.ie(1,jj) = i_tmp(5,jj)/Nbob(1);
        SOL.fa(1,jj) = f(1)*Nbob(1);
        SOL.fb(1,jj) = f(2)*Nbob(1);
        SOL.fc(1,jj) = f(3)*Nbob(1); 
        SOL.fd(1,jj) = f(4)*Nbob(1); 
        SOL.fe(1,jj) = f(5)*Nbob(1); 
    else
        SOL.id(jj) = mean(mean(id)./Nbob); % Divide by Ns (simulation done with one turn per coil)
        SOL.iq(jj) = mean(mean(iq)./Nbob);
        SOL.io(jj) = mean(mean(io)./Nbob);
        SOL.fd(jj) = mean(fd.*Nbob); % Times Ns
        SOL.fq(jj) = mean(fq.*Nbob);
        SOL.fo(jj) = mean(fo.*Nbob);

        for ff=1:n3phase
            SOL.ia(ff,jj) = i_tmp(1+3*(ff-1),jj)/Nbob(ff);
            SOL.ib(ff,jj) = i_tmp(2+3*(ff-1),jj)/Nbob(ff);
            SOL.ic(ff,jj) = i_tmp(3+3*(ff-1),jj)/Nbob(ff);
            SOL.fa(ff,jj) = f(1+3*(ff-1))*Nbob(ff);
            SOL.fb(ff,jj) = f(2+3*(ff-1))*Nbob(ff);
            SOL.fc(ff,jj) = f(3+3*(ff-1))*Nbob(ff);
        end
    end

    if strcmp(geo.RotType,'EESM')
        SOL.if(jj) = If;
        if ~flag_xfemm
            temp_out = mo_getcircuitproperties('field');
        else
            temp_out = myfpproc.getcircuitprops('field');
        end
         
        SOL.ff(jj) = temp_out(3) * 2 * p/ps*Nf; %ps number of poles in FEMM
    end

    switch eval_type
        case 'flxdn'
            if ~flag_xfemm
                for ff=1:angRes-2
                    tmp = mo_getgapb('AGap',angVect(ff));
                    %Bg(ff,jj+1) = (tmp(1)^2+tmp(2)^2)^0.5;
                    Bg(ff,jj+1) = tmp(1);
                    tmp = mo_getb(xTooth(ff),yTooth(ff));
                    Bt(ff,jj+1) = (tmp(1)^2+tmp(2)^2)^0.5;
                    tmp = mo_getb(xYoke(ff),yYoke(ff));
                    By(ff,jj+1) = (tmp(1)^2+tmp(2)^2)^0.5;
                end
            else
                for ff=1:angRes-2
                    tmp = myfpproc.getgapb('AGap',angVect(ff));
                    Bg(ff,jj+1) = tmp(1);
                    tmp = myfpproc.getb(xTooth(ff),yTooth(ff));
                    Bt(ff,jj+1) = (tmp(1)^2+tmp(2)^2)^0.5;
                    tmp = myfpproc.getb(xYoke(ff),yYoke(ff));
                    By(ff,jj+1) = (tmp(1)^2+tmp(2)^2)^0.5;
                end
            end
            SOL.Bg = Bg;
            SOL.Bt = Bt;
            SOL.By = By;

        case 'force'
            for ff=1:length(angRef)
                aS = angRef(ff)-angStp/2;
                aE = angRef(ff)+angStp/2;
                xS = rF*cosd(aS);
                yS = rF*sind(aS);
                xE = rF*cosd(aE);
                yE = rF*sind(aE);
                if ~flag_xfemm
                    mo_clearcontour();
                    mo_addcontour(xS,yS);
                    mo_addcontour(xE,yE);
                    mo_bendcontour(angStp,0.5)
                    tmp = mo_lineintegral(3);
                    mo_clearcontour();
                else
                    myfpproc.clearcontour();
                    myfpproc.addcontour(xS,yS);
                    myfpproc.addcontour(xE,yE);
                    %mo_bendcontour(angStp,0.5)
                    tmp = myfpproc.lineintegral(3);
                end
                fx = tmp(1);
                fy = tmp(2);
                fm = (fx^2+fy^2)^0.5;
                fa = atan2(fy,fx);
                Fr(ff,jj+1) = fm*cos(fa-angRef(ff)*pi/180)/(rF/1000*angStp*pi/180*l/1000);
                Ft(ff,jj+1) = fm*sin(fa-angRef(ff)*pi/180)/(rF/1000*angStp*pi/180*l/1000);
                Fx(ff,jj+1) = fx/(rF/1000*angStp*pi/180*l/1000);
                Fy(ff,jj+1) = fy/(rF/1000*angStp*pi/180*l/1000);
                
            end

            SOL.Fr = Fr;
            SOL.Ft = Ft;
            SOL.Fx = Fx;
            SOL.Fy = Fy;

        case {'idemag','idemagmap','demagArea'}
            Br = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,per.tempPP);
            Bd = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Bd,per.tempPP);
            if ~flag_xfemm
                mo_showdensityplot(1,0,Bd,Br,'bmag');
            end

            Bmin = Br;

            % get the PMs flux density
            if ~flag_xfemm
                EleNo = mo_numelements;               % Number of mesh elements
            else
                EleNo = myfpproc.numelements;
            end
            pos = zeros(EleNo,1);                 % Mesh elements centroid coordinates as complex number
            groNo = zeros(EleNo,1);               % group number
            area = zeros(EleNo,1);
            vert = zeros(EleNo,3);
            for ee = 1:EleNo
                if ~flag_xfemm
                    elm = mo_getelement(ee);
                else
                    elm = myfpproc.getelements(ee);
                end
                pos(ee) = elm(4)+j*elm(5);
                groNo(ee) = elm(7);
                area(ee) = elm(6);
                if ~flag_xfemm
                    for vv=1:3
                        tmpV=mo_getnode(elm(vv));
                        vert(ee,vv)=tmpV(1)+j*tmpV(2);
                    end
                else
                    tmp = myfpproc.getvertices(ee);
                    vert(ee,1) = tmp(1)+j*tmp(2);
                    vert(ee,2) = tmp(3)+j*tmp(4);
                    vert(ee,3) = tmp(5)+j*tmp(6);
                end
            end
            EleIn=1:1:EleNo;
            EleOK=EleIn(groNo>=min(BrGro)&groNo<=max(BrGro));
            groNo=groNo(groNo>=min(BrGro)&groNo<=max(BrGro));
            num=0;
            den=0;
            xyDemagTmpC=[];
            xyDemagTmpV=[];
            xyDemagTmpB=[];
            for ee=1:length(EleOK)
                if ~flag_xfemm
                    tmp=mo_getpointvalues(real(pos(EleOK(ee))),imag(pos(EleOK(ee))));
                else
                    tmp = myfpproc.getpointvalues(real(pos(EleOK(ee))),imag(pos(EleOK(ee))));
                end
                Btmp=abs(tmp(2)+j*tmp(3))*cos(angle(tmp(2)+j*tmp(3))-BrDir(groNo(ee)-200));

                if Btmp<Bmin
                    Bmin=Btmp;
                end
                if Btmp<Bd
                    num=num+area(EleOK(ee));
                    xyDemagTmpC=[xyDemagTmpC, pos(EleOK(ee))];
                    xyDemagTmpV=[xyDemagTmpV, vert(EleOK(ee),:).'];
                    xyDemagTmpB=[xyDemagTmpB, Btmp];
                end
                den=den+area(EleOK(ee));
            end
            VolDem=0;
            VolTot=den;
            if num>0
                dPM=num/den;
                VolDem=num;
                VolTot=den;
            else
                dPM=0;
                VolDem=0;
                VolTot=den;
            end

            SOL.xyDemagTmpC{jj} = xyDemagTmpC;
            SOL.xyDemagTmpV{jj} = xyDemagTmpV;
            SOL.xyDemagTmpB{jj} = xyDemagTmpB;

            SOL.Bmin(jj)      = Bmin;
            SOL.dPM(jj)       = dPM;
            SOL.VolDem(jj)    = VolDem;
            SOL.VolTot(jj)    = VolTot;
        case {'singtIron','singmIron'} % simulation with iron loss
            if jj==1 % store the mesh information
                if ~flag_xfemm
                    EleNo = mo_numelements;               % Number of mesh elements
                else
                    EleNo = myfpproc.numelements;
                end
                pos   = zeros(1,EleNo);               % Mesh elements centroid coordinates as complex number
                groNo = zeros(1,EleNo);               % group number
                area  = zeros(1,EleNo);
                %                 vert  = zeros(EleNo,3);
                for ee = 1:EleNo
                    if ~flag_xfemm
                        elm = mo_getelement(ee);
                    else
                        elm = myfpproc.getelements(ee);
                    end
                    pos(ee) = elm(4)+j*elm(5);
                    groNo(ee) = elm(7);
                    area(ee) = elm(6);
                end
                EleIn = 1:1:EleNo;
                EleOK = [EleIn(groNo==12) EleIn(groNo==22) EleIn(groNo>200)]; % select just the element of stator iron (group=12), rotor iron (group=22) and PMs (group>200)
                EleNo = length(EleOK);
                pos   = pos(EleOK);
                groNo = groNo(EleOK);
                area  = area(EleOK);
                bs    = zeros(nsim,EleNo);
                br    = zeros(nsim,EleNo);
                am    = zeros(nsim,EleNo);
            end

            % download from FEMM the flux density data for iron loss
            % computation (for each mesh element)
            for ee=1:EleNo
                if ~flag_xfemm
                    tmp = mo_getpointvalues(real(pos(ee)),imag(pos(ee)));
                else
                    tmp = myfpproc.getpointvalues(real(pos(ee)),imag(pos(ee)));
                end
                switch groNo(ee)
                    case 12  % stator iron
                        bs(jj,ee) = tmp(2)+j*tmp(3);
                    case 22  % rotor iron
                        br(jj,ee) = tmp(2)+j*tmp(3);
                    otherwise % rotor PM
                        am(jj,ee) = tmp(1);
                end
            end
    end

    if strcmp(eval_type,'singtIron')
        disp(['Single point Iron loss evaluation - evaluated ' int2str(jj) ' of ' int2str(nsim)])
    end

    % flag_xfemm = 1;
    % closefemm;

end

if (strcmp(eval_type,'singtIron')||strcmp(eval_type,'singmIron'))
    % computation of the elements volume
    vol = (area/1e6)*(l/1000);

    SOL.bs    = bs;
    SOL.br    = br;
    SOL.am    = am;
    SOL.pos   = pos;
    SOL.vol   = vol;
    SOL.groNo = groNo;
    if isfield(per,'flag_singtpp')
        flag_singtpp = per.flag_singtpp;
    else
        flag_singtpp = 0;
    end
    if (strcmp(eval_type,'singmIron')||strcmp(eval_type,'singtIron')&&(flag_singtpp==0))
        [SOL] = evalIronLossFEMM(geo,per,mat,SOL,2);

        % the last input of evalIronLossFEMM is the method of the iron loss
        % computation:
        % 0 --> use the standard method introduced by FEMM: harmonic
        %       decomposition both for hysteresis and eddy current
        % 1 --> use the iGSE along the main flux density axis for hysteresis
        %       and harmonic decomposition for eddy current
        % 2 --> use iGSE along the main and the quadrature flux density axis
        %       for hysteresis (to account for the rotational loss) and harmonic
        %       decomposition for eddy current.
        % 3 --> use iGSE along the main and the quadrature flux density axis
        %       for hysteresis (to account for the rotational loss) and harmonic
        %       decomposition for eddy current. For the hysteresis, correction
        %       factor are adopted.

    end
end

if ~flag_xfemm
    mo_close, mi_close
    closefemm
else

end


