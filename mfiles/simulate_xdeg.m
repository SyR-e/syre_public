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

gamma   = per.gamma;
th0     = geo.th0;
p       = geo.p;
ps      = geo.ps;
n3phase = geo.win.n3phase; %AS number of 3-phase circuits
Nbob    = geo.win.Nbob;
l       = geo.l;
if isfield(geo,'slidingGap')
    flagSG=1;
else
    flagSG=0;
end

switch eval_type
    case 'MO_OA' % optimization
        nsim = per.nsim_MOOA;
        xdeg = per.delta_sim_MOOA;
        sim_step=xdeg/(nsim+0.5);
        offset=sim_step*rand;   % during optimization, random position offset
        theta=offset:sim_step:xdeg+offset;
        thetaPark=th0(1)+[theta(1:nsim) theta(1)]; % disregard the last position
        iOffsetPU = 0;
    case 'MO_GA' % optimization
        nsim = per.nsim_MOOA;
        xdeg = per.delta_sim_MOOA;
        sim_step=xdeg/(nsim+0.5);
        offset=sim_step*rand;   % during optimization, random position offset
        theta=offset:sim_step:xdeg+offset;
        thetaPark=th0(1)+[theta(1:nsim) theta(1)]; % disregard the last position
        iOffsetPU = 0;
    case 'singt' % single working point (id,iq)
        xdeg=per.delta_sim_singt;
        nsim = round(per.nsim_singt*xdeg/per.delta_sim_singt);
        sim_step=xdeg/(nsim);               % during re-evaluation, regular position steps
        offset=0;
        theta=offset:sim_step:xdeg+offset;
        thetaPark=th0(1)+[theta(1:nsim) theta(1)]; % disregard the last position
        iOffsetPU = 0;
    case 'singm' % flux map
        xdeg = per.delta_sim_singt;
        nsim = round(per.nsim_singt*xdeg/per.delta_sim_singt);
        sim_step=xdeg/(nsim);               % during re-evaluation, regular position steps
        offset=0;
        theta=offset:sim_step:xdeg+offset;
        thetaPark=th0(1)+[theta(1:nsim) theta(1)]; % disregard the last position
        iOffsetPU = 0;
    case {'idemag','idemagmap','demagArea'} % demagnetization
        nsim = 1;
        theta=per.delta_sim_singt;
        thetaPark=th0(1)+theta;
        iOffsetPU = 0;
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
        iOffsetPU = 0;
        
        % initialize the additional variables for flux densities
        angIni  = 0;
        angFin  = 360/(2*p)*ps;
        angRes  = 1002;
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
    case 'force'
        xdeg = per.delta_sim_singt;
        nsim = round(per.nsim_singt*xdeg/per.delta_sim_singt);
        sim_step=xdeg/(nsim);               % during re-evaluation, regular position steps
        offset=0;
        theta=offset:sim_step:xdeg+offset;
        thetaPark=th0(1)+[theta(1:nsim) theta(1)]; % disregard the last position
        iOffsetPU = 0;
        
        % initialize the additional variables for flux densities
        angIni  = 0;
        angFin  = 360/(2*p)*ps;
        angRes  = 1002;
        angVect = linspace(angIni,angFin,angRes);
        %angVect = angVect(2:end-1); % the first and the last point must be avoided because they are on the boundary
        %angRef = cumsum(diff(angVect));
        angStp = angVect(2)-angVect(1);
        angRef = cumsum(diff(angVect))-angStp/2;
        
        rF = geo.r+geo.g/6;
        xF = rF*cosd(angVect);
        yF = rF*sind(angVect);
        
        Fr = zeros(angRes-1,nsim+1);
        Ft = zeros(angRes-1,nsim+1);
        Fr(:,1) = angRef';
        Ft(:,1) = angRef';
    case {'singtIron','singmIron'} % simulation with iron loss
        xdeg = per.delta_sim_singt;
        nsim = round(per.nsim_singt*xdeg/per.delta_sim_singt);
        sim_step = xdeg/(nsim);               % during re-evaluation, regular position steps
        offset = 0;
        theta = offset:sim_step:xdeg+offset;
        thetaPark = th0(1)+[theta(1:nsim) theta(1)]; % disregard the last position
        iOffsetPU = 0;
end

% evaluation of the phase current values for all positions to be simulated
% iAmp = per.overload*calc_io(geo,per);
iAmp = per.overload*per.i0;
if iAmp==0
    iOff = iOffsetPU*per.i0;
else
    iOff = iOffsetPU*iAmp;
end

iAmpCoil = iAmp*Nbob;
iOffCoil = iOff*Nbob;

id = iAmpCoil * cos(gamma * pi/180);
iq = iAmpCoil * sin(gamma * pi/180);

i_tmp = zeros(3*n3phase,nsim);   %matrix containing all phase current values for the simulated rotor position

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% ciclo for %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open and draw motor once, rotate and simulate nsim positions
openfemm(1);
opendocument([pathname,filename]);

SOL.th = zeros(1,nsim);         % electrical angle in degree
SOL.id = zeros(1,nsim);         % d-axis current
SOL.iq = zeros(1,nsim);         % q-axis current
SOL.fd = zeros(1,nsim);         % d-axis flux linkage
SOL.fq = zeros(1,nsim);         % q-axis flux linkage
SOL.T  = zeros(1,nsim);         % Torque from block integral (suggested by David Meeker)
SOL.ia = zeros(n3phase,nsim);   % phase a current
SOL.ib = zeros(n3phase,nsim);   % phase b current
SOL.ic = zeros(n3phase,nsim);   % phase c current
SOL.fa = zeros(n3phase,nsim);   % phase a flux linkage
SOL.fb = zeros(n3phase,nsim);   % phase b flux linkage
SOL.fc = zeros(n3phase,nsim);   % phase c flux linkage

phase_name = cell(n3phase*3,1);
phase_name_neg = cell(n3phase*3,1);

for jj = 1:nsim    
    % assign the phase current values to the FEMM circuits
    for ik=0:(n3phase-1)
        if geo.win.avv_flag((3*ik)+1)==1 && geo.win.avv_flag((3*ik)+2)==1 && geo.win.avv_flag((3*ik)+3)==1
            % healthy winding set
            %             i123 = dq2abc(id,iq,thetaPark(jj)*pi/180,n3phase,ik);
            %             i123 = dq2abc(id,iq,(thetaPark(jj)-ik*60/n3phase)*pi/180); % AS version
            i123 = dq2abc(id,iq,(thetaPark(jj)+(th0(ik+1)-th0(1)))*pi/180);             % each 3phase set has its own offset angle
            i_tmp((3*ik)+1,jj) = i123(1)+iOffCoil;
            i_tmp((3*ik)+2,jj) = i123(2)+iOffCoil;
            i_tmp((3*ik)+3,jj) = i123(3)+iOffCoil;
        else
            % open-circuit winding set
            if geo.win.avv_flag((3*ik)+1)==0 && geo.win.avv_flag((3*ik)+2)==0 && geo.win.avv_flag((3*ik)+3)==0
                i_tmp((3*ik)+1,jj) = 0;
                i_tmp((3*ik)+2,jj) = 0;
                i_tmp((3*ik)+3,jj) = 0;
            end
        end
        
        phase_name{3*ik+1}=strcat('fase',num2str(3*ik+1));
        phase_name{3*ik+2}=strcat('fase',num2str(3*ik+2));
        phase_name{3*ik+3}=strcat('fase',num2str(3*ik+3));
        phase_name_neg{3*ik+1}=strcat('fase',num2str(3*ik+1),'n');
        phase_name_neg{3*ik+2}=strcat('fase',num2str(3*ik+2),'n');
        phase_name_neg{3*ik+3}=strcat('fase',num2str(3*ik+3),'n');
        
        % change current value in FEMM
        mi_modifycircprop(phase_name{3*ik+1}, 1,i_tmp((3*ik)+1,jj));
        mi_modifycircprop(phase_name{3*ik+2}, 1,i_tmp((3*ik)+2,jj));
        mi_modifycircprop(phase_name{3*ik+3}, 1,i_tmp((3*ik)+3,jj));
        mi_modifycircprop(phase_name_neg{3*ik+1}, 1,-i_tmp((3*ik)+1,jj));
        mi_modifycircprop(phase_name_neg{3*ik+2}, 1,-i_tmp((3*ik)+2,jj));
        mi_modifycircprop(phase_name_neg{3*ik+3}, 1,-i_tmp((3*ik)+3,jj));
        
    end
    
    % assign the Hc property to each of the bonded magnets
    %     if strcmp(geo.RotType,'SPM')
    %         mi_modifymaterial(mat.LayerMag.MatName,3,mat.LayerMag.Hc);
    %     else
    if length(mat.LayerMag.Hc)==1
        Hc_vect = mat.LayerMag.Hc*ones(1,length(geo.BLKLABELS.rotore.BarName));
    else
        %             Hc_vect=[Hc Hc];  // DUBBIO -- a cosa serve? 2018 07 26
    end
    for ii = 1:length(Hc_vect)
        mi_modifymaterial([mat.LayerMag.MatName '_' num2str(ii)],3,Hc_vect(ii));
    end
    %     end
    
    theta_r = theta(jj)/p;
    if flagSG
        % sliding gap (since 2019)
        mi_modifyboundprop('AGap',10,theta_r);  % modify inner boundary angle
    else
        % before sliding gap was introduced - delete the airgap arc prior to moving the rotor
        mi_selectgroup(20), mi_deleteselectedarcsegments;
        % rotate the rotor
        mi_selectgroup(22), mi_selectgroup(2), mi_selectgroup(200),
        if jj > 1
            mi_moverotate(0,0,(thetaPark(jj) - thetaPark(jj-1))/p);
        else
            mi_moverotate(0,0,theta_r);
        end
        % redraw the airgap arc
        if (ps<2*p)
            draw_airgap_arc_with_mesh(geo,theta_r,geo.mesh_res)
        else
            draw_airgap_arc_with_mesh_fullMachine(geo,theta_r,geo.mesh_res)
        end
    end
    
    mi_analyze(1);
    mi_loadsolution;
    
    % load phase flux linkages
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
    
    for ik=0:(n3phase-1) %AS
        fdq = abc2dq(f(3*ik+1),f(3*ik+2),f(3*ik+3),(thetaPark(jj)+(th0(ik+1)-th0(1)))*pi/180);
        %         fdq = abc2dq(f(3*ik+1),f(3*ik+2),f(3*ik+3),(thetaPark(jj)-ik*60/n3phase)*pi/180);
        %         fdq = abc2dq(f(3*ik+1),f(3*ik+2),f(3*ik+3),thetaPark(jj)*pi/180,n3phase,ik);
        fd_temp(ik+1,jj)=fdq(1);
        fq_temp(ik+1,jj)=fdq(2);
        fa_temp(ik+1,jj)=f(3*ik+1);
        fb_temp(ik+1,jj)=f(3*ik+2);
        fc_temp(ik+1,jj)=f(3*ik+3);
    end
    
    fd=mean(fd_temp(:,jj));
    fq=mean(fq_temp(:,jj));
    
    % Torque computation. For old model, the rotor blocks are selected and
    % torque is computed as block integral. For new models (with
    % slidingGap), torque is directly computed from the airgap boundary
    if flagSG
        %mo_groupselectblock(2), mo_groupselectblock(22), mo_groupselectblock(200)
        T = mo_gapintegral('AGap',0);
    else
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
    end
    
    SOL.th(jj) = thetaPark(jj);
    SOL.id(jj) = id/Nbob; % Divide by Ns (simulation done with one turn per coil)
    SOL.iq(jj) = iq/Nbob;
    SOL.fd(jj) = fd*Nbob; % Times Ns
    SOL.fq(jj) = fq*Nbob;
    SOL.T(jj)  = T;
    
    for ff=1:n3phase
        SOL.ia(ff,jj) = i_tmp(1+3*(ff-1),jj)/Nbob;
        SOL.ib(ff,jj) = i_tmp(2+3*(ff-1),jj)/Nbob;
        SOL.ic(ff,jj) = i_tmp(3+3*(ff-1),jj)/Nbob;
        SOL.fa(ff,jj) = f(1+3*(ff-1))*Nbob;
        SOL.fb(ff,jj) = f(2+3*(ff-1))*Nbob;
        SOL.fc(ff,jj) = f(3+3*(ff-1))*Nbob;
    end
    
    switch eval_type
        case 'flxdn'
            for ff=1:angRes-2
                tmp = mo_getgapb('AGap',angVect(ff));
                %Bg(ff,jj+1) = (tmp(1)^2+tmp(2)^2)^0.5;
                Bg(ff,jj+1) = tmp(1);
                tmp = mo_getb(xTooth(ff),yTooth(ff));
                Bt(ff,jj+1) = (tmp(1)^2+tmp(2)^2)^0.5;
                tmp = mo_getb(xYoke(ff),yYoke(ff));
                By(ff,jj+1) = (tmp(1)^2+tmp(2)^2)^0.5;
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
                mo_clearcontour();
                mo_addcontour(xS,yS);
                mo_addcontour(xE,yE);
                mo_bendcontour(angStp,0.5)
                tmp = mo_lineintegral(3);
                fx = tmp(1);
                fy = tmp(2);
                fm = (fx^2+fy^2);
                fa = atan2(fy,fx);
                Fr(ff,jj+1) = fm*cos(fa-angRef(ff)*pi/180);
                Ft(ff,jj+1) = fm*sin(fa-angRef(ff)*pi/180);
                mo_clearcontour();
            end
            
            SOL.Fr = Fr;
            SOL.Ft = Ft;
            
        case {'idemag','idemagmap','demagArea'}
            if (jj == nsim)
                Br = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Br,per.tempPP);
                Bd = interp1(mat.LayerMag.temp.temp,mat.LayerMag.temp.Bd,per.tempPP);
                mo_showdensityplot(1,0,Bd,Br,'bmag');
                
                Bmin = Br;
                
                % get the PMs flux density
                EleNo = mo_numelements;               % Number of mesh elements
                pos = zeros(EleNo,1);                 % Mesh elements centroid coordinates as complex number
                groNo = zeros(EleNo,1);               % group number
                area = zeros(EleNo,1);
                vert = zeros(EleNo,3);
                for ee = 1:EleNo
                    elm = mo_getelement(ee);
                    pos(ee) = elm(4)+j*elm(5);
                    groNo(ee) = elm(7);
                    area(ee) = elm(6);
                    for vv=1:3
                        tmpV=mo_getnode(elm(vv));
                        vert(ee,vv)=tmpV(1)+j*tmpV(2);
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
                    tmp=mo_getpointvalues(real(pos(EleOK(ee))),imag(pos(EleOK(ee))));
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
                    
                    motFolder=[pathname 'critical machines\'];
                    mkdir(motFolder);
                    
                    caseName=[motFolder 'pos_' int2str(theta/(0.5*360/(6*geo.q))) '_cur_' int2str(iAmp) '_temp_' int2str(per.tempPP)];
                    
                    mi_saveas([caseName '.fem']);
                else
                    dPM=0;
                    VolDem=0;
                    VolTot=den;
                end
                if jj>1
                    xyDemagTmpC=xyDemagTmpC.*exp(-j*(thetaPark(jj) - thetaPark(jj-1))/p*pi/180); % rotation in zero position
                    xyDemagTmpV=xyDemagTmpV.*exp(-j*(thetaPark(jj) - thetaPark(jj-1))/p*pi/180); % rotation in zero position
                end
                SOL.xyDemagTmpC=xyDemagTmpC;
                SOL.xyDemagTmpV=xyDemagTmpV;
                SOL.xyDemagTmpB=xyDemagTmpB;
                
                SOL.Bmin      = Bmin;
                SOL.dPM       = dPM;
                SOL.VolDem    = VolDem;
                SOL.VolTot    = VolTot;
            end
        case {'singtIron','singmIron'} % simulation with iron loss
            if jj==1 % store the mesh information
                EleNo = mo_numelements;               % Number of mesh elements
                pos   = zeros(1,EleNo);               % Mesh elements centroid coordinates as complex number
                groNo = zeros(1,EleNo);               % group number
                area  = zeros(1,EleNo);
                %                 vert  = zeros(EleNo,3);
                for ee = 1:EleNo
                    elm = mo_getelement(ee);
                    pos(ee) = elm(4)+j*elm(5);
                    groNo(ee) = elm(7);
                    area(ee) = elm(6);
                end
                EleIn = 1:1:EleNo;
                EleOK = [EleIn(groNo==12) EleIn(groNo==22) EleIn(groNo>200)]; % select just the element of stator iron (group=12), rotor iron (group=22) and PMs (group=200)
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
                tmp = mo_getpointvalues(real(pos(ee)),imag(pos(ee)));
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

mo_close, mi_close
closefemm



