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

function [cost,geo,mat,out,pathname] = FEMMfitness(RQ,geo,per,mat,eval_type,filenameIn)

% FEMMfitness
% runs FEMM simulation from file (existing machine) or file+RQ (MODE, RQ is the set of inputs)
% - creates a temp dir
% - moves or draw the machine into the temp dir
% - simulate_xdeg(eval_type)
% - calc out from SOL
% - saves the .mat file into the temp dir

[~,filename,ext] = fileparts(filenameIn);
filename = [filename ext]; % fem file name

[~,pathname]=createTempDir();


if ~isempty(RQ)

    % MODE optimization (RQ geometry)
    RQ=RQ';
    geo.pathname=pwd();

    %     options.iteration=0;
    %     options.currentgen=1;
    %     options.PopulationSize=1;

    if strcmp(eval_type,'MO_OA')
%         RQ % debug .. when syre crashes it is useful to have visibility of last RQ
    end

    [geo,gamma,mat] = interpretRQ(RQ,geo,mat);
    per.gamma=gamma;

    [geo,mat] = draw_motor_in_FEMM(geo,mat,pathname,filename);
    if geo.XFEMMsimulation
        closefemm;
    end

    [~,geo] = calc_endTurnLength(geo);
    [~,geo] = calc_endTurnFieldLength(geo);
 
    %     flag_OptCurrConst = 1;
    switch per.flag_OptCurrConst
        case 0 % constant thermal loading
            per.loss = NaN;
            per.J    = NaN;
        case 1 % constant current density
            per.kj   = NaN;
            per.Loss = NaN;
        case 2 % constant current
            per.kj   = NaN;
            per.Loss = NaN;
            per.J    = NaN;
    end
    per = calc_i0(geo,per,mat);
    % warning('Define the ratio between stator and rotor current density')
    per.Jf = per.J*per.JfPU;
    per = calc_if(geo,per,mat);
    per.if = per.if0;

    %     if any(strcmp(geo.OBJnames,'Fdq0'))
    %         per0 = per;
    %         per0.overload = 0;
    %         per0.gamma = 0;
    %         per0.nsim_singt = 1;
    %     end

    save([pathname filename(1:end-4) '.mat'],'geo','per','mat')

else
    % post proc or FEMM simulation (existing geometry)
    copyfile(checkPathSyntax(filenameIn),checkPathSyntax([pathname filename])); % copy .fem in the temporary folder
end

mat.LayerMag.Br = per.BrPP;
mat.LayerMag.Hc = per.BrPP/(4e-7*pi*mat.LayerMag.mu);

flagSim = 1;
if ~isempty(RQ)
%     if per.MechStressOptCheck
%         simSetup.evalSpeed = geo.nmax;
%         simSetup.meshSize  = 'coarse';
%         simSetup.flagFull  = 0;
%         simSetup.shaftBC   = 1;
%         warning('off')
%         %     [structModel] = syre2pde(geo,mat,simSetup);
%         simSetup.filename = filename;
%         simSetup.pathname = pathname;
%         [out.structModel] = femm2pde(geo,mat,simSetup);
%         [out.sVonMises,R,out.structModel] = calcVonMisesStress(out.structModel);
%         [outMech] = eval_maxStress(out.structModel,out.sVonMises,R,geo,mat);
%         figure
%         figSetting
%         pdeplot(out.structModel)
%         saveas(gcf,[pathname 'mechMesh.fig']);
%         close
%         warning('on')
% 
%         if sum([outMech.nPointOverRad outMech.nPointOverTan])>100
%             flagSim=0;
%         end
%         %     if any(outMech.stress_T>mat.Rotor.sigma_max*10^6)  || any(outMech.stress_R>mat.Rotor.sigma_max*10^6)
%         %         flagSim = 0;
%         %     end
%         %         if any(outMech.stress_R>mat.Rotor.sigma_max*10^6)
%         %             flagSim = 0;
%         %         end
% %         if outMech.MaxStress>mat.Rotor.sigma_max*10^6
% %             flagSim=0;
% %         end
%     end
% 
    if geo.statorYokeDivision
        % GalFerContest - structural simulation
        simSetup.evalSpeed = geo.nmax;
        simSetup.meshSize  = 'FEMM original';
        simSetup.flagFull  = 0;
        simSetup.shaftBC   = 1;
        simSetup.meshShaft = 0;
        warning('off')
        simSetup.filename = filename;
        simSetup.pathname = pathname;
        warning('on')
        [structModel] = femm2pde(geo,mat,simSetup);
        [sVonMises,R,structModel] = calcVonMisesStress(structModel,0);
        [outStruct] = eval_maxStress(structModel,sVonMises,R,geo,mat);
        % GalFerContest - 3D thermal simulation
        if exist([pathname filename(1:end-4) '.fem'],'file')&&~exist([pathname filename(1:end-4) '.ans'],'file')
            openfemm(1);
            opendocument([pathname filename(1:end-4) '.fem']);
            mi_createmesh;
            mi_analyze(1);
            mi_loadsolution;
            closefemm();
        end
        [tempCuAvg,~,tempCuMax] = solve_therm(geo,per,mat,[pathname filename(1:end-4)]);
    end

end

if (strcmp(eval_type,'idemag')||strcmp(eval_type,'idemagmap')||strcmp(eval_type,'demagArea'))
    flagSim = 0;
    flagDemag = 1;
else
    flagDemag = 0;
end

if ~isfield(per,'offset')
    per.offset = 0;
end

% tic
if flagSim
    if strcmp(geo.RotType,'IM')
        [SOL] = simulate_FOC_IM(geo,per,mat,eval_type,pathname,filename);
    else
        runSimFlag=1;
        if geo.XFEMMsimulation
            if strcmp(geo.RotType,'EESM')
                if per.overload==0 && per.if==0
                    runSimFlag=0;
                    SOL.th = per.offset:per.delta_sim_singt/per.nsim_singt:per.delta_sim_singt+per.offset;
                    SOL.ff = zeros(1,per.nsim_singt);
                    SOL.T  = zeros(1,per.nsim_singt);
                    SOL.we = zeros(1,per.nsim_singt);
                    SOL.wc = zeros(1,per.nsim_singt);
                    if(isnan(geo.win.n3phase))
                        SOL.id1 = zeros(1,per.nsim_singt);
                        SOL.iq1 = zeros(1,per.nsim_singt);
                        SOL.id3 = zeros(1,per.nsim_singt);
                        SOL.iq3 = zeros(1,per.nsim_singt);
                        SOL.io = zeros(1,per.nsim_singt);
                        SOL.if = zeros(1,per.nsim_singt);
                        SOL.fd1 = zeros(1,per.nsim_singt);
                        SOL.fq1 = zeros(1,per.nsim_singt);
                        SOL.fd3 = zeros(1,per.nsim_singt);
                        SOL.fq3 = zeros(1,per.nsim_singt);
                        SOL.fo = zeros(1,per.nsim_singt);

                        SOL.ia = zeros(1,per.nsim_singt);
                        SOL.ib = zeros(1,per.nsim_singt);
                        SOL.ic = zeros(1,per.nsim_singt);
                        SOL.id = zeros(1,per.nsim_singt);
                        SOL.ie = zeros(1,per.nsim_singt);
                        SOL.fa = zeros(1,per.nsim_singt);
                        SOL.fb = zeros(1,per.nsim_singt);
                        SOL.fc = zeros(1,per.nsim_singt);
                        SOL.fd = zeros(1,per.nsim_singt);
                        SOL.fe = zeros(1,per.nsim_singt);
                    else
                        SOL.id = zeros(1,per.nsim_singt);
                        SOL.iq = zeros(1,per.nsim_singt);
                        SOL.io = zeros(1,per.nsim_singt);
                        SOL.if = zeros(1,per.nsim_singt);
                        SOL.fd = zeros(1,per.nsim_singt);
                        SOL.fq = zeros(1,per.nsim_singt);
                        SOL.fo = zeros(1,per.nsim_singt);

                        SOL.ia = zeros(geo.win.n3phase,per.nsim_singt);
                        SOL.ib = zeros(geo.win.n3phase,per.nsim_singt);
                        SOL.ic = zeros(geo.win.n3phase,per.nsim_singt);
                        SOL.fa = zeros(geo.win.n3phase,per.nsim_singt);
                        SOL.fb = zeros(geo.win.n3phase,per.nsim_singt);
                        SOL.fc = zeros(geo.win.n3phase,per.nsim_singt);
                    end
                    if contains(eval_type,'Iron')
                        SOL.psh          = 0;
                        SOL.psc          = 0;
                        SOL.prh          = 0;
                        SOL.prc          = 0;
                        SOL.ppm          = 0;
                        SOL.ppm_no3D     = 0;
                        SOL.ppm_noRFno3D = 0;
                        SOL.ppm_PM       = 0;
                        SOL.ppm_PM       = 0;
                        SOL.freq         = 0;
                        SOL.bs           = 0;
                        SOL.br           = 0;
                        SOL.am           = 0;
                        SOL.Jm           = 0;
                        SOL.pos          = 0;
                        SOL.vol          = 0;
                        SOL.groNo        = 0;
                    end
                end
            elseif ~(strcmp(geo.RotType,'SPM')||strcmp(geo.RotType,'SPM-Hallbach'))&&(strcmp(mat.LayerMag.MatName,'Air')||max(geo.PMdim(:))==0)
                if per.overload==0
                    runSimFlag=0;
                    SOL.th = per.offset:per.delta_sim_singt/per.nsim_singt:per.delta_sim_singt+per.offset;
                    SOL.T  = zeros(1,per.nsim_singt);
                    SOL.we = zeros(1,per.nsim_singt);
                    SOL.wc = zeros(1,per.nsim_singt);
                    if(isnan(geo.win.n3phase))
                        SOL.id1 = zeros(1,per.nsim_singt);
                        SOL.iq1 = zeros(1,per.nsim_singt);
                        SOL.id3 = zeros(1,per.nsim_singt);
                        SOL.iq3 = zeros(1,per.nsim_singt);
                        SOL.io = zeros(1,per.nsim_singt);
                        SOL.fd1 = zeros(1,per.nsim_singt);
                        SOL.fq1 = zeros(1,per.nsim_singt);
                        SOL.fd3 = zeros(1,per.nsim_singt);
                        SOL.fq3= zeros(1,per.nsim_singt);
                        SOL.fo = zeros(1,per.nsim_singt);

                        SOL.ia = zeros(1,per.nsim_singt);
                        SOL.ib = zeros(1,per.nsim_singt);
                        SOL.ic = zeros(1,per.nsim_singt);
                        SOL.id = zeros(1,per.nsim_singt);
                        SOL.ie = zeros(1,per.nsim_singt);
                        SOL.fa = zeros(1,per.nsim_singt);
                        SOL.fb = zeros(1,per.nsim_singt);
                        SOL.fc = zeros(1,per.nsim_singt);
                        SOL.fd = zeros(1,per.nsim_singt);
                        SOL.fe = zeros(1,per.nsim_singt);
                    else
                        SOL.id = zeros(1,per.nsim_singt);
                        SOL.iq = zeros(1,per.nsim_singt);
                        SOL.io = zeros(1,per.nsim_singt);
                        SOL.fd = zeros(1,per.nsim_singt);
                        SOL.fq = zeros(1,per.nsim_singt);
                        SOL.fo = zeros(1,per.nsim_singt);

                        SOL.ia = zeros(geo.win.n3phase,per.nsim_singt);
                        SOL.ib = zeros(geo.win.n3phase,per.nsim_singt);
                        SOL.ic = zeros(geo.win.n3phase,per.nsim_singt);
                        SOL.fa = zeros(geo.win.n3phase,per.nsim_singt);
                        SOL.fb = zeros(geo.win.n3phase,per.nsim_singt);
                        SOL.fc = zeros(geo.win.n3phase,per.nsim_singt);
                    end
                    if contains(eval_type,'Iron')
                        SOL.psh          = 0;
                        SOL.psc          = 0;
                        SOL.prh          = 0;
                        SOL.prc          = 0;
                        SOL.ppm          = 0;
                        SOL.ppm_no3D     = 0;
                        SOL.ppm_noRFno3D = 0;
                        SOL.ppm_PM       = 0;
                        SOL.ppm_PM       = 0;
                        SOL.freq         = 0;
                        SOL.bs           = 0;
                        SOL.br           = 0;
                        SOL.am           = 0;
                        SOL.Jm           = 0;
                        SOL.pos          = 0;
                        SOL.vol          = 0;
                        SOL.groNo        = 0;
                    end
                end
            end
        end
        if runSimFlag
            [SOL] = simulate_xdeg(geo,per,mat,eval_type,pathname,filename);
        end
    end
    % standard results
    
    out.T      = mean(SOL.T);                                       % [Nm]
    out.dT     = std(SOL.T);                                        % [Nm]
    out.dTpu   = std(SOL.T)/out.T;                                  % [pu]
    out.dTpp   = max(SOL.T)-min(SOL.T);                             % [Nm]
    
    if(isnan(geo.win.n3phase))
        if(per.custom_act)
            out.T = SOL.T;
            out.ia = SOL.ia;
            out.ib = SOL.ib;
            out.ic = SOL.ic;
            out.id = SOL.id;
            out.ie = SOL.ie;
        end
    end
    
    out.We     = mean(SOL.we);                                      % [J]
    out.Wc     = mean(SOL.wc);                                      % [J]
    out.SOL    = SOL;

    if(isnan(geo.win.n3phase))

        out.id1     = mean(SOL.id1(:));                                   % [A]
        out.iq1     = mean(SOL.iq1(:));                                   % [A]
        out.id3     = mean(SOL.id3(:));                                   % [A]
        out.iq3     = mean(SOL.iq3(:));                                   % [A
        out.io     = mean(SOL.io(:));                                   % [A]
        out.fd1     = mean(SOL.fd1);                                      % [Vs]
        out.fq1     = mean(SOL.fq1);                                      % [Vs]
        out.fd3     = mean(SOL.fd3);                                      % [Vs]
        out.fq3     = mean(SOL.fq3);                                      % [Vs]
        out.fo     = mean(SOL.fo);                                      % [Vs]

        out.IPF    = sin(atan2(out.iq1,out.id1)-atan2(out.fq1,out.fd1));

    else
        out.id     = mean(SOL.id(:));                                   % [A]
        out.iq     = mean(SOL.iq(:));                                   % [A]
        out.io     = mean(SOL.io(:));                                   % [A]
        out.fd     = mean(SOL.fd);                                      % [Vs]
        out.fq     = mean(SOL.fq);                                      % [Vs]
        out.fo     = mean(SOL.fo);                                      % [Vs]

        out.IPF    = sin(atan2(out.iq,out.id)-atan2(out.fq,out.fd));
    end

    if isfield(out.SOL,'if')
       out.if = mean(SOL.if(:));
       out.ff = mean(SOL.ff(:));
    % else
    %     out.if = NaN;
    %     out.ff = NaN;
    end
    
    if(~isnan(geo.win.n3phase))
        if geo.win.n3phase>1
            for ff=1:geo.win.n3phase
                th = out.SOL.th;
    
                fa = out.SOL.fa(ff,:);
                fb = out.SOL.fb(ff,:);
                fc = out.SOL.fc(ff,:);
                % fdq = abc2dq(fa,fb,fc,(th+geo.th0(ff)-geo.th0(1))*pi/180);
                fdq0 = abc2dq0([fa;fb;fc],(th+geo.th0(ff)-geo.th0(1))*pi/180);
                out.SOL.sets(ff).fa = fa;
                out.SOL.sets(ff).fb = fb;
                out.SOL.sets(ff).fc = fc;
                out.SOL.sets(ff).fd = fdq0(1,:);
                out.SOL.sets(ff).fq = fdq0(2,:);
                out.SOL.sets(ff).f0 = fdq0(3,:);
    
                ia = out.SOL.ia(ff,:);
                ib = out.SOL.ib(ff,:);
                ic = out.SOL.ic(ff,:);
                % idq = abc2dq(ia,ib,ic,(th+geo.th0(ff)-geo.th0(1))*pi/180);
                idq0 = abc2dq0([ia;ib;ic],(th+geo.th0(ff)-geo.th0(1))*pi/180);
                out.SOL.sets(ff).ia = ia;
                out.SOL.sets(ff).ib = ib;
                out.SOL.sets(ff).ic = ic;
                out.SOL.sets(ff).id = idq0(1,:);
                out.SOL.sets(ff).iq = idq0(2,:);
                out.SOL.sets(ff).i0 = idq0(3,:);
    
                out.sets(ff).id = mean(out.SOL.sets(ff).id);
                out.sets(ff).iq = mean(out.SOL.sets(ff).iq);
                out.sets(ff).i0 = mean(out.SOL.sets(ff).i0);
                out.sets(ff).fd = mean(out.SOL.sets(ff).fd);
                out.sets(ff).fq = mean(out.SOL.sets(ff).fq);
                out.sets(ff).f0 = mean(out.SOL.sets(ff).f0);
            end
        end
    end

    % check Torque sign
    if(isnan(geo.win.n3phase))
        if sign(out.T)~=sign(out.fd1*out.iq1-out.fq1*out.id1)
            out.T = -out.T;
            out.SOL.T = -out.SOL.T;
            % warning('Torque sign correction')
        end
    else
        if sign(out.T)~=sign(out.fd*out.iq-out.fq*out.id)
            out.T = -out.T;
            out.SOL.T = -out.SOL.T;
            % warning('Torque sign correction')
        end
    end

    if isfield(SOL,'F')
        out.F=mean(SOL.F);
    end

    if isfield(SOL,'psh')
        out.Pfes_h        = sum(sum(SOL.psh))*(2*geo.p/geo.ps);
        out.Pfes_c        = sum(sum(SOL.psc))*(2*geo.p/geo.ps);
        out.Pfer_h        = sum(sum(SOL.prh))*(2*geo.p/geo.ps);
        out.Pfer_c        = sum(sum(SOL.prc))*(2*geo.p/geo.ps);
        out.Ppm           = sum(sum(SOL.ppm))*(2*geo.p/geo.ps);
        out.ppm_no3D      = sum(sum(SOL.ppm_no3D))*(2*geo.p/geo.ps);
        out.ppm_noRFno3D  = sum(sum(SOL.ppm_noRFno3D))*(2*geo.p/geo.ps);
        out.Ppm_breakdown = SOL.ppm_PM*(2*geo.p/geo.ps);
        out.Pfe           = out.Pfes_h + out.Pfes_c + out.Pfer_h + out.Pfer_c;
        out.velDim        = per.EvalSpeed;

        if strcmp(eval_type,'singmIron')
            % remove all the debug data from SOL, to avoid excessive data size
            SOL = rmfield(SOL,'psh');
            SOL = rmfield(SOL,'psc');
            SOL = rmfield(SOL,'prh');
            SOL = rmfield(SOL,'prc');
            SOL = rmfield(SOL,'ppm');
            %SOL = rmfield(SOL,'ppm_RF');
            %SOL = rmfield(SOL,'ppm_noRF');
            SOL = rmfield(SOL,'ppm_PM');
            SOL = rmfield(SOL,'freq');
            SOL = rmfield(SOL,'bs');
            SOL = rmfield(SOL,'br');
            SOL = rmfield(SOL,'am');
            SOL = rmfield(SOL,'Jm');
            SOL = rmfield(SOL,'pos');
            SOL = rmfield(SOL,'vol');
            SOL = rmfield(SOL,'groNo');
            out.SOL = SOL;
        end
    end

    if strcmp(geo.RotType,'IM')
        out.IM = SOL.IM;
    end

else
    out.T    = 0;
    out.dTpp = 10^50;
    out.IPF  = -10^50;
end
% toc

if flagDemag
    SOL = simulate_xdeg(geo,per,mat,eval_type,pathname,filename);
    
    if(isnan(geo.win.n3phase))
        out.id1   = mean(SOL.id1);
        out.iq1   = mean(SOL.iq1);
        out.id3   = mean(SOL.id3);
        out.iq3   = mean(SOL.iq3);
    else
        out.id   = mean(SOL.id);
        out.iq   = mean(SOL.iq);
    end
    
    out.SOL  = SOL;
    out.Bmin = min(SOL.Bmin);
    out.dPM  = max(SOL.dPM);
end

% Variables necessary for MTPA calculation
flagMTPA = 0;

maxIter = 4;
gammaStep = 2;
direction = 0;

if isfield(per,'if0')
    per.if = per.if0;
else
    per.if = 0;
end
if ~isempty(RQ)
    RQ(end) = 90;
end

if any(strcmp ('gamma', geo.RQnames))
    flagMTPA = 0;
else
    flagMTPA = 1;
end


if ~isempty(RQ)     % MODE optimization (RQ geometry)

    % Cost functions
    cost = nan(1,length(geo.OBJnames));
    %temp1 = 1;
    % Torque
    %if strcmp(geo.OBJnames{temp1},'Torque')
    if per.objs_check(1)
        if ~flagMTPA
            cost(1) = -out.T;
        else
            % aggiungere ricerca MTPA
            gamma0 = RQ(end);
            jj = 1;
            done = 0;
            TVect    = nan(1,maxIter);
            gVect    = nan(1,maxIter);
            dTppVect = nan(1,maxIter);
            idqVect  = nan(1,maxIter);
            fdqVect  = nan(1,maxIter);

            perTmp = per;
            TmpSOL_old = [];
            TmpSOL_new = [];


             while ~done
                    if jj==1
                        gammaSim = gamma0;
                    elseif jj==2
                        gammaSim = gamma0+gammaStep;
                    elseif jj==3
                        gammaSim = gamma0-gammaStep;
                    else
                        gammaSim = gammaSim+direction*gammaStep;
                    end

                     RQ(end) = gammaSim;
                     TmpSOL = simulate_xdeg(geo,perTmp,mat,eval_type,pathname,filename);
                     TVect(jj)    = mean(TmpSOL.T);
                     gVect(jj)    = gammaSim;
                     dTppVect(jj) = max(TmpSOL.T) - min(TmpSOL.T);
                     idqVect(jj)  = mean(TmpSOL.id)+j*mean(TmpSOL.iq);
                     fdqVect(jj)  = mean(TmpSOL.fd)+j*mean(TmpSOL.fq);
                     
                    % To save the last 2 results of the iteractive process
                     TmpSOL_old = TmpSOL_new;
                     TmpSOL_new = TmpSOL;


                     if jj==3
                            [~,index] = max(TVect,[],'omitnan');
                            if index==1
                               done=1;
                            elseif index==2
                                   direction=+1;
                            else
                                   direction=-1;
                            end
                      elseif jj>3
                          if TVect(jj)<TVect(jj-1)
                             done=1;
                          end
                     end

                     if jj==maxIter
                          done=1;
                     end

                     jj = jj+1;
                     disp(['Simulation ' int2str(jj-1) ' done'])
             end

             [~,index] = max(TVect,[],'omitnan');

            

             % Output data
              OUT.geo   = geo;
              OUT.per   = perTmp;
              OUT.mat   = mat;
              OUT.T     = TVect(index);
              OUT.dTpp  = dTppVect(index);
              OUT.gamma = gVect(index);
              OUT.idq   = idqVect(index);
              OUT.fdq   = fdqVect(index);
              OUT.RQ    = RQ;
              %OUT.nFEMM = nFEMM;
              OUT.Pjs   = 3/2*per.Rs*abs(OUT.idq)^2;

              % Identify Rf value depending on rotor geometry
              if strcmp(geo.RotType, 'EESM')
                 Rf=per.Rf;
                 OUT.Pjf   = Rf*out.if;
              else
                  Rf=0;
                  OUT.Pjf = nan;
              end

              if mean(TmpSOL_new.T) > mean(TmpSOL_old.T)
                  OUT.SOL = TmpSOL_new;
              else
                  OUT.SOL = TmpSOL_old;
              end

             out.SOL = OUT.SOL;
             out.id     = real(OUT.idq);                                   % [A]
             out.iq     = imag(OUT.idq);                                   % [A]
             out.fd     = real(OUT.fdq);                                   % [Vs]
             out.fq     = imag(OUT.fdq);                                   % [Vs]
             out.T      = OUT.T;                                           % [Nm]
             out.dTpp   = OUT.dTpp;                                        % [Nm]
             out.gamma = OUT.gamma;                                        % [degrees] MTPA angle

             cost(1) = -mean(OUT.T);
             %temp1 = temp1+1;
   
        end
        
    end

    % Torque Ripple
    %if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'TorRip')
    if per.objs_check(2)
        %         cost(temp1) = out.dTpu*100;
        cost(2) = out.dTpp;
        %temp1 = temp1+1;
    end
    % Copper Mass
    %if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'MassCu')
    if per.objs_check(3)
        if flagSim
            cost(3) = calcMassCu(geo,mat);
            %temp1=temp1+1;
        else
            cost(3) = 10^50;
        end
    end
    % PM Mass
    %if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'MassPM')
    if per.objs_check(4) 
        if flagSim
            cost(4) = calcMassPM(geo,mat);
            %temp1=temp1+1;
        else
            cost(4) = 10^50;
        end
    end

    % Power Factor
    %if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'PF')
    if per.objs_check(5)
        cost(5) = -out.IPF;
        %temp1=temp1+1;
    end

    % No Load flux
    %if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'Fdq0')
    if per.objs_check(6)
        if flagSim
            per0 = per;
            per0.overload = 0;
            per0.gamma = 0;
            per0.nsim_singt = 1;
            per0.nsim_MOOA = 1;
            [SOL0] = simulate_xdeg(geo,per0,mat,eval_type,pathname,filename);
            cost(6) = abs(SOL0.fd+j*SOL0.fq);
        else
            cost(6) = 10^50;
        end
        %temp1 = temp1+1;
    end

    %  % Structural properties
    %if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'MechStress')
    if per.objs_check(7)   

        [geo,~,mat] = interpretRQ(RQ,geo,mat);
        % Set up info for solving structural model
        simSetup0.evalSpeed = geo.nmax;
        simSetup0.meshSize  = 'PDE fine';
        simSetup0.meshShaft = 0;
        simSetup0.flagFull  = 0;
        simSetup0.shaftBC   = 1;
        geo.custom = 0;
        warning('off')
        filemot = [pathname filename];
        [~, filename, ~] = fileparts(filemot);

        % Geometry/points generation
        [geo,mat] = draw_motor_in_FEMM(geo,mat,pathname,filename);
        simSetup0.pathname  = pathname;
        simSetup0.filename  = [filename '.mat'];
    
        % Convertion FEMM - PDE
        [out.structModel] = femm2pde(geo,mat,simSetup0);
        %Solving the PDE model
        [out.sVonMises,R,out.structModel] = calcVonMisesStress(out.structModel);
        % Deformation, clerance and stress calculations
        [outMech] = eval_maxStress(out.structModel,out.sVonMises,R,geo,mat);
    
         if (strcmp(geo.RotType,'Circular')||strcmp(geo.RotType,'Seg'))
             cost(7) = max([outMech.Rad_sigmaTotPrc outMech.Tan_sigmaTotPrc])/1e6; 
        % Deformation magnitude
        %     OUT.MaxDef            = outMech.MaxDef;
        %     % Air-gap clearance
        %     OUT.agclear           = outMech.agclear;
        %     % Maximum stress values
        %     OUT.MaxStress         = outMech.MaxStress/1e6;
        %     OUT.TanRibsStress     = (max(outMech.sigmaTanMax))/1e6;
        %     OUT.RadRibsStress     = (max(outMech.sigmaRadMax))/1e6;
        %     % Percentile stress values
        %     %OUT.MaxStress_prc     = out.Stress.sigmaTotPrc/1e6;
        %     OUT.TanRibsStress_prc = outMech.Tan_sigmaTotPrc/1e6;
        %     OUT.RadRibsStress_prc = outMech.Rad_sigmaTotPrc/1e6;
         else
             cost(7) = outMech.MaxStress/1e6; 
        %     OUT.MaxDef            = outMech.MaxDef;
        %     OUT.agclear           = outMech.agclear;
        %     OUT.MaxStress         = outMech.MaxStress/1e6;
        %     % OUT.MaxStress_prc     = out.Stress.sigmaTotPrc/1e6;
         end
        
        %temp1=temp1+1;
    end
    
    % Stator Joule loss
    %if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'Pjs')
    if per.objs_check(8)
        if flagSim
            cost(8) = 3/2*per.Rs*abs(out.id+j*out.iq)^2;
        else
            cost(8) = 10^50;
        end
        %temp1 = temp1+1;
    end

    % Field Joule loss
    %if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'Pjf')
    if per.objs_check(9)
        if flagSim
            cost(9) = per.Rf*out.if^2;
        else
            cost(9) = 10^50;
        end
        %temp1 = temp1+1;
    end
    

    % Structural GalFerContest
    if geo.statorYokeDivision
        cost = [cost outStruct.sigmaTotPrc/1e6 tempCuMax tempCuAvg];
    else
        % penalize weak solutions
        for ii = 1:length(cost)
            if ~isnan(cost(ii))
                if ii ~= 7 && cost(ii)>per.objs(ii,1) && per.objs(ii,3)==0
                    if per.objs(ii,1)>0
                        cost(ii)=cost(ii)*10;  % minimization problem
                    else
                        cost(ii)=cost(ii)*0.1; % maximization problem
                    end

                elseif ii == 5 && ((cost(ii)<per.objs(ii,1)-per.objs(ii,3)) || (cost(ii)>per.objs(ii,1)+per.objs(ii,3))) %&& per.objs(ii,3) 
                    if per.objs(ii,1)>0
                        cost(ii)=cost(ii)*10;  % minimization problem
                    else
                        cost(ii)=cost(ii)*0.1; % maximization problem
                    end

                elseif ii == 7 % Mechanical penalization. 
                    if cost(ii)>per.objs(ii,1) %&& per.objs(ii,3)==0
                        cost(ii)=cost(ii)*10;
                    end
                end
            end    
        end
    end

    %Rimuovo da cost gli obiettivi non selezionati
    cost = cost(~isnan(cost));

    %     dataSetPath = strcat(thisfilepath,'\dataSet.mat');    %OCT
    load('dataSet.mat');
    geo.RQ = RQ;

    [dataSet] = SaveInformation(geo,mat,dataSet);
    if isoctave()            %OCT
        save('-mat7-binary', strrep(filename,'.fem','.mat'),'geo','cost','per','dataSet','mat');
    else
        save([pathname strrep(filename,'.fem','.mat')],'geo','cost','per','dataSet','mat','out');
    end
else
    cost = [];
    save([pathname strrep(filename,'.fem','.mat')],'geo','out','mat','per');   % save geo and out
end

