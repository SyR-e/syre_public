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

    [~,geo] = calc_endTurnLength(geo);
    per = calc_i0(geo,per,mat);

    %     if any(strcmp(geo.OBJnames,'Fdq0'))
    %         per0 = per;
    %         per0.overload = 0;
    %         per0.gamma = 0;
    %         per0.nsim_singt = 1;
    %     end

else
    % post proc or FEMM simulation (existing geometry)
    copyfile(filenameIn,[pathname filename]); % copy .fem in the temporary folder
end

mat.LayerMag.Br = per.BrPP;
mat.LayerMag.Hc = per.BrPP/(4e-7*pi*mat.LayerMag.mu);

flagSim = 1;
if ~isempty(RQ)
    if per.MechStressOptCheck
        simSetup.evalSpeed = geo.nmax;
        simSetup.meshSize  = 'coarse';
        simSetup.flagFull  = 0;
        simSetup.shaftBC   = 1;
        warning('off')
        %     [structModel] = syre2pde(geo,mat,simSetup);
        simSetup.filename = filename;
        simSetup.pathname = pathname;
        [out.structModel] = femm2pde(geo,mat,simSetup);
        [out.sVonMises,~,out.structModel] = calcVonMisesStress(out.structModel);
        [outMech] = eval_maxStress(out.structModel,out.sVonMises,geo,mat);
        figure
        figSetting
        pdeplot(out.structModel)
        saveas(gcf,[pathname 'mechMesh.fig']);
        close
        warning('on')
        
        if sum([outMech.nPointOverRad outMech.nPointOverTan])>100
            flagSim=0;
        end
        %     if any(outMech.stress_T>mat.Rotor.sigma_max*10^6)  || any(outMech.stress_R>mat.Rotor.sigma_max*10^6)
        %         flagSim = 0;
        %     end
        %         if any(outMech.stress_R>mat.Rotor.sigma_max*10^6)
        %             flagSim = 0;
        %         end
%         if outMech.MaxStress>mat.Rotor.sigma_max*10^6
%             flagSim=0;
%         end
    end
end

if (strcmp(eval_type,'idemag')||strcmp(eval_type,'idemagmap')||strcmp(eval_type,'demagArea'))
    flagSim = 0;
    flagDemag = 1;
else
    flagDemag = 0;
end

if flagSim
    if strcmp(geo.RotType,'IM')
        [SOL] = simulate_FOC_IM(geo,per,mat,eval_type,pathname,filename);
    else
        [SOL] = simulate_xdeg(geo,per,mat,eval_type,pathname,filename);
    end
    % standard results
    out.id     = mean(SOL.id);                                      % [A]
    out.iq     = mean(SOL.iq);                                      % [A]
    out.fd     = mean(SOL.fd);                                      % [Vs]
    out.fq     = mean(SOL.fq);                                      % [Vs]
    out.T      = mean(SOL.T);                                       % [Nm]
    out.dT     = std(SOL.T);                                        % [Nm]
    out.dTpu   = std(SOL.T)/out.T;                                  % [pu]
    out.dTpp   = max(SOL.T)-min(SOL.T);                             % [Nm]
    out.IPF    = sin(atan2(out.iq,out.id)-atan2(out.fq,out.fd));
    out.We     = mean(SOL.we);                                      % [J]
    out.Wc     = mean(SOL.wc);                                      % [J]
    out.SOL    = SOL;

    % check Torque sign
    if sign(out.T)~=sign(out.fd*out.iq-out.fq*out.id)
        out.T = -out.T;
        out.SOL.T = -out.SOL.T;
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
        out.ppm_RF        = sum(sum(SOL.ppm_RF))*(2*geo.p/geo.ps);
        out.ppm_noRF      = sum(sum(SOL.ppm_noRF))*(2*geo.p/geo.ps);
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
            SOL = rmfield(SOL,'ppm_RF');
            SOL = rmfield(SOL,'ppm_noRF');
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

if flagDemag
    SOL = simulate_xdeg(geo,per,mat,eval_type,pathname,filename);
    
    out.id   = mean(SOL.id);
    out.iq   = mean(SOL.iq);
    out.SOL  = SOL;
    out.Bmin = min(SOL.Bmin);
    out.dPM  = max(SOL.dPM);
end


if ~isempty(RQ)     % MODE optimization (RQ geometry)

    % Cost functions
    cost = zeros(1,length(geo.OBJnames));
    temp1 = 1;
    % Torque
    if strcmp(geo.OBJnames{temp1},'Torque')
        cost(temp1) = -out.T;
        temp1 = temp1+1;
    end
    % Torque Ripple
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'TorRip')
        %         cost(temp1) = out.dTpu*100;
        cost(temp1) = out.dTpp;
        temp1 = temp1+1;
    end
    % Copper Mass
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'MassCu')
        if flagSim
            cost(temp1) = calcMassCu(geo,mat);
            temp1=temp1+1;
        else
            cost(temp1) = 10^50;
        end
    end
    % PM Mass
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'MassPM')
        if flagSim
            cost(temp1) = calcMassPM(geo,mat);
            temp1=temp1+1;
        else
            cost(temp1) = 10^50;
        end
    end

    % Power Factor
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'PF')
        cost(temp1) = -out.IPF;
        temp1=temp1+1;
    end

    % No Load flux
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'Fdq0')
        if flagSim
            per0 = per;
            per0.overload = 0;
            per0.gamma = 0;
            per0.nsim_singt = 1;
            per0.nsim_MOOA = 1;
            [SOL0] = simulate_xdeg(geo,per0,mat,eval_type,pathname,filename);
            cost(temp1) = abs(SOL0.fd+j*SOL0.fq);
        else
            cost(temp1) = 10^50;
        end
        temp1 = temp1+1;
    end

     % Structural properties
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'MechStress')
        cost(temp1) = max([outMech.sigmaRadMax outMech.sigmaTanMax])/1e6;
        temp1=temp1+1;
    end

    % penalize weak solutions
    for ii = 1:length(cost)
        if cost(ii)>per.objs(ii,1) && per.objs(ii,3)==0
            if per.objs(ii,1)>0
                cost(ii)=cost(ii)*10;  % minimization problem
            else
                cost(ii)=cost(ii)*0.1; % maximization problem
            end
        end

        if ((cost(ii)<per.objs(ii,1)-per.objs(ii,3)) || (cost(ii)>per.objs(ii,1)+per.objs(ii,3))) && per.objs(ii,3)
            if per.objs(ii,1)>0
                cost(ii)=cost(ii)*10;  % minimization problem
            else
                cost(ii)=cost(ii)*0.1; % maximization problem
            end
        end
    end

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

