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

[thisfilepath,pathname]=createTempDir();

if ~isempty(RQ)
    
    % MODE optimization (RQ geometry)
    RQ=RQ';
    geo.pathname=pwd();
    
    options.iteration=0;
    options.currentgen=1;
    options.PopulationSize=1;
    
    if strcmp(eval_type,'MO_OA')
        options.iteration=options.iteration+1;
        iteration=options.iteration;
        populationSize=options.PopulationSize;
        generation=floor(iteration/populationSize)+1;
        options.currentgen=generation;
        RQ % debug .. when syre crashes it is useful to have visibility of last RQ
    end
    
    [geo,gamma,mat] = interpretRQ(RQ,geo,mat);
    per.gamma=gamma;
    
    [geo,mat] = draw_motor_in_FEMM(geo,mat,pathname,filename);
    
    per.i0 = calc_io(geo,per);
    
else
    
    % post proc or FEMM simulation (existing geometry)
    copyfile(filenameIn,[pathname filename]); % copy .fem in the temporary folder
end

% evaluates the candidate machine (depending on eval_type)
% fem = dimMesh(geo,eval_type);

mat.LayerMag.Br = per.BrPP;
mat.LayerMag.Hc = per.BrPP/(4e-7*pi*mat.LayerMag.mu);

[SOL] = simulate_xdeg(geo,per,mat,eval_type,pathname,filename);

% standard results
out.id     = mean(SOL.id);                                      % [A]
out.iq     = mean(SOL.iq);                                      % [A]
out.fd     = mean(SOL.fd);                                      % [Vs]
out.fq     = mean(SOL.fq);                                      % [Vs]
out.T      = mean(SOL.T);                                       % [Nm]
out.dT     = std(SOL.T);                                        % [Nm]
out.dTpu   = std(SOL.T)/out.T;                                  % [pu]
out.dTpp   = max(SOL.T)-min(SOL.T);                             % [Nm]
out.IPF    = sin(atan(out.iq./out.id)-atan(out.fq./out.fd));
%out.MassPM = mean(SOL.VolPM)*mat.LayerMag.kgm3;                 % [kg]
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
    out.Pfes_h = sum(sum(SOL.psh))*(2*geo.p/geo.ps);
    out.Pfes_c = sum(sum(SOL.psc))*(2*geo.p/geo.ps);
    out.Pfer_h = sum(sum(SOL.prh))*(2*geo.p/geo.ps);
    out.Pfer_c = sum(sum(SOL.prc))*(2*geo.p/geo.ps);
    out.Ppm    = sum(sum(SOL.ppm))*(2*geo.p/geo.ps);
    out.Pfe    = out.Pfes_h + out.Pfes_c + out.Pfer_h + out.Pfer_c;
    out.velDim = per.EvalSpeed;
    
    if strcmp(eval_type,'singmIron')
        % remove all the debug data from SOL, to avoid excessive data size
        SOL = rmfield(SOL,'psh');
        SOL = rmfield(SOL,'psc');
        SOL = rmfield(SOL,'prh');
        SOL = rmfield(SOL,'prc');
        SOL = rmfield(SOL,'ppm');
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
        cost(temp1) = calcMassCu(geo,mat);
        temp1=temp1+1;
    end
    % PM Mass
    if temp1<=length(geo.OBJnames) && strcmp(geo.OBJnames{temp1},'MassPM')
        cost(temp1) = calcMassPM(geo,mat);
    end
    
    % penalize weak solutions
    for j = 1:length(cost)
        if cost(j)>per.objs(j,1)
            if per.objs(j,1)>0
                cost(j)=cost(j)*10;  % minimization problem
            else
                cost(j)=cost(j)*0.1; % maximization problem
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
        save([pathname strrep(filename,'.fem','.mat')],'geo','cost','per','dataSet','mat');
    end
else
    cost = [];
    save([pathname strrep(filename,'.fem','.mat')],'geo','out','mat');   % save geo and out
end

