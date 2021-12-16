
% Copyright 2019
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

function [dataSet] = PMdesign(dataSet)
% 
% [dataSet] = PMdesign(dataSet)
%

%% 1) save the initial machine with all the selected barriers filled with PM

dataSet0 = dataSet;

dataSet.PMdim(dataSet.PMdim>0) = inf;

pathname = dataSet.currentpathname;
filename = [dataSet.currentfilename(1:end-4) '_PMdesign.mat'];

dataSet = DrawAndSaveMachine(dataSet,filename,pathname);

load([pathname filename]);

%% 2)Simulations #1 and #2 (only if the machine is different from the previous)

Mat0 = [dataSet0.PMdesign.geo.stator;dataSet0.PMdesign.geo.rotor];
Mat  = [geo.stator;geo.rotor];

Mat0(isnan(Mat0)) = 0;
Mat(isnan(Mat))   = 0;

disp('PM design simulations...')

filemot = [pathname filename(1:end-4) '.fem'];
per.nsim_singt      = 10;
per.delta_sim_singt = 60;

if prod(size(Mat0)==size(Mat))
    if prod(prod(Mat0==Mat))
        geoCheck = 1;
    else
        geoCheck = 0;
    end
else
    geoCheck = 0;
end

if geoCheck
    disp('FEA simulations #1 and #2 already done for this geometry')
    Br = dataSet0.PMdesign.Br12;
    fM = dataSet0.PMdesign.fM12;
    
    per = calc_i0(geo,per);
    i0 = per.i0;
    if dataSet0.PMdesign.iq0 == dataSet.CurrPM*i0
        disp('FEA simulation #3 already done for this current')
        iq0 = dataSet.PMdesign.iq0;
        fq0 = dataSet.PMdesign.fq0;
    else
        if strcmp(dataSet.TypeOfRotor,'Vtype')
            per.gamma = 180;
        else
            per.gamma = 90;
        end
        per.overload  = dataSet0.CurrPM;
        per.BrPP      = 0;
        [~,~,~,out,~] = FEMMfitness([],geo,per,mat,'singt',filemot);
        if strcmp(dataSet.TypeOfRotor,'Vtype')
            fq0           = abs(out.fd);
            iq0           = abs(out.id);
        else
            fq0           = abs(out.fq);
            iq0           = abs(out.iq);
        end
        disp('Simulation #3 done')
    end
else
    Br = [0.3 0.6];
    fM = zeros(size(Br));

    for ii=1:length(Br)
        per.gamma     = 0;
        per.overload  = 0;
        per.BrPP      = Br(ii);
        [~,~,~,out,~] = FEMMfitness([],geo,per,mat,'singt',filemot);
        if strcmp(dataSet.TypeOfRotor,'Vtype')
            fM(ii)        = abs(out.fd);            
        else
            fM(ii)        = abs(out.fq);
        end
        disp(['Simulation #' int2str(ii) ' done']);
    end

    if strcmp(dataSet.TypeOfRotor,'Vtype')
        per.gamma = 180;
    else
        per.gamma = 90;
    end
    per.overload  = dataSet0.CurrPM;
    per.BrPP      = 0;
    [~,~,~,out,~] = FEMMfitness([],geo,per,mat,'singt',filemot);
    if strcmp(dataSet.TypeOfRotor,'Vtype')
        fq0           = abs(out.fd);
        iq0           = abs(out.id);
    else
        fq0           = abs(out.fq);
        iq0           = abs(out.iq);
    end
    disp('Simulation #3 done')
end

%% 3) Compute Br' (BrDes)

% figure()
% figSetting()
% xlabel('$B_r$ [$T$]')
% ylabel('$\lambda_{PM}$ [$Vs$]')
% plot(Br,fM,'kx','DisplayName','Sim 1 \& 2')
% xPlot = [Br(1)-(Br(2)-Br(1))/2 Br(2)+(Br(2)-Br(1))/2];
% m     = (fM(2)-fM(1))/(Br(2)-Br(1));
% q     = fM(1)-m*Br(1);
% plot(xPlot,m*xPlot+q,'-b','HandleVisibility','off')

BrDes = (Br(2)-Br(1))/(fM(2)-fM(1))*fq0;

disp(['Br during design = ' num2str(BrDes,4) ' T'])

BrPM = dataSet.Br;

%% 4) Compute the PMs volume

AreaCDes = geo.AreaC;
AreaCDes(isnan(AreaCDes))=0;
AreaEDes = geo.AreaE;
AreaEDes(isnan(AreaEDes))=0;

AreaTotDes = AreaCDes+AreaEDes;       % total area, virtual PM

AreaTot = AreaTotDes*BrDes/BrPM;        % total area, real PM

% To decide where to put the real PM, start to fill the 
AreaC = AreaTot;
AreaC(AreaC>AreaCDes) = AreaCDes(AreaC>AreaCDes);
AreaE = AreaTot-AreaC;
AreaE(AreaE>AreaEDes) = AreaEDes(AreaE>AreaEDes);

PMdim = [AreaC;AreaE]./[geo.hc;geo.hc];

if BrDes>BrPM
    warning('off','backtrace')
    warning('The selected PM grade is too low for the target characteristic current')
    warning('on','backtrace')
end


%% 5) save the new machine

dataSet = dataSet0;
dataSet.PMdim = PMdim;

dataSet.PMdesign.Br12  = Br;
dataSet.PMdesign.fM12  = fM;
dataSet.PMdesign.iq0   = iq0;
dataSet.PMdesign.fq0   = fq0;
dataSet.PMdesign.BrPM  = BrPM;
dataSet.PMdesign.BrDes = BrDes;
dataSet.PMdesign.geo   = geo;

%% 6) remove the temporary machine

delete([pathname filename(1:end-4) '.mat']);
delete([pathname filename(1:end-4) '.fem']);

