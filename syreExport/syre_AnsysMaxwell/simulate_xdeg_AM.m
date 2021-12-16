% Copyright 2021
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


function [SOL] = simulate_xdeg_AM(geo,per,mat,eval_type,filepath,filename,ipypath)


th0 = geo.th0(1); %initial angle between CS dq and abc
p   = geo.p;
corelossflag = per.corelossflag;
% ns  = geo.ns;
N_parallel = 1;
N_cond = geo.win.Nbob; %Ncond per slot geo.win.Ns/geo.p/(geo.q)/size(geo.win.avv,1);
q = geo.q;
gamma   = per.gamma;
gamma_  = gamma*pi/180; %current angle dq [rad]
win_avv = geo.win.avv; %layer win
nlay = geo.nlay; %layaer of rotor barriers
radial_ribs_split = nnz(geo.radial_ribs_split*geo.pontR'); %half added barrier area using split
% n_PM = geo.n_PM;                                           %number of PM's
n_PM = 0;                                                  %number of PM's
tempPP   = per.tempPP;                                       %temperature of PMs
CurrLoPP = per.overload;                                   %Current [pu]
iAmp = CurrLoPP*per.i0;


Ia   = per.custom_ia;
Ib   = per.custom_ib;
Ic   = per.custom_ic;
time = per.custom_time;
ansysCount = per.custom_ansyscount;

if per.EvalSpeed==0
    prompt    = {'The Evaluation Speed is 0! Please insert a different value in [rpm]:'};
    data      = inputdlg(prompt,'EvalSpeed',[1 35],{'2000'});
    EvalSpeed = str2double(data{1});
else
    EvalSpeed = per.EvalSpeed;
end

%Solution Type
switch eval_type
    case 'singt'
        nsim = per.nsim_singt; %Simulation Points
        xdeg = per.delta_sim_singt; %Angular Excursion
    case 'singtIron'
        nsim = per.nsim_singt; %Simulation Points
        xdeg = 420; %Angular Excursion
        %xdeg = 360; %Angular Excursion
    otherwise
        error('Ansys Maxwell Simulations not available for the selected evaluation type!');
end

resFolder = [filepath,filename(1:end-4),'_results\FEA results\Ansys_', eval_type , '_', num2str(round(iAmp,2)), 'A'];
mkdir(resFolder);

% evaluation of the phase current values for all positions to be simulated
% id = iAmp * cos(gamma_);
% iq = iAmp * sin(gamma_);
%Imod = abs(id + 1i* iq);

phi_init = th0*pi/180+gamma_; %initial current phase angle

iAmp = round(iAmp,2);

pwd1 = [pwd '\syreExport\syre_AnsysMaxwell'];

%% Custom Current
Tw = 1/(EvalSpeed*p/60);

if per.custom_act
    if length(Ia(:,1))==1
        Ia = Ia';
        Ib = Ib';
        Ic = Ic';
        time = time';
    end
    time = time-time(1);
    if time(end)<1.5*Tw
        time = [time; time(end)+time];
        Ia = [Ia; Ia];
        Ib = [Ib; Ib];
        Ic = [Ic; Ic];
    end
 
    Ia_time = [time Ia];
    Ib_time = [time Ib];
    Ic_time = [time Ic];
    
    % Current txt
    currtxtA = fopen(strcat(pwd1,'\temp\Ia_time.tab'), 'wt' );
    currtxtB = fopen(strcat(pwd1,'\temp\Ib_time.tab'), 'wt' );
    currtxtC = fopen(strcat(pwd1,'\temp\Ic_time.tab'), 'wt' );
    
    fprintf(currtxtA,'%f   %f\n',Ia_time');
    fprintf(currtxtB,'%f   %f\n',Ib_time');
    fprintf(currtxtC,'%f   %f\n',Ic_time');
    
%     ansysCount = ansysCount + 17;
end

%% Export simulation data on python
% pwd1 = [pwd '\syreExport\syre_AnsysMaxwell'];
filepath = strrep(filepath,'\','/');
pyfile = fopen( strcat(pwd1,'\temp\simdata.txt'), 'wt' );
export={
    'import pickle'
    'import sys'
    'simdata={'
    '"filename":"%s", #motfilename'
    '"filepath":"%s", #motfilepath'
    '"eval_type":"%s", #sim type'
    '"iAmp":%f, #Currents Amplitude'
    '%s #matrix with windings index'
    '"N_parallel" :%d,   #'
    '"N_cond" :%d,   #N cond per slot'
    '"EvalSpeed" :%f,   #layer of rotor barrier'
    '"phi_init":%f, #initial current phase'
    '"p":%d,#pole pair'
    '"q":%d,#slot/(phase*p)'
    '"radial_ribs_split":%d, #area added by split barrier'
    '"n_PM":%d, #number of permanent magnet (number of magnetic segments)'
    '"nlay" :%d,   #layer of rotor barrier'
    '"tempPP" :%f, #PMs temp'
    '"gamma":%f,#current angle dq [rad]'
    '"nsim":%d,#simulation points'
    '"xdeg":%f,#angular excursion [deg]'
    '"corelossflag":%d,#flag that activate core loss'
    '"save_fields":%d,#save fields every # steps'
    '"dataset_counter":%d,#count for dataset current'
    '}'
    'with open("%s/temp/temp.pkl", "wb") as export:'
    '    pickle.dump(simdata,export,protocol=2)'
    '###'
    'sys.exit()'
    };
str=sprintf('%s\n',export{:});

win_avvpy='"win_avv":[';
for ii=1:size(win_avv,1)
    win_avvpy=strcat(win_avvpy,'[');
    for jj=1:size(win_avv,2)
        win_avvpy=strcat(win_avvpy,num2str(win_avv(ii,jj)),',');
    end
    win_avvpy=strcat(win_avvpy,'],');
end
win_avvpy=strcat(win_avvpy,'],');


fprintf(pyfile,str,filename,strrep(filepath,'\','/'),eval_type,iAmp,win_avvpy,N_parallel,N_cond,EvalSpeed,...
    phi_init,p,q,radial_ribs_split,n_PM,nlay,tempPP,gamma_,nsim,xdeg,per.corelossflag,per.save_fields,ansysCount,strrep(pwd1,'\','/'));

fclose(pyfile);

%from txt to py
%movefile(strcat(pwd1,'\temp\simdata.txt'),strcat(pwd1,'\temp\simdata.py'),'f');
copyfile ([pwd1 '\temp\simdata.txt'],[pwd1 '\temp\simdata.py'])

command = strcat(ipypath,pwd1,'\temp\simdata.py"');
% dos(command);

[status,cmdout] = dos(command);

if status ~= 0
    fprintf('an error occurred - program aborted - \n %s',cmdout)
end

%%
if per.custom_act
    command = strcat(ipypath,pwd1,'\sim_xdeg_AM_customI.py"');
    [status,cmdout] = dos(command);
else
    command = strcat(ipypath,pwd1,'\sim_xdeg_AM.py"');
    [status,cmdout] = dos(command);
end

if status ~= 0
    fprintf('an error occurred - program aborted - \n %s \n',cmdout)
end

outdatafolder = strcat(filepath,filename(1:end-4),'_results\FEA results\Ansys_', eval_type , '_', num2str(round(iAmp,2)), 'A');
addpath(outdatafolder);
opts = detectImportOptions('FluxesPlotData.csv','NumHeaderLines',1);
opts.VariableNamesLine = 1;
opts.PreserveVariableNames = 1;

warning off
tmp   = readtable('TorquePlotData.csv',opts);   
T     = table2array(tmp(:,2));                                        %[t(ms),T(Nm)]
tmp   = readtable('FluxesPlotData.csv',opts);
F     = table2array(tmp(:,2:4));                                        %[t(ms),f1(Wb),f2(Wb),f3(Wb)]
tmp   = readtable('PositionData.csv',opts);
Theta = table2array(tmp(:,2));                                        %[t(ms),Theta(deg)]
tmp   = readtable('CurrentsPlotData.csv',opts);
I     = table2array(tmp(:,2:4));                                       %[t(ms),i1(A),i2(A),i3(A)]
tmp   = readtable('VoltagesPlotData.csv',opts);
V     = table2array(tmp(:,2:4));                                        %[t(ms),v1(V),v2(V),v3(V)]

SOL.T     = T;
SOL.F     = F;
SOL.Theta = Theta;
SOL.I     = I;
SOL.V     = V;
SOL.xdeg  = xdeg;

if corelossflag==1
    tmp = readtable('CoreLossData.csv',opts);
    CoreLoss   = table2array(tmp(:,2)); 
    CoreLoss_s = table2array(tmp(:,3));
    CoreLoss_r = table2array(tmp(:,4));
    %CoreLoss =  readtable('CoreLossData.txt'); %[t(ms),CoreLoss(W)]
    SOL.CoreLoss   = CoreLoss;
    SOL.CoreLoss_s = CoreLoss_s;
    SOL.CoreLoss_r = CoreLoss_r;
     
%     opts = detectImportOptions('CoreLossAvg.txt','NumHeaderLines',1);
%     opts.VariableNamesLine     = 1;
%     opts.PreserveVariableNames = 1;
%     CoreLossAvg       = readtable('CoreLossAvg.txt',opts);
%     %CoreLossAvg       = xlsread('CoreLossAvg.csv'); %[Tot;Stat;Rot]
%     HysteresisLossAvg = readtable('HysteresisLossAvg.txt',opts); %[Tot;Stat;Rot]
%     EddyLossAvg       = readtable('EddyCurrentsLossAvg.txt',opts); %[Tot;Stat;Rot]
%     ExcessLossAvg     = readtable('ExcessLossAvg.txt',opts); %[Tot;Stat;Rot]
    
%     SOL.CoreLossAvg       = CoreLossAvg;
%     SOL.HysteresisLossAvg = HysteresisLossAvg;
%     SOL.EddyLossAvg       = EddyLossAvg;
%     SOL.ExcessLossAvg     = ExcessLossAvg;
end

end