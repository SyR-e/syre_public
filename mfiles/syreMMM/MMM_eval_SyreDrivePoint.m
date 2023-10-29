% Copyright 2022
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [iabc2] = MMM_eval_SyreDrivePoint(motorModel,syreDriveSingt)

fPWM = motorModel.SyreDrive.Converter.fPWM;
p    = motorModel.geo.p;

% cd([motorModel.data.pathname motorModel.data.motorName '_ctrl_INST'])

open(motorModel.SyreDrive.SIM_path)
init_sim

n = syreDriveSingt.n;  

id_MTPA = syreDriveSingt.Id;
iq_MTPA = syreDriveSingt.Iq;
T0      = syreDriveSingt.T;
n0      = 2*n;   
tSim    = syreDriveSingt.tSim;

flagwp = 1;

options = simset('SrcWorkspace','current', 'DstWorkspace', 'current');
out = sim(motorModel.SyreDrive.SIM_path,tSim,options);   % Starts simulation

% without movemean and percentage gap computed as
% delta% = (|a-b|/(a+b)/2)*100

Time_dq = out.Out_M.Idq_m.Time;
Id = out.Out_M.Idq_m.Data(:,1);
Iq = out.Out_M.Idq_m.Data(:,2);
T_m = out.Out_M.T_m.Data;
% Id_ref = out.Outputs.id_ref.Data(:);
% Iq_ref = out.Outputs.iq_ref.Data(:);

t = 60/(out.Out_M.n_m.Data(length(out.Out_M.n_m.Data(:)))*motorModel.data.p);

%% Method 1
% ind_st = find(Id_ref == Id_ref(length(Id_ref)),1);
% Id_tmp = Id(ind_st:length(Id));
% Iq_tmp = Iq(ind_st:length(Iq));
% 
% t_ind = floor(length(Time_dq)/(Time_dq(length(Time_dq))/t));
% md = numel(Id_tmp);
% mq = numel(Iq_tmp);
% Id_m = mean(reshape( [Id_tmp(:);nan(mod(-md,t_ind),1)],t_ind,[]),'omitnan');
% Iq_m = mean(reshape( [Iq_tmp(:);nan(mod(-mq,t_ind),1)],t_ind,[]),'omitnan');
% Id_ref1 = ones(1,length(Id_m))*Id_ref(length(Id_ref));
% Iq_ref1 = ones(1,length(Id_m))*Iq_ref(length(Iq_ref));
% 
% Id_diff = abs(Id_ref1-Id_m);
% Id_avg = 0.5*(Id_ref1+Id_m);
% Id_per1 = (Id_diff./Id_avg)*100;
% Iq_diff = abs(Iq_ref1-Iq_m);
% Iq_avg = 0.5*(Iq_ref1+Iq_m);
% Iq_per1 = (Iq_diff./Iq_avg)*100;
% 
% ind_exp_per = find(Id_per1 == min(Id_per1));
% ind_exp_in = ind_st+t_ind*ind_exp_per;
% ind_exp_fin = ind_exp_in+t_ind;
% 
% Id_exp1 = Id(ind_exp_in:ind_exp_fin);
% Iq_exp1 = Iq(ind_exp_in:ind_exp_fin);

%% Method 2

[~,ind_st2] = min(abs((Time_dq - (Time_dq(end)-t))));
Id_exp2 = Id(ind_st2:end);
Iq_exp2 = Iq(ind_st2:end);
T_mexp2    = T_m(ind_st2:end);

%% Conversion in abc
NUMSTEP = 2*fPWM/(n*p/60)*100;
th0 = motorModel.geo.th0;
th_vect = linspace(th0*pi/180,(th0+360)*pi/180,NUMSTEP);

% Id1 = interp1(linspace(0,1,length(Id_exp1)),Id_exp1,linspace(0,1,NUMSTEP));
% Iq1 = interp1(linspace(0,1,length(Iq_exp1)),Iq_exp1,linspace(0,1,NUMSTEP));

Id2  = interp1(linspace(0,1,length(Id_exp2)),Id_exp2,linspace(0,1,NUMSTEP));
Iq2  = interp1(linspace(0,1,length(Iq_exp2)),Iq_exp2,linspace(0,1,NUMSTEP));
T_m2 = interp1(linspace(0,1,length(T_mexp2)),T_mexp2,linspace(0,1,NUMSTEP));

% iabc1 = dq2abc(Id1,Iq1,th_vect);
iabc2 = dq2abc(Id2,Iq2,th_vect);

% figure
% figSetting
% plot(iabc1(1,:),'b','LineWidth',2)
% plot(iabc1(2,:),'g','LineWidth',2)
% plot(iabc1(3,:),'r','LineWidth',2)

time = linspace(0,60/n/p,length(iabc2(1,:)));

iabc2 = [iabc2; time];

figure
figSetting(14,10,12)
% title(['$I_d$ = ' num2str(round(id_MTPA,2)) 'A - $I_q$ = ' num2str(round(iq_MTPA,2)) ' A - n = ' num2str(round(n)) ' rpm'])
title(['$T$ = ' num2str(round(mean(T_m2),2)) ' Nm - n = ' num2str(round(n)) ' rpm'])
plot(time*1000,iabc2(1,:),'b','LineWidth',2)
plot(time*1000,iabc2(2,:),'g','LineWidth',2)
plot(time*1000,iabc2(3,:),'r','LineWidth',2)
xlim([0 t(end)*1000])
xlabel('$t$ [ms]')
ylabel('$I$ [A]')
