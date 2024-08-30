function [geo,mat,opt,T,P,M] = preprocess_therm(geo,per,mat,filename)
filename_FEMM = sprintf('%s.ans',filename);
[~,~,~,P,T,M]=readFemm(filename_FEMM);
% filename_mat = sprintf('%s.mat',filename);
% load(filename_mat)
% addpath("readfemm\")
%% FEMM mesh update (air gap meshing, material code definition) 
tol=1e-5; % in meters
mat_code_shaft = max(M)+1; % Repetto dixit
mat_code_gap = max(M)+2; % Repetto dixit
mat_code_liner = max(M)+3; % Repetto dixit
mat_code_air = 7;
% 02/07/24
r_in_sb=(geo.r+geo.g/3)*1e-3; % inner radius of the sliding band
r_out_sb=(geo.r+geo.g*2/3)*1e-3; % outer radius of the sliding band
r_sh=geo.Ar*1e-3; % shaft radius
r_in_stat=(geo.r + geo.g)*1e-3; % inner radius of stator
r_out_rot=(geo.r)*1e-3; % outer radius of rotor
%
% n_nodes=length(P(:,1));
r=sqrt(P(:,1).^2 + P(:,2).^2);
node_aux1=find(abs(r-r_in_sb)<tol);
theta_in=atan2(P(node_aux1,2),P(node_aux1,1));
node_aux2=find(abs(r-r_out_sb)<tol);
theta_out=atan2(P(node_aux2,2),P(node_aux2,1));
n_gap1 = length(node_aux1);
n_gap2 = length(node_aux2);
[theta_in_s, idx_in]=sort(theta_in);
[theta_out_s, idx_out]=sort(theta_out);
T_pos=length(T(:,1))+1;
for j=1:n_gap1-1
    T(T_pos,:)=[node_aux1(idx_in(j)), node_aux2(idx_out(j)), node_aux1(idx_in(j+1))];
    M(T_pos)=mat_code_gap;
    T_pos=T_pos+1;
    T(T_pos,:)=[node_aux1(idx_in(j+1)), node_aux2(idx_out(j)), node_aux2(idx_out(j+1))];
    M(T_pos)=mat_code_gap;
    T_pos=T_pos+1;
end
% 02/07/24 Compute element center
for k=1:length(T(:,1))
    T_bar(k,1)=sum(P(T(k,:),1))/3;
    T_bar(k,2)=sum(P(T(k,:),2))/3;
end
R_bar=sqrt(T_bar(:,1).^2+T_bar(:,2).^2);

for k=1:length(T(:,1))
    if R_bar(k)<r_sh
        M(k)=mat_code_shaft;
    elseif R_bar(k)<r_in_stat && R_bar(k)>r_out_rot
        M(k)=mat_code_gap;
    elseif R_bar(k)>r_in_stat
           if M(k)==mat_code_air
               M(k)=mat_code_liner;
           end
    end
end

%%
geo.lfact=0.1; % relative length of heads

%%
opt.n_repetition_stack=3; % how many mesh elements for the stack
opt.n_repetition_head=1; % how many mesh elements for the head

% opt.lstack=geo.l*1e-3;
% opt.lhead=geo.l*1e-3*0.1;

opt.r_min_head=(geo.r+geo.g+geo.ttd)*1e-3; % min radius of the head
opt.r_max_head=(geo.R-geo.ly+geo.pont0)*1e-3;  % max radius of the head

% % Per validazione UNIPD
% opt.r_min_head=9.375e-3; 
% opt.r_max_head=14.6e-3;

% 09/07/2024 Abbiamo deciso di lasciare le caratteristiche termiche dei
% materiali pluggate dentro 

opt.Tamb=60;   %[°C]
opt.Time=30;   %[s]
opt.Nt=10;     %number of time step
opt.plot=false; 
opt.plot_final=false; % plot della soluzione a t=30s

mat.LayerAir.CondTerm=29e-3; % https://www.engineeringtoolbox.com/air-properties-viscosity-conductivity-heat-capacity-d_1509.html
mat.LayerAir.HeatCap=700;
mat.LayerAir.kgm3=1.2;

mat.SlotAir.CondTerm=1.9; 
mat.SlotAir.HeatCap=733;
mat.SlotAir.kgm3=1400;

mat.Stator.CondTerm=30;
mat.Stator.HeatCap=460;

% SF: non è formalmente corretto, ma per le proprietà del conduttore ho
% mediato quelle del Cu con una resina in base al fattore di riempimento
% cava...da verificare
mat.SlotCond.CondTerm = 400*geo.win.kcu+1.9*(1-geo.win.kcu);
mat.SlotCond.HeatCap  = 385*geo.win.kcu+733*(1-geo.win.kcu);
mat.SlotCond.kgm3     = 8900*geo.win.kcu+1400*(1-geo.win.kcu);

mat.LayerMag.CondTerm=9;
mat.LayerMag.HeatCap=460;

mat.Shaft.CondTerm=25;
mat.Shaft.HeatCap=460;
mat.Shaft.kgm3=7770;

opt.q=per.Loss/(2*geo.p*2);  % potenza dissipata (W) /2 per simmetria /6 
opt.h=2404;%5.6761e+03; %1e3;  % coefficiente di scambio termico convettivo. 
%%
mat.CondTerm(14)=mat.Stator.CondTerm;  % extra iron head
mat.kgm3(14)    =mat.Stator.kgm3;
mat.HeatCap(14) =mat.Stator.HeatCap;

mat.CondTerm(1)=mat.Stator.CondTerm;  % iron
mat.kgm3(1)    =mat.Stator.kgm3;
mat.HeatCap(1) =mat.Stator.HeatCap;
 
mat.CondTerm(13)=mat.SlotCond.CondTerm; % winding head 
mat.kgm3(13)    =mat.SlotCond.kgm3;
mat.HeatCap(13) =mat.SlotCond.HeatCap;

mat.CondTerm(3) =mat.SlotCond.CondTerm; % windings
mat.kgm3(3)     =mat.SlotCond.kgm3;
mat.HeatCap(3)  =mat.SlotCond.HeatCap;

mat.CondTerm(9) =mat.LayerMag.CondTerm; % magnet?
mat.kgm3(9)     =mat.LayerMag.kgm3;
mat.HeatCap(9)  =mat.LayerMag.HeatCap;

mat.CondTerm(8) =mat.LayerMag.CondTerm; % magnet?
mat.kgm3(8)     =mat.LayerMag.kgm3;
mat.HeatCap(8)  =mat.LayerMag.HeatCap;

mat.CondTerm(7) =mat.LayerAir.CondTerm; % air in magnet slot
mat.kgm3(7)     =mat.LayerAir.kgm3;
mat.HeatCap(7)  =mat.LayerAir.HeatCap;

mat.CondTerm(10) =mat.Shaft.CondTerm; % shaft 
mat.kgm3(10)     =mat.Shaft.kgm3;
mat.HeatCap(10)  =mat.Shaft.HeatCap;

mat.CondTerm(11)=mat.LayerAir.CondTerm;  % Air gap (Da modificare in funzione della velocità)  
mat.kgm3(11)    =mat.LayerAir.kgm3;
mat.HeatCap(11) =mat.LayerAir.HeatCap;

mat.CondTerm(12) =mat.SlotAir.CondTerm; % slot liner
mat.kgm3(12)     =mat.SlotAir.kgm3;
mat.HeatCap(12)  =mat.SlotAir.HeatCap;


mat.q = [3 13]; % where to place power

mat.Tout=[3 13]; % mean temperature to extract, in which domain. 
end
