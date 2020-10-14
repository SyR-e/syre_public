
% build_fdfq_idiq.m
% interp maps on a fine grid of nxn points

% input:
% F_map -- machine maps (FEA or Exp)
% pathname -- final path
% n -- number of points

% output:
% results file fdfq_idiq_nxxx.mat saved in pathname

% Content of output file:
% 1) nxn meshgrid Id Iq
% 2) nxn matrixes Fd Fq organised according to meshgrid
% 3) (optional) other quantities, such as T,dT,loss ..

function build_fdfq_idiq(n,F_map,Kr,Kl,Lld,Llq,pathname,motor_name)

% fine id iq grid
i_d=linspace(round(min(min(F_map.Id))),round(max(max(F_map.Id))),n);
i_q=linspace(round(min(min(F_map.Iq))),round(max(max(F_map.Iq))),n);
[Id,Iq]=meshgrid(i_d,i_q);
% interp flux linkages
Fd = interp2(F_map.Id(1,:),F_map.Iq(:,1),F_map.Fd,Id,Iq,'spline');
Fq = interp2(F_map.Id(1,:),F_map.Iq(:,1),F_map.Fq,Id,Iq,'spline');
% change number of turns
Id=Id/Kr;
Iq=Iq/Kr;
Fd=Fd*Kr;
Fq=Fq*Kr;
% change stack length
Fd=Fd*Kl;
Fq=Fq*Kl;
% add end-connections term
Fd = Fd + Lld * Id;
Fq = Fq + Llq * Iq;

save ([pathname 'fdfq_idiq_n' num2str(n) '.mat'],'Fd','Fq','Id','Iq');

% interp torque and torque ripple
if isfield(F_map,'T')
    T = interp2(F_map.Id(1,:),F_map.Iq(:,1),F_map.T,Id,Iq,'spline');
    T = T*Kl;
    save ([pathname 'fdfq_idiq_n' num2str(n) '.mat'],'T','-append');
end
if isfield(F_map,'dT')
    dT = interp2(F_map.Id(1,:),F_map.Iq(:,1),F_map.dT,Id,Iq,'spline');
    dT = dT*Kl;
    save ([pathname 'fdfq_idiq_n' num2str(n) '.mat'],'dT','-append');
end


% pm loss
if isfield(F_map,'Ppm')
    Ppm = F_map.Ppm;
    Ppm = interp2(F_map.Id(1,:),F_map.Iq(:,1),Ppm,Id,Iq,'spline');
    save ([pathname 'fdfq_idiq_n' num2str(n) '.mat'],'Ppm','-append');
end

% core loss (MAGNET)
if isfield(F_map,'Pfes_h')
    
    Pfes_h = F_map.Pfes_h;
    Pfes_c = F_map.Pfes_c;
    Pfer_h = F_map.Pfer_h;
    Pfer_c = F_map.Pfer_c;
    
    Pfes_h = interp2(F_map.Id(1,:),F_map.Iq(:,1),Pfes_h,Id,Iq,'spline');
    Pfes_c = interp2(F_map.Id(1,:),F_map.Iq(:,1),Pfes_c,Id,Iq,'spline');
    Pfer_h = interp2(F_map.Id(1,:),F_map.Iq(:,1),Pfer_h,Id,Iq,'spline');
    Pfer_c = interp2(F_map.Id(1,:),F_map.Iq(:,1),Pfer_c,Id,Iq,'spline');
    
    velDim = F_map.speed;
    
    save ([pathname 'fdfq_idiq_n' num2str(n) '.mat'],'Pfes_h','Pfes_c','Pfer_h','Pfer_c','velDim','-append');
    
end

% IM motor only -- calc w_slip map @100C
if and(exist('Rr','var'),exist('F_IM','var'))
    
    LR = F_IM.LR;
    
    id = F_map.Id;
    iq = F_map.Iq;
    
    w_slip = Rr./LR .* iq./id;
    
    Wslip = interp2(id,iq,w_slip,Id,Iq);
    
    Ir = interp2(id,iq,F_IM.IR,Id,Iq);
    
    IM.Ir = Ir;
    IM.Wslip = Wslip;
    IM.Rr = Rr;
    IM.Rr_temp = Rr_temp;
    IM.Pr = 1.5 * Rr * Ir.^2;
    
    % %     debug
    %     figure
    %     meshc(id,iq,w_slip*30/pi/p), hold on
    %     figure
    %     meshc(Id,Iq,Wslip*30/pi/p)
    
    %     save ([pathname 'fdfq_idiq_n' num2str(n) '.mat'],'Fd','Fq','Id','Iq','Rr','Wslip');
    save ([pathname 'fdfq_idiq_n' num2str(n) '.mat'],'IM','-append');
end


