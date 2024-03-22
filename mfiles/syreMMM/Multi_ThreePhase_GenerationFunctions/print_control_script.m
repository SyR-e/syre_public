
function print_control_script(ctrlFolder_path,n_set,Quad_Maps)

%Initialization 

line.inputs = 36;
line.outputs = 37;
line.ia = 94;

Motor_ctrl_0_path = [ctrlFolder_path '\Motor_ctrl_0.c'];

fid = fopen(Motor_ctrl_0_path,'r');
i = 1;
tline = fgetl(fid);
readData{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    readData{i} = tline;
end
fclose(fid);

Motor_ctrl_0 = string(readData)';

%% -------------------------Aggiunta di ingressi-------------------------%%

Motor_ctrl = Motor_ctrl_0(1:107);
Motor_ctrl(line.ia) = sprintf(blanks(4)+"isabc1.a   			= U(%d);",0);
Motor_ctrl(line.ia+1) = sprintf(blanks(4)+"isabc1.b   			= U(%d);",1);
Motor_ctrl(line.ia+2) = sprintf(blanks(4)+"isabc1.c   			= U(%d);",2);
x=11;
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"Ctrl_type			= U(%d);",x)];
x=x+1;
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"accel				= U(%d);",x)];
x=x+1;
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"Quad_Maps			= U(%d);",x)];
x=x+1;
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"lambda_M			= U(%d);",x)];
x=x+1;
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"th0			        = U(%d);",x)];
x=x+1;       
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"L_sigma   			= U(%d);",x)];
x = x+1;

for i=2:n_set
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"isabc%d.a   			= U(%d);",i,x)];
    x = x+1;
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"isabc%d.b   			= U(%d);",i,x)];
    x = x+1;
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"isabc%d.c   			= U(%d);",i,x)];
    x = x+1;
end


Motor_ctrl(line.inputs) = sprintf('#define NINPUTS  %d', 17+(n_set-1)*3); 
Motor_ctrl(line.outputs) = sprintf('#define NOUTPUTS   %d', 34+4*(n_set-1)); 


%% -------------------------ERROR STATE----------------------------------%%


Motor_ctrl = [Motor_ctrl; strings(2,1)];
Motor_ctrl = [Motor_ctrl; Motor_ctrl_0(121:148)];
Motor_ctrl = [Motor_ctrl; strings(1,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"isdq_refcm.d        = 0.0f;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"isdq_refcm.q        = 0.0f;")];

switch(Quad_Maps)
    case 0 %SyR Motor
        Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_dq.d         = 0.0f;")];
        Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_dq.q         = 0.0f;")];
        Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_CM_dq.d      = 0.0f;")];
        Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_CM_dq.q      = 0.0f;")];
    case 1 %PM-SyR
        Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_dq.d         = 0.0f;")];
        Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_dq.q         = -lambda_M;")];
        Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_CM_dq.d      = 0.0f;")];
        Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_CM_dq.q      = -lambda_M;")];
    case 2 % PM
        Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_dq.d         = lambda_M;")];
        Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_dq.q         = 0.0f;")];
        Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_CM_dq.d      = lambda_M;")];
        Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_CM_dq.q      = 0.0f;")];
end


Motor_ctrl = [Motor_ctrl; strings(2,1)];

index_end = length(Motor_ctrl);


for i=1:n_set
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"isdq_refdm%d.d        = 0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"isdq_refdm%d.q        = 0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"isdq%d.d              = 0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"isdq%d.q              = 0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"vsdq%d_ref.d          = 0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"vsdq%d_ref.q          = 0.0f;",i)];
    switch(Quad_Maps)
        case 0 %SyR
            Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_OBS%d.alpha    = 0.0f;",i)];
            Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_OBS%d.beta     = 0.0f;",i)];
        case 1 %PM-SyR
            Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_OBS%d.alpha    = lambda_M*cos(PI*0.5-PP*th0);",i)];
            Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_OBS%d.beta     = lambda_M*sin(PI*0.5-PP*th0);",i)];
        case 2 % PM    
            Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_OBS%d.alpha    = lambda_M*cos(PP*th0);",i)];
            Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"lambda_OBS%d.beta     = lambda_M*sin(PP*th0);",i)];
    end       
    Motor_ctrl = [Motor_ctrl; strings(1,1)];
end

for i=1:n_set
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"isabc%d.a        = 0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"isabc%d.b        = 0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"isabc%d.c        = 0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"duty_abc%d.a     = 0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"duty_abc%d.b     = 0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"duty_abc%d.c     = 0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; strings(1,1)];
end

Motor_ctrl = [Motor_ctrl; strings(1,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"if(Go) State= WAKE_UP;")];
Motor_ctrl = [Motor_ctrl; strings(1,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"break;")];
Motor_ctrl = [Motor_ctrl; strings(2,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(8)+"case WAKE_UP:")];

%% ----------------------------WAKE UP ----------------------------------%%

for i=1:n_set
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"duty_abc%d.a     = 0.5f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"duty_abc%d.b     = 0.5f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"duty_abc%d.c     = 0.5f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"pwm_stop%d       = 0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; strings(1,1)];
end

Motor_ctrl = [Motor_ctrl; Motor_ctrl_0(238:241)];
Motor_ctrl = [Motor_ctrl;sprintf(blanks(12)+"counter++;")];
Motor_ctrl = [Motor_ctrl; strings(1,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"break;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(8)+"case READY:")];

%% ----------------------------READY ------------------------------------%%

for i=1:n_set
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"duty_abc%d.a     = 0.5f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"duty_abc%d.b     = 0.5f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"duty_abc%d.c     = 0.5f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"pwm_stop%d       = 0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; strings(1,1)];
end
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"counter        = 0.0f;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"PLL_var.intg   = omega_elt_meas_f;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"theta_PLL      = theta_elt_meas;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"if(Go) State = START;")];

Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"break;")];


Motor_ctrl = [Motor_ctrl; sprintf(blanks(8)+"case START:")];


%% ---------------------------------START--------------------------------%%

%--------------------------------Speed Computation-----------------------%%
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"//-------------------Speed Compute----------------------------------//")];
Motor_ctrl = [Motor_ctrl; strings(1,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"theta_elt_meas      = PP * theta_mec_meas;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"while(theta_elt_meas > PI)")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"theta_elt_meas -= TWOPI;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"while(theta_elt_meas < -PI)")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"theta_elt_meas += TWOPI;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"SinCos_elt_meas.sin = sin(theta_elt_meas);")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"SinCos_elt_meas.cos = cos(theta_elt_meas);")];
Motor_ctrl = [Motor_ctrl; strings(1,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"speed_compute_sc(SinCos_elt_meas, &SinCos_elt_meas_old, &omega_elt_meas);")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"_Filter(omega_elt_meas, omega_elt_meas_f, Ts*TWOPI*50);")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"omega_mec_meas      = omega_elt_meas/PP;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"omega_mec_meas_f    = omega_elt_meas_f/PP;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"omega_mec_meas_rpm  = omega_mec_meas_f*30/PI;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"SinCos_elt_meas_old.sin = SinCos_elt_meas.sin;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"	SinCos_elt_meas_old.cos = SinCos_elt_meas.cos;")];


Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"//---------------------------------PLL-------------------------------//")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"if(SS_on) {	")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"PLL_par.kp   = 2*OMEGA_B_PLL;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"PLL_par.ki   = pow(OMEGA_B_PLL,2)*Ts;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"PLL_par.lim  = RPM2RAD * nmax_mot * PP;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"PLL_var.ref  = pos_err;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"PLL_var.fbk  = 0;")];

Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"PIReg(&PLL_par, &PLL_var);")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"omega_PLL  = PLL_var.out;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"theta_PLL += Ts*PLL_var.out; ;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"_Filter(omega_PLL, omega_elt_meas_f, Ts*TWOPI*25);")];

Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"if(theta_PLL >= TWOPI)")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"theta_PLL -= TWOPI;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"if(theta_PLL <0)")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"theta_PLL += TWOPI;")];

Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"SinCos_elt.sin = sin(theta_PLL);")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"SinCos_elt.cos = cos(theta_PLL);")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"position_error_real = asin(sin(theta_elt_meas - theta_PLL));")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"}")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"else { //Encorder	")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"SinCos_elt.sin = sin(theta_elt_meas);")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"SinCos_elt.cos = cos(theta_elt_meas);")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"omega_elt = omega_elt_meas_f;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"position_error_real = 0;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"}")];



%-------------------------------Control trajectory-----------------------%%

Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"//-------------------------------Control Type ------------------------------//")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"switch (Ctrl_type){")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"case 0: //CurrentControl")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"isdq_refcm.d = isdq_ext.d;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"isdq_refcm.q = isdq_ext.q;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"break;")];
Motor_ctrl = [Motor_ctrl; strings(1,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"case 2: //TorqueControl")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"ReadLut(&ID_REF[0], fabs(T_ext), TMAX, TMIN, DT, INV_DT, &isdq_refcm.d);")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"ReadLut(&IQ_REF[0], fabs(T_ext), TMAX, TMIN, DT, INV_DT, &isdq_refcm.q);")];
for i=1:(n_set-1)
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"isdq_refdm%d.d=0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"isdq_refdm%d.q=0.0f;",i)];
end

Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"switch(Quad_Maps){")]; 
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+4)+"case 0:")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+8)+"isdq_refcm.d = sgn(T_ext)*isdq_refcm.d;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+4)+"break;")];       
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+4)+"case 1:")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+8)+"isdq_refcm.d = sgn(T_ext)*isdq_refcm.d;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+4)+"break;")];                
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+4)+"case 2:")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+8)+"isdq_refcm.q = sgn(T_ext)*isdq_refcm.q;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+4)+"break;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"}")];  
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"break;")];  
Motor_ctrl = [Motor_ctrl; strings(1,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"case 3: //Speed Control")];  
Motor_ctrl = [Motor_ctrl; Motor_ctrl_0(358:368)]; 
Motor_ctrl = [Motor_ctrl; strings(1,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"ReadLut(&ID_REF[0], fabs(T_ext), TMAX, TMIN, DT, INV_DT, &isdq_refcm.d);")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"ReadLut(&IQ_REF[0], fabs(T_ext), TMAX, TMIN, DT, INV_DT, &isdq_refcm.q);")];
for i=1:(n_set-1)
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"isdq_refdm%d.d=0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"isdq_refdm%d.q=0.0f;",i)];
end

Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"switch (Quad_Maps){")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+4)+"case 0:")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+8)+"isdq_refcm.d = sgn(T_ext)*isdq_refcm.d;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+4)+"break;")];       
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+4)+"case 1:")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+8)+"isdq_refcm.d = sgn(T_ext)*isdq_refcm.d;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+4)+"break;")];                
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+4)+"case 2:")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+8)+"isdq_refcm.q = sgn(T_ext)*isdq_refcm.q;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20+4)+"break;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"}")];  
Motor_ctrl = [Motor_ctrl; sprintf(blanks(20)+"break;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"}")];
Motor_ctrl = [Motor_ctrl; strings(2,1)];

% ---------------------------------rotational transformation-------------%%
%Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"_clarke(isabc1,isab1);")];
for i=1:n_set
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"_clarke%d(isabc%d,isab%d);",i,i,i)];
end

for i=1:n_set
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"_rot(isab%d,SinCos_elt,isdq%d);",i,i)];
end


Motor_ctrl = [Motor_ctrl; strings(1,1)];


%% ---------------------------------Current Decoupling--------------------%



tmp_s = sprintf('_Decoupling(');

for i=1:n_set
tmp_s = [tmp_s sprintf('isdq%d,',i)];
end

tmp_s = [tmp_s sprintf('isdq_cm')];

for i=1:(n_set-1)
tmp_s = [tmp_s sprintf(',isdq_dm%d',i)];
end

Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"%s);",tmp_s)];
Motor_ctrl = [Motor_ctrl; strings(1,1)];


%% --------------------------------Dead-Time Compensation-----------------%


for i=1:n_set
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"vsab%d_km1=vsab%d_k0;",i,i)];
end
Motor_ctrl = [Motor_ctrl; strings(1,1)];
for i=1:n_set
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"DTComp(duty_abc%d,duty_abc%d_km1,isabc%d,vdc,deadtime,&vsabc%d_k0);",i,i,i,i)];
end
Motor_ctrl = [Motor_ctrl; strings(1,1)];
for i=1:n_set
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"_clarke%d(vsabc%d_k0,vsab%d_k0);",i,i,i)];
end
Motor_ctrl = [Motor_ctrl; strings(1,1)];

Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"FluxObserver();")];

Motor_ctrl = [Motor_ctrl; strings(1,1)];

%% ---------------------------Current Vector Control---------------------%%

Motor_ctrl = [Motor_ctrl; strings(1,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"//--------------------Current Vector Control--------------------//")];
Motor_ctrl = [Motor_ctrl; strings(1,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"kp_cmd=OMEGA_BI*Ld_inic;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"ki_cmd=OMEGA_BI*RS;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"kp_dmd=OMEGA_BI*L_sigma;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"ki_dmd=OMEGA_BI*RS;")];
Motor_ctrl = [Motor_ctrl; strings(1,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"kp_cmq=OMEGA_BI*Lq_inic;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"ki_cmq=OMEGA_BI*RS;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"kp_dmq=OMEGA_BI*L_sigma;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(8+4)+"ki_dmq=OMEGA_BI*RS;")];

Motor_ctrl = [Motor_ctrl; strings(1,1)];

Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"id_par.kp     = kp_cmd;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"id_par.ki     = ki_cmd*Ts;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"iq_par.kp     = kp_cmq;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"iq_par.ki     = ki_cmq*Ts;")];
Motor_ctrl = [Motor_ctrl; strings(1,1)];

% define PI regulators gains
for i=2:n_set
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"id_par%d.kp     = kp_dmd;",i-1)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"id_par%d.ki     = ki_dmd*Ts;",i-1)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"iq_par%d.kp     = kp_dmq;",i-1)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"iq_par%d.ki     = ki_dmq*Ts;",i-1)];
    Motor_ctrl = [Motor_ctrl; strings(1,1)];
end


% define current loops
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"Current_loop(vdc,Imax_mot,isdq_refcm,isdq_cm,&id_par,&id_var,&iq_par,&iq_var,&vsdq_cm_ref);")];
for i=2:n_set
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"Current_loop(vdc,Imax_mot,isdq_refdm%d,isdq_dm%d,&id_par%d,&id_var%d,&iq_par%d,&iq_var%d,&vsdq_dm%d_ref);",i-1,i-1,i-1,i-1,i-1,i-1,i-1)];
end


%% -----------------------------Feedforward-------------------------------%


Motor_ctrl = [Motor_ctrl; strings(2,1)];

Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"vsdq_cm_ref.d += RS*isdq_cm.d-omega_elt*lambda_dq.q;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"vsdq_cm_ref.q += RS*isdq_cm.q+omega_elt*lambda_dq.d;")];
Motor_ctrl = [Motor_ctrl; strings(1,1)];

for i=2:n_set    
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"vsdq_dm%d_ref.d += RS*isdq_dm%d.d-omega_elt*L_sigma*isdq_dm%d.q;",i-1,i-1,i-1)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"vsdq_dm%d_ref.q += RS*isdq_dm%d.q+omega_elt*L_sigma*isdq_dm%d.d;",i-1,i-1,i-1)];
    Motor_ctrl = [Motor_ctrl; strings(1,1)];
end

%% --------------------Voltage Decoupling---------------------------------%

Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"//Voltage Decoupling")];


tmp_s = sprintf('_InvDecoupling(');
tmp_s = [tmp_s sprintf('vsdq_cm_ref')];

for i=1:(n_set-1)
    tmp_s = [tmp_s sprintf(',vsdq_dm%d_ref',i)];
end

for i=1:n_set
    tmp_s = [tmp_s sprintf(',vsdq%d_ref',i)];
end

Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"%s);",tmp_s)];
Motor_ctrl = [Motor_ctrl; strings(1,1)];



%% ---------------------------------PWM GENERATION------------------------%

for i=1:n_set 
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"_invrot(vsdq%d_ref,SinCos_elt,vsab%d_ref);",i,i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"_invclarke%d(vsab%d_ref,vsabc%d_ref);",i,i,i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"PWMduty(vsabc%d_ref,vdc,&duty_abc%d);",i,i)];
    Motor_ctrl = [Motor_ctrl; strings(1,1)];
end

Motor_ctrl = [Motor_ctrl; sprintf(blanks(5)+"}")];
Motor_ctrl = [Motor_ctrl; strings(2,1)];

% duty cycle saturation
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"//Duty cycle saturation")];
for i=1:n_set
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(5)+"if (duty_abc%d.a>0.99f) duty_abc%d.a = 0.99f;",i,i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(5)+"if (duty_abc%d.b>0.99f) duty_abc%d.b = 0.99f;",i,i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(5)+"if (duty_abc%d.c>0.99f) duty_abc%d.c = 0.99f;",i,i)];
    Motor_ctrl = [Motor_ctrl; strings(1,1)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(5)+"if (duty_abc%d.a<0.01f) duty_abc%d.a = 0.01f;",i,i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(5)+"if (duty_abc%d.b<0.01f) duty_abc%d.b = 0.01f;",i,i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(5)+"if (duty_abc%d.c<0.01f) duty_abc%d.c = 0.01f;",i,i)];
    Motor_ctrl = [Motor_ctrl; strings(1,1)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(5)+"if (pwm_stop%d){",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(9)+"duty_abc%d.a = 0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(9)+"duty_abc%d.b = 0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(9)+"duty_abc%d.c = 0.0f;",i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(9)+"}")];
    Motor_ctrl = [Motor_ctrl; strings(1,1)];
end

%% ---------------------------S function Outputs--------------------------% 

Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[0] = duty_abc1.a;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[1] = duty_abc1.b;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[2] = duty_abc1.c;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[3] = pwm_stop1;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[4]  = omega_ref_ramp*60/TWOPI;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[5]  = omega_mec_meas_f*60/TWOPI;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[6]  = vsdq_cm_ref.d;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[7]  = vsdq_cm_ref.q;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[8]  = isab1.alpha;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[9]  = isab1.beta;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[10] = isdq_refcm.d;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[11] = isdq_cm.d;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[12] = isdq_refcm.q;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[13] = isdq_cm.q;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[14] = f_omega;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[15] = T_elt;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[16] = lambda_dq.d;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[17] = lambda_dq.q;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[18] = lambda_CM_dq.d;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[19] = lambda_CM_dq.q;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[20] = T_ext;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[21] = theta_PLL*180/PI;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[22] = omega_elt/PP*60/TWOPI;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[23] = ld;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[24] = lq;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[25] = ldq;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[26] = isabc1.a;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[27] = isabc1.b;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[28] = isabc1.c;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[29] = theta_elt_meas*180/PI;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[30] = pos_err*180/PI;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[31] = position_error_real*180/PI;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[32] = pos_err_LS*180/PI;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[33] = pos_err_HS*180/PI;")];
k = 34;
for i=2:n_set
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[%d] = duty_abc%d.a;",k,i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[%d] = duty_abc%d.b;",k+1,i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[%d] = duty_abc%d.c;",k+2,i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"y[%d] = pwm_stop%d;",k+3,i)];
    k=k+4;
end
			
Motor_ctrl = [Motor_ctrl; strings(2,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(1)+"}")];
Motor_ctrl = [Motor_ctrl; strings(2,1)];

Motor_ctrl = [Motor_ctrl; Motor_ctrl_0(707:719)]; 

%% --------------------------------Flux Observer--------------------------%


Motor_ctrl = [Motor_ctrl; sprintf(blanks(1)+"void FluxObserver(void){")];
Motor_ctrl = [Motor_ctrl; strings(2,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(9)+"lambda_CM_ab_km1 = lambda_CM_ab;")];
Motor_ctrl = [Motor_ctrl; strings(2,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(9)+"switch (Quad_Maps){")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(9)+"case 0:")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"ReadLut2d(&FD_LUT[0][0], fabs(isdq_cm.d), fabs(isdq_cm.q), DIDD, INV_DIDD, DIQD, INV_DIQD , ID_TAB_MAX, ID_TAB_MIN, IQ_TAB_MAX , IQ_TAB_MIN, n_size, &lambda_CM_dq.d);")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"ReadLut2d(&FQ_LUT[0][0], fabs(isdq_cm.q), fabs(isdq_cm.d), DIQQ, INV_DIQQ, DIQD, INV_DIQD , IQ_TAB_MAX, IQ_TAB_MIN, ID_TAB_MAX , ID_TAB_MIN, n_size, &lambda_CM_dq.q);")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"if(isdq_cm.d<0)")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"lambda_CM_dq.d = -lambda_CM_dq.d;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"if(isdq_cm.q<0)")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"lambda_CM_dq.q = -lambda_CM_dq.q;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(9)+"break;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(9)+"case 1:")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"ReadLut2d(&FD_LUT[0][0], fabs(isdq_cm.d), isdq_cm.q, DIDD, INV_DIDD, DIQD, INV_DIQD , ID_TAB_MAX, ID_TAB_MIN, IQ_TAB_MAX , IQ_TAB_MIN, n_size, &lambda_CM_dq.d);")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"ReadLut2d(&FQ_LUT[0][0], isdq_cm.q, fabs(isdq_cm.d), DIQQ, INV_DIQQ, DIDQ, INV_DIDQ , IQ_TAB_MAX, IQ_TAB_MIN, ID_TAB_MAX , ID_TAB_MIN, n_size, &lambda_CM_dq.q);")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"if(isdq_cm.d<0)")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"lambda_CM_dq.d = -lambda_CM_dq.d;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(9)+"break;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(9)+"case 2:")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"ReadLut2d(&FD_LUT[0][0], isdq_cm.d, fabs(isdq_cm.q), DIDD, INV_DIDD, DIQD, INV_DIQD , ID_TAB_MAX, ID_TAB_MIN, IQ_TAB_MAX , IQ_TAB_MIN, n_size, &lambda_CM_dq.d);")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"ReadLut2d(&FQ_LUT[0][0], fabs(isdq_cm.q), isdq_cm.d, DIQQ, INV_DIQQ, DIDQ, INV_DIDQ , IQ_TAB_MAX, IQ_TAB_MIN, ID_TAB_MAX , ID_TAB_MIN, n_size, &lambda_CM_dq.q);")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(12)+"if (isdq_cm.q < 0)")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(16)+"lambda_CM_dq.q = -lambda_CM_dq.q;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(9)+"break;")];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(9)+"}")];
Motor_ctrl = [Motor_ctrl; strings(1,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"_invrot(lambda_CM_dq, SinCos_elt, lambda_CM_ab);")];


for i=1:n_set
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"feedback_OBS%d.alpha = lambda_CM_ab_km1.alpha - lambda_OBS%d.alpha;",i,i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"feedback_OBS%d.beta  = lambda_CM_ab_km1.beta - lambda_OBS%d.beta;",i,i)];
end

Motor_ctrl = [Motor_ctrl; strings(1,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(5)+"//Integration")];


for i=1:n_set
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"lambda_OBS%d.alpha += Ts*(vsab%d_km1.alpha - RS*isab%d.alpha + KOBS*feedback_OBS%d.alpha);",i,i,i,i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"lambda_OBS%d.beta  += Ts*(vsab%d_km1.beta - RS*isab%d.beta + KOBS*feedback_OBS%d.beta);",i,i,i,i)];
    Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"lambda_OBS%d.amp   += sqrtf(pow(lambda_OBS%d.alpha,2)+pow(lambda_OBS%d.beta,2));",i,i,i)];
end


tmp_a = sprintf('aux_ab.alpha=(lambda_OBS1.alpha');
tmp_b = sprintf('aux_ab.beta =(lambda_OBS1.beta');
for i=2:n_set
    tmp_a = [tmp_a sprintf('+lambda_OBS%d.alpha',i)];
    tmp_b = [tmp_b sprintf('+lambda_OBS%d.beta',i)];
end

Motor_ctrl = [Motor_ctrl; strings(2,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"%s)/%d;",tmp_a,n_set)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"%s)/%d;",tmp_b,n_set)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"_rot(aux_ab,SinCos_elt,lambda_dq);")];

Motor_ctrl = [Motor_ctrl; strings(2,1)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"T_elt = %d.0f*1.5f*PP*(lambda_dq.d*isdq_cm.q-lambda_dq.q*isdq_cm.d);",n_set)];
Motor_ctrl = [Motor_ctrl; sprintf(blanks(4)+"delta = atan(lambda_dq.q/lambda_dq.d);")];

Motor_ctrl = [Motor_ctrl; sprintf(blanks(1)+"}")];



%% -------------------------Stampa del nuovo file------------------------%%

Motor_ctrl_path = [ctrlFolder_path '\Motor_ctrl.c'];

fid = fopen(Motor_ctrl_path, 'w');
for i = 1:numel(Motor_ctrl)
    fprintf(fid,'%s\n', Motor_ctrl{i});
end
fclose(fid);

delete(Motor_ctrl_0_path);

end









