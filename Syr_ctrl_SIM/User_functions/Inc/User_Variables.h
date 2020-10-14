#include "User_data_types.h"
#include "User_Constants.h"

//General control variables
float Go;
int	Go_flag;
int State;
int counter=0;
int Ctrl_type;

//Acquisition channels and feedback
Xabc isabc;
Xalphabeta isab;
Xin input,offset_in;
float vdc;
float offset_current_a, offset_current_b, offset_current_c;

//Protection variables
int pwm_stop;

//Reference and open-loop variables
Xabc vsabc_ref, vsabc_k0, duty_abc, duty_abc_km1;
Xalphabeta vsab_ref, vsab_k0, vsab_km1;
Xdq vsdq_ref;
float v0, flux_nom, volt_nom;
XPIRegPars  id_par, iq_par, sp_par;
XPIRegVars  id_var, iq_var, sp_var;
float kp_w, ki_w, kp_id, ki_id, kp_iq, ki_iq;

//Mechanical variables
float omega_ref_in, omega_ref_ramp, accel, theta_ref;
float omega_elt_meas, theta_elt_meas, omega_mec_meas, theta_mec_meas, omega_elt_meas_f, omega_mec_meas_f, omega_mec_meas_rpm;
Xsc SinCos_ref, SinCos_elt_meas, SinCos_elt_meas_old;

// Flux observer
Xdq lambda_dq, lambda_CM_dq; 
Xalphabeta lambda_CM_ab, lambda_CM_ab_km1, feedback_OBS, lambda_OBS;
float delta;
// Compute inductance
float Ld, Ld_unfilt, Lq, Lq_unfilt, ld, ld_unfilt, lq, lq_unfilt, ldq, ldq_unfilt, ldq_d, ldq_q;
float deadtime = 1e-6;

float tmp1, tmp2, tmp3;
float n_ref_in, omega_r_filt, omega_r;

#ifdef SIM
float Reset,  BlueBTN;
#endif

Xdq isdq_ref, isdq, isdq_ext;

float T_ext, T_elt, Tref, Tlim, theta_r;

