#include "User_data_types.h"
#include "User_Constants.h"

//General control variables
float Go;
int	Go_flag;
int State;
int counter=0;
int Ctrl_type;
int Quad_Maps;
int Ctrl_strategy;

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
Xdq vsdq_ref, vsdq_k0;
float v0, flux_nom, volt_nom;
XPIRegPars  id_par, iq_par, sp_par;
XPIRegVars  id_var, iq_var, sp_var;
float kp_w, ki_w, kp_id, ki_id, kp_iq, ki_iq;

//Mechanical variables
float omega_ref_in, omega_ref_ramp, accel, theta_ref, dTheta;
float omega_elt_meas, theta_elt_meas, theta_elt_meas_enc, omega_mec_meas, theta_mec_meas, omega_elt_meas_f, omega_mec_meas_f, omega_mec_meas_rpm;

Xsc SinCos_ref; 
Xsc SinCos_elt_meas;
Xsc SinCos_elt_meas_old;
Xsc SinCos_elt;
Xsc SinCos_elt_dTheta;
Xsc SinCos_thetaS_dTheta;
Xsc SinCos_thetaS;
float omega_elt;

// Flux observer
Xdq lambda_dq;
Xdq lambda_CM_dq; 
Xalphabeta lambda_CM_ab;
Xalphabeta lambda_CM_ab_km1;
Xalphabeta feedback_OBS;
Xalphabeta_amp lambda_OBS;
float delta;
float lambda_M,th0;

// Feed-forward
Xdq vffw_dq;
Xdq vffw_dsqs;



// Compute inductance
float Ld, Ld_unfilt, Lq, Lq_unfilt, ld, ld_unfilt, lq, lq_unfilt, ldq, ldq_unfilt, ldq_d, ldq_q, ldm, Ldm;

float tmp1, tmp2, tmp3;
float n_ref_in, omega_r_filt, omega_r;

// HF injection
int counterHF, ns = 20; // ns = fs/f_inj_puls
float V_inj = 100, v_inj, f_inj_puls = 500, theta_inj = 0, f_inj_sq = fs/2, vhf = 125;
float acc_id, buffer_id[20], ishf, ishf_mag;
float ishf, lambdahf, ishf_mag, ishf_mag_old, lambdahf_mag, pos_err_LS;
float isq_old, lambdaq_old;

//Direct Flux Vector Control
XPIRegPars lambda_par,iqs_par,delta_par;
XPIRegVars lambda_var,iqs_var,delta_var;
Xsc SinCos_delta;
Xdq vdsqs_ref;
Xdq idsqs;
float lambda_MTPA,Flux_Lim,lambda_ref,iqs_ref;
Xdq isdq_aux,lambda_aux;
float delta_mtpv,delta_max;
float a,b;
float tmp_1;
float delta_ref;
float iqs_Lim;
float T_ref;
float i_MTPV;

//Flux Control
float flux_ref,flux_ext_ref;
float T_flux_ref;
float T_Ref;


// High speed control
float pos_err_HS;
int HS_ctrl; // AF or APP

// AF
float k_AF;
Xalphabeta_amp lambda_AF;
Xsc SinCos_AF;

// APP
Xdq lambda_aux, lambda_diff;
float k_APP;

// PLL
XPIRegPars PLL_par;
XPIRegVars PLL_var;
float OMEGA_B_PLL;
float omega_PLL, omega_PLL_filtered, theta_PLL, k_e, f_omega = 0;
float pos_err, position_error_real;
int SS_on, cnt = 0;
int dem; 	 // CRNTDEM: current demodulation, FLUXDEM: flux demodulation
int inj_waveform; // SINUS: pulsating, SQUARE: squarewave


// Mechanical PLL
Xsc SinCos_mec_PLL,SinCos_elt_meas,SinCos_mec_meas;
XPIRegPars PLL_mec_par;
XPIRegVars PLL_mec_var;
float omega_mec_PLL,omega_mec,theta_mec_PLL;


#ifdef SIM
float Reset,  BlueBTN;
#endif

Xdq isdq_ref, isdq, isdq_ext;

float T_ext, T_elt, Tref, Tlim, theta_r;