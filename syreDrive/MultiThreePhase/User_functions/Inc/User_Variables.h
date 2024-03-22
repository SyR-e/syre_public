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
float L_sigma;

//Currents in abc alpha-beta and dq axis
Xabc isabc, isabc1, isabc2, isabc3,isabc4;
Xalphabeta isab, isab1, isab2, isab3,isab4;
Xdq isdq_cm,isdq_dm1,isdq_dm2,isdq_dm3;
Xdq isdq_refcm,isdq_ref,isdq_refdm, isdq_refdm1, isdq_refdm2, isdq_refdm3,isdq_refdm4;
Xdq isdq1, isdq2, isdq3,isdq4;
Xdq isdq_ext;
float offset_current_a,offset_current_b,offset_current_c;

//Voltages in abc alpha-beta and dq axis
Xabc vsabc1_ref, vsabc2_ref, vsabc3_ref, vsabc4_ref;
Xabc vsabc1_k0,vsabc2_k0,vsabc3_k0,vsabc4_k0;
Xalphabeta vsab1_ref,vsab2_ref,vsab3_ref,vsab4_ref;
Xalphabeta vsab1_k0,vsab2_k0,vsab3_k0,vsab4_k0;
Xalphabeta vsab1_km1,vsab2_km1,vsab3_km1,vsab4_km1;
Xdq vsdq1_ref, vsdq2_ref, vsdq3_ref,vsdq4_ref;
Xdq vsdq_cm_ref,vsdq_dm1_ref,vsdq_dm2_ref,vsdq_dm3_ref,vsdq_dm4_ref;


// duty cycles
Xabc duty_abc1, duty_abc2, duty_abc3,duty_abc4;
Xabc duty_abc1_km1,duty_abc2_km1,duty_abc3_km1,duty_abc4_km1;
float vdc;

//Protection variables
int pwm_stop, pwm_stop1, pwm_stop2, pwm_stop3,pwm_stop4;
int pwm_enable1,pwm_enable2,pwm_enable3,pwm_enable4;

//PI regulators
XPIRegPars  id_par, iq_par, id_par1, iq_par1, id_par2, iq_par2, id_par3, iq_par3,id_par4, iq_par4,sp_par;
XPIRegVars  id_var, iq_var, id_var1, iq_var1, id_var2, iq_var2, id_var3, iq_var3,id_var4, iq_var4, sp_var;
float kp_w, ki_w, kp_cmd, ki_cmd, kp_cmq, ki_cmq, kp_dmd, ki_dmd, kp_dmq, ki_dmq;

//Mechanical variables
float omega_ref_in, omega_ref_ramp, accel, theta_ref, dTheta;
float omega_elt_meas, theta_elt_meas, theta_elt_meas_enc, omega_mec_meas, theta_mec_meas, omega_elt_meas_f, omega_mec_meas_f, omega_mec_meas_rpm;
Xsc SinCos_ref, SinCos_elt_meas, SinCos_elt_meas_old, SinCos_elt, SinCos_elt_dTheta;
float omega_elt;

// Flux observer
Xdq lambda_dq, lambda_CM_dq; 
Xalphabeta lambda_CM_ab, lambda_CM_ab_km1, feedback_OBS, feedback_OBS1, feedback_OBS2, feedback_OBS3,feedback_OBS4;
Xalphabeta_amp lambda_OBS, lambda_OBS1, lambda_OBS2, lambda_OBS3,lambda_OBS4,aux_ab;
float delta;
float lambda_M,th0;
float flux_nom,v0;

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

#ifdef SIM
float Reset,  BlueBTN;
#endif


float T_ext, T_elt, Tref, Tlim, theta_r;