#ifndef USER_MOTORCONTROL_FUNCTIONS_H
#define USER_MOTORCONTROL_FUNCTIONS_H

#include "math.h"
#include "User_Constants.h"
#include "User_data_types.h"
#include "User_Macros.h"
//#include "../Inc/User_Variables.h"
//#include "../Inc/MotorData.h"

void PIReg(XPIRegPars *par,XPIRegVars *var);
void PIRegAsy(XPIRegPars *par,XPIRegVars *var, float lim_max, float lim_min);
void PWMduty(Xabc vsabc_ref, float vdc, Xabc* duty_abc);
void speed_compute_sc(Xsc sincos, Xsc *sincos_old, float* omega);
void DTComp1(Xabc isabc, float amp_dt, Xabc *duty_abc);
void DTComp2(Xabc isabc, float amp_dt_V, Xabc *vsabc);
void DTComp(Xabc duty, Xabc duty_km1, Xabc isabc,float vdc,float dt, Xabc *vsabc);
void Current_loop(float vdc, float Imax, Xdq isdq_ref, Xdq  isdq, Xdq vffw_dq,XPIRegPars* id_par, XPIRegVars* id_var, XPIRegPars* iq_par, XPIRegVars* iq_var,Xdq* vsdq_ref);
void Gen_theta_ref(float omega_ref_ramp, float* theta_ref, Xsc* SinCos_ref );
void CurrentProtection(Xabc isabc, int* State, int* pwm_stop);
void ramp(float target, float delta, float *output);
void HF_position_detect2(Xalphabeta isab, int commiss_counter,int* counter ,float* theta_hf,Xalphabeta* vsab_ref);
void ReadLut(float *tab0, float Xin, float Xmax,float Xmin, float DX, float inv_DX, float* Yout);
void interp2d(float *P, float x, float y, float Dx, float invDx, float Dy, float invDy, float Xmax, float Xmin, float Ymax , float Ymin, int Npointx, float* V);
void ReadLut2d(float *tab0,float x,float y,float Dx,float invDx,float Dy,float invDy,float Xmax,float Xmin,float Ymax,float Ymin,int Npointx,float* V);
float sgn(float i);
void FluxObserver(void);
void Compute_Inductance(void);
void HF_dem(float *sig, float *sig_HF, float *acc, float *buffer, int counterHF, int ns);

#endif

