#include "../Inc/User_MotorControl_Functions.h"

void PIReg(XPIRegPars *par,XPIRegVars *var)
{
    float int_lim;
    
    //Error computation
    var->err  = var->ref-var->fbk;    //Proportional regulator
    var->prop = par->kp*var->err;
    //Saturation of the proportional part to lim
    if (var->prop>par->lim)
        var->prop=par->lim;
    if (var->prop<(-par->lim))
        var->prop=-par->lim;
    
    //Limit of the intg part
    int_lim = par->lim-fabs(var->prop);
    
    //Integral part
    var->intg+=par->ki*var->err;
    
    //Saturation of the intg part
    if (var->intg > int_lim)
        var->intg = int_lim;
    if (var->intg < (-int_lim))
        var->intg = -int_lim;
    
    //Output computation
    var->out=var->prop+var->intg;
}

void PIRegAsy(XPIRegPars *par,XPIRegVars *var, float lim_max, float lim_min)
{
    float int_lim_max;
    float int_lim_min;
    
    //Error computation
    var->err  = var->fbk-var->ref;
    //var->err  = var->ref-var->fbk;  //Proportional regulator
    var->prop = par->kp*var->err;
    //Saturation of the proportional part to lim
    if (var->prop>lim_max)
        var->prop=lim_max;
    if (var->prop<lim_min)
        var->prop=lim_min;
    
    //Limit of the intg part
    int_lim_max = lim_max-fabs(var->prop);
    int_lim_min = lim_min-fabs(var->prop);
    //Integral part
    var->intg+=par->ki*var->err;
    
    //Saturation of the intg part
    if (var->intg > int_lim_max)
        var->intg = int_lim_max;
    if (var->intg < int_lim_min)
        var->intg = int_lim_min;
    
    //Output computation
    var->out=var->prop+var->intg;
}

//----------------------------------------------------------------------------------------------------//

void PWMduty(Xabc vsabc_ref, float vdc, Xabc* duty_abc) {
    
		float tmp1, tmp2, tmp3, vdc_inv, vzs;
		
		if (vdc>0)	vdc_inv = 1.0f/vdc;	// avoid division by zero
		
		// zero-sequence detection:
		// tmp1 = max,  tmp2 = min
		// tmp3 is the phase in between: tmp3 = -(max + min)
		
		tmp1 = vsabc_ref.a;	// max
		tmp2 = vsabc_ref.b; // min
		
		if (vsabc_ref.a < vsabc_ref.b) {
				tmp1 = vsabc_ref.b;		// max
				tmp2 = vsabc_ref.a;		// min
		}
		
		if (tmp1<vsabc_ref.c) tmp3 = tmp1;
		else {
				if (tmp2>vsabc_ref.c) tmp3 = tmp2;
				else 	tmp3 = vsabc_ref.c;
		}
		
		vzs=tmp3*0.5f;	// zero sequence voltage
				
		// Duty-cycles, with zero-sequence
		duty_abc->a = 0.5f + (vsabc_ref.a + vzs)*vdc_inv;;
		duty_abc->b = 0.5f + (vsabc_ref.b + vzs)*vdc_inv;;
		duty_abc->c = 0.5f + (vsabc_ref.c + vzs)*vdc_inv;;  
        
//         duty_abc->a += deadtime/Ts*sgn(isabc.a);
//         duty_abc->b += deadtime/Ts*sgn(isabc.b);
//         duty_abc->c += deadtime/Ts*sgn(isabc.c);
}

//----------------------------------------------------------------------------------------------------//

void speed_compute_sc(Xsc sincos, Xsc *sincos_old, float* omega) { 
    
	  float tmp1, tmp2;
    
    tmp1 = sincos.sin * sincos_old -> cos;
    tmp2 = sincos.cos * sincos_old -> sin;
    
    tmp1 = tmp1-tmp2;
    *omega = fs * tmp1;
    
    sincos_old -> sin = sincos.sin;
    sincos_old -> cos = sincos.cos;
   
}

//----------------------------------------------------------------------------------------------------//

void DTComp1(Xabc isabc, float amp_dt, Xabc *duty_abc) {
    
    if (isabc.a>0) duty_abc->a += amp_dt;
    else duty_abc->a += -amp_dt;
    
    if (isabc.b>0) duty_abc->b += amp_dt;
    else duty_abc->b += -amp_dt;
    
    if (isabc.c>0) duty_abc->c += amp_dt;
    else duty_abc->c += -amp_dt;
    
}

//----------------------------------------------------------------------------------------------------//

void DTComp2(Xabc isabc, float amp_dt_V, Xabc *vsabc) {
    
    if (isabc.a >0) vsabc->a += -amp_dt_V;
    else vsabc->a += amp_dt_V;
    
    if (isabc.b>0) 	vsabc->b += -amp_dt_V;
    else vsabc->b += amp_dt_V;
    
    if (isabc.c>0) vsabc->c += -amp_dt_V;
    else vsabc->c += amp_dt_V;
    
}

//----------------------------------------------------------------------------------------------------//

void DTComp(Xabc duty, Xabc duty_km1, Xabc isabc,float vdc,float dt, Xabc *vsabc) 
{
  float v_error_A, v_error_B, v_error_C, vA, vB, vC;
  float Ir,slopeA,slopeB,slopeC,max,min;
  Ir = 15;
  max = Ir*0.01;
  min = -Ir*0.01;
	
    if(duty.a>1) duty.a = 1.;
	if(duty.a<0) duty.a = 0.;
	if(duty.b>1) duty.b = 1.;
	if(duty.b<0) duty.b = 0.;
	if(duty.c>1) duty.c = 1.;
	if(duty.c<0) duty.c = 0.;

//   v_error_A = vdc*dt/Ts;
// 	if (isabc.a < 0)
// 		v_error_A = -v_error_A;
//   v_error_B = vdc*dt/Ts;
// 	if (isabc.b < 0)
// 		v_error_B = -v_error_B;
//   v_error_C = vdc*dt/Ts;
// 	if (isabc.c < 0)
// 		v_error_C = -v_error_C;
  
    v_error_A = vdc*dt/Ts;
    slopeA = v_error_A/max;
    if(isabc.a<max && isabc.a>min)
        v_error_A = slopeA*(isabc.a-max)+v_error_A;
    else{
        if(isabc.a<min)
            v_error_A = -v_error_A;
    }
  
        
    v_error_B = vdc*dt/Ts;
    slopeB = v_error_B/max;
    if(isabc.b<max && isabc.b>min)
        v_error_B = slopeB*(isabc.b-max)+v_error_B;
    else{
        if(isabc.b<min)
            v_error_B = -v_error_B;
    }
  
    
    
    v_error_C = vdc*dt/Ts;
    slopeC = v_error_C/max;
    if(isabc.c<max && isabc.c>min)
        v_error_C = slopeC*(isabc.c-max)+v_error_C;
    else{
        if(isabc.c<min)
            v_error_C = -v_error_C;
    }
  
	
	if ((duty.a  == 0 && duty_km1.a == 0) || (duty.a  == 1 && duty_km1.a == 1))
       v_error_A = 0;
    if ((duty.b  == 0 && duty_km1.b == 0) || (duty.b  == 1 && duty_km1.b == 1))
       v_error_B = 0;
    if ((duty.c  == 0 && duty_km1.c == 0) || (duty.c  == 1 && duty_km1.c == 1))
       v_error_C = 0;
    
    if ((duty.a == 1 && duty_km1.a == 0 && isabc.a<0) || (duty.a == 0 && duty_km1.a == 1 && isabc.a>0))
       v_error_A = 0;
    if ((duty.b == 1 && duty_km1.b == 0 && isabc.b<0) || (duty.b == 0 && duty_km1.b == 1 && isabc.b>0))
       v_error_B = 0;
    if ((duty.c == 1 && duty_km1.c == 0 && isabc.c<0) || (duty.c == 0 && duty_km1.c == 1 && isabc.c>0))
       v_error_C = 0;
    
    vA	= vdc*(duty.a-0.5)- v_error_A;
    vB	= vdc*(duty.b-0.5)- v_error_B;
    vC	= vdc*(duty.c-0.5)- v_error_C;
    
    vsabc->a = 2.0*vA/3.0-(vB+vC)/3.0;
    vsabc->b = 2.0*vB/3.0-(vA+vC)/3.0;
    vsabc->c = 2.0*vC/3.0-(vA+vB)/3.0;
		
		duty_km1.a = duty.a;
		duty_km1.b = duty.b;
		duty_km1.c = duty.c;
}

//----------------------------------------------------------------------------------------------------//

//Current control loops
void Current_loop(float vdc, float Imax, Xdq isdq_ref, Xdq  isdq, Xdq vffw_dq ,XPIRegPars* id_par, XPIRegVars* id_var, XPIRegPars* iq_par, XPIRegVars* iq_var,Xdq* vsdq_ref) {
	
    float vs_max;
		float tmp1;
	
	
    vs_max = vdc * SQRT1OVER3;
    
    //d-axis control loop
    id_par->lim = vs_max-vffw_dq.d;
    id_var->ref = isdq_ref.d;
    id_var->fbk = isdq.d;
    PIReg(id_par, id_var);
    vsdq_ref->d=id_var->out;
    
    // saturation to Imax
    tmp1 = sqrtf(Imax*Imax-isdq_ref.d*isdq_ref.d);
    if (isdq_ref.q>tmp1)	isdq_ref.q = tmp1;
    if (isdq_ref.q<-tmp1)	isdq_ref.q =-tmp1;
    
    //q-axis control loop
    iq_par->lim   =sqrtf(vs_max*vs_max-vsdq_ref->d*vsdq_ref->d)-vffw_dq.q;
    iq_var->ref   =isdq_ref.q;
    iq_var->fbk   =isdq.q;
    PIReg(iq_par, iq_var);
    vsdq_ref->q=iq_var->out;
    
}

//----------------------------------------------------------------------------------------------------//

// General formulation for a low pass filter with only one coefficient
// sarebbe da sostituire con una macro
//float Filter(float xk,float xf_k, float kfilt) 
//{
//    
//    xf_k+= kfilt * (xk-xf_k);
//    
//    return xf_k;
//}

//----------------------------------------------------------------------------------------------------//

//Reference angle generation (V/Hz, I-Hz)
void Gen_theta_ref(float omega_ref_ramp, float* theta_ref, Xsc* SinCos_ref )
{
    *theta_ref += Ts * omega_ref_ramp;
    if (*theta_ref > TWOPI)	*theta_ref -= TWOPI;
    if (*theta_ref < 0.0f)	*theta_ref += TWOPI;
    
    SinCos_ref->sin = sinf(*theta_ref);
    SinCos_ref->cos = cosf(*theta_ref);
}

//----------------------------------------------------------------------------------------------------//

void CurrentProtection(Xabc isabc, int* State, int* pwm_stop) {
    
    
    if (fabs(isabc.a)>CRT_PROT) {
        *pwm_stop = 1;
        *State=ERROR;
    }
    if (fabs(isabc.b)>CRT_PROT) {
        *pwm_stop = 1;
        *State=ERROR;
    }
    if (fabs(isabc.c)>CRT_PROT) {
        *pwm_stop = 1;
        *State=ERROR;
    }
    
}

//----------------------------------------------------------------------------------------------------//

void ramp(float target, float delta, float *output) {
    
    if (*output <= target) {
        *output = (*output + delta);
        if (*output>target)
            *output=target;
    }
    else {
        *output = (*output - delta);
        if (*output<target)
            *output=target;
    }
}

//----------------------------------------------------------------------------------------------------//

void HF_position_detect2(Xalphabeta isab, int commiss_counter,int* counter ,float* theta_hf,Xalphabeta* vsab_ref) {
    //==========================================================================
    //	High Frequency Position Detection
    //	Rotating voltage --> Imax detection
    //==========================================================================
    
	float IHFamp, IHFmax, theta0;
	
	// HF_position_detect	
   #define	nsamples_HF  20.025f	// fPWM/fHF
	//#define	20.0f	// fPWM/fHF
    
	#ifdef SPM
		#define V_HF  20.0f			// HF voltage amplitude for SPM motor
	#else
		#define V_HF 120.0f			// HF voltage amplitude for PMA motor
	#endif
	
    IHFamp=isab.alpha*isab.alpha+isab.beta*isab.beta;	// Current (amplitude)^2
    
    // find maximum amplitude position (after 400 samples)
    if((IHFamp>=IHFmax)&&(commiss_counter>400))	{
        IHFmax = IHFamp;
#if defined(SPM)
        theta0 = *theta_hf - HALFPI;
#else
        theta0 = *theta_hf - PI;
#endif
    }
    
    *theta_hf += TWOPI/nsamples_HF;	    // HF phase angle on nsamples_HF positions
    if (*theta_hf>TWOPI) *theta_hf -= TWOPI;
    
    if (commiss_counter<1000) {
        // Voltage signals
        *counter=*counter+1;
        if(*counter>200) *counter = 200;
        vsab_ref->alpha = (*counter/200.0f)*V_HF*cosf(*theta_hf);
        vsab_ref->beta  = (*counter/200.0f)*V_HF*sinf(*theta_hf);
    }
    else {
        // stop injection
        *counter = 0;
        vsab_ref->alpha = 0.0f;
        vsab_ref->beta  = 0.0f;
    }
}

//----------------------------------------------------------------------------------------------------//

void  ReadLut(float *tab0, float Xin, float Xmax,float Xmin, float DX, float inv_DX, float* Yout)
{

    int X0;
    
    //Limit input range
    if (Xin > Xmax) Xin = Xmax;
    if (Xin < Xmin) Xin = Xmin;
    
    //Select element X0
    X0=floor((Xin-Xmin)*inv_DX);
    //address of X0
    tab0 = tab0 + X0;
    //Interpolate between Y(tab0) and Y(tab0+1)
    *Yout=*tab0 +(*(tab0+1)-*tab0)*inv_DX*(Xin-Xmin-DX*X0);
    
}

//----------------------------------------------------------------------------------------------------//

// Read a 2 dimensional look-up-table (matrix)
void interp2d(float *P, float x, float y, float Dx, float invDx, float Dy, float invDy, float Xmax, float Xmin, float Ymax , float Ymin, int Npointx, float*  V) {
    float ix, iy;
    float rx, ry;
    float V1, V2;
    int delta;
    // Funzione per l'interpolazione di Tabelle 2D
    // Limitazione degli ingressi
    if (x > Xmax) {x=Xmax-0.1f*Dx;}
    if (x < Xmin) {x=Xmin+0.1f*Dx;;}
    if (y > Ymax) {y=Ymax-0.1f*Dy;}
    if (y < Ymin) {y=Ymin+0.1f*Dy;}
    //calcolo delle posizione inziale per l'interpolazione
    rx=(x-Xmin)*invDx;
    ry=(y-Ymin)*invDy+1;
    // indici interi
    ix=floor(rx);
    iy=floor(ry);
    // resti per l'interpolazione
    rx=rx-ix;
    ry=ry-iy;
    //primo punto
    delta=(ix)+(iy-1)*Npointx;
    P=P+delta;
    //interpolazione quadrata
    V1=*P+(*(P+1)-*P)*rx;
    P=P+Npointx;
    V2=*P+(*(P+1)-*P)*rx;
    *V=V1+(V2-V1)*ry;
    
}

//----------------------------------------------------------------------------------------------------//

// Read a 2 dimensional look-up-table (matrix)
void ReadLut2d(float *tab0,float x,float y,float Dx,float invDx,float Dy,float invDy,float Xmax,float Xmin,float Ymax,float Ymin,int Npointx,float* V)
{
    float ix, iy;
    float rx, ry;
    float V1, V2;
    int delta;
   
	//Limit input range
	if (x > Xmax) {x=Xmax-0.1f*Dx;}
    if (x < Xmin) {x=Xmin+0.1f*Dx;}
    if (y > Ymax) {y=Ymax-0.1f*Dy;}
    if (y < Ymin) {y=Ymin+0.1f*Dy;}
    //Select element 0 for x and y
    rx=(x-Xmin)*invDx;
	ry=(y-Ymin)*invDy;

    // indexes of element 0
    ix=floor(rx);
    iy=floor(ry);
    // distance from element 0
    rx=rx-ix;
    ry=ry-iy;
    //address of element 0
	delta=(ix)+(iy)*Npointx;
   
    tab0=tab0+delta;
    //interpolate on x
	//V1 is exact on x and FLOOR along y
    V1=*tab0+(*(tab0+1)-*tab0)*rx;
	//V2 is exact on x and CEIL along y
    tab0=tab0+Npointx;
    V2=*tab0+(*(tab0+1)-*tab0)*rx;
	// interpolate on y
    *V=V1+(V2-V1)*ry;
}

//----------------------------------------------------------------------------------------------------//

float sgn(float i)
{
 float sgn;
    if (i>0.000001)
        sgn=1.0;
    else
     if(i<-0.000001)
		sgn=-1.0;
     else
         sgn=0;
    

 return (sgn);
}

//----------------------------------------------------------------------------------------------------//

// Bandpass Fitler

void HF_dem(float *sig, float *sig_HF, float *acc, float *buffer, int counterHF, int ns)
{
	*acc += (*sig - *(buffer+counterHF))/ns;
	*(buffer+counterHF) = *sig;
	*sig_HF = *sig - *acc;
}
