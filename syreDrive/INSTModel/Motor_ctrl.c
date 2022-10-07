// // /*----------------------------------------------------------------------*/
/* File: Syr_ctrl.c                                               	*/
/* the parameter CTRL_TYPE decides the type of control		*/
/*  CTRL_TYPE = 0 - Current control												*/
/*  CTRL_TYPE = 1 - Flux control											  	*/
/*  CTRL_TYPE = 2 - Torque control												*/
/*  CTRL_TYPE = 3 - Speed control													*/
/*----------------------------------------------------------------------*/
/* Motor type: 	                                   */
/*----------------------------------------------------------------------*/

#define S_FUNCTION_NAME Motor_ctrl
#define S_FUNCTION_LEVEL 2

#define U(element) (*uPtrs[element])

#include "simstruc.h"

// USER CODE BEGINS HERE

#define SIM

#include "math.h"
#include "User_functions\Inc\User_data_types.h"
#include "User_functions\Inc\User_Variables.h"
#include "User_functions\Inc\User_Constants.h"
#include "User_functions\Inc\User_Macros.h"
#include "User_functions\Inc\User_MotorControl_Functions.h"
#include "User_functions\Inc\MotorData.h"

//#include "User_functions\Inc\VS_SimpleP36.h"
//#include "User_functions\Inc\stm32f303xe.h"	// valid for F303 MCU

// USER CODE ENDS HERE

#define NINPUTS 	20
#define NOUTPUTS 	34
#define NPARAMS 	0


static void mdlInitializeSizes(SimStruct *S) {
    ssSetNumSFcnParams(S, NPARAMS);  /* Number of expected parameters */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        return; /* Parameter mismatch will be reported by Simulink */
    }
    
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);
    
    if (!ssSetNumInputPorts(S, 1)) return;
    ssSetInputPortWidth(S, 0, NINPUTS);
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    
    if (!ssSetNumOutputPorts(S, 1)) return;
    ssSetOutputPortWidth(S, 0, NOUTPUTS);
    
    ssSetNumSampleTimes(S, 1);
    ssSetNumRWork(S, 0);
    ssSetNumIWork(S, 0);
    ssSetNumPWork(S, 0);
    ssSetNumModes(S, 0);
    ssSetNumNonsampledZCs(S, 0);
    
    /* Take care when specifying exception free code - see sfuntmpl_doc.c */
    ssSetOptions(S, SS_OPTION_EXCEPTION_FREE_CODE);
}

static void mdlInitializeSampleTimes(SimStruct *S) {
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);/*inherited*/
    ssSetOffsetTime(S, 0, 0.);
}

#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================*/
static void mdlInitializeConditions(SimStruct *S)
/*this routine is called at the start of simulation */
{
    
}

// controllo inizia qui
static void mdlOutputs(SimStruct *S, int_T tid)
/*this routine compute the outputs of S-Function block */
{
    
    double            *y      = ssGetOutputPortRealSignal(S, 0);
    InputRealPtrsType uPtrs   = ssGetInputPortRealSignalPtrs(S, 0);
    
    
    // ISR - USER CODE BEGINS HERE
    

    
    isabc.a   			= U(0);
    isabc.b   			= U(1);
    isabc.c   			= U(2);
    vdc  	  	  		= U(3);
    
    theta_mec_meas  	= U(4);			// encoder
    n_ref_in   			= U(5);			// rpm
    T_ext				= U(6);			// reference torque
    isdq_ext.d 			= U(7);  		// Current reference
    isdq_ext.q 			= U(8);  		// Current reference
	
    Reset 				= U(9);			// Black button
    Go    				= U(10);	  // Blue  button
		
	Ctrl_type			= U(11);			
    inj_waveform		= U(12);		// Injected waveform (sinusoidal or squarewave)
    dem 			  	= U(13);		// Current or Flux demoduation
	HS_ctrl		  		= U(14);		// High speed position error estimation technique (AF or APP)
	SS_on				= U(15);		// Sensorless ON or OFF
	accel				= U(16); 		// speed acceleration rpm/s
	Quad_Maps			= U(17);
    lambda_M            = U(18);
    th0                 = U(19);

    
  
    
   
	switch(State){
        
        case ERROR:
            // Variables Init
            pwm_stop 					= 1;
            counter 					= 0;
            
            offset_current_a			= 2120;
            offset_current_b			= 2120;
            offset_current_c			= 2120;
            
            n_ref_in 					= 0.0f;
            omega_ref_in				= 0.0f;
            omega_ref_ramp				= 0.0f;
            accel               		= 1000;	// rpm/s
			omega_elt_meas           	= 0.0f;
			omega_elt_meas_f         	= 0.0f;
            SinCos_elt_meas.sin      	= 0.0f;
            SinCos_elt_meas.cos      	= 1.0f;
            SinCos_elt_meas_old.sin  	= 0.0f;
            SinCos_elt_meas_old.cos  	= 1.0f;

			Ld 							= Ld_inic;
			Lq							= Lq_inic;
			ld 							= Ld_inic;
			lq 							= Lq_inic;
			flux_nom 					= 0.0f;	// rated flux (Vs)
            v0 							= 0.0f;			// phase dc voltage (V)
			isdq_ref.d 					= 0.0f;
            isdq_ref.q 					= 0.0f;
            isdq.d                      = 0.0f;
            isdq.q                      = 0.0f;
            vsdq_ref.d                  = 0.0f;
            vsdq_ref.q                  = 0.0f;
            
            lambda_CM_dq.d              = 0.0f;
            lambda_CM_dq.q              = 0.0f;
            
            switch (Quad_Maps){
				case 0:
                    lambda_dq.d         = 0.0f;
                    lambda_dq.q         = 0.0f;
                    lambda_CM_dq.d      = 0.0f;
                    lambda_CM_dq.q      = 0.0f;
				    lambda_OBS.alpha    = 0.0f;
                    lambda_OBS.beta     = 0.0f;
				break;
				case 1:
                    lambda_dq.d         = 0.0f;
                    lambda_dq.q         = -lambda_M;
                    lambda_CM_dq.d      = 0.0f;
                    lambda_CM_dq.q      = -lambda_M;
                    lambda_OBS.alpha    = lambda_M*sin(PP*th0);
                    lambda_OBS.beta     = -lambda_M*cos(PP*th0);
				break;
				case 2:
			        lambda_dq.d         = lambda_M;
                    lambda_dq.q         = 0.0f;
                    lambda_CM_dq.d      = lambda_M;
                    lambda_CM_dq.q      = 0.0f;
                    lambda_OBS.alpha    = lambda_M*cos(PP*th0);
                    lambda_OBS.beta     = lambda_M*sin(PP*th0);
				break;
            }
            
            isabc.a                     = 0.0f;
            isabc.b                     = 0.0f;
            isabc.c                     = 0.0f;
            
            duty_abc.a					= 0.0f;
            duty_abc.b					= 0.0f;
            duty_abc.c					= 0.0f;
            
            offset_in.ch0				= 0.0f;
            offset_in.ch1				= 0.0f;
            offset_in.ch2				= 0.0f;
            offset_in.ch3				= 0.0f;
						
            if(Go) State = WAKE_UP;	// Wait for Go state
            
            break;
            
        case WAKE_UP:
            
            duty_abc.a = 0.5f;
            duty_abc.b = 0.5f;
            duty_abc.c = 0.5f;
            
            pwm_stop = 0;
            
            // Current offset accumulation
            
            if (counter>100){
                offset_in.ch0+=input.ch0;
                offset_in.ch1+=input.ch1;
                offset_in.ch2+=input.ch2;
            }
            
            // Offset computation and bootstrap
                     
            if (counter==(100+200)){ 
                offset_current_a=(float)(offset_in.ch0/200.0f);
                offset_current_b=(float)(offset_in.ch1/200.0f);
                offset_current_c=(float)(offset_in.ch2/200.0f);
            }
            						
			switch (inj_waveform) {
				case 0: // sinusoidal 
					OMEGA_B_PLL = WB_PLL_PULS; 
					break;
				case 1: // squarewave
					OMEGA_B_PLL = WB_PLL_SQUARE; 
					break;
			}
						
            // Switch to READY state
            
            if (counter > 300){   
                counter=0;
                State = READY;
            }
            
            counter++;
            
            break;
            
        case READY:
            
            duty_abc.a=0.5f;
            duty_abc.b=0.5f;
            duty_abc.c=0.5f;
            counter = 0;
			PLL_var.intg = omega_elt_meas_f;
			theta_PLL = theta_elt_meas;
            if (Go)	State = START;
      
            break;
            
        case START:
			
			//-------------------Speed Compute----------------------------------//
			
			theta_elt_meas      = PP * theta_mec_meas;
			while(theta_elt_meas > PI)
				theta_elt_meas -= TWOPI;
			while(theta_elt_meas < -PI)
				theta_elt_meas += TWOPI;
			SinCos_elt_meas.sin = sin(theta_elt_meas);
			SinCos_elt_meas.cos = cos(theta_elt_meas);
				
			speed_compute_sc(SinCos_elt_meas, &SinCos_elt_meas_old, &omega_elt_meas);
			_Filter(omega_elt_meas, omega_elt_meas_f, Ts*TWOPI*10);
			omega_mec_meas      = omega_elt_meas/PP;
			omega_mec_meas_f    = omega_elt_meas_f/PP;
			omega_mec_meas_rpm  = omega_mec_meas_f*30/PI;
				
			SinCos_elt_meas_old.sin = SinCos_elt_meas.sin;
			SinCos_elt_meas_old.cos = SinCos_elt_meas.cos;
								
			//---------------------------------PLL-------------------------------//
			
			if(SS_on) {	
				PLL_par.kp   = 2*OMEGA_B_PLL;
				PLL_par.ki   = pow(OMEGA_B_PLL,2)*Ts;
				PLL_par.lim  = RPM2RAD * nmax_mot * PP;
				PLL_var.ref  = pos_err;
				PLL_var.fbk  = 0;
				
				PIReg(&PLL_par, &PLL_var);
				omega_PLL  = PLL_var.out;
				theta_PLL += Ts*PLL_var.out; 
				_Filter(omega_PLL, omega_elt, Ts*TWOPI*25);
				
				if(theta_PLL >= TWOPI)
					theta_PLL -= TWOPI;
				if (theta_PLL < 0)
					theta_PLL +=TWOPI;
				
				SinCos_elt.sin = sin(theta_PLL);
				SinCos_elt.cos = cos(theta_PLL);
				position_error_real = asin(sin(theta_elt_meas - theta_PLL));

			}		
			
			else { // Encoder
				SinCos_elt.sin = sin(theta_elt_meas);
				SinCos_elt.cos = cos(theta_elt_meas);
				omega_elt = omega_elt_meas_f;
				position_error_real = 0;
			}
			
			//-------------------------------Control Type ------------------------------//										
							
			switch (Ctrl_type){
			
				case 0: //CurrentControl
					isdq_ref.d = isdq_ext.d;
					isdq_ref.q = isdq_ext.q;
				break;
															
				case 2: //TorqueControl
					ReadLut(&ID_REF[0], fabs(T_ext), TMAX, TMIN, DT, INV_DT, &isdq_ref.d); 
					ReadLut(&IQ_REF[0], fabs(T_ext), TMAX, TMIN, DT, INV_DT, &isdq_ref.q);
					switch (Quad_Maps){
						case 0:
							isdq_ref.d = sgn(T_ext)*isdq_ref.d; 
						break;
						
						case 1:
							isdq_ref.d = sgn(T_ext)*isdq_ref.d; 
						break;
						
						case 2:
							isdq_ref.q = sgn(T_ext)*isdq_ref.q;
						break;
					}																
				break;

				case 3: //SpeedControl
					omega_ref_in = n_ref_in * RPM2RAD;
					ramp(omega_ref_in, accel * RPM2RAD*Ts, &omega_ref_ramp);
					sp_var.ref  = omega_ref_ramp;
					sp_var.fbk  = omega_elt/PP;
					sp_par.lim  = T_rated;
					kp_w        = 2*OMEGA_BW*J;
					ki_w        = pow(OMEGA_BW,2)*J;
					sp_par.ki   = ki_w*Ts;
					sp_par.kp   = kp_w;
					PIReg(&sp_par, &sp_var);
					T_ext       = sp_var.out;
																	
					ReadLut(&ID_REF[0], fabs(T_ext), TMAX, TMIN, DT, INV_DT, &isdq_ref.d); 
					ReadLut(&IQ_REF[0], fabs(T_ext), TMAX, TMIN, DT, INV_DT, &isdq_ref.q);
					switch (Quad_Maps){
						case 0:
							isdq_ref.d = sgn(T_ext)*isdq_ref.d; 
						break;
						
						case 1:
							isdq_ref.d = sgn(T_ext)*isdq_ref.d; 
						break;
						
						case 2:
							isdq_ref.q = sgn(T_ext)*isdq_ref.q;  
						break;
					}							
				break;
			}
														
			_clarke(isabc, isab);
			_rot(isab, SinCos_elt, isdq);
						
			vsab_km1 = vsab_k0;
			DTComp(duty_abc, duty_abc_km1, isabc, vdc, deadtime, &vsabc_k0);
			_clarke(vsabc_k0, vsab_k0)
			_rot(vsab_k0,SinCos_elt,vsdq_k0);
			
			FluxObserver();
			//Compute_Inductance();
			
			//---------------------------Sensorless Position Estimation--------------------------//
			
			pos_err_LS = 0;
			pos_err_HS = 0;

			// HF demodulation
			if (abs(omega_elt) < KOBS + 1.2*OMEGA_G) {
				switch (dem) {
						case 0: // current dem
							switch(inj_waveform) {
								
								case 0: // sinusoidal
									HF_dem(&isdq.q, &ishf, &acc_id, &buffer_id[0], counterHF, ns);
									tmp3 = 2*ishf*sin(theta_inj - TWOPI*f_inj_puls*Ts); 
									_Filter(tmp3, ishf_mag, Ts*OMEGA_0_INJ);         // k = Ts*TWOPI*f0
									k_e  = - V_inj/(TWOPI*f_inj_puls) * (ld-lq)/(2*ld*lq);
									pos_err_LS = ishf_mag / k_e;
									break;
									
								case 1: // squarewave
									ishf = isdq.q - isq_old;
									ishf_mag = ishf*cos(theta_inj - TWOPI*f_inj_sq*Ts) + ishf_mag_old; 
									ishf_mag_old = ishf*cos(theta_inj - TWOPI*f_inj_sq*Ts);
									k_e  = V_inj/f_inj_sq * (ld-lq)/(2*ld*lq);
									pos_err_LS = ishf_mag / k_e;
									isq_old = isdq.q;
									break;
							}
							break;
							
						case 1: // flux dem
							switch(inj_waveform) {
								
								case 0: // sinusoidal
									HF_dem(&lambda_CM_dq.q, &lambdahf, &acc_id, &buffer_id[0], counterHF, ns);
									tmp3 = 2*lambdahf*sin(theta_inj - TWOPI*f_inj_puls*Ts); 
									_Filter(tmp3, lambdahf_mag, Ts*OMEGA_0_INJ);         // k = Ts*TWOPI*f0
									k_e  = - V_inj/(TWOPI*f_inj_puls) * (lq*ldm - ldq*ldq)/(ld*lq - ldq*ldq);
									pos_err_LS = lambdahf_mag / k_e;	
									break;
									
								case 1: // squarewave
									lambdahf = lambda_CM_dq.q - lambdaq_old;
									lambdahf_mag = lambdahf*cos(theta_inj - TWOPI*f_inj_sq*Ts); 
									k_e  = V_inj/f_inj_sq * (lq*ldm - ldq*ldq)/(ld*lq - ldq*ldq);
									pos_err_LS = lambdahf_mag / k_e;	
									lambdaq_old = lambda_CM_dq.q;
									break;
							}								
							break;
				}
			}

			lambda_aux.d = (ld - Lq)*isdq.q - ldq*isdq.d;
			lambda_aux.q = (Ld - lq)*isdq.d + ldq*isdq.q;
		  
			// High speed control
			switch (HS_ctrl) {
				case AF: // Active Flux
					if (f_omega > 0.01) {
						if (abs(Ldm) < 1e-3)
							Ldm = 1e-3;
						if (abs(isdq.d) < 1e-3)
							isdq.d = 1e-3*sgn(isdq.d);
						k_AF = 1/(2*Ldm*isdq.d);
						pos_err_HS = k_AF * (lambda_dq.q - Lq*isdq.q);
					}
				break;
				
				case APP:
					if (abs(omega_elt) > KOBS - 2*OMEGA_G) {
						lambda_diff.d = lambda_dq.d - lambda_CM_dq.d;
						lambda_diff.q = lambda_dq.q - lambda_CM_dq.q;
						k_APP = 1 / (omega_elt * (pow(lambda_aux.d,2) + pow(lambda_aux.q,2)));
						pos_err_HS = k_APP * (KOBS * (lambda_aux.d*lambda_diff.q - lambda_aux.q*lambda_diff.d) + omega_elt_meas_f * (lambda_aux.d*lambda_diff.d + lambda_aux.q*lambda_diff.q));
					}
				break;
			} 
			if (abs(omega_elt) < KOBS - OMEGA_G) 
				f_omega = 0;
			else if (abs(omega_elt) > KOBS + OMEGA_G)
				f_omega = 1;	
			else
				f_omega = (OMEGA_G+abs(omega_elt)-KOBS)/(2*OMEGA_G);
			
			pos_err = f_omega * pos_err_HS + (1-f_omega)*pos_err_LS;
			
			counterHF++;
			if (counterHF >= ns)
				counterHF = 0;
			
			//-----------------------Current Vector Control--------------------//
			
// 			kp_id = OMEGA_BI*ld;
// 			ki_id = kp_id*OMEGA_BI/10;
// 			kp_iq = OMEGA_BI*lq;
// 			ki_iq = kp_iq*OMEGA_BI/10;
//             
            kp_id = OMEGA_BI*Ld_inic;
			ki_id = kp_id*OMEGA_BI/10;
			kp_iq = OMEGA_BI*Lq_inic;
			ki_iq = kp_iq*OMEGA_BI/10;

			id_par.kp   = kp_id;
			id_par.ki   = ki_id*Ts;
			iq_par.kp   = kp_iq;
			iq_par.ki   = ki_iq*Ts;

			//Current loop
			Current_loop(vdc, Imax_mot, isdq_ref, isdq, &id_par, &id_var, &iq_par, &iq_var, &vsdq_ref);
			vsdq_ref.d += RS*isdq.d - omega_elt*lambda_dq.q;
			vsdq_ref.q += RS*isdq.q + omega_elt*lambda_dq.d;


            //Phase advance

            dTheta = 1.5f*omega_elt*Ts;
            SinCos_elt_dTheta.sin = SinCos_elt.sin*cosf(dTheta) +SinCos_elt.cos*sinf(dTheta);
            SinCos_elt_dTheta.cos = SinCos_elt.cos*cosf(dTheta) -SinCos_elt.sin*sinf(dTheta);
			
			// HF injection
			if (abs(omega_elt) < KOBS + 1.2*OMEGA_G && SS_on == 1) {
				switch(inj_waveform) {				
					case 0: // pulsating
						theta_inj += TWOPI*f_inj_puls*Ts;
						if (theta_inj > TWOPI)
							theta_inj -= TWOPI;
						v_inj = V_inj * cos(theta_inj);
						vsdq_ref.d += v_inj;
						break;
						
					case 1: // squarewave
						theta_inj += TWOPI*f_inj_sq*Ts;
						if (theta_inj > TWOPI)
							theta_inj -= TWOPI;  
						vsdq_ref.d += vhf;
						vhf = vhf * (-1);
						break;
				}
			}	
			
		_invrot(vsdq_ref, SinCos_elt, vsab_ref);
		_invclarke(vsab_ref, vsabc_ref);
		PWMduty(vsabc_ref, vdc,&duty_abc);								
	
	}
    
    //duty cycles saturation
    if (duty_abc.a > 0.99f) duty_abc.a=0.99f;
    if (duty_abc.b > 0.99f) duty_abc.b=0.99f;
    if (duty_abc.c > 0.99f) duty_abc.c=0.99f;
    
    if (duty_abc.a < 0.01f) duty_abc.a=0.01f;
    if (duty_abc.b < 0.01f) duty_abc.b=0.01f;
    if (duty_abc.c < 0.01f) duty_abc.c=0.01f;
    
    if(pwm_stop){
        duty_abc.a=0.0f;
        duty_abc.b=0.0f;
        duty_abc.c=0.0f;
    }
    
#ifndef SIM
    TIM1->CCR1=TIM1_cycle*duty_abc.a;
    TIM1->CCR2=TIM1_cycle*duty_abc.b;
    TIM1->CCR3=TIM1_cycle*duty_abc.c;
    
    //Virtual Scope
    VS_simpleP36_1f(duty_abc.a);
    //    VS_simpleP36_2f(duty_abc.a,duty_abc.b);
    //    VS_simpleP36_4f(duty_abc.a,duty_abc.b,ia_var.ref,isabc.a);
#endif
    
    /* Outputs */
    y[0]  = duty_abc.a;
    y[1]  = duty_abc.b;
    y[2]  = duty_abc.c;
    y[3]  = pwm_stop;
    y[4]  = omega_ref_ramp*60/TWOPI;
    y[5]  = omega_mec_meas_f*60/TWOPI;
    y[6]  = vsdq_ref.d;
    y[7]  = vsdq_ref.q;
    y[8]  = isab.alpha;
    y[9]  = isab.beta;
    y[10] = isdq_ref.d;
    y[11] = isdq.d;
    y[12] = isdq_ref.q;
    y[13] = isdq.q;
	y[14] = f_omega;
	y[15] = T_elt;
	y[16] = lambda_dq.d;
	y[17] = lambda_dq.q;
	y[18] = lambda_CM_dq.d;
	y[19] = lambda_CM_dq.q;
	y[20] = T_ext;
	y[21] = theta_PLL*180/PI;
	y[22] = omega_elt/PP*60/TWOPI;
	y[23] = ld;
	y[24] = lq;
	y[25] = ldq;
	y[26] = isabc.a;
	y[27] = isabc.b;
	y[28] = isabc.c;
	y[29] = theta_elt_meas*180/PI;
	y[30] = pos_err*180/PI;
	y[31] = position_error_real*180/PI;
	y[32] = pos_err_LS*180/PI;
	y[33] = pos_err_HS*180/PI;
   
}

/* Function: mdlTerminate ===================================================== */
static void mdlTerminate(SimStruct *S) {
    UNUSED_ARG(S); /* unused input argument */
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif

//----------------------------------------------------------------------------------------------------//

void FluxObserver(void) {
	
	lambda_CM_ab_km1 = lambda_CM_ab;
	
						
		switch (Quad_Maps){
		case 0:
			ReadLut2d(&FD_LUT[0][0], fabs(isdq.d), fabs(isdq.q), DIDD, INV_DIDD, DIQD, INV_DIQD , ID_TAB_MAX, ID_TAB_MIN, IQ_TAB_MAX , IQ_TAB_MIN, n_size, &lambda_CM_dq.d);
			ReadLut2d(&FQ_LUT[0][0], fabs(isdq.q), fabs(isdq.d), DIQQ, INV_DIQQ, DIQD, INV_DIQD , IQ_TAB_MAX, IQ_TAB_MIN, ID_TAB_MAX , ID_TAB_MIN, n_size, &lambda_CM_dq.q);
			if (isdq.d < 0)
				lambda_CM_dq.d = -lambda_CM_dq.d;
			if (isdq.q < 0)
				lambda_CM_dq.q = -lambda_CM_dq.q;
		break;
		
		case 1:
			ReadLut2d(&FD_LUT[0][0], fabs(isdq.d), isdq.q, DIDD, INV_DIDD, DIQD, INV_DIQD , ID_TAB_MAX, ID_TAB_MIN, IQ_TAB_MAX , IQ_TAB_MIN, n_size, &lambda_CM_dq.d);
			ReadLut2d(&FQ_LUT[0][0], isdq.q, fabs(isdq.d), DIQQ, INV_DIQQ, DIDQ, INV_DIDQ , IQ_TAB_MAX, IQ_TAB_MIN, ID_TAB_MAX , ID_TAB_MIN, n_size, &lambda_CM_dq.q);
			if (isdq.d < 0)
				lambda_CM_dq.d = -lambda_CM_dq.d;
		break;
		
		case 2:
			ReadLut2d(&FD_LUT[0][0], isdq.d, fabs(isdq.q), DIDD, INV_DIDD, DIQD, INV_DIQD , ID_TAB_MAX, ID_TAB_MIN, IQ_TAB_MAX , IQ_TAB_MIN, n_size, &lambda_CM_dq.d);
			ReadLut2d(&FQ_LUT[0][0], fabs(isdq.q), isdq.d, DIQQ, INV_DIQQ, DIDQ, INV_DIDQ , IQ_TAB_MAX, IQ_TAB_MIN, ID_TAB_MAX , ID_TAB_MIN, n_size, &lambda_CM_dq.q);
			if (isdq.q < 0)
				lambda_CM_dq.q = -lambda_CM_dq.q;
		break;
	}

	_invrot(lambda_CM_dq, SinCos_elt, lambda_CM_ab);
	
	feedback_OBS.alpha	=	lambda_CM_ab_km1.alpha - lambda_OBS.alpha;
	feedback_OBS.beta		=	lambda_CM_ab_km1.beta - lambda_OBS.beta;
    
	// Integration
	lambda_OBS.alpha	+= Ts*(vsab_km1.alpha - RS*isab.alpha + KOBS*feedback_OBS.alpha);
	lambda_OBS.beta 	+= Ts*(vsab_km1.beta  - RS*isab.beta  + KOBS*feedback_OBS.beta);
	lambda_OBS.amp     = sqrt(pow(lambda_OBS.alpha,2) + pow(lambda_OBS.beta,2));
	
	_rot(lambda_OBS, SinCos_elt, lambda_dq);
	T_elt = 1.5f*PP*(lambda_OBS.alpha*isab.beta - lambda_OBS.beta*isab.alpha);	
	delta = atan(lambda_dq.q/lambda_dq.d);
	
}

//----------------------------------------------------------------------------------------------------//

void Compute_Inductance(void) {
	
	float ldq_d, ldq_q;
	
	    //---------Apparent Inductances-----------//
    if (fabs(isdq.d) > 0.05) Ld_unfilt = fabs(lambda_CM_dq.d/isdq.d);
    else {
        interp2d(&FD_LUT[0][0], 0.05, fabs(isdq.q), DIDD, INV_DIDD, DIQD, INV_DIQD , ID_TAB_MAX, ID_TAB_MIN, IQ_TAB_MAX , IQ_TAB_MIN, n_size, &tmp1);
        Ld_unfilt = fabs(tmp1/0.05);
    }
	_Filter(Ld_unfilt, Ld, Ts*TWOPI*100);
	// Ld = Ld_unfilt;
    
		if (fabs(isdq.q) > 0.05) Lq_unfilt = fabs(lambda_CM_dq.q/isdq.q);
    else {
		interp2d(&FQ_LUT[0][0], 0.05, fabs(isdq.d), DIQQ, INV_DIQQ, DIQD, INV_DIQD , IQ_TAB_MAX, IQ_TAB_MIN, ID_TAB_MAX , ID_TAB_MIN, n_size, &tmp1);
        Lq_unfilt = fabs(tmp1/0.05);
    }
	_Filter(Lq_unfilt, Lq, Ts*TWOPI*100);	
	// Lq = Lq_unfilt;
    	
    //--------Incremental Inductances---------//
    
	switch (Quad_Maps){
        case 0:
            isdq.d = fabs(isdq.d);
            isdq.q = fabs(isdq.q);
        case 1:
			isdq.d = fabs(isdq.d);
        case 2:
            isdq.q = fabs(isdq.q);
    }
	
		interp2d(&FD_LUT[0][0], isdq.d+0.1, isdq.q, DIDD, INV_DIDD, DIQD, INV_DIQD , ID_TAB_MAX, ID_TAB_MIN, IQ_TAB_MAX , IQ_TAB_MIN, n_size, &tmp1);
    if(Quad_Maps == 2)
        tmp1 = fabs(tmp1);
    ld_unfilt = (tmp1 - fabs(lambda_CM_dq.d))/0.1;
    
    interp2d(&FQ_LUT[0][0], isdq.q+0.1, isdq.d, DIQQ, INV_DIQQ, DIQD, INV_DIQD , IQ_TAB_MAX, IQ_TAB_MIN, ID_TAB_MAX , ID_TAB_MIN, n_size, &tmp1);
    if(Quad_Maps == 1)
        tmp1 = fabs(tmp1);
    lq_unfilt = (tmp1 - fabs(lambda_CM_dq.q))/0.1;

    interp2d(&FD_LUT[0][0], isdq.d, isdq.q+0.1, DIDD, INV_DIDD, DIQD, INV_DIQD , ID_TAB_MAX, ID_TAB_MIN, IQ_TAB_MAX , IQ_TAB_MIN, n_size, &tmp1);
    if(Quad_Maps == 2)
        tmp1 = fabs(tmp1);
    ldq_d = (tmp1 - fabs(lambda_CM_dq.d))/0.1;
    
    interp2d(&FQ_LUT[0][0], isdq.q, isdq.d+0.1, DIQQ, INV_DIQQ, DIQD, INV_DIQD , IQ_TAB_MAX, IQ_TAB_MIN, ID_TAB_MAX , ID_TAB_MIN, n_size, &tmp1);
    if(Quad_Maps == 1)
        tmp1 = fabs(tmp1);
    ldq_q = (tmp1 - fabs(lambda_CM_dq.q))/0.1;
    ldq_unfilt = (ldq_d + ldq_q)/2*sgn(isdq.d)*sgn(isdq.q);
	
	_Filter(ld_unfilt,  ld,  Ts*TWOPI*100);		
	_Filter(lq_unfilt,  lq,  Ts*TWOPI*100);		
	_Filter(ldq_unfilt, ldq, Ts*TWOPI*100);		
	
	ldm = 0.5*(ld - lq);
	Ldm = 0.5*(Ld - Lq);
	
}