/*----------------------------------------------------------------------*/
/* File: Syr_ctrl.c                                               	*/
/* the parameter CTRL_TYPE decides the type of control		*/
/*  CTRL_TYPE = 0 - Current control												*/
/*  CTRL_TYPE = 1 - Flux control											  	*/
/*  CTRL_TYPE = 2 - Torque control												*/
/*  CTRL_TYPE = 3 - Speed control													*/
/*----------------------------------------------------------------------*/
/* Motor type: 	BariSyr                                   */
/*----------------------------------------------------------------------*/

#define S_FUNCTION_NAME Syr_ctrl
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

#define NINPUTS 	12
#define NOUTPUTS 	18
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
    
#ifdef SIM
    
    isabc.a   			= U(0);
    isabc.b   			= U(1);
    isabc.c   			= U(2);
    vdc  	  	  		= U(3);
    
    theta_mec_meas  = U(4);			// encoder
    n_ref_in   			= U(5);			// rpm
    T_ext				    = U(6);			// reference torque
    isdq_ext.d 			= U(7);  		// Current reference
    isdq_ext.q 			= U(8);  		// Current reference
	
    Reset 				  = U(9);			// Black button
    Go    				  = U(10);	  // Blue  button
    Ctrl_type   		= U(11);    
#else
    
    /* USER CODE BEGIN TIM1_UP_TIM16_IRQn 0 */
    GPIOC->ODR^=(1<<3);
    
    IWDG->KR = 0xAAAA;		// reload watchdog TIMER (IDWG RLR)
    while(((ADC1->ISR)&(1<<6))==0){} 	//wait end of ADC conversions
    ADC1->ISR|=((1<<6));	//clear flag JEOC (the flag JEOC is cleared by setting the correspondent bit to 1)
    
    // feedback
    input.ch0=(float)ADC1->JDR1;
    input.ch1=(float)ADC1->JDR2;
    input.ch2=(float)ADC1->JDR3;
    
    isabc.a=-(input.ch0-offset_current_a)*scala_current; // current phase a
    isabc.b=-(input.ch1-offset_current_b)*scala_current; // current phase b
    isabc.c=-(input.ch2-offset_current_c)*scala_current; // current phase c
    vdc = ((int)ADC1->JDR4)*scala_voltage;               // vdc-link voltage
    
    // Read Button State
    if(((GPIOC->IDR)& (1<<13))==0 && Go_flag==1){
        Go = 1.0f;	// Go Button
        Go_flag=0;
    }
    else	Go = 0.0f;
    if ((GPIOC->IDR)& (1<<13)) Go_flag=1;
    
#endif
    
    //Over current protection
    CurrentProtection(isabc,&State,&pwm_stop);
    
    //Clarke transformation (a,b,c)--> (alpha,beta)
    _clarke(isabc, isab);
    
	switch(State){
        
        case ERROR:
            // Variables Init
            pwm_stop 									= 1;
            counter 									= 0;
            
            offset_current_a					=2120;
            offset_current_b					=2120;
            offset_current_c					=2120;
            
            n_ref_in 									= 0.0f;
            omega_ref_in							= 0.0f;
            omega_ref_ramp						= 0.0f;
            accel               			= 1000;	// rpm/s
						omega_elt_meas           	= 0.0f;
						omega_elt_meas_f         	= 0.0f;
            SinCos_elt_meas.sin      	= 0.0f;
            SinCos_elt_meas.cos      	= 1.0f;
            SinCos_elt_meas_old.sin  	= 0.0f;
            SinCos_elt_meas_old.cos  	= 1.0f;

						Ld 												= Ld_inic;
						Lq												= Lq_inic;
						ld 												= 0;
						lq 												= 0;
						
            flux_nom 									= 0.0f;	// rated flux (Vs)
            v0 												= 0.0f;			// phase dc voltage (V)
            isdq_ref.d 								= 0.0f;
            
            duty_abc.a								= 0.0f;
            duty_abc.b								= 0.0f;
            duty_abc.c								= 0.0f;
            
            offset_in.ch0							= 0.0f;
            offset_in.ch1							= 0.0f;
            offset_in.ch2							= 0.0f;
            offset_in.ch3							= 0.0f;
            
            
            if(Go) State = WAKE_UP;	// Wait for Go state
            
            break;
            
        case WAKE_UP:
            
            duty_abc.a=0.5f;
            duty_abc.b=0.5f;
            duty_abc.c=0.5f;
            
            pwm_stop = 0;
            
            // Current offset accumulation
            if (counter>100){
                offset_in.ch0+=input.ch0;
                offset_in.ch1+=input.ch1;
                offset_in.ch2+=input.ch2;
            }
            
            // Offset computation and bootstrap
            if (counter==(100+200)) {
                offset_current_a=(float)(offset_in.ch0/200.0f);
                offset_current_b=(float)(offset_in.ch1/200.0f);
                offset_current_c=(float)(offset_in.ch2/200.0f);
            }
            
            // Switch to READY state
            if (counter > 300) {
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
            if (Go)	State = START;
            
            break;
            
        case START:
            
            omega_ref_in        = n_ref_in * RPM2RAD * PP;
            ramp(omega_ref_in, accel * RPM2RAD * PP * Ts, &omega_ref_ramp);
            
						// - CALCULATION OF SPEED
						theta_elt_meas      = PP * theta_mec_meas;
						SinCos_elt_meas.sin = sin(theta_elt_meas);
						SinCos_elt_meas.cos = cos(theta_elt_meas);
						speed_compute_sc(SinCos_elt_meas, &SinCos_elt_meas_old, &omega_elt_meas);
						_Filter(omega_elt_meas, omega_elt_meas_f, Ts*TWOPI*100);
						omega_mec_meas      = omega_elt_meas/PP;
						omega_mec_meas_f    = omega_elt_meas_f/PP;
						omega_mec_meas_rpm  = omega_mec_meas_f*30/PI;
							
						SinCos_elt_meas_old.sin = SinCos_elt_meas.sin;
						SinCos_elt_meas_old.cos = SinCos_elt_meas.cos;
				
            // Ignore reference speed during start-up
//            if (counter<BLANKTIME) {
//                omega_ref_ramp = 0.0f;
//                isdq_ref.q = 0.0f;
//                counter++;
//            }               										
							
			switch (Ctrl_type){
			
				case 0: //CurrentControl
					isdq_ref.d = isdq_ext.d;
					isdq_ref.q = isdq_ext.q;
				break;
															
				case 2: //TorqueControl
					ReadLut(&ID_REF[0], fabs(T_ext), TMAX, TMIN, DT, INV_DT, &isdq_ref.d); 
					ReadLut(&IQ_REF[0], fabs(T_ext), TMAX, TMIN, DT, INV_DT, &isdq_ref.q);
					isdq_ref.q = sgn(T_ext)*isdq_ref.q;																	
				break;

				case 3: //SpeedControl
					sp_var.ref  = omega_ref_in/PP;
					sp_var.fbk  = omega_mec_meas_f;
					sp_par.lim  = Tmax_mot;
					kp_w        = 2*OMEGA_BW*J;
					ki_w        = pow(OMEGA_BW,2)*J;
					sp_par.ki   = ki_w*Ts;
					sp_par.kp   = kp_w;
					PIReg(&sp_par, &sp_var);
					T_ext       = sp_var.out;
																	
					ReadLut(&ID_REF[0], fabs(T_ext), TMAX, TMIN, DT, INV_DT, &isdq_ref.d); 
					ReadLut(&IQ_REF[0], fabs(T_ext), TMAX, TMIN, DT, INV_DT, &isdq_ref.q);
					isdq_ref.q = sgn(T_ext)*isdq_ref.q;								
				break;
			}
														
			_clarke(isabc, isab);
			_rot(isab, SinCos_elt_meas, isdq);
			
			vsab_km1 = vsab_k0;
			DTComp(duty_abc, duty_abc_km1, isabc, vdc, deadtime, &vsabc_k0);
			_clarke(vsabc_k0, vsab_k0)
			
			FluxObserver();
			Compute_Inductance();
				
			kp_id = OMEGA_BI*ld;
			ki_id = kp_id*OMEGA_BI/10;
			kp_iq = OMEGA_BI*lq;
			ki_iq = kp_iq*OMEGA_BI/10;

			id_par.kp   = kp_id;
			id_par.ki   = ki_id*Ts;
			iq_par.kp   = kp_iq;
			iq_par.ki   = ki_iq*Ts;

			//Current loop
			Current_loop(vdc, Imax_mot, isdq_ref, isdq, &id_par, &id_var, &iq_par, &iq_var, &vsdq_ref);
			vsdq_ref.d += RS*isdq.d - omega_elt_meas_f*lambda_dq.q;
			vsdq_ref.q += RS*isdq.q + omega_elt_meas_f*lambda_dq.d;

			_invrot(vsdq_ref, SinCos_elt_meas, vsab_ref);
			_invclarke(vsab_ref, vsabc_ref);
			PWMduty(vsabc_ref, vdc, &duty_abc);								
		}
    
    //duty cycles saturation
    if (duty_abc.a > 0.95f) duty_abc.a=0.95f;
    if (duty_abc.b > 0.95f) duty_abc.b=0.95f;
    if (duty_abc.c > 0.95f) duty_abc.c=0.95f;
    
    if (duty_abc.a < 0.05f) duty_abc.a=0.05f;
    if (duty_abc.b < 0.05f) duty_abc.b=0.05f;
    if (duty_abc.c < 0.05f) duty_abc.c=0.05f;
    
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
    y[4]  = omega_ref_in/PP;
    y[5]  = omega_mec_meas_f;
    y[6]  = vsab_ref.alpha;
    y[7]  = vsab_ref.beta;
    y[8]  = isab.alpha;
    y[9]  = isab.beta;
    y[10] = isdq_ref.d;
    y[11] = isdq.d;
    y[12] = isdq_ref.q;
    y[13] = isdq.q;
		y[14] = State;
		y[15] = T_elt;
		y[16] = lambda_dq.d;
		y[17] = lambda_dq.q;
    
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

void FluxObserver(void) {
	
	lambda_CM_ab_km1 = lambda_CM_ab;
	
  interp2d(&FD_LUT[0][0], fabs(isdq.d), fabs(isdq.q), DIDD, INV_DIDD, DIQD, INV_DIQD , ID_TAB_MAX, ID_TAB_MIN, IQ_TAB_MAX , IQ_TAB_MIN, n_size, &lambda_CM_dq.d);
  interp2d(&FQ_LUT[0][0], fabs(isdq.q), fabs(isdq.d), DIQQ, INV_DIQQ, DIQD, INV_DIQD , IQ_TAB_MAX, IQ_TAB_MIN, ID_TAB_MAX , ID_TAB_MIN, n_size, &lambda_CM_dq.q);
	if (isdq.d < 0)
		lambda_CM_dq.d = -lambda_CM_dq.d;
	if (isdq.q < 0)
		lambda_CM_dq.q = -lambda_CM_dq.q;

	_invrot(lambda_CM_dq, SinCos_elt_meas, lambda_CM_ab);
	
	feedback_OBS.alpha	=	lambda_CM_ab_km1.alpha - lambda_OBS.alpha;
  feedback_OBS.beta		=	lambda_CM_ab_km1.beta - lambda_OBS.beta;
    
	// Integration
  lambda_OBS.alpha	+= Ts*(vsab_km1.alpha - RS*isab.alpha + KOBS*feedback_OBS.alpha);
  lambda_OBS.beta 	+= Ts*(vsab_km1.beta  - RS*isab.beta  + KOBS*feedback_OBS.beta);
	
  _rot(lambda_OBS, SinCos_elt_meas, lambda_dq);
	
	T_elt = 1.5f * PP*(lambda_OBS.alpha*isab.beta - lambda_OBS.beta*isab.alpha);	
	
	// lambda = sqrt(pow(lambda_OBS.alpha,2) + pow(lambda_OBS.beta,2));
	delta = atan(lambda_dq.q/lambda_dq.d);
		
}

void Compute_Inductance(void) {
	
//    float ldq_d, ldq_q;
	
	    //---------Apparent Inductances-----------//
    if (fabs(isdq.d) > 0.05) Ld_unfilt = fabs(lambda_CM_dq.d/isdq.d);
    else {
        interp2d(&FD_LUT[0][0], 0.05, fabs(isdq.q), DIDD, INV_DIDD, DIQD, INV_DIQD , ID_TAB_MAX, ID_TAB_MIN, IQ_TAB_MAX , IQ_TAB_MIN, n_size, &tmp1);
        Ld_unfilt = fabs(tmp1/0.05);
    }
//	_Filter(Ld_unfilt, Ld, 150);
	Ld = Ld_unfilt;
    
	if (fabs(isdq.q) > 0.05) Lq_unfilt = fabs(lambda_CM_dq.q/isdq.q);
    else {
		interp2d(&FQ_LUT[0][0], 0.05, fabs(isdq.d), DIQQ, INV_DIQQ, DIQD, INV_DIQD , IQ_TAB_MAX, IQ_TAB_MIN, ID_TAB_MAX , ID_TAB_MIN, n_size, &tmp1);
        Lq_unfilt = fabs(tmp1/0.05);
    }
//	_Filter(Lq_unfilt, Lq, 150);	
	Lq = Lq_unfilt;
    	
    //--------Incremental Inductances---------//
    
    interp2d(&FD_LUT[0][0], fabs(isdq.d)+0.1, fabs(isdq.q), DIDD, INV_DIDD, DIQD, INV_DIQD , ID_TAB_MAX, ID_TAB_MIN, IQ_TAB_MAX , IQ_TAB_MIN, n_size, &tmp1);
    ld_unfilt = (tmp1 - fabs(lambda_CM_dq.d))/0.1;
    interp2d(&FQ_LUT[0][0], fabs(isdq.q)+0.1, fabs(isdq.d), DIQQ, INV_DIQQ, DIQD, INV_DIQD , IQ_TAB_MAX, IQ_TAB_MIN, ID_TAB_MAX , ID_TAB_MIN, n_size, &tmp1);
    lq_unfilt = (tmp1 - fabs(lambda_CM_dq.q))/0.1;
    
    interp2d(&FD_LUT[0][0], fabs(isdq.d), fabs(isdq.q)+0.1, DIDD, INV_DIDD, DIQD, INV_DIQD , ID_TAB_MAX, ID_TAB_MIN, IQ_TAB_MAX , IQ_TAB_MIN, n_size, &tmp1);
    ldq_d = (tmp1 - fabs(lambda_CM_dq.d))/0.1;
    interp2d(&FQ_LUT[0][0], fabs(isdq.q), fabs(isdq.d)+0.1, DIQQ, INV_DIQQ, DIQD, INV_DIQD , IQ_TAB_MAX, IQ_TAB_MIN, ID_TAB_MAX , ID_TAB_MIN, n_size, &tmp1);
    ldq_q = (tmp1 - fabs(lambda_CM_dq.q))/0.1;
    ldq_unfilt = (ldq_d + ldq_q)/2*sgn(isdq.d)*sgn(isdq.q);
	
//	_Filter(ld_unfilt, ld, 150);		
//	_Filter(lq_unfilt, lq, 150);		
//	_Filter(ldq_unfilt, ldq, 150);		
		
		ld  =  ld_unfilt;
		lq  =  lq_unfilt;
		ldq = ldq_unfilt;
	}