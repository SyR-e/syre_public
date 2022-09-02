/*----------------------------------------------------------------------*/
/* Filename:    generic inverter                                        */
/* Date:        05.02.2010                                              */
/* Description: Average model of a generic three-phase inverter         */
/*              supplying a three-phase, wye-connected simmetrical load */
/* Status:      OK                                                      */
/*----------------------------------------------------------------------*/

#define S_FUNCTION_NAME inverter
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"
#include "math.h"

#define U(element) (*uPtrs[element])  /* Pointer to Input Port0 */

/* Define input and output widths  */
#define NINPUTS 7
#define NOUTPUTS 6

/* Mask parameters   */
#define fsw				(mxGetPr(ssGetSFcnParam(S,0))[0])	/* switching frequency (Hz)			 */
#define dead_time		(mxGetPr(ssGetSFcnParam(S,1))[0])	/* dead time duration (us)			 */
#define enable_dt		(mxGetPr(ssGetSFcnParam(S,2))[0])	/* enable dead-time					 */
#define enable_vdrop	(mxGetPr(ssGetSFcnParam(S,3))[0])	/* enable on-state drops			 */
#define vON1_TR			(mxGetPr(ssGetSFcnParam(S,4))[0])	/* vON IGBT, low current range (V)	 */
#define Rd1_TR			(mxGetPr(ssGetSFcnParam(S,5))[0])	/* Rd  IGBT, low current range (Ohm) */
#define Iknee_TR		(mxGetPr(ssGetSFcnParam(S,6))[0])	/* IGBT low current upper limit (A)	 */
#define vON2_TR			(mxGetPr(ssGetSFcnParam(S,7))[0])	/* vON IGBT, high current range (V)	 */
#define Rd2_TR        	(mxGetPr(ssGetSFcnParam(S,8))[0])	/* Rd  IGBT, high current range (Ohm)*/
#define vON1_D        	(mxGetPr(ssGetSFcnParam(S,9))[0])	/* vON diode, low current range (V)	 */
#define Rd1_D         	(mxGetPr(ssGetSFcnParam(S,10))[0])	/* Rd  diode, low current range (Ohm)*/
#define Iknee_D       	(mxGetPr(ssGetSFcnParam(S,11))[0])	/* diode low current upper limit (A) */
#define vON2_D          (mxGetPr(ssGetSFcnParam(S,12))[0])	/* vON diode, high current range (V) */
#define Rd2_D           (mxGetPr(ssGetSFcnParam(S,13))[0])	/* Rd  diode, high current range (Ohm)*/
#define TR_mdl          (mxGetPr(ssGetSFcnParam(S,14))[0])	/* IGBT model (2 or 5 parameters)     */
#define D_mdl           (mxGetPr(ssGetSFcnParam(S,15))[0])	/* diode model (2 or 5 parameters)     */
#define NPARAMS 16

/* Threshold values to define the signum function   */
#define  crt_max 1.e-3
#define  crt_min -1.e-3

//Component parameters structure
typedef struct { double vON1;
				 double Rd1;
				 double Iknee;
				 double vON2;
				 double Rd2;
                 double curr;
                 double mdl;
               } XCompPar;
               
XCompPar TR_par, D_par;

/* Global functions   */
double sgn(double i)
{
 double sgn;
    if (i>crt_max)
        sgn=1.0;
    else
        {
            if (i<crt_min)
                sgn=-1.0;
            else
                sgn=0.0;
        }
    return (sgn);
}

/* Voltage drop computation for a generic component */
double deltaV(XCompPar *CompPar)
{
	double voltage_drop;
    if (CompPar->mdl < 1.5)
        voltage_drop=CompPar->vON1+CompPar->Rd1*fabs(CompPar->curr);
    else
    {
	if(fabs(CompPar->curr)<=CompPar->Iknee)
		voltage_drop=CompPar->vON1+CompPar->Rd1*fabs(CompPar->curr);
	else
		voltage_drop=CompPar->vON2+CompPar->Rd2*(fabs(CompPar->curr)-CompPar->Iknee);
    }
	return(voltage_drop);
}

/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ============================================== */
static void mdlInitializeSizes(SimStruct *S)
{
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



/* Function: mdlInitializeSampleTimes =========================================*/
static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
}

#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================*/
static void mdlInitializeConditions(SimStruct *S)
{
}

/* Function: mdlOutputs =======================================================*/
static void mdlOutputs(SimStruct *S, int_T tid)
{
    real_T            *y    = ssGetOutputPortRealSignal(S,0);
    real_T            *x    = ssGetRealDiscStates(S);
    InputRealPtrsType uPtrs = ssGetInputPortRealSignalPtrs(S,0);
  
/* Local variables  */
    real_T dutyA,dutyB,dutyC;
    real_T dutyA_km1,dutyB_km1,dutyC_km1;
    real_T vdc;
    real_T iA,iB,iC;
    real_T vA,vB,vC;
    real_T v_error_dt_A,v_error_dt_B,v_error_dt_C;
    real_T v_error_vdrop_A,v_error_vdrop_B,v_error_vdrop_C;
    real_T v_errorA,v_errorB,v_errorC;
    
/* Assign inputs to local variables   */
    dutyA=U(0);
    dutyB=U(1);
    dutyC=U(2);
    vdc  =U(3);
    iA   =U(4);
    iB   =U(5);
    iC   =U(6);
    
	if(dutyA>1) dutyA = 1.;
	if(dutyA<0) dutyA = 0.;
	if(dutyB>1) dutyB = 1.;
	if(dutyB<0) dutyB = 0.;
	if(dutyC>1) dutyC = 1.;
	if(dutyC<0) dutyC = 0.;
	
	
/* Assign components parameters */
    TR_par.vON1 = vON1_TR;
    TR_par.Rd1 = Rd1_TR;
    TR_par.Iknee = Iknee_TR;
    TR_par.vON2 = vON2_TR;
    TR_par.Rd2 = Rd2_TR;
    TR_par.mdl = TR_mdl;
    D_par.vON1 = vON1_D;
    D_par.Rd1 = Rd1_D;
    D_par.Iknee = Iknee_D;
    D_par.vON2 = vON2_D;
    D_par.Rd2 = Rd2_D;
    D_par.mdl = D_mdl;
    
/* Inverter phase voltage errors due to the dead-time effects    */
    v_error_dt_A	= vdc*dead_time*fsw*sgn(iA);
    v_error_dt_B	= vdc*dead_time*fsw*sgn(iB);
    v_error_dt_C	= vdc*dead_time*fsw*sgn(iC);
    
    // if ((dutyA  == 0 && dutyA_km1 == 0) || (dutyA  == 1 && dutyA_km1 == 1))
        // v_error_dt_A = 0;
    // if ((dutyB  == 0 && dutyB_km1 == 0) || (dutyB  == 1 && dutyB_km1 == 1))
        // v_error_dt_B = 0;
    // if ((dutyC  == 0 && dutyC_km1 == 0) || (dutyC  == 1 && dutyC_km1 == 1))
        // v_error_dt_C = 0;
    
    // if ((dutyA == 1 && dutyA_km1 == 0 && sgn(iA)<0) || (dutyA == 0 && dutyA_km1 == 1 && sgn(iA)>0))
        // v_error_dt_A = 0;
    // if ((dutyB == 1 && dutyB_km1 == 0 && sgn(iB)<0) || (dutyB == 0 && dutyB_km1 == 1 && sgn(iB)>0))
        // v_error_dt_B = 0;
    // if ((dutyC == 1 && dutyC_km1 == 0 && sgn(iC)<0) || (dutyC == 0 && dutyC_km1 == 1 && sgn(iC)>0))
        // v_error_dt_C = 0;
    
    dutyA_km1 = dutyA; dutyB_km1 = dutyB; dutyC_km1 = dutyC;
    
/* Voltage errors due to the on-state voltage drops	*/
    TR_par.curr = iA;
    D_par.curr = iA;
    if (iA>crt_max)
        {
        v_error_vdrop_A= deltaV(&D_par)+dutyA*(deltaV(&TR_par)-deltaV(&D_par));
        }
         else
            {
    		if (iA<crt_min)
                v_error_vdrop_A=-deltaV(&TR_par)+dutyA*(deltaV(&TR_par)-deltaV(&D_par));
        	else 
        	    v_error_vdrop_A=0.0;   
            } 
       
    TR_par.curr = iB;
    D_par.curr = iB; 
    if (iB>crt_max)
        {
        v_error_vdrop_B= deltaV(&D_par)+dutyB*(deltaV(&TR_par)-deltaV(&D_par));
        }
         else
            {
    		if (iB<crt_min)
                v_error_vdrop_B=-deltaV(&TR_par)+dutyB*(deltaV(&TR_par)-deltaV(&D_par));
        	else 
        	    v_error_vdrop_B=0.0;   
            }
    
    TR_par.curr = iC;
    D_par.curr = iC;
    if (iC>crt_max)
        {
        v_error_vdrop_C= deltaV(&D_par)+dutyC*(deltaV(&TR_par)-deltaV(&D_par));
        }
         else
            {
    		if (iC<crt_min)
                v_error_vdrop_C=-deltaV(&TR_par)+dutyC*(deltaV(&TR_par)-deltaV(&D_par));
        	else 
        	    v_error_vdrop_C=0.0;   
            } 
        
   
    /*Errors in the inverters' output voltages respect to its neutral point*/ 
   
   v_errorA=v_error_dt_A*enable_dt+v_error_vdrop_A*enable_vdrop;
   v_errorB=v_error_dt_B*enable_dt+v_error_vdrop_B*enable_vdrop;
   v_errorC=v_error_dt_C*enable_dt+v_error_vdrop_C*enable_vdrop;

   
   /* Inverter output voltages respect to the inverter's neutral point	*/
   vA=vdc*(dutyA-0.5)-v_errorA;
   vB=vdc*(dutyB-0.5)-v_errorB;
   vC=vdc*(dutyC-0.5)-v_errorC;

   v_errorA=v_error_vdrop_A*enable_vdrop;
   v_errorB=v_error_vdrop_B*enable_vdrop;
   v_errorC=v_error_vdrop_C*enable_vdrop;

   
   /* Load phase voltages for wye connection*/
    y[0]=2.0*vA/3.0-(vB+vC)/3.0;
    y[1]=2.0*vB/3.0-(vA+vC)/3.0;
    y[2]=2.0*vC/3.0-(vA+vB)/3.0;
    
	// y[0]=vdc*dutyA-v_errorA;
    // y[1]=vdc*dutyB-v_errorB;
    // y[2]=vdc*dutyC-v_errorC;

	
   /* Generic outputs */ 
    y[3]=v_errorA;
    y[4]=v_errorB;
    y[5]=v_errorC;
    
}

#define MDL_UPDATE
/* Function: mdlUpdate ====================================================== */
static void mdlUpdate(SimStruct *S, int_T tid)
{
    real_T            *x       = ssGetRealDiscStates(S);
    InputRealPtrsType uPtrs    = ssGetInputPortRealSignalPtrs(S,0);

    UNUSED_ARG(tid); /* not used in single tasking mode */
}

/* Function: mdlTerminate ===================================================== */
static void mdlTerminate(SimStruct *S)
{
    UNUSED_ARG(S); /* unused input argument */
}

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif
