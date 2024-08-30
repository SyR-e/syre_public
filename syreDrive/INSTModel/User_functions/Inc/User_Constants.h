#ifndef USER_CONSTANTS_H
#define USER_CONSTANTS_H

// numbers
#define SQRT3OVER2		0.8660254f
#define SQRT1OVER3 		0.57735027f
#define PI				3.1415926f
#define TWOPI			6.2831853f

#define ONE_HALF   		0.5f
#define TWO_THIRDS      0.666666666666f
#define RPM2RAD         0.1046667f	// rpm to mechanical rad/s
#define RAD2RPM         9.5493f

//States
#define START      			3
#define READY      			2
#define WAKE_UP    			1
#define ERROR 	   			0

//Ctrl_type
#define CurrentControl 	    0
#define FluxControl 	    1
#define TorqueControl 	    2
#define SpeedControl 	    3

//Ctrl strategy
#define FOC                 0
#define DFVC                1


#define CRT_PROT 		30.0f
				
#define fs 				10.0e3f
#define Ts 				(1.0f/(float)fs)

#define scala_current   0.01567f
#define scala_voltage   0.0160f
#define COEFF_DAC_SPI   2047
#define ZERO_OFFSET     2048

#define CRT_MAX 		13.0f

#define BLANKTIME		1750

//Encoder 
#define ENC_OFF      0
#define ENC_FAC		 0.00076699039f

// Flux Observer
#define KOBS 		    TWOPI*10.0f  // 300 rpm

// PI Reg
#define OMEGA_BW    	TWOPI*5.0f
#define OMEGA_BI    	TWOPI*700.0f
#define OMEGA_0_INJ 	TWOPI*50.0f 
#define WB_PLL_PULS 	TWOPI*10.00f  
#define OMEGA_PLL       TWOPI*25.00f
#define OMEGA_DELTA     TWOPI*5.00f
#define WB_PLL_SQUARE   TWOPI*25.0f

// Control fusion
#define OMEGA_G				TWOPI*4.0f   // 120 rpm

// Voltage injection
#define SINUS 	0
#define SQUARE 	1
// Demodulation
#define CRNTDEM 0
#define FLUXDEM 1

// High speed control
#define AF      0
#define APP 		1
#endif
