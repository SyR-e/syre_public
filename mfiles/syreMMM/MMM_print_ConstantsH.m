% Copyright 2020
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function MMM_print_ConstantsH(motorModel)

ConstantsH_path = checkPathSyntax([motorModel.data.pathname motorModel.data.motorName '_ctrl\User_functions\Inc\User_Constants.h']);

fs       = motorModel.SyreDrive.Converter.fPWM;
CRT_PROT = motorModel.data.Imax;

MTPA = motorModel.controlTrajectories.MTPA;
MTPA.T(isnan(MTPA.T)) = 0;
i0       = motorModel.data.i0;
T_rated(~isnan(motorModel.data.T0)) = motorModel.data.T0;
T_rated(isempty(T_rated)) = interp1(abs(MTPA.id+1i*MTPA.iq),MTPA.T,i0);

fd_MTPA  = interp1(MTPA.T,MTPA.fd,T_rated);
fq_MTPA  = interp1(MTPA.T,MTPA.fq,T_rated);
f_MTPA   = abs(fd_MTPA+1j*fq_MTPA);
w_base   = motorModel.data.Vdc/(sqrt(3)*f_MTPA);
KOBS     = w_base/3;

fid = fopen(ConstantsH_path,'w');
fprintf(fid,'#ifndef USER_CONSTANTS_H\n');
fprintf(fid,'#define USER_CONSTANTS_H\n\n');
fprintf(fid,'// numbers \n');
fprintf(fid,'#define SQRT3OVER2		0.8660254f\n');
fprintf(fid,'#define SQRT1OVER3 		0.57735027f\n');
fprintf(fid,'#define PI				3.1415926f\n');
fprintf(fid,'#define TWOPI			6.2831853f\n');
fprintf(fid,'#define ONE_HALF   		0.5f\n');
fprintf(fid,'#define TWO_THIRDS      0.666666666666f\n');
fprintf(fid,'#define RPM2RAD         0.1046667f	// rpm to mechanical rad/s\n');
fprintf(fid,'#define RAD2RPM         9.5493f\n\n');

fprintf(fid,'// States \n');
fprintf(fid,'#define START      3\n');
fprintf(fid,'#define READY      2\n');
fprintf(fid,'#define WAKE_UP    1\n');
fprintf(fid,'#define ERROR      0\n\n');

fprintf(fid,'// Ctrl_type \n');
fprintf(fid,'#define CurrentControl    0\n');
fprintf(fid,'#define FluxControl       1\n');
fprintf(fid,'#define TorqueControl     2\n');
fprintf(fid,'#define SpeedControl      3\n\n');

fprintf(fid,'// Ctrl strategy \n');
fprintf(fid,'#define FOC    0\n');
fprintf(fid,'#define DFVC   1\n\n');

fprintf(fid,'// Inverter settings \n');
fprintf(fid,'#define CRT_PROT   %4.2f\n',CRT_PROT);
fprintf(fid,'#define fs 	    %4.2f\n',fs);
fprintf(fid,'#define Ts     	(1.0f/(float)fs)\n');

fprintf(fid,' \n');
fprintf(fid,'// Flux Observer \n');
fprintf(fid,'#define KOBS   %4.2f\n\n',KOBS);

fprintf(fid,'// PI reg \n');
fprintf(fid,'#define OMEGA_BW    	TWOPI*5.0f\n');
fprintf(fid,'#define OMEGA_BI    	TWOPI*700.0f\n');
fprintf(fid,'#define OMEGA_0_INJ   	TWOPI*50.0f\n');
fprintf(fid,'#define WB_PLL_PULS   	TWOPI*10.0f\n');
fprintf(fid,'#define OMEGA_PLL    	TWOPI*25.0f\n');
fprintf(fid,'#define OMEGA_DELTA   	TWOPI*5.0f\n');
fprintf(fid,'#define WB_PLL_SQUARE 	TWOPI*25.0f\n\n');

fprintf(fid,'// SENSORLESS\n');
fprintf(fid,'// Control fusion \n');
fprintf(fid,'#define OMEGA_G   %4.2f\n',w_base/4);
fprintf(fid,'// Voltage injection \n');
fprintf(fid,'#define SINUS    	0\n');
fprintf(fid,'#define SQUARE    	1\n');
fprintf(fid,'// Demodulation \n');
fprintf(fid,'#define CRNTDEM   	0\n');
fprintf(fid,'#define FLUXDEM   	1\n');
fprintf(fid,'// High speed control \n');
fprintf(fid,'#define AF   	0\n');
fprintf(fid,'#define APP   	1\n');

fprintf(fid,' \n');
fprintf(fid,'#endif \n');

fclose(fid);

