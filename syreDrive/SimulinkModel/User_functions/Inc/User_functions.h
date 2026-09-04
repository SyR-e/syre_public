
extern char data_output_vector[]; //vettore uscite

void command_int (void) {
extern char data_input_vector[]; //vettore uscite
variabile_int=(int*)((data_input_vector[0]<<24)+(data_input_vector[1]<<16)+(data_input_vector[2]<<8)+(data_input_vector[3]<<0));
valore_int=((data_input_vector[4]<<24)+(data_input_vector[5]<<16)+(data_input_vector[6]<<8)+(data_input_vector[7]<<0));
	valore=*(float*)&valore_int;
}

void char_to_float (void) {
extern char data_input_vector[]; //vettore uscite
add1=(float*)((data_input_vector[0]<<24)+(data_input_vector[1]<<16)+(data_input_vector[2]<<8)+(data_input_vector[3]<<0));
add2=(float*)((data_input_vector[4]<<24)+(data_input_vector[5]<<16)+(data_input_vector[6]<<8)+(data_input_vector[7]<<0));
add3=(float*)((data_input_vector[8]<<24)+(data_input_vector[9]<<16)+(data_input_vector[10]<<8)+(data_input_vector[11]<<0));
add4=(float*)((data_input_vector[12]<<24)+(data_input_vector[13]<<16)+(data_input_vector[14]<<8)+(data_input_vector[15]<<0));
	
}

void float_two_char( float num_float, char *x1, char *x2, char *x3 ,char *x4)
{

  int num_int;
	
	num_int=*(int*)&num_float;
	
	*x1=((num_int>>24)&	0xff); // ultimi 8 bit MSB
	*x2=((num_int>>16)&	0xff);
	*x3=((num_int>>8)	&	0xff); 
  *x4=((num_int)		&	0xff); //primi otto LSB
	
}
void write_data_output_vector(float num_float){
	
	extern int posizione_vettore_usart;//posizione nel vettore
	//la funzione prende il float a 32bit e lo spezza in 3 char a 8 bit e impila tutto nel vettore
	float_two_char(  num_float, &data_output_vector[posizione_vettore_usart+3], &data_output_vector[posizione_vettore_usart+2], &data_output_vector[posizione_vettore_usart+1] ,&data_output_vector[posizione_vettore_usart]);
	
	posizione_vettore_usart=posizione_vettore_usart+4;//4=32bit/8bit;
}

void command_float(void) {
extern char data_input_vector[]; //vettore uscite
variabile=(float*)((data_input_vector[0]<<24)+(data_input_vector[1]<<16)+(data_input_vector[2]<<8)+(data_input_vector[3]<<0));
valore_int=((data_input_vector[4]<<24)+(data_input_vector[5]<<16)+(data_input_vector[6]<<8)+(data_input_vector[7]<<0));
	valore=*(float*)&valore_int;
}

