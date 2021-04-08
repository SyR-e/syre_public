//#include "User_data_types.h"

// saturate
#define _sat(A, Pos, Neg)	if (A>Pos)	\
							A = Pos;	\
							if (A<Neg)	\
							A = Neg;
					 			
//Direct Clarke transformation (a,b,c)--> (alpha,beta)
#define _clarke(abc,ab)	ab.alpha = TWO_THIRDS*(abc.a-0.5f * abc.b-0.5f * abc.c);		\
							ab.beta = SQRT1OVER3 * (abc.b - abc.c);					

//Inverse Clarke transformation (alpha,beta)-->(a,b,c) 		 
#define _invclarke(ab,abc) abc.a=ab.alpha;			\
		tmp1=-(ab.alpha/2);							\
		tmp2 = ab.beta * SQRT3OVER2;				\
		abc.b = tmp1 + tmp2;						\
		abc.c = tmp1 - tmp2;
		
//Direct rotational transformation (alpha,beta)-->(d,q)
#define _rot(ab,sc,dq) tmp1 = ab.alpha * sc.cos;	\
		tmp2 = ab.beta * sc.sin;					\
		dq.d = tmp1 + tmp2;							\
		tmp1 = ab.alpha * sc.sin;					\
		tmp2 = ab.beta * sc.cos;					\
		dq.q = -tmp1 + tmp2;
		
//Inverse rotational transformation (alpha,beta)-->(d,q)
#define _invrot(dq,sc,ab)	tmp1 = dq.d * sc.cos;		\
		tmp2 = dq.q * sc.sin;							\
		ab.alpha = tmp1 - tmp2;							\
		tmp1 = dq.d * sc.sin;							\
		tmp2 = dq.q * sc.cos;							\
		ab.beta = tmp1 + tmp2;		
	
		
#define _Filter( x, x_filt, k_filt) x_filt= x_filt + k_filt*(x-x_filt) ;



