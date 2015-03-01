/*************************************************************************/
/*									 */
/*		Definitions used in C4.5				 */
/*              ------------------------				 */
/*									 */
/*************************************************************************/


#include <stdio.h>
#include <math.h>

#define	 Eof			EOF             /*char read on end of file*/
#define	 Nil			0               /*null pointer*/
#define	 false			0 
#define	 true			1 
#define	 None			-1
#define	 Epsilon                1E-3

long	 random();
#define	 Random			((random()&2147483647) / 2147483648.0)

#define	 Max(a,b)               ((a)>(b) ? a : b) 
#define	 Min(a,b)               ((a)<(b) ? a : b) 
#define	 Round(x)		((int) (x+0.5))
#define	 Log2			0.69314718055994530942
#define	 Log(x)			((x) <= 0 ? 0.0 : log((float)x) / Log2)

#define	 Bit(b)			(1 << (b))
#define	 In(b,s)		((s[(b) >> 3]) & Bit((b) & 07))
#define	 ClearBits(n,s)		memset(s,0,n)
#define	 CopyBits(n,f,t)	memcpy(t,f,n)
#define	 SetBit(b,s)		(s[(b) >> 3] |= Bit((b) & 07))

#define	 ForEach(v,f,l)		for(v=f ; v<=l ; ++v) 

#define	 Verbosity(d)		if(VERBOSITY >= d)

#define	 Check(v,l,h)\
	     if ( v<l||v>h ) {printf("\t** illegal value **\n"); exit(1);}
