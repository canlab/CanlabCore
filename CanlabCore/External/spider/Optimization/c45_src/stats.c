/*************************************************************************/
/*									 */
/*  Statistical routines for C4.5					 */
/*  -----------------------------					 */
/*									 */
/*************************************************************************/


#include "defns.i"
#include "types.i"
#include "extern.i"

									
/*************************************************************************/
/*									 */
/*  Compute the additional errors if the error rate increases to the	 */
/*  upper limit of the confidence level.  The coefficient is the	 */
/*  square of the number of standard deviations corresponding to the	 */
/*  selected confidence level.  (Taken from Documenta Geigy Scientific	 */
/*  Tables (Sixth Edition), p185 (with modifications).)			 */
/*									 */
/*************************************************************************/


float Val[] = {  0,  0.001, 0.005, 0.01, 0.05, 0.10, 0.20, 0.40, 1.00},
      Dev[] = {4.0,  3.09,  2.58,  2.33, 1.65, 1.28, 0.84, 0.25, 0.00};


float AddErrs(N, e)
/*    -------  */
    ItemCount N, e;
{
    static float Coeff=0;
    float Val0, Pr;

    if ( ! Coeff )
    {
	/*  Compute and retain the coefficient value, interpolating from
	    the values in Val and Dev  */

	int i;

	i = 0;
	while ( CF > Val[i] ) i++;

	Coeff = Dev[i-1] +
		  (Dev[i] - Dev[i-1]) * (CF - Val[i-1]) /(Val[i] - Val[i-1]);
	Coeff = Coeff * Coeff;
    }

    if ( e < 1E-6 )
    {
	return N * (1 - exp(log(CF) / N));
    }
    else
    if ( e < 0.9999 )
    {
	Val0 = N * (1 - exp(log(CF) / N));
	return Val0 + e * (AddErrs(N, 1.0) - Val0);
    }
    else
    if ( e + 0.5 >= N )
    {
	return 0.67 * (N - e);
    }
    else
    {
	Pr = (e + 0.5 + Coeff/2
	        + sqrt(Coeff * ((e + 0.5) * (1 - (e + 0.5)/N) + Coeff/4)) )
             / (N + Coeff);
	return (N * Pr - e);
    }
}
