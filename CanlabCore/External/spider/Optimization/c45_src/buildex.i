/*************************************************************************/
/*									 */
/*	  Global data for C4.5 used for building decision trees		 */
/*	  -----------------------------------------------------		 */
/*									 */
/*************************************************************************/

#include "defns.i"
#include "types.i"
#include "extern.i"


extern ItemCount
	*Weight,	/* Weight[i]  = current fraction of item i */
	**Freq,		/* Freq[x][c] = no. items of class c with outcome x */
	*ValFreq;	/* ValFreq[x] = no. items with att value v */

extern float
	*Gain,		/* Gain[a] = info gain by split on att a */
	*Info,		/* Info[a] = potential info from split on att a */
	*Bar,		/* Bar[a]  = best threshold for contin att a */
	*UnknownRate;	/* UnknownRate[a] = current unknown rate for att a */

extern char
	*Tested;	/* Tested[a] = true if att a already tested */
