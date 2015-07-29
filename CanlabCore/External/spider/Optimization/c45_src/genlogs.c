/*************************************************************************/
/*									 */
/*	Tabluate logs and log factorials (to improve speed)		 */
/*	--------------------------------				 */
/*									 */
/*************************************************************************/


#include "defns.i"
#include "types.i"
#include "extern.i"


float	*LogItemNo;
double	*LogFact;


/*************************************************************************/
/*									 */
/*  Set up the array LogItemNo to contain the logs of integers and	 */
/*  the array LogFact to contain logs of factorials (all to base 2)	 */
/*									 */
/*************************************************************************/


    GenerateLogs()
/*  ------------  */
{
    ItemNo i;

    LogItemNo = (float *) malloc((MaxItem+100) * sizeof(float));
    LogFact = (double *) malloc((MaxItem+100) * sizeof(double));

    LogItemNo[0] = -1E38;
    LogItemNo[1] = 0;
    LogFact[0] = LogFact[1] = 0;

    ForEach(i, 2, MaxItem+99)
    {
	LogItemNo[i] = log((float) i) / Log2;
	LogFact[i] = LogFact[i-1] + LogItemNo[i];
    }
}
