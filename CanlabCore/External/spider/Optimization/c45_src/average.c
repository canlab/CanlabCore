/*************************************************************************/
/*									 */
/*	Average results for training and test sets			 */
/*	------------------------------------------			 */
/*									 */
/*	This is a generic program that averages any numbers found on	 */
/*	a set of lines of the same pattern.				 */
/*									 */
/*************************************************************************/


#include <stdio.h>

#define	MAXLINE		200	/* max line length */
#define MAXVALS		 10	/* max values to be averaged */


main()
{
    char Line[MAXLINE], *p1, *p2;
    int Numbers=0, Lines=0, i, TrainTest;
    float Val, Sum[2][MAXVALS];
    double strtod();

    for ( i = 0 ; i < MAXVALS ; i++ )
    {
	Sum[0][i] = Sum[1][i] = 0;
    }

    while ( fgets(Line, MAXLINE, stdin) )
    {
	i = 0;
	TrainTest = Lines % 2;
	printf("%s", Line);

	/*  Count the numbers appearing on the line  */

	for ( p1 = Line ; *p1 != '\n' ; p1++ )
	{
	    if ( *p1 < '0' || *p1 > '9' ) continue;

	    Val = strtod(p1, &p2);
	    Sum[TrainTest][i++] += Val;
	    p1 = p2-1;
	}

	/*  The number of numbers must match any previous lines  */

	if ( Lines )
	{
	    if ( i != Numbers ) exit();
	}
	else
	{
	    Numbers = i;
	}

	Lines++;
    }

    putchar('\n');
    for ( TrainTest = 0 ; TrainTest <= 1 ; TrainTest++ )
    {
	i = 0;
	printf("%s:\t", TrainTest ? "test" : "train");

	for ( p1 = Line ; *p1 != '\n' ; p1++ )
	{
	    if ( *p1 < '0' || *p1 > '9' )
	    {
		putchar(*p1);
	    }
	    else
	    {
		printf("%.1f", Sum[TrainTest][i++] / (0.5 * Lines));
		strtod(p1, &p2);
		p1 = p2-1;
	    }
	}
	putchar('\n');
    }
}
