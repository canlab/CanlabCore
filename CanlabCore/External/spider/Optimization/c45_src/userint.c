/*************************************************************************/
/*								 	 */
/*	User interface for consulting trees and rulesets	 	 */
/*	------------------------------------------------	 	 */
/*								 	 */
/*************************************************************************/


#include "defns.i"
#include "types.i"
#include "extern.i"


typedef	struct	ValRange	*RangeDescRec;
struct	ValRange
	{
	    Boolean		Known, Asked;
	    float		LowerBound, UpperBound, *Probability;
	};

extern	RangeDescRec		RangeDesc;


#define Fuzz 0.01


/*************************************************************************/
/*									 */
/*	Ask for the value of attribute Att if necessary			 */
/*									 */
/*************************************************************************/


    CheckValue(Att, T)
/*  ----------  */
    Attribute Att;
    Tree T;
{
    if ( RangeDesc[Att].Asked ) return;

    printf("%s", AttName[Att]);
    if ( RangeDesc[Att].Known )
    {
	printf(" [ ");
	PrintRange(Att);
	printf(" ]");
    }
    printf(": ");

    ReadRange(Att, T);
}


 
/*************************************************************************/
/*									 */
/*	Print the range of values for attribute Att			 */
/*									 */
/*************************************************************************/


    PrintRange(Att)
/*  -----------  */
    Attribute Att;
{
    DiscrValue dv;
    Boolean First=true;
    float p;

    if ( MaxAttVal[Att] )  /*  discrete attribute  */
    {
	ForEach(dv, 1, MaxAttVal[Att] )
	{
	    if ( (p = RangeDesc[Att].Probability[dv]) > Fuzz )
	    {
		if ( ! First )
		{
		    printf(", ");
		}
		First = false;

		printf("%s", AttValName[Att][dv]);
		if ( p < 1-Fuzz )
		{
		    printf(": %.2f", p);
		}
	    }
	}
    }
    else  /*  continuous attribute  */
    {
	printf("%g", RangeDesc[Att].LowerBound);
	if ( RangeDesc[Att].UpperBound > RangeDesc[Att].LowerBound + Fuzz )
	{
	    printf(" - %g", RangeDesc[Att].UpperBound);
	}
    }
}



extern	char		Delimiter;
#define	SkipSpace	while ( (c = getchar()) == ' ' || c == '\t' )


/*************************************************************************/
/*									 */
/*	Read a range of values for attribute Att or <cr>		 */
/*									 */
/*************************************************************************/


    ReadRange(Att, T)
/*  ----------  */
    Attribute Att;
    Tree T;
{
    char c;

    RangeDesc[Att].Asked=true;

    SkipSpace;

    if ( c == '\n' )
    {
	return;
    }
    if ( c == '?' )
    {
	if ( (c = getchar()) == 't' )
	{
	    if ( T ) PrintTree(T);
	    SkipLine(c);
	    RangeDesc[Att].Asked = false;
	    CheckValue(Att, T);
	}
	else
	{
	    RangeDesc[Att].Known = false;
	    SkipLine(c);
	}
	return;
    }

    ungetc(c, stdin);
    RangeDesc[Att].Known = true;

    if ( MaxAttVal[Att] )
    {
	ReadDiscr(Att, T);
    }
    else
    {
	ReadContin(Att, T);
    }
}



/*************************************************************************/
/*									 */
/*	Read a discrete attribute value or range			 */
/*									 */
/*************************************************************************/

    ReadDiscr(Att, T)
/*  ---------  */
    Attribute Att;
    Tree T;
{
    char Name[500];
    Boolean ReadName();
    DiscrValue dv, PNo;
    float P, PSum;

    ForEach(dv, 1, MaxAttVal[Att])
    {
	RangeDesc[Att].Probability[dv] = 0.0;
    }

    do
    {
	ReadName(stdin, Name);

	dv = Which(Name, AttValName[Att], 1, MaxAttVal[Att]);
	if ( ! dv )
	{
	    printf("\tPermissible values are %s", AttValName[Att][1]);
	    ForEach(dv, 2, MaxAttVal[Att])
	    {
		printf(", %s", AttValName[Att][dv]);
	    }
	    printf("\n");

	    SkipLine(Delimiter);
	    Retry(Att, T);
	    return;
	}

	if ( Delimiter == ':' )
	{
	    ReadName(stdin, Name);
	    sscanf(Name, "%f", &P);	/* get probability */
	}
	else
	{
	    P = 1.0;		/*  only one attribute value  */
	}

	RangeDesc[Att].Probability[dv] = P;
    }
    while ( Delimiter == ',' );

    /*  Check that sum of probabilities is not > 1  */

    PNo = MaxAttVal[Att];
    PSum = 1.0;
    ForEach(dv, 1, MaxAttVal[Att])
    {
	if ( RangeDesc[Att].Probability[dv] > Fuzz )
	{
	    PSum -= RangeDesc[Att].Probability[dv];
	    PNo--;
	}
    }

    if ( PSum < 0 || ! PNo && PSum > Fuzz )
    {
	printf("Probability values must sum to 1\n");
	SkipLine(Delimiter);
	Retry(Att, T);
	return;
    }

    /*  Distribute the remaining probability equally among
	the unspecified attribute values  */

    PSum /= PNo;
    ForEach(dv, 1, MaxAttVal[Att])
    {
	if ( RangeDesc[Att].Probability[dv] < Fuzz )
	{
	    RangeDesc[Att].Probability[dv] = PSum;
	}
    }
}



/*************************************************************************/
/*									 */
/*	Read a continuous attribute value or range			 */
/*									 */
/*************************************************************************/


    ReadContin(Att, T)
/*  ----------  */
    Attribute Att;
    Tree T;
{
    char c;

    scanf("%f", &RangeDesc[Att].LowerBound);
    SkipSpace;

    if ( c == '-'  )
    {
	scanf("%f", &RangeDesc[Att].UpperBound);
	SkipSpace;
    }
    else
    {
	RangeDesc[Att].UpperBound = RangeDesc[Att].LowerBound;
    }

    if ( c != '\n' )
    {
	printf("Must be a continuous value or range\n");
	SkipLine(c);
	Retry(Att, T);
    }
}



/*************************************************************************/
/*									 */
/*	Try again to obtain a value for attribute Att			 */
/*									 */
/*************************************************************************/

    Retry(Att, T)
/*  -----  */
    Attribute Att;
    Tree T;
{
    RangeDesc[Att].Asked = false;
    RangeDesc[Att].Known = false;
    CheckValue(Att, T);
}



/*************************************************************************/
/*									 */
/*	Skip to the end of the line of input				 */
/*									 */
/*************************************************************************/


    SkipLine(c)
/*  --------  */
    char c;
{
    while ( c != '\n' ) c = getchar();
}



/*************************************************************************/
/*									 */
/*		Clear the range description				 */
/*									 */
/*************************************************************************/

    Clear()
/*  -----  */
{
    Attribute Att;

    ForEach(Att, 0, MaxAtt)
    {
	RangeDesc[Att].Known = false;
    }
}
