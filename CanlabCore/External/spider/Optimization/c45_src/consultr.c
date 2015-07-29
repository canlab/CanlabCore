/*************************************************************************/
/*								   	 */
/*	Classify items interactively using a set of rules		 */
/*	-------------------------------------------------		 */
/*								   	 */
/*************************************************************************/


#include "defns.i"
#include "types.i"


		/*  External data  */

short		MaxAtt, MaxClass, MaxDiscrVal;

ItemNo		MaxItem;

Description	*Item;

DiscrValue	*MaxAttVal;

String		*ClassName,
		*AttName,
		**AttValName,
		FileName = "DF";

char		*SpecialStatus;

short		VERBOSITY = 0,
		TRACE     = 0;


Boolean		FirstTime = true;

PR		*Rule;

RuleNo		NRules = 0,
		*RuleIndex;
 
short		RuleSpace = 0;
 
ClassNo		DefaultClass;


	/*  The interview module uses a more complex description of an
	    case called a "Range Description" with lower and upper bounds.
	    For further information, see consult.c  */

typedef	struct ValRange *RangeDescRec;

struct ValRange
	{
	    Boolean Known, Asked;
	    float LowerBound, UpperBound, *Probability;
	};

RangeDescRec RangeDesc;

float Confidence;

#define	MINCF	0.50		/* minimum cf for useable rule */


/*************************************************************************/
/*								  	 */
/*  Find the best rule for the current case.			  	 */
/*  Note:  leave probability in Confidence				 */
/*								  	 */
/*************************************************************************/


RuleNo BestRule()
/*     -------- 	 */
{
    RuleNo r;
    float cf, RuleStrength();

    Confidence = 0.0;
    
    ForEach(r, 1, NRules)
    {
	cf = RuleStrength(Rule[r]);

	if ( cf > 0.3 )
	{
	    Confidence = cf;
	    return r;
	}
    }

    return 0;
}



/*************************************************************************/
/*								    	 */
/*  Given a rule, determine the strength with which we can conclude	 */
/*  that the current case belongs to the class specified by the RHS	 */
/*  of the rule								 */
/*								    	 */
/*************************************************************************/


float RuleStrength(Rule)
/*    ------------ 	 */
    PR Rule;
{
    short d;
    float RuleProb=1.0, ProbSatisfied();

    ForEach(d, 1, Rule.Size)
    {
	RuleProb *= ProbSatisfied(Rule.Lhs[d]);
	if ( RuleProb < MINCF )
	{
	    return 0.0;
	}
    }

    return ( (1 - Rule.Error) * RuleProb );
}



/*************************************************************************/
/*							   		 */
/*  Determine the probability of the current case description		 */
/*  satisfying the given condition					 */
/*							   		 */
/*************************************************************************/


float ProbSatisfied(c)
/*    ------------- 	 */
    Condition c;
{
    Attribute a;
    char v;
    float AddProb=0.0;
    Test t;
    DiscrValue i;

    t = c->CondTest;
    a = t->Tested;
    v = c->TestValue;

    CheckValue(a, Nil);

    if ( ! RangeDesc[a].Known )
    {
	return 0.0;
    }

    switch ( t->NodeType )
    {
	case BrDiscr:  /* test of discrete attribute */

	    return RangeDesc[a].Probability[v];

	case ThreshContin:  /* test of continuous attribute */

	    if ( RangeDesc[a].UpperBound <= t->Cut )
	    {
		return ( v == 1 ? 1.0 : 0.0 );
	    }
	    else
	    if ( RangeDesc[a].LowerBound > t->Cut )
	    {
		return ( v == 2 ? 1.0 : 0.0 );
	    }
	    else
	    if ( v == 1 )
	    {
		return (t->Cut - RangeDesc[a].LowerBound) /
		       (RangeDesc[a].UpperBound - RangeDesc[a].LowerBound);
	    }
	    else
	    {
		return (RangeDesc[a].UpperBound - t->Cut) /
		       (RangeDesc[a].UpperBound - RangeDesc[a].LowerBound);
	    }

	case BrSubset:  /* subset test on discrete attribute  */

	    ForEach(i, 1, MaxAttVal[a])
	    {
		if ( In(i, t->Subset[v]) )
		{
		    AddProb += RangeDesc[a].Probability[i];
		}
	    }
	    return AddProb;

    } 
    return 0.0;
}



/*************************************************************************/
/*								  	 */
/*		Process a single case				  	 */
/*								  	 */
/*************************************************************************/


    InterpretRuleset()
/*  ---------------- 	 */
{ 
    char Reply;
    Attribute a;
    RuleNo r;

    /*  Initialise  */

    ForEach(a, 0, MaxAtt)
    {
	RangeDesc[a].Asked = false;
    }

    if ( FirstTime )
    {
	FirstTime = false;
	printf("\n");
    }
    else
    {
	printf("\n-------------------------------------------\n\n");
    }

    /*  Find the first rule that fires on the item  */

    if ( r = BestRule() )
    {
	printf("\nDecision:\n");
	printf("\t%s", ClassName[Rule[r].Rhs]);

	if ( Confidence < 1.0 )
	{
	    printf("  CF = %.2f", Confidence);
	}

	printf("\n");
    }
    else
    {
	printf("\nDecision:\n");
	printf("\t%s (default class)\n", ClassName[DefaultClass]);
    }

    /*  Prompt for what to do next  */

    while ( true )
    {
	printf("\nRetry, new case or quit [r,n,q]: ");
	Reply = getchar();
	SkipLine(Reply);
	switch ( Reply )
	{
	  case 'r':  return;
	  case 'n':  Clear(); return;
	  case 'q':  exit(0);
	  default:   printf("Please enter 'r', 'n' or 'q'");
	}
    }
}
    


/*************************************************************************/
/*								  	 */
/*  Main routine for classifying items using a set of rules	  	 */
/*								  	 */
/*************************************************************************/

    main(Argc, Argv)
/*  ---- 	 */
    int Argc;
    char *Argv[];
{
    int o;
    extern char *optarg;
    extern int optind;
    Attribute a;
    RuleNo r;

    PrintHeader("production rule interpreter");

    /*  Process options  */

    while ( (o = getopt(Argc, Argv, "tvf:")) != EOF )
    {
	switch (o)
	{
	    case 't':	TRACE = 1;
			break;
	    case 'v':	VERBOSITY = 1;
			break;
	    case 'f':	FileName = optarg;
			break;
	    case '?':	printf("unrecognised option\n");
			exit(1);
	}
    }

    /*  Initialise  */

    GetNames();
    GetRules();

    if ( TRACE )
    {
	ForEach(r, 1, NRules)
	{
	    PrintRule(r);
	}
	printf("\nDefault class: %s\n", ClassName[DefaultClass]);
    }

    /*  Allocate value ranges  */

    RangeDesc = (struct ValRange *) calloc(MaxAtt+1, sizeof(struct ValRange));

    ForEach(a, 0, MaxAtt)
    {
	if ( MaxAttVal[a] )
	{
	    RangeDesc[a].Probability =
		(float *) calloc(MaxAttVal[a]+1, sizeof(float));
	}
    }

    /*  Consult  */

    Clear();
    while ( true )
    {
	InterpretRuleset();
    }
}
