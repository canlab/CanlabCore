/*************************************************************************/
/*									 */
/*	Generate all rulesets from the decision trees		  	 */
/*	---------------------------------------------		  	 */
/*								  	 */
/*************************************************************************/


#include "defns.i"
#include "types.i"
#include "extern.i"
#include "rulex.i"


/*************************************************************************/
/*								  	 */
/*  For each tree, form a set of rules and process them, then form a	 */
/*  composite set of rules from all of these sets.		  	 */
/*  If there is only one tree, then no composite set is formed.	  	 */
/*								  	 */
/*  Rulesets are stored in  PRSet[0]  to  PRSet[TRIALS], where    	 */
/*  PRSet[TRIALS] contains the composite ruleset.		  	 */
/*								  	 */
/*  On completion, the current ruleset is the composite ruleset (if one	 */
/*  has been made), otherwise the ruleset from the single tree. 	 */
/*								  	 */
/*************************************************************************/


    GenerateRules()
/*  -------------  */
{
    Tree DecisionTree, GetTree();
    short t=0, RuleSetSpace=0, r;

    /*  Find bits to encode attributes and branches  */

    FindTestCodes();

    /*  Now process each decision tree  */

    while ( DecisionTree = GetTree(".unpruned") )
    {
	printf("\n------------------\n");
	printf("Processing tree %d\n", t);

	/*  Form a set of rules from the next tree  */

	FormRules(DecisionTree);

	/*  Process the set of rules for this trial  */

	ConstructRuleset();

	printf("\nFinal rules from tree %d:\n", t);
	PrintIndexedRules();
	    
	/*  Make sure there is enough room for the new ruleset  */

	if ( t + 1 >= RuleSetSpace )
	{
	    RuleSetSpace += 10;

	    if ( RuleSetSpace > 10 )
	    {
		PRSet = (RuleSet *) realloc(PRSet, RuleSetSpace * sizeof(RuleSet));
	    }
	    else
	    {
		PRSet = (RuleSet *) malloc(RuleSetSpace * sizeof(RuleSet));
	    }

	}

	PRSet[t].SNRules = NRules;
	PRSet[t].SRule = Rule;
	PRSet[t].SRuleIndex = RuleIndex;
	PRSet[t].SDefaultClass = DefaultClass;

	++t;
    }

    if ( ! t )
    {
	printf("\nERROR:  can't find any decision trees\n");
	exit(1);
    }

    TRIALS = t;

    /*  If there is more than one tree in the trees file,
	make a composite ruleset of the rules from all trees  */

    if ( TRIALS > 1 )
    {
	CompositeRuleset();
    }
}



/*************************************************************************/
/*								  	 */
/*	Determine code lengths for attributes and branches		 */
/*								  	 */
/*************************************************************************/


    FindTestCodes()
/*  -------------  */
{
    Attribute Att;
    DiscrValue v, V;
    ItemNo i, *ValFreq;
    int PossibleCuts;
    float Sum, SumBranches=0, p;
    void SwapUnweighted();

    BranchBits  = (float *) malloc((MaxAtt+1) * sizeof(float));

    ForEach(Att, 0, MaxAtt)
    {
	if ( (V = MaxAttVal[Att]) )
	{
	    ValFreq = (ItemNo *) calloc(V+1, sizeof(ItemNo));

	    ForEach(i, 0, MaxItem)
	    {
		ValFreq[DVal(Item[i],Att)]++;
	    }

	    Sum = 0;
	    ForEach(v, 1, V)
	    {
		if ( ValFreq[v] )
		{
		    Sum += (ValFreq[v] / (MaxItem+1.0)) *
			   (LogItemNo[MaxItem+1] - LogItemNo[ValFreq[v]]);
		}
	    }
	    free(ValFreq);

	    BranchBits[Att] = Sum;
	}
	else
	{
	    Quicksort(0, MaxItem, Att, SwapUnweighted);

	    PossibleCuts = 1;
	    ForEach(i, 1, MaxItem)
	    {
		if ( CVal(Item[i],Att) > CVal(Item[i-1],Att) )
		{
		    PossibleCuts++;
		}
	    }

	    BranchBits[Att] = PossibleCuts > 1 ?
			      1 + LogItemNo[PossibleCuts] / 2 : 0 ;
	}

	SumBranches += BranchBits[Att];
    }

    AttTestBits = 0;
    ForEach(Att, 0, MaxAtt)
    {
	if ( (p = BranchBits[Att] / SumBranches) > 0 )
	{
	    AttTestBits -= p * log(p) / log(2.0);
	}
    }
}



/*************************************************************************/
/*                                                                	 */
/*  Exchange items at a and b.  Note:  unlike the similar routine in	 */
/*  buildtree, this does not assume that items have a Weight to be	 */
/*  swapped as well!							 */
/*                                                                	 */
/*************************************************************************/


void SwapUnweighted(a, b)
/*   --------------  */
    ItemNo a, b;
{
    Description Hold;

    Hold = Item[a];
    Item[a] = Item[b];
    Item[b] = Hold;
}



/*************************************************************************/
/*								  	 */
/*	Form composite ruleset for all trials			  	 */
/*								  	 */
/*************************************************************************/


    CompositeRuleset()
/*  ----------------  */
{
    RuleNo r;
    short t, ri;
    Boolean NewRule();
    
    InitialiseRules();

    /*  Lump together all the rules from each ruleset  */

    ForEach(t, 0, TRIALS-1)
    {
	ForEach(ri, 1, PRSet[t].SNRules)
	{
	    r = PRSet[t].SRuleIndex[ri];
	    NewRule(PRSet[t].SRule[r].Lhs, PRSet[t].SRule[r].Size,
		     PRSet[t].SRule[r].Rhs, PRSet[t].SRule[r].Error);
	}
    }

    /*  ... and select a subset in the usual way  */

    ConstructRuleset();

    printf("\nComposite ruleset:\n");
    PrintIndexedRules();

    PRSet[TRIALS].SNRules    = NRules;
    PRSet[TRIALS].SRule      = Rule;
    PRSet[TRIALS].SRuleIndex = RuleIndex;
    PRSet[TRIALS].SDefaultClass = DefaultClass;
}
