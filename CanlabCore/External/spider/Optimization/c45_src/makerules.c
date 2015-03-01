/*************************************************************************/
/*								  	 */
/*	Form a set of rules from a decision tree			 */
/*	----------------------------------------			 */
/*								  	 */
/*************************************************************************/


#include "defns.i"
#include "types.i"
#include "extern.i"
#include "rulex.i"


ItemNo	*TargetClassFreq,	/* [Boolean] */
	*Errors,		/* [Condition] */
	*Total;			/* [Condition] */

float	*Pessimistic,		/* [Condition] */
	*Actual,		/* [Condition] */
	*CondSigLevel;		/* [Condition] */

Boolean	**CondSatisfiedBy,	/* [Condition][ItemNo] */
	*Deleted;		/* [Condition] */

DiscrValue *SingleValue;	/* [Attribute] */

Condition *Stack;

short	MaxDisjuncts,
	MaxDepth;



/*************************************************************************/
/*								  	 */
/*	Form a ruleset from decision tree t			  	 */
/*								  	 */
/*************************************************************************/


    FormRules(t)
/*  ----------  */
    Tree t;
{
    short i;

    /*  Find essential parameters and allocate storage  */

    MaxDepth = 0;
    MaxDisjuncts = 0;

    TreeParameters(t, 0);

    Actual = (float *) calloc(MaxDepth+2, sizeof(float));
    Total = (ItemNo *) calloc(MaxDepth+2, sizeof(ItemNo));
    Errors = (ItemNo *) calloc(MaxDepth+2, sizeof(ItemNo));
    Pessimistic = (float *) calloc(MaxDepth+2, sizeof(float));

    CondSigLevel = (float *) calloc(MaxDepth+2, sizeof(float));

    TargetClassFreq = (ItemNo *) calloc(2, sizeof(ItemNo));

    Deleted = (Boolean *) calloc(MaxDepth+2, sizeof(Boolean));
    CondSatisfiedBy = (char **) calloc(MaxDepth+2, sizeof(char *));
    Stack = (Condition *) calloc(MaxDepth+2, sizeof(Condition));

    ForEach(i, 0, MaxDepth+1)
    {
	CondSatisfiedBy[i] = (char *) calloc(MaxItem+1, sizeof(char));
	Stack[i] = (Condition) malloc(sizeof(struct CondRec));
    }

    SingleValue = (DiscrValue *) calloc(MaxAtt+1, sizeof(DiscrValue));

    InitialiseRules();

    /*  Extract and prune disjuncts  */

    Scan(t, 0);

    /*  Deallocate storage  */

    ForEach(i, 0, MaxDepth+1)
    {
	cfree(CondSatisfiedBy[i]);
	cfree(Stack[i]);
    }
    cfree(Deleted);
    cfree(CondSatisfiedBy);
    cfree(Stack);

    cfree(Actual);
    cfree(Total);
    cfree(Errors);
    cfree(Pessimistic);

    cfree(CondSigLevel);

    cfree(TargetClassFreq);
}



/*************************************************************************/
/*                                                                	 */
/*  Find the maximum depth and the number of leaves in tree t	  	 */
/*  with initial depth d					  	 */
/*                                                                	 */
/*************************************************************************/


    TreeParameters(t, d)
/*  ---------------  */
    Tree t;
    short d;
{
    DiscrValue v;

    if ( t->NodeType )
    {
	ForEach(v, 1, t->Forks)
	{
	    TreeParameters(t->Branch[v], d+1);
	}
    }
    else
    {
	/*  This is a leaf  */

	if ( d > MaxDepth ) MaxDepth = d;
	MaxDisjuncts++;
    }
}



/*************************************************************************/
/*								  	 */
/*  Extract disjuncts from tree t at depth d, and process them	  	 */
/*								  	 */
/*************************************************************************/


    Scan(t, d)
/*  ----  */
    Tree t;
    short d;
{
    DiscrValue v;
    short i;
    Condition *Term;
    Test x, FindTest();

    if ( t->NodeType )
    {
	d++;

	x = (Test) malloc(sizeof(struct TestRec));
	x->NodeType = t->NodeType;
	x->Tested = t->Tested;
	x->Forks = t->Forks;
	x->Cut = ( t->NodeType == ThreshContin ? t->Cut : 0 );
	if ( t->NodeType == BrSubset )
	{
	    x->Subset = (Set *) calloc(t->Forks + 1, sizeof(Set));
	    ForEach(v, 1, t->Forks)
	    {
		x->Subset[v] = t->Subset[v];
	    }
	}

	Stack[d]->CondTest = FindTest(x);

	ForEach(v, 1, t->Forks)
	{
	    Stack[d]->TestValue = v;
	    Scan(t->Branch[v], d);
	}
    }
    else
    if ( t->Items >= 1 )
    {
	/*  Leaf of decision tree - construct the set of
 	    conditions associated with this leaf and prune  */

	Term = (Condition *) calloc(d+1, sizeof(Condition));
	ForEach(i, 1, d)
	{
	    Term[i] = (Condition) malloc(sizeof(struct CondRec));
	    Term[i]->CondTest = Stack[i]->CondTest;
	    Term[i]->TestValue = Stack[i]->TestValue;
	}

	PruneRule(Term, d, t->Leaf);

	cfree(Term);
    }
}
