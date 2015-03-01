/*************************************************************************/
/*								   	 */
/*	Pruning single rules					   	 */
/*	--------------------				   	      	 */
/*								   	 */
/*************************************************************************/


#include "defns.i"
#include "types.i"
#include "extern.i"
#include "rulex.i"


/*  External data structures used in building rules  */

extern ItemNo	*TargetClassFreq,	/* [Boolean] */
		*Errors,		/* [Condition] */
		*Total;			/* [Condition] */

extern float	*Pessimistic,		/* [Condition] */
		*Actual,		/* [Condition] */
		*CondSigLevel;		/* [Condition] */

extern Boolean	**CondSatisfiedBy,	/* [Condition][ItemNo] */
		*Deleted;

#define Before(n1,n2)  (n1->Tested < n2->Tested ||\
			n1->NodeType < n2->NodeType ||\
		        n1->Tested == n2->Tested && n1->Cut < n2->Cut)

#define IsTarget(case)  (Class(case) == TargetClass ? 1 : 0)



/*************************************************************************/
/*									 */
/*  Prune the rule given by the conditions Cond, and the number of	 */
/*  conditions NCond, and add the resulting rule to the current		 */
/*  ruleset if it is sufficiently accurate				 */
/*									 */
/*************************************************************************/


    PruneRule(Cond, NCond, TargetClass)
/*  ---------  */
    Condition Cond[];
    short NCond;
    ClassNo TargetClass;
{
    short d, dd, id, Bestd, Bestid, Remaining=NCond;
    float DefaultError, Extra, AddErrs(), TableProb();
    Boolean Alter, Satisfies(), NewRule(), Redundant();
    Condition Hold;
    ItemNo i;

    ForEach(d, 0, NCond)
    {
	Deleted[d] = false;
    }

    /*  Evaluate the satisfaction matrix  */

    TargetClassFreq[0] = TargetClassFreq[1] = 0;

    ForEach(i, 0, MaxItem)
    {
	ForEach(d, 1, NCond)
	{
	    CondSatisfiedBy[d][i] = Satisfies(Item[i], Cond[d]);
	}
	TargetClassFreq[IsTarget(Item[i])]++;
    }

    DefaultError = 1.0 - (TargetClassFreq[true] + 1.0) / (MaxItem + 3.0);

    /*  Find conditions to delete  */

    Verbosity(1) printf("\n  pruning rule for %s", ClassName[TargetClass]);

    do
    {
	Alter = false;

	FindTables(NCond, TargetClass);

	/*  Find the condition, deleting which would most improve
	    the accuracy of the rule.
	    Notes: The pessimistic error rate, and not the actual
		   error rate, is currently used.
	    	   When d is 0, we are dealing with all conditions.  */

	Bestd = id = 0;

	Verbosity(1)
	    printf("\n     Err Used   Pess\tAbsent condition\n");

	ForEach(d, 0, NCond)
	{
	    if ( Deleted[d] ) continue;

	    if ( Total[d] )
	    {
		Actual[d] = Errors[d] / (float) Total[d];
		Extra = AddErrs((float) Total[d], (float) Errors[d]);
		Pessimistic[d] = (Errors[d] + Extra) / Total[d];
	    }
	    else
	    {
		Actual[d] = 0;
		Pessimistic[d] = DefaultError;
	    }

	    Verbosity(1)
		printf("   %5d%5d  %4.1f%%",
		       Errors[d], Total[d], 100 * Pessimistic[d]);

	    if ( ! d )
	    {
		Verbosity(1) printf("\t<base rule>\n");
	    }
	    else
	    {
		id++;

		/*  If significance testing option used, invoke Fisher's
		    exact test here to assess probability that division
		    by d arises from chance.  */

		if ( SIGTEST )
		{
		    CondSigLevel[d] =
			TableProb(Errors[0],
				   Errors[d]-Errors[0],
				   Total[0]-Errors[0],
			           Total[d]-Total[0]-Errors[d]+Errors[0]);

		    Verbosity(1) printf("  Sig=%.3f", CondSigLevel[d]);
		}

		Verbosity(1) PrintCondition(Cond[d]);

		/*  Bestd identifies the condition with lowest pessimistic
		    error  estimate  */

		if ( ! Bestd || Pessimistic[d] <= Pessimistic[Bestd] )
		{
		    Bestd = d;
		    Bestid = id;
		}

		/*  Alter is set true if we are going to drop a condition
		    (either because we get lower pessimistic est, or
		    because one of the conditions fails a significance test)  */

		if ( Pessimistic[d] <= Pessimistic[0] ||
		     Actual[d] <= Actual[0]  ||
		     SIGTEST && CondSigLevel[d] > SIGTHRESH )
		{
		    Alter = true;
		}
	    }
	}

	if ( Alter )
	{
	    Verbosity(1) printf("\teliminate test %d\n", Bestid);

	    Deleted[Bestd] = true;
	    Remaining--;
	}

    } while ( Alter && Remaining );

    if ( ! Remaining || ! Total[0] )
    {
	return;
    }

    if ( Pessimistic[0] >= DefaultError )
    {
	Verbosity(1) printf("\ttoo inaccurate\n");
	return;
    }

    /*  Sort the conditions  */

    ForEach(d, 1, Remaining)
    {
	dd =  0;
	ForEach(id, d, NCond)
	{
	    if ( ! Deleted[id] &&
		 ( ! dd ||
		   Before(Cond[id]->CondTest, Cond[dd]->CondTest) ) )
	    {
		dd = id;
	    }
	}
	if ( dd != d )
	{
	    Hold    = Cond[d];
	    Cond[d] = Cond[dd];
	    Cond[dd] = Hold;
	    Deleted[dd] = Deleted[d];
	}
	Deleted[d] = true;
    }

    NewRule(Cond, Remaining, TargetClass, Pessimistic[0]);
}



/*************************************************************************/
/*									 */
/*  See whether condition R is redundant                           	 */
/*									 */
/*************************************************************************/


Boolean Redundant(R, Cond, NCond)
/*      ---------  */
    Condition Cond[];
    short R, NCond;
{
    short d, v, vv;
    Test t, Rt;
    Boolean IsSubset();

    Rt = Cond[R]->CondTest;
    v =  Cond[R]->TestValue;

    ForEach(d, 1, NCond)
    {
	if ( Deleted[d] || d == R ) continue;

	t = Cond[d]->CondTest;
	vv = Cond[d]->TestValue;

	if ( t->Tested != Rt->Tested ) continue;

	switch ( t->NodeType )
	{
	    case BrDiscr:  /* test of discrete attribute */

		return false;

	    case ThreshContin:  /* test of continuous attribute */

		if ( vv == v &&
		     ( v == 1 ? t->Cut < Rt->Cut : t->Cut > Rt->Cut ) )
		{
		    return true;
		}

		break;

	    case BrSubset:  /* subset test on discrete attribute  */

		if ( IsSubset(t->Subset[vv], Rt->Subset[v], Rt->Tested) )
		{
		    return true;
		}
	}
    }

    return false;
}



/*************************************************************************/
/*									 */
/*  Decide whether subset S1 of values is contained in subset S2	 */
/*									 */
/*************************************************************************/


Boolean IsSubset(S1, S2, Att)
/*      --------  */
    Set S1, S2;
    Attribute Att;
{
    DiscrValue v;

    ForEach(v, 1, MaxAttVal[Att])
    {
	if ( In(v, S1) && ! In(v, S2) ) return false;
    }

    return true;
}



/*************************************************************************/
/*									 */
/*  Find the frequency distribution tables for the current conditions: 	 */
/*									 */
/*	Total[0] = items matching all conditions		   	 */
/*	Total[d] = items matching all except condition d	   	 */
/*									 */
/*	Errors[0] = wrong-class items matching all conditions	   	 */
/*	Errors[d] = wrong-class items matching all but cond d	   	 */
/*									 */
/*  This routine is critical to the efficiency of rule pruning. It	 */
/*  computes the information above in one pass through the data,	 */
/*  looking at cases that fail to satisfy at most one of the		 */
/*  non-deleted conditions						 */
/*									 */
/*************************************************************************/


    FindTables(NCond, TargetClass)
/*  -----------  */
    short NCond;
    ClassNo TargetClass;
{
    ItemNo i;
    short Misses, Missed[2], d;
    Boolean CorrectClass;

    /*  Clear distributions  */

    ForEach(d, 0, NCond)
    {
	Total[d] = Errors[d] = 0;
    }

    /*  Set distributions  */

    ForEach(i, 0, MaxItem)
    {
	Misses = 0;
	CorrectClass = IsTarget(Item[i]);

	for ( d = 1 ; d <= NCond && Misses <= 1 ; d++ )
	{
	    if ( ! Deleted[d] && ! CondSatisfiedBy[d][i] )
	    {
		Missed[Misses++] = d;
	    }
	}

	if ( ! Misses )
	{
	    UpdateCount(Total, Errors, 0, CorrectClass);
	}
	else
	if ( Misses == 1 )
	{
	    UpdateCount(Total, Errors, Missed[0], CorrectClass);
	}
    }

    /*  Adjust counts to reflect cases that met all conditions  */

    ForEach(d, 1, NCond)
    {
	if ( ! Deleted[d] )
	{
	    Total[d] += Total[0];
	    Errors[d] += Errors[0];
	}
    }
}



/*************************************************************************/
/*									 */
/*  Increment the counts Total[d] and Errors[d]				 */
/*									 */
/*************************************************************************/


    UpdateCount(T, E, d, OK)
/*  -----------  */
    ItemNo T[], E[];
    short d;
    Boolean OK;
{
    T[d]++;
    if ( ! OK ) E[d]++;
}



/*************************************************************************/
/*									 */
/*  Determine whether the given case description satisfies the given	 */
/*  condition.								 */
/*									 */
/*************************************************************************/


Boolean Satisfies(CaseDesc, OneCond)
/*      ---------  */
    Description CaseDesc; 
    Condition OneCond;
{
    DiscrValue v;
    float cv;
    Test t;
    short s;
    Boolean Outcome;

    t = OneCond->CondTest;

    /*  Determine the outcome of this test on this item  */

    switch ( t->NodeType )
    {
	case BrDiscr:  /* test of discrete attribute */

	    v = DVal(CaseDesc, t->Tested);
	    Outcome = ( v == 0 ? -1 : v );
	    break;

	case ThreshContin:  /* test of continuous attribute */

	    cv = CVal(CaseDesc, t->Tested);
	    Outcome = ( cv == Unknown ? -1 : cv <= t->Cut ? 1 : 2 );
	    break;

	case BrSubset:  /* subset test on discrete attribute  */

	    v = DVal(CaseDesc, t->Tested);
	    Outcome = -1;
	    ForEach(s, 1, t->Forks)
	    {
		if ( In(v, t->Subset[s]) )
		{
		    Outcome = s;
		    break;
		}
	    }
    }

    return ( Outcome == OneCond->TestValue );
}



/*************************************************************************/
/*									 */
/*  Hypergeometric distribution	(uses tabulated log factorials)		 */
/*									 */
/*************************************************************************/


double Hypergeom(a, r, A, B)
/*               ---------  */
    int a, r, A, B;
{
    return exp( LogFact[A] + LogFact[B] + LogFact[r] + LogFact[A+B-r] -
	        ( LogFact[a] + LogFact[r-a] + LogFact[A-a]
		  + LogFact[B-(r-a)] + LogFact[A+B]) );
}



/*************************************************************************/
/*									 */
/*  TableProb examines the 2x2 contingency table t and computes the      */
/*  probability that a random division could have produced a split at	 */
/*  least as extreme as this.  Also known as "Fisher's Exact Test",	 */
/*  after its inventor, R.A. Fisher.					 */
/*									 */
/*************************************************************************/


float TableProb(t11, t12, t21, t22)
/*    ---------  */
    int t11, t12, t21, t22;
{
    double Sum=0.0;
    int A, B, r, a, k, a0;

    /*  First, rearrange the rows and columns of the table to get it into
	canonical form  */

    if ( t11 + t12 > t21 + t22 )
    {
	A = t11 + t12;
	B = t21 + t22;

	if ( t11 * (t21 + t22) > t21 * (t11 + t12) )
	{
	    a0 = t11;
	    r  = t11 + t21;
	}
	else
	{
	    a0 = t12;
	    r  = t12 + t22;
	}
    }
    else
    {
	A = t21 + t22;
	B = t11 + t12;
	if ( t21 * (t11 + t12) > t11 * (t21 + t22) )
	{
	    a0 = t21;
	    r  = t21 + t11;
	}
	else
	{
	    a0 = t22;
	    r  = t22 + t12;
	}
    }

    /*  Now compute the probability  */

    k = Min(r, A);
    ForEach(a, a0, k)
    {
	Sum += Hypergeom(a, r, A, B);
    }

    return Sum;
}
