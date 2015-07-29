/*************************************************************************/
/*								  	 */
/*	Miscellaneous routines for rule handling		  	 */
/*	----------------------------------------		  	 */
/*								  	 */
/*************************************************************************/


#include "defns.i"
#include "types.i"
#include "extern.i"
#include "rulex.i"

extern  FILE	*TRf;		/* rules file */


Test	*TestVec;
short	NTests = 0;

FILE	*fopen();
extern char	Fn[500];	/* file name */


/*************************************************************************/
/*								  	 */
/*  Save the current ruleset in rules file in order of the index  	 */
/*								  	 */
/*************************************************************************/


    SaveRules()
/*  ----------  */
{
    short ri, d, v, Bytes;
    RuleNo r;
    Test Tst;

    if ( TRf ) fclose(TRf);

    strcpy(Fn, FileName);
    strcat(Fn, ".rules");
    if ( ! ( TRf = fopen(Fn, "w") ) ) Error(0, Fn, " for writing");
    
    StreamOut((char *) &NRules, sizeof(RuleNo));
    StreamOut((char *) &DefaultClass, sizeof(ClassNo));

    ForEach(ri, 1, NRules)
    {
	r = RuleIndex[ri];
        StreamOut((char *) &Rule[r].Size, sizeof(short));
	ForEach(d, 1, Rule[r].Size)
	{
	    Tst = Rule[r].Lhs[d]->CondTest;

	    StreamOut((char *) &Tst->NodeType, sizeof(short));
	    StreamOut((char *) &Tst->Tested, sizeof(Attribute));
	    StreamOut((char *) &Tst->Forks, sizeof(short));
	    StreamOut((char *) &Tst->Cut, sizeof(float));
	    if ( Tst->NodeType == BrSubset )
	    {
		Bytes = (MaxAttVal[Tst->Tested]>>3) + 1;
		ForEach(v, 1, Tst->Forks)
		{
		    StreamOut((char *) Tst->Subset[v], Bytes);
		}
	    }
	    StreamOut((char *) &Rule[r].Lhs[d]->TestValue, sizeof(short));
	}
	StreamOut((char *) &Rule[r].Rhs, sizeof(ClassNo));
	StreamOut((char *) &Rule[r].Error, sizeof(float));
    }

    SaveDiscreteNames();
}



/*************************************************************************/
/*                                                                	 */
/*	Get a new ruleset from rules file			  	 */
/*                                                                	 */
/*************************************************************************/


    GetRules()
/*  ---------  */
{
    RuleNo nr, r;
    short n, d, v, Bytes;
    Condition *Cond;
    Test Tst, FindTest();
    ClassNo c;
    float e;
    Boolean NewRule();

    if ( TRf ) fclose(TRf);

    strcpy(Fn, FileName);
    strcat(Fn, ".rules");
    if ( ! ( TRf = fopen(Fn, "r") ) ) Error(0, Fn, "");
    
    StreamIn((char *) &nr, sizeof(RuleNo));
    StreamIn((char *) &DefaultClass, sizeof(ClassNo));

    ForEach(r, 1, nr)
    {
        StreamIn((char *) &n, sizeof(short));
	Cond = (Condition *) calloc(n+1, sizeof(Condition));
	ForEach(d, 1, n)
	{
	    Tst = (Test) malloc(sizeof(struct TestRec));

	    StreamIn((char *) &Tst->NodeType, sizeof(short));
	    StreamIn((char *) &Tst->Tested, sizeof(Attribute));
	    StreamIn((char *) &Tst->Forks, sizeof(short));
	    StreamIn((char *) &Tst->Cut, sizeof(float));
	    if ( Tst->NodeType == BrSubset )
	    {
		Tst->Subset = (Set *) calloc(Tst->Forks + 1, sizeof(Set));

		Bytes = (MaxAttVal[Tst->Tested]>>3) + 1;
		ForEach(v, 1, Tst->Forks)
		{
		    Tst->Subset[v] = (Set) malloc(Bytes);
		    StreamIn((char *) Tst->Subset[v], Bytes);
		}
	    }

	    Cond[d] = (Condition) malloc(sizeof(struct CondRec));
	    Cond[d]->CondTest = FindTest(Tst);
	    StreamIn((char *) &Cond[d]->TestValue, sizeof(short));
	}
	StreamIn((char *) &c, sizeof(ClassNo));
	StreamIn((char *) &e, sizeof(float));
	NewRule(Cond, n, c, e);
	cfree(Cond);
    }

    RecoverDiscreteNames();
}



/*************************************************************************/
/*								  	 */
/*  Find a test in the test vector; if it's not there already, add it	 */
/*								  	 */
/*************************************************************************/


Test FindTest(Newtest)
/*   ---------  */
    Test Newtest;
{
    static short TestSpace=0;
    short i;
    Boolean SameTest();

    ForEach(i, 1, NTests)
    {
	if ( SameTest(Newtest, TestVec[i]) )
	{
	    cfree(Newtest);
	    return TestVec[i];
	}
    }

    NTests++;
    if ( NTests >= TestSpace )
    {
	TestSpace += 1000;
	if ( TestSpace > 1000 )
	{
	    TestVec = (Test *) realloc(TestVec, TestSpace * sizeof(Test));
	}
	else
	{
	    TestVec = (Test *) malloc(TestSpace * sizeof(Test));
	}
    }

    TestVec[NTests] = Newtest;

    return TestVec[NTests];
}



/*************************************************************************/
/*								  	 */
/*	See if test t1 is the same test as test t2		  	 */
/*								  	 */
/*************************************************************************/


Boolean SameTest(t1, t2)
/*      ---------  */
    Test t1, t2;
{
    short i;

    if ( t1->NodeType != t2->NodeType ||
	t1->Tested != t2->Tested )
    {
	return false;
    }

    switch ( t1->NodeType )
    {
	case BrDiscr:       return true;
	case ThreshContin:  return  t1->Cut == t2->Cut;
	case BrSubset:      ForEach(i, 1, t1->Forks)
			    {
		 		if ( t1->Subset[i] != t2->Subset[i] )
				{
				    return false;
				}
			    }
    }
    return true;
}



/*************************************************************************/
/*								  	 */
/*		Clear for new set of rules			  	 */
/*								  	 */
/*************************************************************************/


    InitialiseRules()
/*  ----------------  */
{
    NRules = 0;
    Rule = 0;
    RuleSpace = 0;
}



/*************************************************************************/
/*								  	 */
/*  Add a new rule to the current ruleset, by updating Rule[],	  	 */
/*  NRules and, if necessary, RuleSpace				  	 */
/*								  	 */
/*************************************************************************/


Boolean NewRule(Cond, NConds, TargetClass, Err)
/*       -------  */
    Condition Cond[];
    short NConds;
    ClassNo TargetClass;
    float Err;
{
    short d, r;
    Boolean SameRule();

    /*  See if rule already exists  */

    ForEach(r, 1, NRules)
    {
	if ( SameRule(r, Cond, NConds, TargetClass) )
	{
	    Verbosity(1) printf("\tduplicates rule %d\n", r);

	    /*  Keep the most pessimistic error estimate  */

	    if ( Err > Rule[r].Error )
	    {
		Rule[r].Error = Err;
	    }

	    return false;
	}
    }

    /*  Make sure there is enough room for the new rule  */

    NRules++;
    if ( NRules >= RuleSpace )
    {
	RuleSpace += 100;
	if ( RuleSpace > 100 )
	{
	    Rule = (PR *) realloc(Rule, RuleSpace * sizeof(PR));
	}
	else
	{
	    Rule = (PR *) malloc(RuleSpace * sizeof(PR));
	}
    }

    /*  Form the new rule  */

    Rule[NRules].Size = NConds;
    Rule[NRules].Lhs = (Condition *) calloc(NConds+1, sizeof(Condition));
    ForEach(d, 1, NConds)
    {
	Rule[NRules].Lhs[d] = (Condition) malloc(sizeof(struct CondRec));

	Rule[NRules].Lhs[d]->CondTest = Cond[d]->CondTest;
	Rule[NRules].Lhs[d]->TestValue = Cond[d]->TestValue;
    }
    Rule[NRules].Rhs = TargetClass;
    Rule[NRules].Error = Err;

    Verbosity(1) PrintRule(NRules);

    return true;
}



/*************************************************************************/
/*								  	 */
/*  Decide whether the given rule duplicates rule r		  	 */
/*								  	 */
/*************************************************************************/


Boolean SameRule(r, Cond, NConds, TargetClass)
/*      --------  */
    RuleNo r;
    Condition Cond[];
    short NConds;
    ClassNo TargetClass;
{
    short d, i;
    Test SubTest1, SubTest2;

    if ( Rule[r].Size != NConds || Rule[r].Rhs != TargetClass )
    {
	return false;
    }

    ForEach(d, 1, NConds)
    {
	if ( Rule[r].Lhs[d]->CondTest->NodeType != Cond[d]->CondTest->NodeType ||
	     Rule[r].Lhs[d]->CondTest->Tested   != Cond[d]->CondTest->Tested )
	{
	    return false;
	}

	switch ( Cond[d]->CondTest->NodeType )
	{
	    case BrDiscr:
		if ( Rule[r].Lhs[d]->TestValue != Cond[d]->TestValue )
		{
		    return false;
		}
		break;

	    case ThreshContin:
		if ( Rule[r].Lhs[d]->CondTest->Cut != Cond[d]->CondTest->Cut )
		{
		    return false;
		}
		break;

	    case BrSubset:
		SubTest1 = Rule[r].Lhs[d]->CondTest;
		SubTest2 = Cond[d]->CondTest;
		ForEach(i, 1, SubTest1->Forks)
		{
		    if ( SubTest1->Subset[i] != SubTest2->Subset[i] )
		    {
			return false;
		    }
		}
	}
    }

    return true;
}



/*************************************************************************/
/*								  	 */
/*		Print the current indexed ruleset		  	 */
/*								  	 */
/*************************************************************************/


    PrintIndexedRules()
/*  -----------------  */
{
    short ri;

    ForEach(ri, 1, NRules )
    {
	PrintRule(RuleIndex[ri]);
    }
    printf("\nDefault class: %s\n", ClassName[DefaultClass]);
}



/*************************************************************************/
/*								  	 */
/*		Print the rule r				  	 */
/*								  	 */
/*************************************************************************/


    PrintRule(r)
/*  ---------  */
    RuleNo r;
{
    short d;

    printf("\nRule %d:\n", r);
    ForEach(d, 1, Rule[r].Size)
    {
        printf("    ");
        PrintCondition(Rule[r].Lhs[d]);
    }
    printf("\t->  class %s  [%.1f%%]\n",
	    ClassName[Rule[r].Rhs], 100 * (1 - Rule[r].Error));
}



/*************************************************************************/
/*								  	 */
/*	Print a condition c of a production rule		  	 */
/*								  	 */
/*************************************************************************/


    PrintCondition(c)
/*  --------------  */
    Condition c;
{
    Test tp;
    DiscrValue v, pv, Last, Values=0;
    Boolean First=true;
    Attribute Att;

    tp = c->CondTest;
    v = c->TestValue;
    Att = tp->Tested;

    printf("\t%s", AttName[Att]);

    if ( v < 0 )
    {
	printf(" is unknown\n");
	return;
    }

    switch ( tp->NodeType )
    {
	case BrDiscr:
	    printf(" = %s\n", AttValName[Att][v]);
	    break;

	case ThreshContin:
	    printf(" %s %g\n", ( v == 1 ? "<=" : ">" ), tp->Cut);
	    break;

	case BrSubset:
	    /*  Count values at this branch  */

	    for ( pv=1 ; Values <= 1 && pv <= MaxAttVal[Att] ; pv++ )
	    {
		if ( In(pv, tp->Subset[v]) )
		{
		    Last = pv;
		    Values++;
		}
	    }

	    if ( Values == 1 )
	    {
		printf(" = %s\n", AttValName[Att][Last]);
		break;
	    }

	    printf(" in ");
	    ForEach(pv, 1, MaxAttVal[Att])
	    {
		if ( In(pv, tp->Subset[v]) )
		{
		    if ( First )
		    {
			printf("{");
			First = false;
		    }
		    else
		    {
			printf(", ");
		    }
		    printf("%s", AttValName[Att][pv]);
		}
	    }
	    printf("}\n");
    }
}
