/*************************************************************************/
/*								 	 */
/*    Central tree-forming algorithm incorporating all criteria  	 */
/*    ---------------------------------------------------------	 	 */
/*								 	 */
/*************************************************************************/


#include "defns.i"
#include "types.i"
#include "extern.i"


ItemCount
	*Weight,	/* Weight[i]  = current fraction of item i */
	**Freq,		/* Freq[x][c] = no. items of class c with outcome x */
	*ValFreq,	/* ValFreq[x]   = no. items with outcome x */
	*ClassFreq;	/* ClassFreq[c] = no. items of class c */

float
	*Gain,		/* Gain[a] = info gain by split on att a */
	*Info,		/* Info[a] = potential info of split on att a */
	*Bar,		/* Bar[a]  = best threshold for contin att a */
	*UnknownRate;	/* UnknownRate[a] = current unknown rate for att a */

Boolean
	*Tested,	/* Tested[a] set if att a has already been tested */
	MultiVal;	/* true when all atts have many values */


	/*  External variables initialised here  */

extern float
	*SplitGain,	/* SplitGain[i] = gain with att value of item i as threshold */
	*SplitInfo;	/* SplitInfo[i] = potential info ditto */

extern ItemCount
	*Slice1,	/* Slice1[c]    = saved values of Freq[x][c] in subset.c */
	*Slice2;	/* Slice2[c]    = saved values of Freq[y][c] */

extern Set
	**Subset;	/* Subset[a][s] = subset s for att a */

extern short
	*Subsets;	/* Subsets[a] = no. subsets for att a */



/*************************************************************************/
/*								 	 */
/*		Allocate space for tree tables			 	 */
/*								 	 */
/*************************************************************************/


    InitialiseTreeData()
/*  ------------------  */
{ 
    DiscrValue v;
    Attribute a;

    Tested	= (char *) calloc(MaxAtt+1, sizeof(char));

    Gain	= (float *) calloc(MaxAtt+1, sizeof(float));
    Info	= (float *) calloc(MaxAtt+1, sizeof(float));
    Bar		= (float *) calloc(MaxAtt+1, sizeof(float));

    Subset = (Set **) calloc(MaxAtt+1, sizeof(Set *));
    ForEach(a, 0, MaxAtt)
    {
	if ( MaxAttVal[a] )
	{
	    Subset[a]  = (Set *) calloc(MaxDiscrVal+1, sizeof(Set));
	    ForEach(v, 0, MaxAttVal[a])
	    {
		Subset[a][v] = (Set) malloc((MaxAttVal[a]>>3) + 1);
	    }
	}
    }
    Subsets = (short *) calloc(MaxAtt+1, sizeof(short));

    SplitGain = (float *) calloc(MaxItem+1, sizeof(float));
    SplitInfo = (float *) calloc(MaxItem+1, sizeof(float));

    Weight = (ItemCount *) calloc(MaxItem+1, sizeof(ItemCount));

    Freq  = (ItemCount **) calloc(MaxDiscrVal+1, sizeof(ItemCount *));
    ForEach(v, 0, MaxDiscrVal)
    {
	Freq[v]  = (ItemCount *) calloc(MaxClass+1, sizeof(ItemCount));
    }

    ValFreq = (ItemCount *) calloc(MaxDiscrVal+1, sizeof(ItemCount));
    ClassFreq = (ItemCount *) calloc(MaxClass+1, sizeof(ItemCount));

    Slice1 = (ItemCount *) calloc(MaxClass+2, sizeof(ItemCount));
    Slice2 = (ItemCount *) calloc(MaxClass+2, sizeof(ItemCount));

    UnknownRate = (float *) calloc(MaxAtt+1, sizeof(float));

    /*  Check whether all attributes have many discrete values  */

    MultiVal = true;
    if ( ! SUBSET )
    {
	for ( a = 0 ; MultiVal && a <= MaxAtt ; a++ )
	{
	    if ( SpecialStatus[a] != IGNORE )
	    {
		MultiVal = MaxAttVal[a] >= 0.3 * (MaxItem + 1);
	    }
	}
    }
}



/*************************************************************************/
/*								 	 */
/*		Initialise the weight of each item		 	 */
/*								 	 */
/*************************************************************************/


    InitialiseWeights()
/*  -----------------  */
{
    ItemNo i;

    ForEach(i, 0, MaxItem)
    {
        Weight[i] = 1.0;
    }
}



/*************************************************************************/
/*								 	 */
/*  Build a decision tree for the cases Fp through Lp:		 	 */
/*								 	 */
/*  - if all cases are of the same class, the tree is a leaf and so	 */
/*      the leaf is returned labelled with this class		 	 */
/*								 	 */
/*  - for each attribute, calculate the potential information provided 	 */
/*	by a test on the attribute (based on the probabilities of each	 */
/*	case having a particular value for the attribute), and the gain	 */
/*	in information that would result from a test on the attribute	 */
/*	(based on the probabilities of each case with a particular	 */
/*	value for the attribute being of a particular class)		 */
/*								 	 */
/*  - on the basis of these figures, and depending on the current	 */
/*	selection criterion, find the best attribute to branch on. 	 */
/*	Note:  this version will not allow a split on an attribute	 */
/*	unless two or more subsets have at least MINOBJS items. 	 */
/*								 	 */
/*  - try branching and test whether better than forming a leaf	 	 */
/*								 	 */
/*************************************************************************/


Tree FormTree(Fp, Lp)
/*   ---------  */
    ItemNo Fp, Lp; 
{ 
    ItemNo i, Kp, Ep, Group();
    ItemCount Cases, NoBestClass, KnownCases, CountItems();
    float Factor, BestVal, Val, AvGain=0, Worth();
    Attribute Att, BestAtt, Possible=0;
    ClassNo c, BestClass;
    Tree Node, Leaf();
    DiscrValue v;
    Boolean PrevAllKnown;

    Cases = CountItems(Fp, Lp);

    /*  Generate the class frequency distribution  */

    ForEach(c, 0, MaxClass)
    {
	ClassFreq[c] = 0;
    }
    ForEach(i, Fp, Lp)
    { 
	ClassFreq[ Class(Item[i]) ] += Weight[i];
    } 

    /*  Find the most frequent class  */

    BestClass = 0;
    ForEach(c, 0, MaxClass)
    {
	if ( ClassFreq[c] > ClassFreq[BestClass] )
	{
	    BestClass = c;
	}
    }
    NoBestClass = ClassFreq[BestClass];

    Node = Leaf(ClassFreq, BestClass, Cases, Cases - NoBestClass);

    /*  If all cases are of the same class or there are not enough
	cases to divide, the tree is a leaf  */

    if ( NoBestClass == Cases  || Cases < 2 * MINOBJS )
    { 
	return Node;
    } 

    Verbosity(1)
    	printf("\n%d items, total weight %.1f\n", Lp - Fp + 1, Cases);

    /*  For each available attribute, find the information and gain  */

    ForEach(Att, 0, MaxAtt) 
    { 
	Gain[Att] = -Epsilon;

	if ( SpecialStatus[Att] == IGNORE ) continue;

        if ( MaxAttVal[Att] )
	{
	    /*  discrete valued attribute  */

	    if ( SUBSET && MaxAttVal[Att] > 2 )
	    {
		EvalSubset(Att, Fp, Lp, Cases);
	    }
	    else
	    if ( ! Tested[Att] )
	    {
	        EvalDiscreteAtt(Att, Fp, Lp, Cases);
	    }
	}
	else
	{ 
	    /*  continuous attribute  */

	    EvalContinuousAtt(Att, Fp, Lp);
	} 

	/*  Update average gain, excluding attributes with very many values  */

	if ( Gain[Att] > -Epsilon &&
	     ( MultiVal || MaxAttVal[Att] < 0.3 * (MaxItem + 1) ) )
	{
	    Possible++;
	    AvGain += Gain[Att];
	}
    } 

    /*  Find the best attribute according to the given criterion  */

    BestVal = -Epsilon;
    BestAtt = None;
    AvGain  = ( Possible ? AvGain / Possible : 1E6 );

    Verbosity(2)
    {
	if ( AvGain < 1E6 ) printf("\taverage gain %.3f\n", AvGain);
    }

    ForEach(Att, 0, MaxAtt) 
    { 
	if ( Gain[Att] > -Epsilon )
	{ 
	    Val = Worth(Info[Att], Gain[Att], AvGain);
	    if ( Val > BestVal ) 
	    { 
	        BestAtt  = Att; 
	        BestVal = Val;
	    } 
	} 
    } 

    /*  Decide whether to branch or not  */ 

    if ( BestAtt != None )
    { 
	Verbosity(1)
	{
	    printf("\tbest attribute %s", AttName[BestAtt]);
	    if ( ! MaxAttVal[BestAtt] )
	    {
		printf(" cut %.3f", Bar[BestAtt]);
	    }
	    printf(" inf %.3f gain %.3f val %.3f\n",
		   Info[BestAtt], Gain[BestAtt], BestVal);
	}	

	/*  Build a node of the selected test  */

        if ( MaxAttVal[BestAtt] )
	{
	    /*  Discrete valued attribute  */

	    if ( SUBSET && MaxAttVal[BestAtt] > 2 )
	    {
	        SubsetTest(Node, BestAtt);
	    }
	    else
	    {
	        DiscreteTest(Node, BestAtt);
	    }
	}
	else
	{ 
	    /*  Continuous attribute  */

	    ContinTest(Node, BestAtt);
	} 

	/*  Remove unknown attribute values  */

	PrevAllKnown = AllKnown;

	Kp = Group(0, Fp, Lp, Node) + 1;
	if ( Kp != Fp ) AllKnown = false;
	KnownCases = Cases - CountItems(Fp, Kp-1);
	UnknownRate[BestAtt] = (Cases - KnownCases) / (Cases + 0.001);

	Verbosity(1)
	{
	    if ( UnknownRate[BestAtt] > 0 )
	    {
		printf("\tunknown rate for %s = %.3f\n",
		       AttName[BestAtt], UnknownRate[BestAtt]);
	    }
	}

	/*  Recursive divide and conquer  */

	++Tested[BestAtt];

	Ep = Kp - 1;
	Node->Errors = 0;

	ForEach(v, 1, Node->Forks)
	{
	    Ep = Group(v, Kp, Lp, Node);

	    if ( Kp <= Ep )
	    {
		Factor = CountItems(Kp, Ep) / KnownCases;

		ForEach(i, Fp, Kp-1)
		{
		    Weight[i] *= Factor;
		}

		Node->Branch[v] = FormTree(Fp, Ep);
		Node->Errors += Node->Branch[v]->Errors;

		Group(0, Fp, Ep, Node);
		ForEach(i, Fp, Kp-1)
		{
		    Weight[i] /= Factor;
		}
	    }
	    else
	    {
		Node->Branch[v] = Leaf(Node->ClassDist, BestClass, 0.0, 0.0);
	    }
	}

	--Tested[BestAtt];
	AllKnown = PrevAllKnown;

	/*  See whether we would have been no worse off with a leaf  */

	if ( Node->Errors >= Cases - NoBestClass - Epsilon )
	{ 
	    Verbosity(1)
		printf("Collapse tree for %d items to leaf %s\n",
	                Lp - Fp + 1, ClassName[BestClass]);

	    Node->NodeType = 0;
	} 
    }
    else
    { 
	Verbosity(1)
	    printf("\tno sensible splits  %.1f/%.1f\n",
		   Cases, Cases - NoBestClass);
    } 

    return Node; 
} 



/*************************************************************************/
/*								 	 */
/*  Group together the items corresponding to branch V of a test 	 */
/*  and return the index of the last such			 	 */
/*								 	 */
/*  Note: if V equals zero, group the unknown values		 	 */
/*								 	 */
/*************************************************************************/


ItemNo Group(V, Fp, Lp, TestNode)
/*     -----  */
    DiscrValue V;
    ItemNo Fp, Lp;
    Tree TestNode;
{
    ItemNo i;
    Attribute Att;
    float Thresh;
    Set SS;
    void Swap();

    Att = TestNode->Tested;

    if ( V )
    {
	/*  Group items on the value of attribute Att, and depending
	    on the type of branch  */

	switch ( TestNode->NodeType )
	{
	    case BrDiscr:

		ForEach(i, Fp, Lp)
		{
		    if ( DVal(Item[i], Att) == V ) Swap(Fp++, i);
		}
		break;

	    case ThreshContin:

		Thresh = TestNode->Cut;
		ForEach(i, Fp, Lp)
		{
		    if ( (CVal(Item[i], Att) <= Thresh) == (V == 1) ) Swap(Fp++, i);
		}
		break;

	    case BrSubset:

		SS = TestNode->Subset[V];
		ForEach(i, Fp, Lp)
		{
		    if ( In(DVal(Item[i], Att), SS) ) Swap(Fp++, i);
		}
		break;
	}
    }
    else
    {
	/*  Group together unknown values  */

	switch ( TestNode->NodeType )
	{
	    case BrDiscr:
	    case BrSubset:

		ForEach(i, Fp, Lp)
		{
		    if ( ! DVal(Item[i], Att) ) Swap(Fp++, i);
		}
		break;

	    case ThreshContin:

		ForEach(i, Fp, Lp)
		{
		    if ( CVal(Item[i], Att) == Unknown ) Swap(Fp++, i);
		}
		break;
	}
    }

    return Fp - 1;
}



/*************************************************************************/
/*								 	 */
/*	Return the total weight of items from Fp to Lp		 	 */
/*								 	 */
/*************************************************************************/


ItemCount CountItems(Fp, Lp)
/*        ----------  */
    ItemNo Fp, Lp;
{
    register ItemCount Sum=0.0, *Wt, *LWt;

    if ( AllKnown ) return Lp - Fp + 1;

    for ( Wt = Weight + Fp, LWt = Weight + Lp ; Wt <= LWt ; )
    {
	Sum += *Wt++;
    }

    return Sum;
}



/*************************************************************************/
/*                                                               	 */
/*		Exchange items at a and b			 	 */
/*									 */
/*************************************************************************/


void Swap(a,b)
/*   ----  */
    ItemNo a, b;
{
    register Description Hold;
    register ItemCount HoldW;

    Hold = Item[a];
    Item[a] = Item[b];
    Item[b] = Hold;

    HoldW = Weight[a];
    Weight[a] = Weight[b];
    Weight[b] = HoldW;
}
