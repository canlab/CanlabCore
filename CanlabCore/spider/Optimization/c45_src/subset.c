/*************************************************************************/
/*									 */
/*      Evaluation of the subsetting of a discrete attribute		 */
/*      ----------------------------------------------------		 */
/*									 */
/*************************************************************************/


#include "buildex.i"


ItemCount
	*Slice1,	/* Slice1[c]    = saved values of Freq[x][c] in subset.c */
	*Slice2;	/* Slice2[c]    = saved values of Freq[y][c] */

Set
	**Subset;	/* Subset[a][s] = subset s for att a */

short
	*Subsets;	/* Subsets[a] = no. subsets for att a */



/*************************************************************************/
/*									 */
/*  Evaluate subsetting a discrete attribute and form the chosen	 */
/*  subsets Subset[Att][], setting Subsets[Att] to the number of	 */
/*  subsets, and the Info[] and Gain[] of a test on the attribute	 */
/*									 */
/*************************************************************************/


    EvalSubset(Att, Fp, Lp, Items)
/*  ----------  */ 
    Attribute Att;
    ItemNo Fp, Lp; 
    ItemCount Items;
{ 
    DiscrValue V1, V2, BestV1, BestV2, Barred;
    ItemCount KnownItems;
    ClassNo c;
    float BaseInfo, MinGain, ThisGain, ThisInfo,
	Val, BestVal, BestGain, BestInfo,
	PrevVal, PrevGain, PrevInfo,
	DiscrKnownBaseInfo(), Worth(), ComputeGain(), TotalInfo();
    short Blocks=0, MissingValues=0, ReasonableSubsets, Bytes, b;
    Boolean MergedSubsets = false;
    int SaveMINOBJS;

    SaveMINOBJS = MINOBJS;
    MINOBJS = 1;

    /*  First compute Freq[][], ValFreq[], base info, and the gain
	and total info of a split on discrete attribute Att  */

    ComputeFrequencies(Att, Fp, Lp);

    KnownItems = Items - ValFreq[0];
    if ( KnownItems < Epsilon )
    {
	Verbosity(2) printf("\tAtt %s: no known values\n", AttName[Att]);

	Gain[Att] = -Epsilon;
	Info[Att] = 0;
	return;
    }

    BaseInfo = DiscrKnownBaseInfo(KnownItems, MaxAttVal[Att]);

    PrevGain = ComputeGain(BaseInfo, UnknownRate[Att], MaxAttVal[Att],KnownItems);
    PrevInfo = TotalInfo(ValFreq, 0, MaxAttVal[Att]) / Items;
    PrevVal = Worth(PrevInfo, PrevGain, Epsilon);

    Verbosity(2)
    {
	printf("\tAtt %s", AttName[Att]);

	Verbosity(3) PrintDistribution(Att, MaxAttVal[Att], true);

	printf("\tinf %.3f, gain %.3f, val=%.3f\n",
		PrevInfo, PrevGain, PrevVal);
    }

    /*  Eliminate unrepresented attribute values from Freq[] and ValFreq[]
	and form a separate subset for each represented attribute value  */

    Bytes = (MaxAttVal[Att]>>3) + 1;
    ClearBits(Bytes, Subset[Att][0]);

    ForEach(V1, 1, MaxAttVal[Att])
    {
	if ( ValFreq[V1] > 0.5 )
	{
	    if ( ++Blocks < V1 )
	    {
		ValFreq[Blocks] = ValFreq[V1];
		ForEach(c, 0, MaxClass)
		{
		    Freq[Blocks][c] = Freq[V1][c];
		}
	    }
	    ClearBits(Bytes, Subset[Att][Blocks]);
	    SetBit(V1, Subset[Att][Blocks]);
	}
	else
	{
	    SetBit(V1, Subset[Att][0]);
	    MissingValues++;
	}
    }

    /*  Merge any single-class subsets with others of the same class  */
    /*  Note: have ValFreq[V] > 0 for all V  */

    ForEach(V1, 1, Blocks-1)
    {
	for ( c = 0 ; Freq[V1][c] < 0.1 ; c++ )
	    ;

	if ( Freq[V1][c] < ValFreq[V1] - 0.1 ) continue;

	/*  Now have a single class -- look for others  */

	for ( V2 = V1+1 ; V2 <= Blocks ; )
	{
	    if ( Freq[V2][c] < ValFreq[V2] - 0.1 )
	    {
		V2++;
	    }
	    else
	    {
		/*  Merge these subsets  */

		Combine(V1, V2, Blocks);

		ForEach(b, 0, Bytes-1)
		{
		    Subset[Att][V1][b] |= Subset[Att][V2][b];
		    Subset[Att][V2][b] = Subset[Att][Blocks][b];
		}

		Blocks--;
		MergedSubsets = true;
	    }
	}
    }

    if ( MergedSubsets )
    {
	PrevGain = ComputeGain(BaseInfo, UnknownRate[Att], Blocks, KnownItems);
	PrevInfo = TotalInfo(ValFreq, 0, Blocks) / Items;
	PrevVal = Worth(PrevInfo, PrevGain, Epsilon);

	Verbosity(2)
	{
	    printf("\tAfter merging single-class subsets:");

	    Verbosity(3) PrintDistribution(Att, Blocks, false);

	    printf("\tinf %.3f, gain %.3f, val=%.3f\n",
		    PrevInfo, PrevGain, PrevVal);
	}
    }

    /*  Examine possible pair mergers and hill-climb  */

    MinGain = PrevGain / 2;

    while ( Blocks > 2 )
    {
	BestVal = BestV1 = 0;
	BestGain = -Epsilon;

	/*  Check reasonable subsets; if less than 3, bar mergers
	    involving the largest block  */

	ReasonableSubsets = 0;
	Barred = 1;

	ForEach(V1, 1, Blocks)
	{
	    if ( ValFreq[V1] >= SaveMINOBJS ) ReasonableSubsets++;

	    if ( ValFreq[V1] > ValFreq[Barred] ) Barred = V1;
	}

	if ( ReasonableSubsets >= 3 ) Barred = 0;

	/*  For each possible pair of values, calculate the gain and
	    total info of a split in which they are treated as one.
	    Keep track of the pair with the best gain.  */

	ForEach(V1, 1, Blocks-1)
	{
	    ForEach(V2, V1+1, Blocks)
	    {
		if ( V1 == Barred || V2 == Barred ) continue;

		Combine(V1, V2, Blocks);

		ThisGain = ComputeGain(BaseInfo, UnknownRate[Att],
					Blocks-1, KnownItems);
		ThisInfo = TotalInfo(ValFreq, 0, Blocks-1) / Items;
		Val      = Worth(ThisInfo, ThisGain, Epsilon);

		Verbosity(4)
		{
		    printf("\tcombine %d %d info %.3f gain %.3f val %.3f",
		           V1, V2, ThisInfo, ThisGain, Val);
		    PrintDistribution(Att, Blocks-1, false);
		}

		/*  Force a split if
			less than two reasonable subsets, or
			using GAIN criterion
		    Prefer this split to the previous one if
			gain >= MinGain (and previous < MinGain), or
			val >= previous best val  */

		if ( ThisGain >= MinGain && BestGain < MinGain ||
		     Val >= BestVal ||
		     ! BestV1 && ( ! GAINRATIO || ReasonableSubsets < 2 ) )
		{
		    BestVal  = Val;
		    BestGain = ThisGain;
		    BestInfo = ThisInfo;
		    BestV1   = V1;
		    BestV2   = V2;
		}

		Uncombine(V1, V2);
	    }
	}

	if ( GAINRATIO &&
	     ReasonableSubsets >= 2 &&
	     ( ! BestV1 ||
	       BestVal < PrevVal + 1E-5 ||
	       BestVal == PrevVal && BestGain < PrevGain ) ) break;

	PrevGain = BestGain;
	PrevInfo = BestInfo;
        PrevVal = BestVal;

	Combine(BestV1, BestV2, Blocks);

	ForEach(b, 0, Bytes-1)
	{
	    Subset[Att][BestV1][b] |= Subset[Att][BestV2][b];
	    Subset[Att][BestV2][b] = Subset[Att][Blocks][b];
	}

	Blocks--;

	Verbosity(2)
	{
	    printf("\t\tform subset ");
	    PrintSubset(Att, Subset[Att][BestV1]);
	    printf(": %d subsets, inf %.3f, gain %.3f, val %.3f\n",
		   Blocks, BestInfo, BestGain, BestVal);
	    Verbosity(3)
	    {
		printf("\t\tcombine %d, %d", BestV1, BestV2);
		PrintDistribution(Att, Blocks, false);
	    }
	}
    }

    MINOBJS = SaveMINOBJS;

    if ( PrevVal <= 0 )
    {
	Gain[Att] = -Epsilon;
	Info[Att] = 0;
    }
    else
    {
	Gain[Att] = ComputeGain(BaseInfo, UnknownRate[Att], Blocks, KnownItems);
	Info[Att] = PrevInfo;

	if ( MissingValues )
	{
	    Blocks++;
	    CopyBits(Bytes, Subset[Att][0], Subset[Att][Blocks]);
	}

	Subsets[Att] = Blocks;

	Verbosity(2) printf("\tFinal subsets:");
	Verbosity(3) PrintDistribution(Att, Blocks, false);
	Verbosity(2)
	    printf("\tinf %.3f gain %.3f val %.3f\n", 
		   Info[Att], Gain[Att], Worth(Info[Att], Gain[Att], Epsilon));
    }
}



/*************************************************************************/
/*									 */
/*  Combine the distribution figures of discrete attribute values	 */
/*  x and y, putting the combined figures in Freq[x][] and		 */
/*  ValFreq[x][], and saving old values in Slice1 and Slice2		 */
/*									 */
/*************************************************************************/


    Combine(x, y, Last)
/*  -------  */
    DiscrValue x, y, Last;
{
    ClassNo c;

    ForEach(c, 0, MaxClass)
    {
	Slice1[c] = Freq[x][c];
	Slice2[c] = Freq[y][c];

	Freq[x][c] += Freq[y][c];
	Freq[y][c]  = Freq[Last][c];
    }

    Slice1[MaxClass+1] = ValFreq[x];
    Slice2[MaxClass+1] = ValFreq[y];

    ValFreq[x] += ValFreq[y];
    ValFreq[y]  = ValFreq[Last];
}



/*************************************************************************/
/*									 */
/*  Restore old class distribution figures of discrete attribute	 */
/*  values x and y from Slice1 and Slice2				 */
/*									 */
/*************************************************************************/


    Uncombine(x, y)
/*  ---------  */
    DiscrValue x, y;
{
    ClassNo c;

    ForEach(c, 0, MaxClass)
    {
	Freq[x][c] = Slice1[c];
	Freq[y][c] = Slice2[c];
    }

    ValFreq[x] = Slice1[MaxClass+1];
    ValFreq[y] = Slice2[MaxClass+1];
}



/*************************************************************************/
/*									 */
/*  Print the values of attribute Att which are in the subset Ss	 */
/*									 */
/*************************************************************************/


    PrintSubset(Att, Ss)
/*  -----------  */
    Attribute Att;
    Set Ss;
{
    DiscrValue V1;
    Boolean First=true;

    ForEach(V1, 1, MaxAttVal[Att])
    {
	if ( In(V1, Ss) )
	{
	    if ( First )
	    {
		First = false;
	    }
	    else
	    {
		printf(", ");
	    }

	    printf("%s", AttValName[Att][V1]);
	}
    }
}



/*************************************************************************/
/*									 */
/*  Construct and return a node for a test on a subset of values	 */
/*									 */
/*************************************************************************/


    SubsetTest(Node, Att)
/*  -----------  */
    Tree Node;
    Attribute Att;
{ 
    ItemCount CountItems();
    short S, Bytes;

    Sprout(Node, Subsets[Att]);

    Node->NodeType	= BrSubset;
    Node->Tested	= Att;
    Node->Errors	= 0;
    
    Bytes = (MaxAttVal[Att]>>3) + 1;
    Node->Subset = (Set *) calloc(Subsets[Att] + 1, sizeof(Set));
    ForEach(S, 1, Node->Forks)
    {
	Node->Subset[S] = (Set) malloc(Bytes);
	CopyBits(Bytes, Subset[Att][S], Node->Subset[S]);
    }
} 
