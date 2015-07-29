/*************************************************************************/
/*									 */
/*		Global data for C4.5					 */
/*		--------------------					 */
/*									 */
/*************************************************************************/


extern  short		MaxAtt,		/* max att number */
			MaxClass,	/* max class number */
			MaxDiscrVal;	/* max discrete values for any att */

extern  ItemNo		MaxItem;	/* max data item number */

extern  Description	*Item;		/* data items */

extern  DiscrValue	*MaxAttVal;	/* number of values for each att */

extern  char		*SpecialStatus;	/* special att treatment */

extern  String		*ClassName,	/* class names */
		  	*AttName,	/* att names */
		  	**AttValName,	/* att value names */
			FileName;	/* family name of files */

extern  Boolean		AllKnown;	/* true if there have been no splits
					   on atts with missing values above
					   the current position in the tree */


/*************************************************************************/
/*									 */
/*		Global parameters for C4.5				 */
/*		--------------------------				 */
/*									 */
/*************************************************************************/


extern  short		VERBOSITY,	/* verbosity level (0 = none) */
			TRIALS;		/* number of trees to be grown */

extern  Boolean		GAINRATIO,	/* true=gain ratio, false=gain */
			SUBSET,		/* true if subset tests allowed */
			BATCH,		/* true if windowing turned off */
			UNSEENS,	/* true if to evaluate on test data */
			PROBTHRESH;	/* true if to use soft thresholds */

extern  ItemNo		MINOBJS,	/* minimum items each side of a cut */
			WINDOW,		/* initial window size */
			INCREMENT;	/* max window increment each iteration */

extern  float		CF;		/* confidence limit for tree pruning */
