/*************************************************************************/
/*									 */
/*	Global data for constructing and applying rules			 */
/*      -----------------------------------------------			 */
/*									 */
/*************************************************************************/


extern PR	*Rule;		/* production rules */

extern RuleNo	NRules,		/* number of production rules */
		*RuleIndex;	/* index to production rules */

extern short	RuleSpace;	/* space currently allocated for rules */

extern RuleSet	*PRSet;		/* set of rulesets */

extern ClassNo	DefaultClass;	/* default class associated with ruleset */

extern Boolean	SIGTEST,	/* use Fisher's test in rule pruning */
		SIMANNEAL;	/* use simulated annealing */

extern float	SIGTHRESH,	/* sig level used in rule pruning */
		REDUNDANCY,	/* factor governing encoding tradeoff
				   between rules and exceptions */
		AttTestBits,	/* average bits to encode tested attribute */
		*BranchBits;	/* ditto attribute value */

extern float	*LogItemNo;	/* LogItemNo[i] = log2(i) */
extern double	*LogFact;	/* LogFact[i] = log2(i!) */
