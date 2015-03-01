/*************************************************************************/
/*									 */
/*		Type definitions for C4.5				 */
/*              -------------------------				 */
/*									 */
/*************************************************************************/


typedef  char	Boolean, *String, *Set;

typedef  int	ItemNo;		/* data item number */
typedef  float	ItemCount;	/* count of (partial) items */

typedef  short	ClassNo,	/* class number, 0..MaxClass */
		DiscrValue;	/* discrete attribute value (0 = ?) */
typedef  short	Attribute;	/* attribute number, 0..MaxAtt */

typedef  union  _attribute_value
	 {
	    DiscrValue	_discr_val;
	    float	_cont_val;
	 }
	 	AttValue, *Description;

#define  CVal(Case,Attribute)   Case[Attribute]._cont_val
#define  DVal(Case,Attribute)   Case[Attribute]._discr_val
#define  Class(Case)		Case[MaxAtt+1]._discr_val

#define  Unknown  -999		/* unknown value for continuous attribute */


#define  BrDiscr	1	/* node types:	branch */
#define  ThreshContin	2	/*		threshold cut */
#define  BrSubset	3	/*		subset test */

typedef  struct _tree_record *Tree;
typedef  struct _tree_record
	 {
	    short	NodeType;	/* 0=leaf 1=branch 2=cut 3=subset */
	    ClassNo	Leaf;		/* most frequent class at this node */
	    ItemCount	Items,		/* no of items at this node */
			*ClassDist,	/* class distribution of items */
	    		Errors;		/* no of errors at this node */
	    Attribute	Tested; 	/* attribute referenced in test */
	    short	Forks;		/* number of branches at this node */
	    float	Cut,		/* threshold for continuous attribute */
		  	Lower,		/* lower limit of soft threshold */
		  	Upper;		/* upper limit ditto */
	    Set         *Subset;	/* subsets of discrete values  */
	    Tree	*Branch;	/* Branch[x] = (sub)tree for outcome x */
	 }
		TreeRec;

#define  IGNORE		1	/* special attribute status: do not use */
#define  DISCRETE	2	/* ditto: collect values as data read */


typedef short	RuleNo;			/* rule number */


typedef struct TestRec *Test;

struct TestRec
       {
	   short	NodeType;	/* test type (see tree nodes) */
	   Attribute	Tested;		/* attribute tested */
	   short	Forks;		/* possible branches */
	   float	Cut;		/* threshold (if relevant) */
	   Set		*Subset;	/* subset (if relevant) */
       };


typedef struct CondRec *Condition;

struct CondRec
       {
	   Test		CondTest;	/* test part of condition */
	   short	TestValue;	/* specified outcome of test */
       };


typedef struct ProdRuleRec PR;

struct ProdRuleRec
       {
	   short	Size;		/* number of conditions */
	   Condition	*Lhs;		/* conditions themselves */
	   ClassNo	Rhs;		/* class given by rule */
	   float	Error,		/* estimated error rate */
			Bits;		/* bits to encode rule */
	   ItemNo	Used,		/* times rule used */
			Incorrect;	/* times rule incorrect */
       };


typedef struct RuleSetRec RuleSet;

struct RuleSetRec
       {
	   PR		*SRule;		/* rules */
	   RuleNo	SNRules,	/* number of rules */
			*SRuleIndex;	/* ranking of rules */
	   ClassNo	SDefaultClass;	/* default class for this ruleset */
    };
