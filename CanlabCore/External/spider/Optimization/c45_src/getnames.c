/*************************************************************************/
/*									 */
/*	Get names of classes, attributes and attribute values		 */
/*	-----------------------------------------------------		 */
/*									 */
/*************************************************************************/


#include "defns.i"
#include "types.i"
#include "extern.i"

#include <sys/types.h>
#include <sys/stat.h>


#define  Space(s)	(s == ' ' || s == '\n' || s == '\t')
#define  SkipComment	while ( ( c = getc(f) ) != '\n' )

char	Delimiter;
String	CopyString();



/*************************************************************************/
/*									 */
/*  Read a name from file f into string s, setting Delimiter.		 */
/*									 */
/*  - Embedded periods are permitted, but periods followed by space	 */
/*    characters act as delimiters.					 */
/*  - Embedded spaces are permitted, but multiple spaces are replaced	 */
/*    by a single space.						 */
/*  - Any character can be escaped by '\'.				 */
/*  - The remainder of a line following '|' is ignored.			 */
/*									 */
/*************************************************************************/


Boolean ReadName(f, s)
/*      ---------  */
    FILE *f;
    String s;
{
    register char *Sp=s;
    register int c;

    /*  Skip to first non-space character  */

    while ( ( c = getc(f) ) == '|' || Space(c) )
    {
	if ( c == '|' ) SkipComment;
    }

    /*  Return false if no names to read  */

    if ( c == EOF )
    {
	Delimiter = EOF;
	return false;
    }

    /*  Read in characters up to the next delimiter  */

    while ( c != ':' && c != ',' && c != '\n' && c != '|' && c != EOF )
    {
	if ( c == '.' )
	{
	    if ( ( c = getc(f) ) == '|' || Space(c) ) break;
	    *Sp++ = '.';
	}

	if ( c == '\\' )
	{
	    c = getc(f);
	}

	*Sp++ = c;

	if ( c == ' ' )
	{
	    while ( ( c = getc(f) ) == ' ' )
		;
	}
	else
	{
	    c = getc(f);
	}
    }

    if ( c == '|' ) SkipComment;
    Delimiter = c;

    /*  Strip trailing spaces  */

    while ( Space(*(Sp-1)) ) Sp--;

    *Sp++ = '\0';
    return true;
}



/*************************************************************************/
/*									 */
/*  Read the names of classes, attributes and legal attribute values.	 */
/*  On completion, these names are stored in:				 */
/*	ClassName	-  class names					 */
/*	AttName		-  attribute names				 */
/*	AttValName	-  attribute value names			 */
/*  with:								 */
/*	MaxAttVal	-  number of values for each attribute		 */
/*									 */
/*  Other global variables set are:					 */
/*	MaxAtt		-  maximum attribute number			 */
/*	MaxClass	-  maximum class number				 */
/*	MaxDiscrVal	-  maximum discrete values for any attribute	 */
/*									 */
/*  Note:  until the number of attributes is known, the name		 */
/*	   information is assembled in local arrays			 */
/*									 */
/*************************************************************************/


    GetNames()
/*  ---------  */
{
    FILE *Nf, *fopen();
    char Fn[100], Buffer[1000];
    DiscrValue v;
    int AttCeiling=100, ClassCeiling=100, ValCeiling;

    /*  Open names file  */

    strcpy(Fn, FileName);
    strcat(Fn, ".names");
    if ( ! ( Nf = fopen(Fn, "r") ) ) Error(0, Fn, "");

    /*  Get class names from names file  */

    ClassName = (String *) calloc(ClassCeiling, sizeof(String));
    MaxClass = -1;
    do
    {
	ReadName(Nf, Buffer);

	if ( ++MaxClass >= ClassCeiling)
	{
	    ClassCeiling += 100;
	    ClassName = (String *) realloc(ClassName, ClassCeiling*sizeof(String));
	}
	ClassName[MaxClass] = CopyString(Buffer);
    }
    while ( Delimiter == ',' );

    /*  Get attribute and attribute value names from names file  */

    AttName = (String *) calloc(AttCeiling, sizeof(String));
    MaxAttVal = (DiscrValue *) calloc(AttCeiling, sizeof(DiscrValue));
    AttValName = (String **) calloc(AttCeiling, sizeof(String *));
    SpecialStatus = (char *) malloc(AttCeiling);

    MaxAtt = -1;
    while ( ReadName(Nf, Buffer) )
    {
	if ( Delimiter != ':' ) Error(1, Buffer, "");

	if ( ++MaxAtt >= AttCeiling )
	{
	    AttCeiling += 100;
	    AttName = (String *) realloc(AttName, AttCeiling*sizeof(String));
	    MaxAttVal = (DiscrValue *) realloc(MaxAttVal, AttCeiling*sizeof(DiscrValue));
	    AttValName = (String **) realloc(AttValName, AttCeiling*sizeof(String *));
	    SpecialStatus = (char *) realloc(SpecialStatus, AttCeiling);
	}

	AttName[MaxAtt] = CopyString(Buffer);
	SpecialStatus[MaxAtt] = Nil;
	MaxAttVal[MaxAtt] = 0;
	ValCeiling = 100;
	AttValName[MaxAtt] = (String *) calloc(ValCeiling, sizeof(String));

	do
	{
	    if ( ! ( ReadName(Nf, Buffer) ) ) Error(2, AttName[MaxAtt], "");

	    if ( ++MaxAttVal[MaxAtt] >= ValCeiling )
	    {
		ValCeiling += 100;
		AttValName[MaxAtt] =
		    (String *) realloc(AttValName[MaxAtt], ValCeiling*sizeof(String));
	    }

	    AttValName[MaxAtt][MaxAttVal[MaxAtt]] = CopyString(Buffer);
	}
	while ( Delimiter == ',' );

	if ( MaxAttVal[MaxAtt] == 1 )
	{
	    /*  Check for special treatment  */

	    if ( ! strcmp(Buffer, "continuous") )
	    {}
	    else
	    if ( ! memcmp(Buffer, "discrete", 8) )
	    {
		SpecialStatus[MaxAtt] = DISCRETE;

		/*  Read max values, reserve space and check MaxDiscrVal  */

		v = atoi(&Buffer[8]);
		if ( v < 2 )
		{
		    printf("** %s: illegal number of discrete values\n",
			   AttName[MaxAtt]);
		    exit(1);
		}

		AttValName[MaxAtt] =
		    (String *) realloc(AttValName[MaxAtt], (v+2)*sizeof(String));
		AttValName[MaxAtt] = (char **) v;
		if ( v > MaxDiscrVal ) MaxDiscrVal = v;
	    }
	    else
	    if ( ! strcmp(Buffer, "ignore") )
	    {
		SpecialStatus[MaxAtt] = IGNORE;
	    }
	    else
	    {
		/*  Cannot have only one discrete value for an attribute  */

		Error(3, AttName[MaxAtt], "");
	    }

	    MaxAttVal[MaxAtt] = 0;
	}
	else
	if ( MaxAttVal[MaxAtt] > MaxDiscrVal ) MaxDiscrVal = MaxAttVal[MaxAtt];
    }

    fclose(Nf);
}



/*************************************************************************/
/*									 */
/*	Locate value Val in List[First] to List[Last]			 */
/*									 */
/*************************************************************************/


int Which(Val, List, First, Last)
/*  -----  */
    String Val, List[];
    short First, Last;
{
    short n=First;

    while ( n <= Last && strcmp(Val, List[n]) ) n++;

    return ( n <= Last ? n : First-1 );
}



/*************************************************************************/
/*									 */
/*	Allocate space then copy string into it				 */
/*									 */
/*************************************************************************/

String CopyString(x)
/*     -----------  */
    String x;
{
    char *s;

    s = (char *) calloc(strlen(x)+1, sizeof(char));
    strcpy(s, x);
    return s;
}



/*************************************************************************/
/*									 */
/*			Error messages					 */
/*									 */
/*************************************************************************/

    Error(n, s1, s2)
/*  -----  */
    short n;
    String s1, s2;
{
    static char Messages=0;

    printf("\nERROR:  ");
    switch(n)
    {
	case 0: printf("cannot open file %s%s\n", s1, s2);
		exit(1);

	case 1:	printf("colon expected after attribute name %s\n", s1);
		break;

	case 2:	printf("unexpected eof while reading attribute %s\n", s1);
		break;

	case 3: printf("attribute %s has only one value\n", s1);
		break;

	case 4: printf("case %d's value of '%s' for attribute %s is illegal\n",
		    MaxItem+1, s2, s1);
		break;

	case 5: printf("case %d's class of '%s' is illegal\n", MaxItem+1, s2);
    }

    if ( ++Messages > 10 )
    {
	printf("Error limit exceeded\n");
	exit(1);
    }
}
