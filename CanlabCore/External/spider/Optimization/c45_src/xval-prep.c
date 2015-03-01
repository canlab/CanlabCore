/*************************************************************************/
/*									 */
/*	Program to prepare data file for cross-validation		 */
/*	-------------------------------------------------		 */
/*									 */
/*	The number of blocks for the cross-validation appears as the	 */
/*	first argument.  The data are shuffled and divided into the	 */
/*	specified number of blocks, with class distributions as even	 */
/*	as possible in each block.					 */
/*									 */
/*************************************************************************/


#include <math.h>
#include <stdio.h>


long	random();
#define	randf			((random()&2147483647) / 2147483648.0)
#define	ForEach(var,F,L)	for(var=F; var<=L; ++var) 
#define	dig(x)			(x >= '0' && x <= '9')

#define MAXLINE 5000		/* maximum line length */

char	**Item;
int	ItemSpace=1000, MaxItem=0;


    main(argc, argv)
/*  ----  */
    int argc;
    char *argv[];
{
    int i, First=0, Last, Length, Splits;
    char Line[MAXLINE], **ClassPtr, *Temp, *BeginClass();

    sscanf(argv[1], "%d", &Splits);

    Item = (char **) malloc(ItemSpace * sizeof(char *));

    while ( fgets(Line, MAXLINE, stdin) )
    {
	if ( MaxItem >= ItemSpace )
	{
	    ItemSpace += 1000;
	    Item = (char **) realloc(Item, ItemSpace * sizeof(char *));
	}

	Length = strlen(Line)+2;
	Item[MaxItem] = (char *) malloc(Length);
	memcpy(Item[MaxItem], Line, Length);
	MaxItem++;
    }

    if ( ! MaxItem-- ) exit(1);

    Shuffle();

    /*  Find classes  */

    ClassPtr = (char **) malloc((MaxItem+1) * sizeof(char *));
    ForEach(i, 0, MaxItem)
    {
	ClassPtr[i] = BeginClass(Item[i]);
    }

    /*  Sort by class  */

    fprintf(stderr, "\nClass frequencies:\n");

    while ( First <= MaxItem )
    {
	Last = First;

	ForEach(i, First+1, MaxItem)
	{
	    if ( ! strcmp(ClassPtr[i], ClassPtr[First]) )
	    {
		Last++;
		Temp = Item[Last];
		Item[Last] = Item[i];
		Item[i] = Temp;

		Temp = ClassPtr[Last];
		ClassPtr[Last] = ClassPtr[i];
		ClassPtr[i] = Temp;
	    }
	}

	fprintf(stderr, "%6d class %s\n", Last-First+1, ClassPtr[First]);

	First = Last+1;
    }

    ForEach(First, 0, Splits-1)
    {
	for ( i = First ; i <= MaxItem ; i += Splits )
	{
	    printf("%s\n", Item[i]);
	}
    }
}



/*************************************************************************/
/*									 */
/*	Find the beginning character of a class name			 */
/*									 */
/*************************************************************************/


char *BeginClass(S)
/*    ----------  */
    char *S;
{
    char *F;

    F = S - 1;
    do
    {
	S = F + 1;
	while ( *S == ' ' || *S == '\t' || *S == '\n' ) S++;
	F = S;
	while ( *F != ',' && (*F != '.' || dig(*(F+1))) && *F != '\n' ) F++;
    } while ( *F == ',' );

    if ( *F != '.' ) *F = '.';
    *(F+1) = '\0';

    return S;
}



/*************************************************************************/
/*									 */
/*	Shuffle the data items						 */
/*									 */
/*************************************************************************/


    Shuffle()
/*  -------  */
{
    int this, alt, left = MaxItem+1;
    char *hold;

    this = 0;
    while ( left )
    {
        alt = this + (left--) * randf;
	if ( alt > MaxItem || alt < this )
	{
	    fprintf(stderr, "ERROR!\n");
	    exit(1);
	}
        hold = Item[this];
        Item[this++] = Item[alt];
        Item[alt] = hold;
    }
}
