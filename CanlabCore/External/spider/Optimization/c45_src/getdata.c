/*************************************************************************/
/*									 */
/*	Get case descriptions from data file				 */
/*	--------------------------------------				 */
/*									 */
/*************************************************************************/


#include "defns.i"
#include "types.i"
#include "extern.i"

#define Inc 2048



/*************************************************************************/
/*									 */
/*  Read raw case descriptions from file with given extension.		 */
/*									 */
/*  On completion, cases are stored in array Item in the form		 */
/*  of Descriptions (i.e. arrays of attribute values), and		 */
/*  MaxItem is set to the number of data items.				 */
/*									 */
/*************************************************************************/

    GetData(Extension)
/*  --------  */
    String Extension;
{
    FILE *Df, *fopen();
    char Fn[100];
    ItemNo i=0, j, ItemSpace=0;
    Description GetDescription();

    /*  Open data file  */

    strcpy(Fn, FileName);
    strcat(Fn, Extension);
    if ( ! ( Df = fopen(Fn, "r") ) ) Error(0, Fn, "");

    do
    {
	MaxItem = i;

	/*  Make sure there is room for another item  */

	if ( i >= ItemSpace )
	{
	    if ( ItemSpace )
	    {
		ItemSpace += Inc;
		Item = (Description *)
			realloc(Item, ItemSpace*sizeof(Description));
	    }
	    else
	    {
		Item = (Description *)
			malloc((ItemSpace=Inc)*sizeof(Description));
	    }
	}

	Item[i] = GetDescription(Df);

    } while ( Item[i] != Nil && ++i );

    fclose(Df);
    MaxItem = i - 1;
}



/*************************************************************************/
/*									 */
/*  Read a raw case description from file Df.				 */
/*									 */
/*  For each attribute, read the attribute value from the file.		 */
/*  If it is a discrete valued attribute, find the associated no.	 */
/*  of this attribute value (if the value is unknown this is 0).	 */
/*									 */
/*  Returns the Description of the case (i.e. the array of		 */
/*  attribute values).							 */
/*									 */
/*************************************************************************/


Description GetDescription(Df)
/*          ---------------  */
    FILE *Df;
{
    Attribute Att;
    char name[500], *endname, *CopyString();
    Boolean ReadName();
    int Dv;
    float Cv;
    Description Dvec;
    double strtod();

    if ( ReadName(Df, name) )
    {
	Dvec = (Description) calloc(MaxAtt+2, sizeof(AttValue));

        ForEach(Att, 0, MaxAtt)
        {
	    if ( SpecialStatus[Att] == IGNORE )
	    {
		/*  Skip this value  */

		DVal(Dvec, Att) = 0;
	    }
	    else
	    if ( MaxAttVal[Att] || SpecialStatus[Att] == DISCRETE )
	    {
		/*  Discrete value  */ 

	        if ( ! ( strcmp(name, "?") ) )
		{
		    Dv = 0;
		}
	        else
	        {
		    Dv = Which(name, AttValName[Att], 1, MaxAttVal[Att]);
		    if ( ! Dv )
		    {
			if ( SpecialStatus[Att] == DISCRETE )
			{
			    /*  Add value to list  */

			    Dv = ++MaxAttVal[Att];
			    if ( Dv > (int) AttValName[Att][0] )
			    {
				printf("\nToo many values for %s (max %d)\n",
					AttName[Att], (int) AttValName[Att][0]);
				exit(1);
			    }

			    AttValName[Att][Dv] = CopyString(name);
			}
			else
			{
			    Error(4, AttName[Att], name);
			}
		    }
	        }
	        DVal(Dvec, Att) = Dv;
	    }
	    else
	    {
		/*  Continuous value  */

	        if ( ! ( strcmp(name, "?") ) )
		{
		    Cv = Unknown;
		}
	        else
		{
		    Cv = strtod(name, &endname);
		    if ( endname == name || *endname != '\0' )
			Error(4, AttName[Att], name);
		}
		CVal(Dvec, Att) = Cv;
	    }

	    ReadName(Df, name);
        }

        if ( (Dv = Which(name, ClassName, 0, MaxClass)) < 0 )
        {
	    Error(5, "", name);
	    Dv = 0;
        }
	Class(Dvec) = Dv;

	return Dvec;
    }
    else
    {
	return Nil;
    }
}
