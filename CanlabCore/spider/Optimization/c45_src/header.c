/*************************************************************************/
/*									 */
/*	Print header for all C4.5 programs				 */
/*	----------------------------------				 */
/*									 */
/*************************************************************************/


#define  RELEASE "8"


    PrintHeader(Title)
/*  -----------  */
    char *Title;
{
    char *ctime(), TitleLine[80];
    long clock, time();
    short Underline;

    clock = time(0);
    sprintf(TitleLine, "C4.5 [release %s] %s", RELEASE, Title);
    printf("\n%s\t%s", TitleLine, ctime(&clock));

    Underline = strlen(TitleLine);
    while ( Underline-- ) putchar('-');
    putchar('\n');
}
