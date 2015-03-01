#*************************************************************************#
#*									 *#
#*		Makefile for the C4.5 induction system			 *#
#*		--------------------------------------			 *#
#*									 *#
#*************************************************************************#
#	Set CFLAGS for your C compiler, e.g.
#
#	  * if you need to debug the programs (sorry about that!)
#	      CFLAGS = -g
#	    (You might also want to reactivate the lint commands that
#	    have been commented out below.)
#
#	  * if your C compiler supports higher optimization levels
#	      CFLAGS = -O3 etc
CFLAGS = -O2
#	Definitions of file sets
core.ln = \
	getnames.ln getdata.ln trees.ln getopt.ln header.ln
trees.ln = \
	besttree.ln build.ln info.ln discr.ln contin.ln subset.ln\
	prune.ln stats.ln st-thresh.ln classify.ln confmat.ln\
	sort.ln $(core.ln)
rules.ln = \
	rules.ln genlogs.ln genrules.ln makerules.ln prunerule.ln\
	siftrules.ln testrules.ln stats.ln confmat.ln sort.ln $(core.ln)
core = \
	getnames.o getdata.o trees.o getopt.o header.o
trees = \
	besttree.o build.o info.o discr.o contin.o subset.o prune.o\
	stats.o st-thresh.o classify.o confmat.o sort.o $(core)
rules = \
	rules.o genlogs.o genrules.o makerules.o prunerule.o\
	siftrules.o testrules.o stats.o confmat.o sort.o $(core)
#  C4.5 decision tree generator
c4.5:	c4.5.o $(trees)
#	lint -x c4.5.ln $(trees.ln) -lm
	lcc -o c4.5 c4.5.o $(trees) -lm
#	(Sample only: for compiler that supports global optimization,
#	 e.g. DECstation 3100)
c4.5gt:
	cat	defns.i types.i\
		c4.5.c\
		besttree.c build.c info.c discr.c contin.c subset.c\
		prune.c stats.c st-thresh.c confmat.c sort.c\
		getnames.c getdata.c classify.c trees.c header.c\
		| egrep -v 'defns.i|types.i|extern.i|buildex.i' >c4.5gt.c
	lcc -O4 -o c4.5gt c4.5gt.c -lm
	rm c4.5gt.c
#  C4.5 production rule generator
c4.5rules: c4.5rules.o $(rules)
#	lint -x c4.5rules.ln $(rules.ln) -lm
	lcc -o c4.5rules c4.5rules.o $(rules) -lm
c4.5rulesgt:
	cat	defns.i types.i\
		c4.5rules.c\
		rules.c genlogs.c genrules.c makerules.c prunerule.c\
		siftrules.c testrules.c stats.c confmat.c sort.c\
		getnames.c getdata.c classify.c trees.c header.c\
		| egrep -v 'defns.i|types.i|extern.i|rulex.i' >c4.5rulesgt.c
	lcc -O4 -o c4.5rulesgt c4.5rulesgt.c -lm
	rm c4.5rulesgt.c
#  C4.5 decision tree interpreter
consult: consult.o userint.o $(core)
#	lint -x consult.ln userint.ln $(core.ln)
	lcc -o consult consult.o userint.o $(core)
#  C4.5 production rule interpreter
consultr: consultr.o rules.o userint.o $(core)
#	lint -x consultr.ln rules.ln userint.ln $(core.ln)
	lcc -o consultr consultr.o rules.o userint.o $(core)
.c.o:
#	lint -c $<
	lcc $(CFLAGS) -c $<
#  Make everything
all:
	make c4.5
	make c4.5rules
	make consult
	make consultr
	lcc -o xval-prep xval-prep.c
	lcc -o average average.c -lm
$(trees): defns.i types.i extern.i
$(rules): defns.i types.i extern.i
