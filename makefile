CC = gcc

SOURCES = cluster.c
LIBS = -lm 
PROGRAM = cl
MAKEDEP = domakedep

OBJECTS= $(SOURCES:.c=.o)
LINTFILES=$(SOURCES:.c=.ln)

CFLAGS =  -O3 

$(PROGRAM): $(OBJECTS) 
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) $(LIBS)

$(LINTFILES):
	$(LINT) $(LINTFILES) 

clean:
	rm -f *.o core

#clall : clean all
clall:
	rm -f *.o core a.out $(PROGRAM)

depend:
	$(MAKEDEP) *.[chC] > deps_top

#include deps_top
# DO NOT DELETE THIS LINE -- ccdep uses it.
# DO NOT PUT ANYTHING AFTER THIS LINE, IT WILL GO AWAY.



# IF YOU PUT ANYTHING HERE IT WILL GO AWAY

