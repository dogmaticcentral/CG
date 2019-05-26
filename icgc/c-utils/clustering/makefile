CFLAGS=     # options for C compiler
FFLAGS =
FCOMPL= f77   # name of Fortran compiler
CCOMPL = gcc -Wall
CPLUSCOMPL = g++

HERE = $(shell pwd)

NAME = pc

OPTM =  -g

VPATH = ./

PLATFORM = $(shell uname)

INCLUDE = 

CLIBPATH =
LIBS =  -lm
 
ifeq (${PLATFORM},Linux)
	PD =  -DLINUX
endif

DEFINES =  $(PD)  -DVERBOSE 
SOURCES.c =  main.c  read_pdb.c  read_selection.c utils.c dist_matrix.c clustercounter.c clustering_score.c





AUX =  pdbclust.h pdb.h utils.h makefile

SOURCES = $(SOURCES.c)  $(SOURCES.F) $(SOURCES.f)
OBJECTS =  $(SOURCES.F:.F=.o)  $(SOURCES.f:.f=.o) $(SOURCES.c:.c=.o)


$(NAME): $(OBJECTS)
	$(CCOMPL) $(FFLAGS) -o $@  $(OBJECTS) $(OPTM)  $(LIBS)    $(INCLUDE)

.c.o:
	$(CCOMPL)  $(CFLAGS)  -c   $(DEFINES) $<  $(OPTM) $(INCLUDE)

.F.o:
	$(FCOMPL)  $(CFLAGS)  -c   $(DEFINES) $<  $(OPTM) $(INCLUDE)

clean:
	rm $(OBJECTS)

tar:
	tar -cf $(NAME).tar $(SOURCES) $(AUX)



ls:
	 ls $(SOURCES) $(AUX) 


