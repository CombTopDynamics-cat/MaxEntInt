# Makefile of maxentint
FONTDIR=$(PWD)
BASEFUNCTIONS=$(PWD)/../BaseFunctions

OBJS=permutations.o maxent.o
A_OBJS=$(OBJS) markovmatrix.o
HEADERS= $(BASEFUNCTIONS)/permutations.h $(BASEFUNCTIONS)/markovmatrix.h $(BASEFUNCTIONS)/maxent.h

#OPT=-O3 -g
OPT=-O3
CC=g++ $(OPT) -Wall
#FC=g77 $(OPT) -Wall

maxent-1ex-crossprod: $(A_OBJS) maxent-1ex-crossprod.o ; $(CC) -o maxent-1ex-crossprod.exe maxent-1ex-crossprod.o $(A_OBJS)

GenerateYZLists1Except: $(OBJS) GenerateYZLists1Except.o ; $(CC) -o GenerateYZLists1Except.exe GenerateYZLists1Except.o $(OBJS)

maxentint1except: $(A_OBJS) maxentint1except.o ; $(CC) -o maxentint1except.exe maxentint1except.o $(A_OBJS)

permutations.o: $(BASEFUNCTIONS)/permutations.cpp $(BASEFUNCTIONS)/permutations.h
	$(CC) -c $(BASEFUNCTIONS)/permutations.cpp

markovmatrix.o: $(BASEFUNCTIONS)/markovmatrix.cpp $(BASEFUNCTIONS)/markovmatrix.h
	$(CC) -c $(BASEFUNCTIONS)/markovmatrix.cpp

maxent.o: $(BASEFUNCTIONS)/maxent.cpp $(BASEFUNCTIONS)/maxent.h
	$(CC) -c $(BASEFUNCTIONS)/maxent.cpp

all: GenerateYZLists1Except maxent-1ex-crossprod maxentint1except
clean:
	rm -f *~
	rm -f *.o
	rm -f $(BASEFUNCTIONS)/*~