# Makefile of maxentint
FONTDIR=$(PWD)
BASEFUNCTIONS=$(PWD)/../BaseFunctions

OBJS=permutations.o maxent.o
HEADERS= $(BASEFUNCTIONS)/permutations.h $(BASEFUNCTIONS)/maxent.h

OPT=-O3
CC=g++ $(OPT) -Wall

GenerateRestrictedAlphaList: $(OBJS) GenerateRestrictedAlphaList.o ; $(CC) -o GenerateRestrictedAlphaList.exe GenerateRestrictedAlphaList.o $(OBJS)

permutations.o: $(BASEFUNCTIONS)/permutations.cpp $(BASEFUNCTIONS)/permutations.h
	$(CC) -c $(BASEFUNCTIONS)/permutations.cpp

maxent.o: $(BASEFUNCTIONS)/maxent.cpp $(BASEFUNCTIONS)/maxent.h
	$(CC) -c $(BASEFUNCTIONS)/maxent.cpp

all: GenerateRestrictedAlphaList
clean:
	rm -f *~
	rm -f *.o
	rm -f $(BASEFUNCTIONS)/*~