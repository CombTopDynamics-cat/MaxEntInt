# Makefile of maxentint
FONTDIR=$(PWD)
BASEFUNCTIONS=$(PWD)/../BaseFunctions

OBJS=permutations.o maxent.o
A_OBJS=$(OBJS) markovmatrix.o
HEADERS= $(BASEFUNCTIONS)/permutations.h $(BASEFUNCTIONS)/markovmatrix.h  $(BASEFUNCTIONS)/maxent.h
#OPT=-O3 -g
OPT=-O3
CC=g++ $(OPT) -Wall
#FC=g77 $(OPT) -Wall

GenerateFreeAList: $(OBJS) GenerateFreeAList.o ; $(CC) -o GenerateFreeAList.exe GenerateFreeAList.o $(OBJS)
ReadAList: $(A_OBJS) ReadAList.o ; $(CC) -o ReadAList.exe ReadAList.o $(A_OBJS)
maxentcrossprod: $(A_OBJS) maxentcrossprod.o ; $(CC) -o maxentcrossprod.exe maxentcrossprod.o $(A_OBJS)
maxentdotprod: $(A_OBJS) maxentdotprod.o ; $(CC) -o maxentdotprod.exe maxentdotprod.o $(A_OBJS)
singleprodinfo: $(OBJS) markovmatrix.o singleprodinfo.o ; $(CC) -o singleprodinfo.exe singleprodinfo.o $(OBJS) markovmatrix.o

permutations.o: $(BASEFUNCTIONS)/permutations.cpp $(BASEFUNCTIONS)/permutations.h
	$(CC) -c $(BASEFUNCTIONS)/permutations.cpp

markovmatrix.o: $(BASEFUNCTIONS)/markovmatrix.cpp $(BASEFUNCTIONS)/markovmatrix.h
	$(CC) -c $(BASEFUNCTIONS)/markovmatrix.cpp

maxent.o: $(BASEFUNCTIONS)/maxent.cpp $(BASEFUNCTIONS)/maxent.h
	$(CC) -c $(BASEFUNCTIONS)/maxent.cpp

all: GenerateFreeAList ReadAList maxentcrossprod maxentdotprod singleprodinfo
clean:
	rm -f *~
	rm -f *.o
	rm -f $(BASEFUNCTIONS)/*~
