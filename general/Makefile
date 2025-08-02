FONTDIR=$(PWD)
BASEFUNCTIONS=$(PWD)/../BaseFunctions

OBJS=permutations.o markovmatrix.o maxent.o
HEADERS= $(BASEFUNCTIONS)/permutations.h $(BASEFUNCTIONS)/markovmatrix.h $(BASEFUNCTIONS)/maxent.h

L_OBJS=longpermutations.o markovmatrix.o maxent.o
L_HEADERS= $(BASEFUNCTIONS)/longpermutations.h $(BASEFUNCTIONS)/markovmatrix.h $(BASEFUNCTIONS)/maxent.h

#OPT=-O3 -g
OPT=-O3
CC=g++ $(OPT) -Wall
#FC=g77 $(OPT) -Wall

maxentallbounds: maxentallbounds.o $(L_OBJS); $(CC) -o maxentallbounds.exe maxentallbounds.o $(L_OBJS)
mapentint: mapentint.o $(L_OBJS); $(CC) -o mapentint.exe mapentint.o $(L_OBJS)
maxentint: $(OBJS) maxentint.o ; $(CC) -o maxentint.exe maxentint.o $(OBJS)
CountJungreisCycles: $(OBJS) CountJungreisCycles.o ; $(CC) -o CountJungreisCycles.exe CountJungreisCycles.o $(OBJS)

permutations.o: $(BASEFUNCTIONS)/permutations.cpp $(BASEFUNCTIONS)/permutations.h
	$(CC) -c $(BASEFUNCTIONS)/permutations.cpp

longpermutations.o: $(BASEFUNCTIONS)/longpermutations.cpp $(BASEFUNCTIONS)/longpermutations.h
	$(CC) -c $(BASEFUNCTIONS)/longpermutations.cpp

markovmatrix.o: $(BASEFUNCTIONS)/markovmatrix.cpp $(BASEFUNCTIONS)/markovmatrix.h
	$(CC) -c $(BASEFUNCTIONS)/markovmatrix.cpp

maxent.o: $(BASEFUNCTIONS)/maxent.cpp $(BASEFUNCTIONS)/maxent.h
	$(CC) -c $(BASEFUNCTIONS)/maxent.cpp

all: maxentallbounds mapentint maxentint CountJungreisCycles
clean:
	rm -f *~
	rm -f *.o
	rm -f $(BASEFUNCTIONS)/*~