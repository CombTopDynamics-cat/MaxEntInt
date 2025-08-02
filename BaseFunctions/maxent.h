/**************************************************************
** General auxiliary functions declarations                  **
** Used in maxent.cpp, permutations.cpp and markovmatrix.cpp **
**************************************************************/
class MaximodalCiclicIntervalPermutation;
class IntervalMarkovMatrix;

/* To decode integer parameters in a predefined interval */
int IntParam(const char *, unsigned long int *, unsigned long int *);

/****************************
** Error message functions **
****************************/
void ErrorMessage(const char *);
void ErrorAborting(const char *);
void Usage(const char *, const char *);
FILE *sfopen(const char *, const char *, const char *);

/********************************************
** Store and print results for maxentprods **
********************************************/
void ResetInfoResults(unsigned long int, unsigned long int, char, double);
void AddInfoToResults(unsigned long int, unsigned long int, char, double);
void PrintInfoResults(FILE *, const char *, unsigned char, unsigned long int, unsigned long int, double);
