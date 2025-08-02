/**************************************************************
** General auxiliary functions declarations                  **
** Used in maxent.cpp, permutations.cpp and markovmatrix.cpp **
**************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../BaseFunctions/maxent.h"


/* To decode integer parameters in a predefined interval */
int IntParam(const char *parm, unsigned long int *n1, unsigned long int *n2)
{  char *p; long int n;

   if((n = atoi(parm)) <= 0) return 0;
   *n1 = n-1; *n2=0; if((p = strstr(parm,":")) == NULL) return 1;
   p++; n = atoi(p); if(n-1 < (signed long int) *n1) return 0; else *n2 = n; return 1;
}

/****************************
** Error message functions **
****************************/
void ErrorMessage(const char *message) { fprintf (stderr, "ERROR: %s", message); }

void ErrorAborting(const char *message)
{
   fprintf (stderr, "\n"); ErrorMessage(message);
   fprintf (stderr, ".\nAborting...\n\n");
   exit(1);
}

void Usage(const char *programname, const char *message)
{
   fprintf (stderr, "\nUSAGE: %s %s\n\n", programname, message); exit(1);
}

FILE *sfopen(const char *path, const char *mode, const char *message)
{
   FILE *f = fopen (path, mode); if (f == NULL) { ErrorAborting(message); return NULL; } return f;
}

/********************************************
** Store and print results for maxentprods **
********************************************/
char type[1000];
int numcycles=-1;
unsigned long int alphaindmax[1000], betaindmax[1000];
double topentcurrcycle[1000];

void ResetInfoResults(unsigned long int alphaind, unsigned long int betaind, char prodt, double entropy)
{ numcycles=0; alphaindmax[0] = alphaind; betaindmax[0] = betaind; type[0] = prodt; topentcurrcycle[0] = entropy; }

void AddInfoToResults(unsigned long int alphaind, unsigned long int betaind, char prodt, double entropy)
{ numcycles++;
  if(numcycles >= 1000) {
     fprintf(stderr, "Warning: Too many equal entropy cycles at entropy level %.16lf\n", entropy);
     return;
  }
  alphaindmax[numcycles] = alphaind; betaindmax[numcycles] = betaind; type[numcycles] = prodt;
  topentcurrcycle[numcycles] = entropy;
}

void PrintInfoResults(FILE *fout, const char *filebasename, unsigned char per, unsigned long int ini, unsigned long int end, double maxentini){
  int i;
  fprintf(fout, "Computing maximal entropy of products for period: %d\n",per);
  fprintf(fout, "Listname: %s\n",filebasename);
  fprintf(fout, "Alpha interval checked: %ld to %ld\n", ini+1, end);
  if(maxentini > 1.0) fprintf(fout, "Filtering level %f\n", maxentini);
  if(numcycles < 0) {
    fprintf(fout, "All entropies of this interval are smaller than the filtering level\n");
    return;
  }
  fprintf(fout, "Maximum of entropy found for the following cycles:\n");
  for(i=0; i <= numcycles ; i++)
    fprintf(fout, "%.16lf (%ld %c %ld)\n", topentcurrcycle[i], alphaindmax[i]+1, type[i], betaindmax[i]+1);
}


