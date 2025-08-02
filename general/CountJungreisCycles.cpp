/*
Computes the number of non-conjugate Jungreis cycles
of a given period and how many of them have zero
exceptions (blong to C^0_n)

It works for period smaller than 200.

Usage: (3 <) period (<= 200)
*/
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "../BaseFunctions/permutations.h"
#include "../BaseFunctions/maxent.h"

int main (int argc, char *argv[])
{
  unsigned int p=4, fmzi=0, fmz=0;
  unsigned long int mzi=0L, mz=0L;

   if (argc != 2 || (p = (unsigned int) atoi (argv[1])) <= 3 || p > 200)
      Usage(argv[0], "(3 <) period (<= 200)");

  MaximodalCiclicIntervalPermutation perm (p);
  do{
/*
      int test; test = perm.ComputeNumberOfExceptionsIfMaximodal();
      if(test <= 1) { printf("%d ",test); perm.Print(stdout); }
*/
      switch (perm.ComputeNumberOfExceptionsIfMaximodal())
      {
	 case 0: mz++; mzi++;
                 if(mzi == ULONG_MAX) { fmzi++; mzi=0L;}
                 if(mz == ULONG_MAX) { fmz++; mz=0L;}
                 break;
	 case 1: mzi++;
                 if(mzi == ULONG_MAX) { fmzi++; mzi=0L;}
                 break;
      }
  } while (perm.NextNonDual ()) ;
  if(!fmzi) printf("Total number of Jungreis cycles %lu\n", mzi);
  else printf("Total number of Jungreis cycles %lu + %u * %lu\n", mzi, fmzi, ULONG_MAX);
  if(!fmz) printf("Total number of Jungreis cycles with zero exceptions %lu\n", mz);
  else printf("Total number of Jungreis cycles with zero exceptions %lu + %u * %lu\n", mz, fmz, ULONG_MAX);
  return 0;
}
