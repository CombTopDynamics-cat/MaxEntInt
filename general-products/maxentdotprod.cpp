/*
Computes the maximum entropy of the dot-product permutations
a <dot> b of the form descibed below WHICH ARE CYCLES (Section 3):

On input a file that must be in a valid A-list format must be
specified. Optionally, a segment of the list can be specified
(the default segment is the whole list).

Then the program generates all products a <dot> b with:
a varying in the chosen segment
b \ge a in the order of the list
AND KEEPS ONLY THOSE WHICH ARE CYCLES.

For these cycles the program computes the associated Markov Matrix
and then its spectral radius. During this process the absolute
maximumis selected. The output gives the entropy and the product cycle
giving it. The cycle is specified by giving the position occupied by
a and b in the alpha list. To know the actual elements a and b, the
cycle in map form and the Markov matrix with all row and colum sums
together with the total sum the program singleprodinfo.exe must be used.

On input the parameter entropy-filter-level must be specified. We only
look for cycles with entropy larger than this specified entropy-filter-level.
This is a devide to increase the speed of the program in critical cases when
we know in advance that the entropy must be larger than an a priori value.

Usage: filename entropy-filter-level [ini:end]
*/
#include <stdio.h>
#include <stdlib.h>
#include "../BaseFunctions/permutations.h"
#include "../BaseFunctions/markovmatrix.h"
#include "../BaseFunctions/maxent.h"
#define TOL 1.e-8

int main (int argc, char *argv[])
{
  unsigned long int ini=0, end=0, endlist;
  register unsigned long int a, b;
  double maxentini=0.0, maxent=0.0, topent;

  // Testing parameters and open output file
  if( argc < 3 || argc > 4 || (maxent = maxentini = atof(argv[2])) < 1.0 \
      || (argc == 4 && !IntParam(argv[3],&ini,&end) ) ) Usage(argv[0], "list-filename entropy-filter-level [ini:end]");

  ProductPermutation alphabeta(argv[1]); endlist = alphabeta.SizeOfAlphaList();
  if(!end) end = endlist;
  IntervalMarkovMatrix MarMat(alphabeta.Period());

  for(a=ini; a < end ; a++) { alphabeta.Readalpha(a);
    for (b=a; b < endlist ; b++){ alphabeta.Readbeta(b);
      if(!alphabeta.Dot()) continue;
      MarMat.Compute((BasePermutation) alphabeta);
      topent = MarMat.ComputeSpectralRadiusIfNotSmallerThan (maxent);
      if (topent < maxent - TOL) MarMat.StopIfErrorInPowerMethod (maxent,(BasePermutation) alphabeta);
      else if (topent > maxent + TOL) { maxent = topent; ResetInfoResults(a, b, '.', topent); }
      else AddInfoToResults(a, b, '.', topent);
    }
  }

  PrintInfoResults(stdout, argv[1], alphabeta.Period(), ini, end, maxentini);
  return 0;
}
