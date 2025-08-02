/*
Computes the maximum entropy of Jungreis cycles by the
method of generating all maximodal cyclic permutations of a given
length and then checking whether they have zero or one exceptions.
Then, computing the associated Markov Matrix and then its spectral
radius. During this process the absolute maximum is selected.

The output gives the entropy, the cycle in map form and the Markov
matrix with all row and colum sums together with the total sum.

It works for period smaller than 14(17?). Then the number of cases is to
big to finish in reasonable time.

Usage: (3 <) period (<= 200) output-filename [entropy-filter-level]
*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "../BaseFunctions/permutations.h"
#include "../BaseFunctions/markovmatrix.h"
#include "../BaseFunctions/maxent.h"

#define TOL 1.e-8

int main (int argc, char *argv[])
{
  unsigned char p=3;
  double maxent=0.0, topent;
//   long unsigned int count = 0;
  FILE *fp, *OpenOutputFile (char *, unsigned char);

// Testing parameters and open output file
  if (argc < 3 || argc > 4 || (p = (unsigned char) atoi (argv[1])) <= 3 || p > 200){
USAGE:;
      Usage(argv[0], "(3 <) period (<= 200) output-filename [entropy-filter-level]");
  } if(argc == 4 && (maxent = atof (argv[3])) <= 0.0) goto USAGE;
  if (!access (argv[2], F_OK)) ErrorAborting("output file already exists");
  fp = OpenOutputFile (argv[2], p);

// Initializations
  MaximodalCiclicIntervalPermutation perm (p);
  IntervalMarkovMatrix MarMat (p);

  do { if(perm.ComputeNumberOfExceptionsIfMaximodal() > 1) continue;
      MarMat.Compute (perm);
      topent = MarMat.ComputeSpectralRadiusIfNotSmallerThan (maxent);
      if (topent < maxent - TOL) { MarMat.StopIfErrorInPowerMethod (maxent, perm); continue; }
      if (topent > maxent + TOL) { maxent = topent; fclose (fp); fp = OpenOutputFile (argv[2], p); }
      MarMat.PrintAll(fp, topent, perm);
  } while (perm.NextNonDual ());
  fclose (fp);
  return 0;
}

// Utility functions for the main program
FILE *OpenOutputFile (char *name, unsigned char period)
{ FILE *f = sfopen(name, "w", "when opening output file");
  fprintf (f, "\n--- Computing maximal entropy for the interval for period %u ---\n\n", period);
  return f;
}
