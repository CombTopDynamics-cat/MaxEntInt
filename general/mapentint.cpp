/*
Computes the entropy of a cycle specified by a permutation in
map form. The output gives the entropy, the cycle
in map form and the Markov matrix with all row and colum sums
together with the total sum.

It works for period smaller than 1000.

Usage: <map vector i_1 i_2 i_3 ... i_n (2 < length)>
*/
#include <stdio.h>
#include <stdlib.h>
#include "../BaseFunctions/longpermutations.h"
#include "../BaseFunctions/markovmatrix.h"
#include "../BaseFunctions/maxent.h"

int main (int argc, char *argv[])
{
  register int i;
  int map[LONGMAXPERM];
  double ent;

  if(argc < 4){
USAGE:; Usage(argv[0], "<map vector i_1 i_2 i_3 ... i_n (2 < length)>");
  } for(i=1; i<argc;i++) { if((map[i-1] = atoi (argv[i])) <= 0) goto USAGE; }

  LongBasePermutation perm;  if(!perm.DefineManual(argc-1, map)) goto USAGE;

  IntervalMarkovMatrix MarMat (argc-1);
  MarMat.MapCompute (perm.perm);
  ent = MarMat.SpectralRadius ();

  perm.Print (stdout);
  fprintf(stdout, "Entropy = log %.16lf\n", ent);
  MarMat.Print (stdout);

  if (ent < -10) ErrorAborting("no convergence in Power Method");
  else return 0;
}
