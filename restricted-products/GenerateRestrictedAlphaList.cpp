/*
Generates an A list of length n (n odd) determined by a pattern
i_1 i_2 i_3 ... i_n. The positions i_j \ne 0 are fixed in the
list and taken from the input. The positions i_j = 0 are free
and filled in a permutation like way.

Set F = {i_j : j=1,2..,n  and i_j = 0 } \subset {1,2,...n}
and k = Card F.

The program first constructs a translation map from {1,2,...,k}
to F and a list of the free positions in the pattern:
             P = {j=1,2..,n : i_j = 0 }.
Then, each permutation of k elements is translated to F
and substituted in the P positions to generate, together with
the non-free positions, a restricted permutation of length n.

Usage pattern i_1 i_2 i_3 ... i_n
NOTES: n must be odd
       0 means a free position
       No index should be repeated except 0
*/
#include <stdio.h>
#include <stdlib.h>
#include "../BaseFunctions/permutations.h"
#include "../BaseFunctions/maxent.h"

int main (int argc, char *argv[])
{
  unsigned char p = argc-1;
  FILE *f; char name[257];
  unsigned long int number = 0L;
  register unsigned char i;
  int map[256];

  // Testing parameters and open output file
  if((argc < 4 || !(p%2) )){
USAGE:;
    Usage(argv[0], "<alpha pattern i_1 i_2 i_3 ... i_n>\n\
NOTES: n must be odd\n\
       0 means a free position\n\
       No index should be repeated except 0");
  }

  RestrictedAlphaPermutation perm (p);
  for(i=0; i<p;i++) map[i] = atoi (argv[i+1]);
  if(!perm.InitializeRestrictedAlphaPermutation(map)) goto USAGE;

  sprintf(name,"restricted-A-list.%d",p);
  f = sfopen(name, "wb", "when opening A list output file");
  fwrite(&p,  sizeof(unsigned char),  1, f);
  fwrite(&number,  sizeof(unsigned long int),  1, f);

  do{ fwrite(perm.perm,  sizeof(unsigned char),  p, f); number ++; } while (perm.Next());
  fseek(f, sizeof(unsigned char) , SEEK_SET); fwrite(&number,  sizeof(unsigned long int),  1, f);
  fclose (f);
  printf("Total number of elements in the restricted A list: %lu\n\n", number);
  return 0;
}
