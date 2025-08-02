/*
Reads a segment (or the whole) A list (Section 3) of the specified
file that must be in a valid A-list format. Optionally a segment
of the list can be specified (the default segment is the whole list).

It can also be used to know the real length of the inyernal permutations
and how may elements has the list.

The output is (stdout):
* the odd-length of the permutations,
* the number of elements in the A list
* and the chosen ([ini:end]) segment of the list

Usage: filename [ini:end]
*/
#include <stdio.h>
#include <stdlib.h>
#include "../BaseFunctions/permutations.h"
#include "../BaseFunctions/maxent.h"

int main (int argc, char *argv[])
{
  unsigned char p;
  FILE *f;
  unsigned long int number, ini=0, end=0;
  register unsigned long int i;

  // Testing parameters and open output file
  if( argc < 2 || argc > 3 || (argc == 3 && !IntParam(argv[2],&ini,&end) ) ) Usage(argv[0], "filename [ini:end]");

  f = sfopen(argv[1], "rb", "when opening list file");
  fread(&p,  sizeof(unsigned char),  1, f);
  fread(&number,  sizeof(unsigned long int),  1, f);

  BasePermutation perm; perm.AllocateBasePermutation(p);

  printf("A list permutation length: %u\n", p);
  printf("Total number of elements in the A list: %lu\n\n", number);
  if(!end) end = number;

  fseek(f, sizeof(unsigned char) + sizeof(unsigned long int) + ini*p*sizeof(unsigned char) , SEEK_SET);
  for(i=ini; i < end ;i++){
     fread(perm.perm,  sizeof(unsigned char),  p, f);
     perm.Print(stdout); // int j; for(j=0; j < p; j++) { printf("%d,", perm.perm[j]+1); }  printf("\n");
  }
  fclose (f);
  return 0;
}
