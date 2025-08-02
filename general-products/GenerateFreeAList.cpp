/*
Generates the A list for cross and dot products
two obtain the maximum entropy cycle in C^0_n (Section 3).

The output is the number of elements in the A list (stdout)
and the file containing the actual list. Its name is
"free-A-list.<odd-length>"

The file is a direct acces file and it is written in an "A-list format".
This format consists in initializing the file with the following information
at the beginning:
odd-length: unsigned char
the number <n> of elements of the list: unsigned long int
<n> permutations in the list: <n> x (odd-length x unsigned char)

Usage: odd-length
*/
#include <stdio.h>
#include <stdlib.h>
#include "../BaseFunctions/permutations.h"
#include "../BaseFunctions/maxent.h"

int main (int argc, char *argv[])
{
  unsigned char p;
  FILE *f; char name[257];
  unsigned long int number = 0L;

  // Testing parameters and open output file
  if( argc != 2 || (p = atoi (argv[1])) <= 0 || !(p%2)) Usage(argv[0], "odd-length");

  AlphaPermutation perm (p);

  sprintf(name,"free-A-list.%d",p);
  f = sfopen(name, "wb", "when opening the A list output file");
  fwrite(&p,  sizeof(unsigned char),  1, f);
  fwrite(&number,  sizeof(unsigned long int),  1, f);

  do{ fwrite(perm.perm,  sizeof(unsigned char),  p, f); number ++; } while (perm.Next());
  fseek(f, sizeof(unsigned char) , SEEK_SET); fwrite(&number,  sizeof(unsigned long int),  1, f);
  fclose (f);
  printf("Total number of elements in the A list: %lu\n\n", number);
  return 0;
}
