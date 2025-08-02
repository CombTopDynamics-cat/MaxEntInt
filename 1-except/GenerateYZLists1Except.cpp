/*
Generates the Y and Z lists two obtain the maximum entropy cycle
in C^1_n (Section 4).

The output is the number of elements in Y and Z lists (stdout)
respectively, and the files containing the actual lists.
Their names are "freealphalist-1ex-Y.<odd-length>" and
"freealphalist-1ex-Z.<odd-length>"

The files are direct acces and are written in an "A-list format".
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
   FILE *f; char filename[257];
   unsigned long int number = 0L;

   // Testing parameters and open output file
   if( argc != 2 || (p = atoi (argv[1])) <= 0 || !(p%2)) Usage(argv[0], "odd-length");

   AlphaLists1Exc perm (p);

   sprintf(filename,"freealphalist-1ex-Y.%d",p);
   f = sfopen(filename, "wb", "when opening binary Y list output file");
   fwrite(&p,  sizeof(unsigned char),  1, f);
   fwrite(&number,  sizeof(unsigned long int),  1, f);

   puts("Generating the Y list ...");
   perm.InitializeY();
   do{ fwrite(perm.perm,  sizeof(unsigned char),  p, f); number ++; } while (perm.NextY());
   printf("Total number of elements in the Y list: %lu\n\n", number);
   fseek(f, sizeof(unsigned char) , SEEK_SET); fwrite(&number,  sizeof(unsigned long int),  1, f);
   fclose (f);

   number = 0L;
   sprintf(filename,"freealphalist-1ex-Z.%d",p);
   f = sfopen(filename, "wb", "when opening binary Z list output file");
   fwrite(&p,  sizeof(unsigned char),  1, f);
   fwrite(&number,  sizeof(unsigned long int),  1, f);

   puts("Generating the Z list ...");
   perm.InitializeZ();
   do{ fwrite(perm.perm,  sizeof(unsigned char),  p, f); number ++; } while (perm.NextZ());
   printf("Total number of elements in the Z list: %lu\n\n", number);
   fseek(f, sizeof(unsigned char) , SEEK_SET); fwrite(&number,  sizeof(unsigned long int),  1, f);
   fclose (f);

   return 0;
}
