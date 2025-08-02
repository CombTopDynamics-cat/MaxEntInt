/*********************************************************
   Class LongBasePermutation declarations
   Used by mapentint.cpp and maxentallbounds.cpp
   This is another instance of class BasePermutations
   that allows lengths up to LONGMAXPERM (=1000).

   For comments we refer to the class BasePermutations
**********************************************************/
#define LONGMAXPERM 1000
#define SHIFT_NOTATION 1

class LongBasePermutation
{ public:

  unsigned int size, sm1, sm2;
  unsigned int *perm;
  LongBasePermutation(void) { perm=NULL; size=0; }
  void AllocateLongBasePermutation(unsigned int);
  int DefineManual(int, int *);
  void SimplePrint (FILE *fp){ register unsigned int i;
    fprintf (fp, "(%d", *perm + SHIFT_NOTATION);
    for (i = 1; i < size; i++) fprintf (fp, ",%u", perm[i] + SHIFT_NOTATION);
    fprintf (fp, ")");
  };
  void Print (FILE *fp){ fprintf (fp, "Permutation: "); SimplePrint(fp); fprintf(fp,"\n"); }
  int IsCycle(void);
  int IsMaximodal(void);
  int ComputeNumberOfExceptionsIfMaximodal(void);
};
