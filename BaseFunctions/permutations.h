/*********************************************************************

  This file contains the declarations of the following classes
  (all of them derived from class BasePermutation (base class)):
    class BasePermutation (base class)
    class CiclicPermutation
    class MaximodalCiclicIntervalPermutation
    class AlphaPermutation
    class RestrictedAlphaPermutation
    class ProductPermutation
    class AlphaLists1Exc

Used by GenerateYZLists1Except.cpp, maxent-1ex-crossprod.cpp,
  mapentint.cpp, maxentint.cpp, maxentallbounds.cpp,
  CountJungreisCycles.cpp, singleprodinfo.cpp, maxentcrossprod.cpp,
  maxentdotprod.cpp, ReadAList.cpp,GenerateFreeAList.cpp,
  GenerateRestrictedAlphaList.cpp
*********************************************************************/

#define MAXPERM 250
#define SHIFT_NOTATION 1
#define SHIFT_CHAR 48

/*********************************************************************
                      class BasePermutation

perm is understood (as in the whole file unless otherwise said) in MAP
FORM. That is perm[i] denotes the image of i by the permutation perm
(see the first two lines of Section 3 of the paper for the notation).

Also for programming reasons i ranges from 0 to size-1 and perm[i]
belongs to {0,1,...,size-1}. This has to be taken into account (and
corrected in the writing routines of permutations --- SHIFT_NOTATION
and SHIFT_CHAR)
*********************************************************************/
class BasePermutation
{ public:
  void __BaseSTDPermutation(unsigned char order){register unsigned char i;
        AllocateBasePermutation(order);
	for (i = 0; i < size; i++) perm[i] = i;
  };

  unsigned char size, sm1, sm2; // sm1 = size -1; sm2 = size -2;
  unsigned char *perm;
  BasePermutation(void) { perm=NULL; size=0; } // constructor
  void AllocateBasePermutation(unsigned char);
  int DefineManual(int, int *);
  int Next ();
  void SimplePrint (FILE *fp){ register unsigned char i;
    fprintf (fp, "(%u", *perm + SHIFT_NOTATION);
    for (i = 1; i < size; i++) fprintf (fp, ",%u", perm[i] + SHIFT_NOTATION);
    fprintf (fp, ")");
   };
  void Print (FILE *fp){ fprintf (fp, "Permutation: "); SimplePrint(fp); fprintf(fp,"\n"); }
  void RawPrint(FILE * fp){ register unsigned char i;
        for (i = 0; i < size; i++) fprintf (fp, "%c", perm[i]+SHIFT_CHAR); fprintf(fp,"\n"); };
  int IsCycle(void);
  int IsMaximodal(void);
  int TestIfDualIsSmaller (void);
  int ComputeNumberOfExceptionsIfMaximodal(void);
};

/*********************************************************************

             class MaximodalCiclicIntervalPermutation

This aim of class is to successively generate in lexicographical
order maximodal cycles of the interval. The crucial function is. of
course, the function Next(). The jump from a maximodal cyclic
permutation to the ext one is done by carefully choosing the place to
do the shift and then checking that the obtained permutation is
maximodal and cyclic. The operation of choosing the place to
do the shift is done by the function ChooseLegalTryForAPosition.
Then the function NexNonDual() in implemented inline as usual.

Used by maxentint.cpp and CountJungreisCycles.cpp.

*********************************************************************/

class MaximodalCiclicIntervalPermutation:public BasePermutation
{
  int ChooseLegalTryForAPosition(unsigned char, unsigned char, unsigned char *);
  int Next ();
 public:
  MaximodalCiclicIntervalPermutation(unsigned char); // Constructor
  int NextNonDual (){ while (Next ()) { if (!TestIfDualIsSmaller ()) return 1; }  return 0; }
};

/*********************************************************************

             class AlphaPermutation

This is the class that deals with the A-lists (and derived ones)
Then, since the period n is of the form 4k+2, the order of the list is
of the form 2k+1 and, hence, it must be odd. half = (order - 1)/2 = k.

The function AllowsCyclesAfterProducts() implements Corollary 3.10 of
the paper. That is, it checks  for permutations alpha that when
participate as a first element of a product (cross or dot; no matter)
will give a cycle. It returns 0 when the map phi has a cycle and thus
no product starting by alpha can be a cycle. Otherwise it returns 1.

Used by AlphaPermutation

*********************************************************************/

class AlphaPermutation:public BasePermutation
{
  unsigned char half; // = (order - 1)/2
  int AllowsCyclesAfterProducts();
  friend class RestrictedAlphaPermutation;
public:
  AlphaPermutation(unsigned char); // Constructor
  int Next(){while (BasePermutation::Next()){if(AllowsCyclesAfterProducts()) return 1;} return 0;}
};

/*********************************************************************

             class RestrictedAlphaPermutation

This is the class that is used to build the restricted A^*-lists
considered at the end of Section 6. Since it is a derived class of
AlphaPermutation can use its features. In particular size must be
odd and the ability of checking Corollary 3.10 through
AllowsCyclesAfterProducts().

This class works as follows: InitializeRestrictedAlphaPermutation
receives a permutation pattern consisting of natural numbers
(elements that we want to fix in the permutation) and zeroes
(free elements in the permutation). Then
InitializeRestrictedAlphaPermutation checks that the length
of the pattern agrees with the order declared in the
constructor, evaluates which are the elements of the permutation that
we want to fix (different from zero and pairwise disjoint) and the
free elements (marked as zero). It also counts the number
freepart.size of zeroes in the pattern. Then a freepart permutation
of order freepart.size is constructed and initialized. The vectors
freepos and transl (of length freepart.size) are assigned memory and
filled. The first one is a listing of the places which are free in the
permutation. The second one is a "dictionary" to translate a
permutation from the "code" \{1,2,....,freepart.size\} to the code
\{1,2,....,order\} \cap \{elements that are left free by the pattern\}

The function Next() (which is the core of this class) works as
follows:
1- Calls the function freepart.Next() to generate a new permutation
   of length freepart.size (stored in freepart.perm).
2- The inline function copyfreepart() (making use of the vectors
   freepos and transl fills the free part of perm by copying
   freepart.perm to it and translating this vector to the code
\{1,2,....,order\} \cap \{elements that are left free by the pattern\}
3- The function AllowsCyclesAfterProducts() is called. If it returns
   a 1 then Next() also returns a 1. Otherwise Next() keeps trying.

Used by GenerateRestrictedAlphaList.cpp

*********************************************************************/

class RestrictedAlphaPermutation:public AlphaPermutation
{
 BasePermutation freepart; // Non fixed part pf the permutation
 unsigned char *freepos, *transl;
 void copyfreepart(void){ register unsigned char i;
  for(i=0; i < freepart.size ; i++) perm[freepos[i]] = transl[freepart.perm[i]];
};
 public:
 // Constructor is the same as AlphaPermutation
 RestrictedAlphaPermutation(unsigned char order):AlphaPermutation(order){};
 // Crucial initialization function
 int InitializeRestrictedAlphaPermutation(int *);
 int Next(){while (freepart.Next()){ copyfreepart(); if(AllowsCyclesAfterProducts()) return 1;} return 0;}
};

/*********************************************************************

             class ProductPermutation

This class implements the products of two permutations read from a
list. Most of the functions are inline and transparent (based in
the elementary functions of the base class BasePermutation).

In this class the crucial element (to understand how it works is the
constructor). The function ReadAlphaList strongly depends on what is
done by the constructor.

The constructor receives the name of the file where the list is as a
paramater and performs the following actions:
1- Tries to open the file for reading as binary file. Upon failure
   issues an error message and stops the program.
2- Tries to read the size hs of the permutations in the file. As
   indicated above, hs must be odd. Upon failure issues an error
   message and stops the program.
   It is set: size = 2*hs;    hsm1 = hs - 1; type = 'N' ;
3- Reads alphalistsize, the number of elements in the list and
   computes tl = alphalistsize*hs;
4- Allocates memory for alphalist to contain the whole list:
alphalist = (unsigned char *) malloc (sizeof (unsigned char) * tl);
   Upon failure issues an error message and stops the program.
   This is to avoid accessing continuously the disc and speeding the
   program.
5- The whole list is then read and loaded into the vector alphalist.
   Upon failure an error message is issued and the program is stopped.
6- Finally it is allocated and initialized the base permutation to
   store the product:
      BasePermutation::AllocateBasePermutation(size);
   and some constants referring to the permutations alpha and beta:
   alpha.size = beta.size = hs; alpha.sm1 = beta.sm1 = hsm1;
   alpha.sm2 = beta.sm2 = hs - 2;

Notice that alpha and beta have not been allocated. Since the list is
stored in a vector The function ReadAlphaList returns a pointer to the
beginning of the desired element of the list, which is copied into
alpha.perm or beta.perm. So, it is not necessary to allocate alpha
and beta. It is only necessary to set appropriately the constants.

The function ReadAlphaList(permnum):
if(permnum >= alphalistsize) this is an error because this
element does not exist in the list (the elements in the list are
numbered from 0 until alphalistsize-1). So it returns a NULL.
Otherwise it computes the position where the desired element starts:
(pointer to beginning of the list) + (desired element) * (length of
an element) = alphalist + permnum * hs

The functions Cross() and Dot() perform the corresponding product of
the current alpha and beta permutations and store it in perm. They
check whether the result is a cycle. In the affirmative they set type
equal to the type of the product ('x'  for cross and '.' for  dot)
and return 1 in the negative type is set to 'N' and return 0.

Used by maxent-1ex-crossprod.cpp, singleprodinfo.cpp,
        maxentcrossprod.cpp, maxentdotprod.cpp

*********************************************************************/

class ProductPermutation:public BasePermutation
{
  unsigned char hs, hsm1; // Size of elements in the list hsm1=hs-1
  unsigned long int alphalistsize; // Length of the list
  BasePermutation alpha, beta;
  char type;
  unsigned char *alphalist; // Stores the whole list
public:
  // We are fixing (and omitting) the initial symbol 0 (or 1 when printing).
  ProductPermutation (char *); // Constructor
  unsigned long int SizeOfAlphaList(void){ return alphalistsize; } ;
  unsigned char Period(void){ return size; } ;
  unsigned char *ReadAlphaList(long permnum)
         { if(permnum >= (long) alphalistsize) return NULL; else return alphalist+(permnum*hs); };
  void Readalpha(long permnum){ alpha.perm = ReadAlphaList(permnum); }
  void Readbeta(long permnum){ beta.perm = ReadAlphaList(permnum); }
  int Cross(void);
  int Dot(void);
  void AlphaPrint (FILE *fp){ fprintf (fp, "Alpha: "); alpha.SimplePrint(fp); }
  void BetaPrint (FILE *fp){ fprintf (fp, "Beta: "); beta.SimplePrint(fp); }
  void PrintAll(FILE *fp){
         alpha.SimplePrint(fp); fprintf (fp," %c ", type);
         beta.SimplePrint(fp);  fprintf (fp," = ");
         SimplePrint(fp); fprintf(fp,"\n");
  }
};

/*********************************************************************

             class AlphaLists1Exc

This is the class that is used to build the X  and Y lists considered
at the end of Section 4. As above, since the period is n=4k+2, the
order of the list is of the form 2k+1 and, hence, it must be odd. half
= (order - 1)/2 = k.

Except the constructor and the function DoesThereExistATCTGCycle(),
which are common to both lists, the rest of procedures need to be
duplicated due to the difference and symmetry of the elements forming
the permutations in the lists.

The algorithm to generate any of the list is as follows (see the
functions  NextY() and NextZ()). It is performed a double "loop". The
initial loop runs on the permutations of length "size" which are
generated by BasePermutation::Next() starting from an appropriate
point (set by the initialization functions InitializeY() and
InitializeZ() --- see the comments to these functions). The base
permutations generated in this way do not belong to A^- and A^+ since
all their elements are either <=p or > p, respectively.

For each of these permutations there is a second "loop" (iterative
process) that runs among all possible positions and tries to insert in
each of these positions all possible elements a > p in the case of A^-
(list Y) and a <= p in the case of A^+ (list Z). This process is here
called the "perturbation process" and is controlled with the help of
the variables
     fixedpos, perturbp, basepermperturbp and  bdry.
It is initialized (it means that these variables are initialized) by
the functions InitializeYPerturbation() and InitializeZPerturbation()
and a step in this process is made by the functions
NextYPerturbation() and NextZPerturbation(). Both of these functions
return a 0 when is not possible to go further in the perturbation
process for the given base permutation and so it is necessary to go
to to the next base permutation. By Lemma 4.5 the elements of the list
Y must verify properties (i-iv.m) of Lemma 4.3. By Lemma 4.9 the
elements of the list Z  must verify properties (I-III.M) of Lemma 4.7.

This description is clearly implemented (and seen) in the inline
functions NextY() and NextZ().


Used by GenerateYZLists1Except.cpp

*********************************************************************/

class AlphaLists1Exc:public BasePermutation
{
  // half=sm1/2; tsm1= 2*size-1; tsm2=tsm1-1; tsm3=tsm2-1;
  unsigned char half, tsm1, tsm2, tsm3;
  // Auxiliary variables for perturbing the exceptions
  unsigned char fixedpos, perturbp, basepermperturbp, bdry;
  /* Implements Lemma 3.7 of the paper.
     Used for both lists to discard the permutations that cannot give
     a cycle for the cross product */
  int DoesThereExistATCTGCycle(void);
  // The initialization functions of the perturbation process
  void InitializeYPerturbation(void);
  void InitializeZPerturbation(void);
  // Do next step in the perturbation process
  int NextYPerturbation(void);
  int NextZPerturbation(void);
public:
  AlphaLists1Exc(unsigned char); // Constructor
  // Initialization functions for the lists: sets the first element
  void InitializeY(void);
  void InitializeZ(void);
  int NextY(void){
    while(!NextYPerturbation()){ if(!BasePermutation::Next()){ return 0; } InitializeYPerturbation(); }
    return 1;
  }
  int NextZ(void){
    while(!NextZPerturbation()){ if(!BasePermutation::Next()){ return 0; } InitializeZPerturbation(); }
    return 1;
  }
};
