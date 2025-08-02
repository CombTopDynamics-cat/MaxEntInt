/************************************************************************
  This file contains the code for the functions of the following classes
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
  maxentdotprod.cpp, ReadAList.cpp, GenerateFreeAList.cpp,
  GenerateRestrictedAlphaList.cpp
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "../BaseFunctions/permutations.h"
#include "../BaseFunctions/maxent.h"
#define DIR(x,y) ((x) < (y) ?  +1 : ((x) > (y) ? -1 : 0))

/****************************************************************************
                      class BasePermutation

perm is understood (as in the whole file unless otherwise said) in MAP FORM.
That is perm[i] denotes the image of i by the permutation perm (see the
first two lines of Section 3 of the paper for the notation).

Also for programming reasons i ranges from 0 to size-1
and perm[i] belongs to {0,1,...,size-1}.
****************************************************************************/

/*
- checks whether the length is allowed.
- sets size, sm1 and sm2
- allocates the corresponding memory to it
- initializes perm to the smallest permutation in the lexicographical ordering:
  (0,0,....,0)
*/
void BasePermutation::AllocateBasePermutation(unsigned char order)
{
  register unsigned char i;
  if (order < 2 || order > MAXPERM) ErrorAborting("Permutation length not allowed!!");
  if(perm) free(perm);
  size = order; sm1 = size - 1; sm2 = size - 2;
  perm = (unsigned char *) malloc (sizeof (unsigned char) * size);
  if (!perm) ErrorAborting("memory problems initializing the permutations vector");
  for(i=0; i<size; i++) perm[i] = 0;
  return;
};

int BasePermutation::DefineManual(int order, int * map)
{  register unsigned char i;

   BasePermutation::AllocateBasePermutation((unsigned char) order);
   for(i=0; i<size; i++) { if( map[i] < 1 || map [i] > order) { return 0; } perm[i] = (unsigned char) (map[i]-1); map[i] = 0; }
   for(i=0; i<size; i++) {if(map[perm[i]]) return 0; else map[perm[i]] = 1; }
   return 1;
}

/*
THIS IS ONE OF THE CRUCIAL FUNCTIONS.

Given the permutation perm it computes the next one in the
lexicographical ordering until the last one is found. In such case a
fail return code of 0 is issued.

It is an implementation of the function
    void PermLexSuccessor(int n ,permutation pi, int * flag)
(Algorithm 2.14)
which replaces a permutation by its successor
from the book
Combinatorial Algorithms:Generation, Enumeration & Search
by Donald L. Kreher and Douglas Stinson
See Chapter 2: Generating Elementary Combinatorial Objects
See Section 2.4: Permutations (2.4.1 Lexicographic ordering)

The changes with respect to the original function are due to
- programming style :-)
- and, mainly, to the fact that, here, the permutation vector
  is written in MAP FORM.
*/
int BasePermutation::Next ()
{
  unsigned char i, s = sm2, a;

  while (perm[s] >= perm[s + 1]) { if (!s) return 0; s--; }
  for (i = s + 2; i < size; i++) { if (perm[i] <= perm[s]) break; }
  i--;
  a = perm[s];  perm[s] = perm[i];  perm[i] = a;
  for (s++, i = sm1; s < i; i--, s++) { a = perm[i]; perm[i] = perm[s]; perm[s] = a; }
  return 1;
}

/* Modified on 20080929
It checks whether the map vector is a cycle by the straightforward method:
that is it must be a cycle of length size

Fact 1: * If 0 id periodic then, after the while loop, r will be equals to it period.
        * If 0 is pre-periodic (that is, perm[i] \ne 0 for every i in the orbit of 0)
          then, after the while loop, r will be equals to size+1.

Fact 2 : For any map $\theta: \{0,1,....,size-1\} \longmapsto \{0,1,....,size-1\}$,
         the orbit of any point is finite. Thus, every element from $\{0,1,....,size-1\}$
         is either periodic or pre-periodic.

Fact 3: If a map \theta is not bijective then
           $\theta(\{0,1,....,size-1\}) \subsetneq \{0,1,....,size-1\}$. Consequently,
        if an element of $\{0,1,....,size-1\}$ is periodic then it has period smaller
        than size.

Summary: If a map $\theta: \{0,1,....,size-1\} \longmapsto \{0,1,....,size-1\}$, is not
         a permutation then either 0 is pre-periodic (and the while loop of the function
         IsCycle ends with r=size+1) or it is periodic (and the while loop of the function
         IsCycle ends with r<size). In both cases the function returns a 0 (IsCycle = fail).

Fact 4: If a map $\theta: \{0,1,....,size-1\} \longmapsto \{0,1,....,size-1\}$, is a
        permutation then every point is periodic. When $\theta$ is cyclic every point
        belongs to the same periodic orbit which has period size. If $\theta$ is not
        cyclic then there is more than one periodic orbit and, consequently, the period
        of each periodic orbit is smaller than size.

Summary: If $\theta$ is not cyclic then 0 is periodic of period smaller than size. Then
         the while loop of the function IsCycle ends with r<size and the function returns
         a 0 (IsCycle = fail). When $\theta$ is cyclic then 0 is periodic of period size.
         Then the while loop of the function IsCycle ends with r=size and the function
         returns a 1 (IsCycle = true).

So, the next function returns true when the map is a cyclic permutation. Otherwise (not
a permutation or not cyclic) returns false.
*/
int BasePermutation::IsCycle (void)
{ unsigned char i=0, r=1;
  while ((i=perm[i]) != 0 && r <= size ) {r++;}
  if( r != size ) return 0; else return 1;
}

/*
It checks whether the permutation is maximodal by the straightforward method:
that is it determines whether it is increasing or decreasing at the beginning
and it checks for a change of slope at each interval.

NOTE: It assumes that perm[i] != perm[i+1] for all i; which definitely happens
      if the permutation is a cycle.
*/
int BasePermutation::IsMaximodal (void)
{
  register unsigned char i;
  int signo = DIR(perm[0],perm[1]), signn;

  for (i = 1; i < sm1; i++)
  {
    signn=DIR(perm[i],perm[i+1]);
    if (signo == signn) return 0;
    signo = signn;
  }
  return 1;
 }

/*
   Smaller refers to the lexicographical ordering.
*/
int BasePermutation::TestIfDualIsSmaller (void)
{ register unsigned char i;
  for(i=0; i < size; i++){ unsigned char c; c = sm1 - perm[sm1-i];
   if(c > perm[i]) return 0;  if(c < perm[i]) return 1;
  } return 0;
}

/*
For a Maximodal cycle computes the number of maximum values which are less than minimum values
The formula depends on the parity of size.
*/
int BasePermutation::ComputeNumberOfExceptionsIfMaximodal(void)
{ register unsigned char i; unsigned char startmin = 0, separation = (sm2)/2, count = 0;

  if(DIR(perm[0],perm[1]) < 0) startmin = 1;
  else if (size%2) separation++;

  for(i=startmin; i < size; i+=2) if(perm[i] <= separation) count++;
  return separation-count+1;
}

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

/*
- Auxiliary vector busy for ChooseLegalTryForAPosition
- Auxiliary function SetBusy to fill the vector busy.
*/
 unsigned char busy[MAXPERM+1];
 void SetBusy(unsigned char val, unsigned char len) {
   unsigned char *bb, *bbf = busy+len;
   for(bb=busy; bb < bbf; bb++) *bb = val;
}

/*
CONSTRUCTOR.
The cycle implemented by this constructor is, of course, maximodal and cyclic and it
is the smallest in the lexicographical ordering. Also, it has one exception.
*/
MaximodalCiclicIntervalPermutation::MaximodalCiclicIntervalPermutation(unsigned char order)
 {  register unsigned char i;
    if (order < 4) ErrorAborting("permutation length not allowed in maximodal!!");
    BasePermutation::AllocateBasePermutation(order);
    perm[0] = 1; perm[1] = 3; perm[2] = 0;
    for (i = 3; i < size; i++) perm[i] = perm[i-2] + 2;
    if(size%2) perm[sm2] = sm1; else perm[sm1] = sm2;
 }

/*
FUNCTION ChooseLegalTryForAPosition
It is the base on which the function Next() is built
*/
int MaximodalCiclicIntervalPermutation::ChooseLegalTryForAPosition(
          unsigned char pos, unsigned char keeppos, unsigned char *busy)
{  register unsigned char i = 0;
   unsigned char ini = perm[pos-1]+1, end=sm1;
   int dir = DIR(perm[pos-2],perm[pos-1]);


   if(pos == 0){ i = perm[0] + 1; if(i < sm1) return i; else return -1;}
/*
   If keeppos = 0 then i=0. It does not coincide with perm[0] since it
   it is initialized to 1 and afterwards incremented.
   Cyclicity: By initialization, when perm[0]= 1; perm[1]>=3. Otherwise,
   perm[0]>=2. So, it cannot happen that 0 -> 1 -> 0. The only cyclicity
   that can happen is 1 -> 1.
   Obviously, perm[0] must be different from perm[1].
   Observe that pos = 1. So, when i = pos = 1, then perm[0]>=2 > i. Otherwise,
   if i = perm[0] then i = perm[0]>=2>1 = pos.
*/
   else if(pos == 1){
     if(keeppos){i = perm[1]+1; if(i == pos) i++; if(i == perm[0]) i++;}
     if(i < size) return i; else return -1;
   } else if (pos == sm1){
     for(i=0; i < size && busy[i]; i++){} if(i == size)
         ErrorAborting("unexpected error in MaximodalCiclicIntervalPermutation::ChooseLegalTryForAPosition");
/*
      This 'not busy' i MUST be the only one since pos = size-1,
      all entries of perm are between 0 and size-1 and appear only once.
      In this case it is not necessary to check the cyclicity.
      Since this  condition holds for all pos from 0 to size-2 it
      must hold also for pos = size-1;
*/
      if(dir > 0 && i < perm[pos-1] || dir < 0 && i > perm[pos-1]) return i; else return -1;
   } else {
       if(dir > 0) {  ini = 0; end = perm[pos-1]-1; }
       if(keeppos) ini = perm[pos] + 1;

      for(i=ini; i <= end; i++){ unsigned char p;
          if(busy[i]) continue; // Not good; point already occupied
          p = i; while(p < pos) p = perm[p];
          if(p == pos) continue; // Not good; the choice perm[pos] = i generates a cycle.
          return i; // GOOD. perm[pos] = i is free and does not generate a cycle
       }
       return -1;
   }
}


int MaximodalCiclicIntervalPermutation::Next ()
{  register unsigned char pos=size - 3;
   unsigned char keppos = 1;
   int t;

   SetBusy(1,size);
   busy[perm[sm1]] = 0; busy[perm[sm2]] = 0; busy[perm[pos]] = 0;// table of assigned points

   for(;;){
      t = ChooseLegalTryForAPosition(pos, keppos, busy);
      if(t >= 0) { perm[pos] = t; pos++; if(pos == size) return 1;// Legal Try found
        busy[t] = 1; keppos = 0;
      } else { if(!pos) return 0;
        pos--; keppos = 1; busy[perm[pos]] = 0;
      }
   }
}

/*********************************************************************

             class AlphaPermutation

This is the class that deals with the A-lists (and derived ones)
Then, since the period n is of the form 4k+2, the order of the list is
of the form 2k+1 and, hence, it must be odd. half = (order - 1)/2 = k.

Used by AlphaPermutation

*********************************************************************/
/*
CONSTRUCTOR
- Checks that order is odd
- allocates the base permutation of "size" length
- sets half=sm1/2;
- initializes the permutation to the smallest possible element that
  verifies Corollary 3.10
- For security checks this fact.
*/
AlphaPermutation::AlphaPermutation(unsigned char order)
 {  register unsigned char i;

    if (!(order%2)) ErrorAborting("AlphaPermutation length must be odd!!");
    BasePermutation::AllocateBasePermutation(order); half = sm1/2;
    perm[0] = 1; perm[1] = 0; for (i = 2; i < size; i++) perm[i] = i;
    if(!AllowsCyclesAfterProducts())
          ErrorAborting("unexpected error in AlphaPermutation::AlphaPermutation initialization!!!");
}

/*
The next function implements Corollary 3.10 of the paper. That is, it
checks  for permutations alpha that when participate as a first
element of a product (cross or dot; no matter) will give a cycle.

It returns 0 when the map phi has a cycle and thus no product
starting by alpha can be a cycle. Otherwise it returns 1.

If we rewrite the permutation i -> perm[i] as
    phi(t) = perm[(t-1)/2] + 1 (t=1,3,5,...,2*size-1)
we want that the map phi has no cycles.

1st Obs.: phi(2*i+1) = perm[i] + 1 <= size < 2*i+1 for i > half.
Thus, there is no cycle of phi for t = 2*i+1 with i > half.

2nd Obs.: If phi has a cycle, then there is t such that
   phi(t) \in \{1,3,5,...,2*size -1}.
This means that if phi(t) is even, then t does not belong to a cycle of
phi. That is, any cycle of phi is contained in
   {2*i+1: 0 <=  i  <= half and phi(2*i+1) is odd}.

3rd Obs.: Consider i \in {0,1,...,half} and let
          psi(i) = 'nothing' if perm[i] is odd, and
          psi(i) = perm[i]/2  if perm[i] is even.

Then, if phi(2*i+1) is odd, phi(2*i+1) = 2*psi(i)+1.
In other words, phi has a cycle (on the set
   {2*i+1: 0 <=  i  <= half and phi(2*i+1) is odd})
if and only if psi has a cycle on the set
  {i: 0 <= i <= half and perm[i] is even (that is psi(i) is defined)}.

4rth Obs.: The map psi on the set
  {i: 0 <= i <= half and perm[i] is even (that is psi(i) is defined)}
has only two types of orbits: cycles and orbits ending at 'nothing'.
Two different orbits are disjoint. That is, there are no pre-periodic
points in psi.

This is the condition that the permutation generates cycles for the
cross product.

Note: when p < i the orbit of p has been controlled before.
      If it were bad then we would not be in this routine anymore.
      Since we are still here the orbit of p was good.
      Thus; the orbit of i is also good.

[ For the dot product ]
The condition for the dot product is that the cycle must not contain
'size'. So, we can exclude i=half and perm[i] + 1 = size.

What a small function  fur such a big explanation!! :-)
*/

int AlphaPermutation::AllowsCyclesAfterProducts()
{  register unsigned char i;

   for(i=0; i< half; i++){ unsigned char p=i;
      while( !(perm[p]%2) && p >= i && perm[p] < sm1) { p = perm[p]/2; if(p==i) return 0; }
   }
   return 1;
}

/*********************************************************************

             class RestrictedAlphaPermutation

This is the class that is used to build the restricted A^*-lists
considered at the end of Section 6. Since it is a derived class of
AlphaPermutation can use its features. In particular the ability of
checking Corollary 3.10 through AllowsCyclesAfterProducts().

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


int RestrictedAlphaPermutation::InitializeRestrictedAlphaPermutation(int * map){
    register unsigned char i,p=0; unsigned char n=size;

    for(i=0; i<size; i++) perm[i] = 0;
    for(i=0; i<size; i++) { if( map[i] < 0 || (unsigned char) map[i] > size) return 0;
      if(map[i]){ if(perm[map[i]-1]) return 0; n--; perm[map[i]-1] = 1; }
    } if(n == size) ErrorAborting("this is a free alpha list. You better use the 'Free' program");

    freepos = (unsigned char *) malloc (sizeof (unsigned char) * n);
    transl = (unsigned char *) malloc (sizeof (unsigned char) * n);
    if (!freepos || !transl) ErrorAborting("memory problems initializing the pattern vector");
    for(i=0, p=0; i<size; i++) if(!perm[i]) { transl[p] = i; p++; }
    for(i=0, p=0; i<size; i++) if (map[i]) { perm[i] = map[i] - 1; } else { freepos[p] = i; p++; }

    freepart.__BaseSTDPermutation(n); copyfreepart();
    if(AllowsCyclesAfterProducts()) return 1;
    return Next();
}

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

 ProductPermutation::ProductPermutation (char *filename){
   unsigned long int tl;

   FILE *f = sfopen(filename, "rb", "when opening associated List");
   fread(&hs,  sizeof(unsigned char),  1, f); size = 2*hs;    hsm1 = hs - 1; type = 'N' ;
   if (!(hs%2)) ErrorAborting("ProductPermutation period must have odd half!!");
   fread(&alphalistsize,  sizeof(unsigned long int),  1,f); tl = alphalistsize*hs;

   alphalist = (unsigned char *) malloc (sizeof (unsigned char) * tl);
   if (!alphalist) ErrorAborting("memory problems initializing the memory alphalist vector");

   if(tl != fread(alphalist,  sizeof(unsigned char),  tl ,f)) ErrorAborting("reading the alphalist");
   fclose(f);

   BasePermutation::AllocateBasePermutation(size);
   alpha.size = beta.size = hs; alpha.sm1 = beta.sm1 = hsm1; alpha.sm2 = beta.sm2 = hs - 2;
}

int ProductPermutation::Cross(void)
{
  register unsigned char i;
  for(i=0; i< hs; i++) { unsigned char ii=2*i;
   perm[ii] = alpha.perm[i]; perm[ii+1] = sm1 - beta.perm[hsm1-i] ;
  }
  if(IsCycle()) {type = 'x'; return 1; } else { type = 'N'; return 0; }
}

int ProductPermutation::Dot(void)
{
  register unsigned char i;
  for(i=0; i< hs; i++) { unsigned char ii=2*i, aux;
   aux = sm2 - alpha.perm[hsm1-i]; perm[ii] = ((aux > hsm1) ? aux : sm1) ;
   aux = beta.perm[i]+1; perm[ii+1] = ((aux <= hsm1) ? aux : 0) ;
  }
  if(IsCycle()) {type = '.'; return 1; } else { type = 'N'; return 0; }
}

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

/*
CONSTRUCTOR
- Checks that order is odd
- allocates the base permutation of "size" length
- sets half=sm1/2; tsm1= 2*size-1; tsm2=tsm1-1; tsm3=tsm2-1;
*/
AlphaLists1Exc::AlphaLists1Exc(unsigned char order)
 {  if (!(order%2)) ErrorAborting("AlphaPermutation length must be odd!!");
    BasePermutation::AllocateBasePermutation(order);
    half = sm1/2; tsm1 = 2*size - 1; tsm2 = tsm1 - 1; tsm3 = tsm2 - 1;
}

/*
INITIALIZATION FUNCTION FOR THE Y LIST: SETS THE FIRST ELEMENT
1.- Initializes the permutation to the smallest possible element
    that has all elements smallest than size. Also, according to
    Lemma 4.3(iv.m) the first element must be different from 0
2.- Initializes the perturbation process calling
    InitializeYPerturbation()
3.- The current permutation is not valid (it has 0 exceptions and
    could not pass the test DoesThereExistATCTGCycle()). So, to obtain
    the first valid permutation in the list one hast to call the
    function NextY(). If it fails this means a consistency error in
    the class.
*/
void AlphaLists1Exc::InitializeY(void)
{ unsigned char i;
  perm[0] = 1; perm[1] = 0; for (i = 2; i < size; i++) { perm[i] = i; } InitializeYPerturbation();
  if(!NextY()) ErrorAborting("unexpected error in AlphaLists1Exc::InitializeY.\nNothing to do!!!");
}

/*
INITIALIZATION FUNCTION FOR THE Z LIST: SETS THE FIRST ELEMENT
1.- Initializes the permutation to the smallest possible element
    that has all elements >= size. Also, according to Lemma
    4.7(II.M) there is one element equals to tsm1=2*size-1
2.- Initializes the perturbation process calling
    InitializeZPerturbation()
3.- The current permutation is not valid (it has 0 exceptions and
    could not pass the test DoesThereExistATCTGCycle()). So, to obtain
    the first valid permutation in the list one hast to call the
    function NextZ(). If it fails this means a consistency error in
    the class.
*/
void AlphaLists1Exc::InitializeZ(void)
{ unsigned char i;
  for (i = 0; i < size; i++) { perm[i] = size+i; } InitializeZPerturbation();
  if(!NextZ()) ErrorAborting("unexpected error in AlphaLists1Exc::InitializeZ.\nNothing to do!!!");
}

/*
The next function implements Lemma 3.7 of the paper. That is, it
checks  for permutations alpha that when participate as a first
element of a cross product will give a cycle.

It returns 1 when the map phi has a cycle and thus no product
starting by alpha can be a cycle. Otherwise it returns 0.

Compare with the above function
   int AlphaPermutation::AllowsCyclesAfterProducts()
(and see the comments to it).
*/
int AlphaLists1Exc::DoesThereExistATCTGCycle(void)
{  register unsigned char i;
   for(i=0; i< size; i++){ unsigned char p=i;
      while(!(perm[p]%2) && p >= i) { p = perm[p]/2; if(p==i) return 1; }
   }
   return 0;
}

/*
INITIALIZATION FUNCTION FOR THE PERTURBATION PROCESS FOR THE Y LIST

 perturbp = 0; is the element that we are currently perturbing
 bdry=tsm2; is the maximum until which we can increase perm[perturbp]
 basepermperturbp=perm[0]; is the original value (to be restored when
            we modify another position of this permutation and hence
            increasing perturbp
perm[0] = sm1; initial (illegal) value of the perturbed position.
               It will be made a legal perturbation when it will be
               incremented by the first call of NextYPerturbation()

fixedpos is the position of 0 (1 in the written notation). According
to Lemma 4.3(ii.m) this position cannot be perturbed.

According to Lemma 4.3(1,3.m) bdry <= tsm3 except when perturb=0 that
bdry can be tsm2 (precisely what is established by this initializer).
*/

void AlphaLists1Exc::InitializeYPerturbation(void)
{
  perturbp = 0; bdry=tsm2; basepermperturbp=perm[0]; perm[0] = sm1;
  for(fixedpos=1; fixedpos < size; fixedpos++) if(!perm[fixedpos]) break;
}

/*
NEXT STEP IN THE PERTURBATION PROCESS FOR THE Y LIST

This function implements the perturbation process taking into account
Lemma 4.3(1-iv.m).

Each time that perm[perturbp] > bdry the original value
(basepermperturbp) is returned to this position. Then perturbp is
increased (unless perturbp == fixedpos that is bi-increased to leave
fixedpos unperturbed), bdry is fixed to tsm3 (Lemma 4.3(1,3.m)).

If perturbp > sm1 this means that there are no more positions to
perturb and the perturbation process is finished. Then a value of
zero is returned and a new base permutation is obtained.

Otherwise, perm[perturbp] is stored in basepermperturbp to be recovered
later and perm[perturbp] is initially perturbed to size. Later on we
will only have to do perm[perturbp]++. In each case we perform the
test DoesThereExistATCTGCycle() to check whether we can discard the
current perturbed permutation. In this case we continue the
perturbation process.

If the current perturbed permutation is accepted a return value 1 is
issued.
*/

int AlphaLists1Exc::NextYPerturbation(void)
{  do { perm[perturbp]++;
     if(perm[perturbp] > bdry) { perm[perturbp] = basepermperturbp;
       perturbp++; bdry=tsm3; if(perturbp == fixedpos) perturbp++;
       if(perturbp > sm1) return 0;
       basepermperturbp=perm[perturbp]; perm[perturbp] = size;
     }
   } while(DoesThereExistATCTGCycle());
   return 1;
}


/*
INITIALIZATION FUNCTION FOR THE PERTURBATION PROCESS FOR THE Z LIST

fixedpos is the position of tsm1=2*size-1 (n in the written notation).
According to Lemma 4.7(II.M) this position cannot be perturbed.

perturbp = (fixedpos ? 0 : 1) fixed according to fixedpos is the
first available position for perturbing

basepermperturbp=perm[perturbp]; is the original value (to be restored
            when we modify another position of this permutation and
            hence increasing perturbp

perm[perturbp] = size; initial (illegal) value of the perturbed
               position.
               It will be made a legal perturbation when it will be
               decremented by the first call of NextZPerturbation()
*/
void AlphaLists1Exc::InitializeZPerturbation(void)
{
  for(fixedpos=0; fixedpos < size; fixedpos++) if(perm[fixedpos] == tsm1) break;
  perturbp = (fixedpos ? 0 : 1);
  basepermperturbp=perm[perturbp]; perm[perturbp] = size;
}

/*
NEXT STEP IN THE PERTURBATION PROCESS FOR THE Z LIST

This function implements the perturbation process taking into account
Lemma 4.7(I-III.M). In particular by Lemma 4.7(I,III.M),
perm[perturbp] >= 2.

Thus, once perm[perturbp] < 2, the original value (basepermperturbp)
is returned to this position. Then perturbp is increased (unless
perturbp == fixedpos that is bi-increased to leave fixedpos
unperturbed).

If perturbp > sm1 this means that there are no more positions to
perturb and the perturbation process is finished. Then a value of
zero is returned and a new base permutation is obtained.

Otherwise, perm[perturbp] is stored in basepermperturbp to be recovered
later and perm[perturbp] is initially perturbed to sm1. Later on we
will only have to do perm[perturbp]--. In each case we perform the
test DoesThereExistATCTGCycle() to check whether we can discard the
current perturbed permutation. In this case we continue the
perturbation process.

If the current perturbed permutation is accepted a return value 1 is
issued.
*/

int AlphaLists1Exc::NextZPerturbation(void)
{  do { perm[perturbp]--;
     if(perm[perturbp] < 2) { perm[perturbp] = basepermperturbp;
       perturbp++; if(perturbp == fixedpos) perturbp++;
       if(perturbp > sm1) return 0;
       basepermperturbp=perm[perturbp]; perm[perturbp] = sm1;
     }
   } while(DoesThereExistATCTGCycle());
   return 1;
}
