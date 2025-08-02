********* GenerateFreeAList.cpp **************
---------------------------------------------------------------------
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
---------------------------------------------------------------------

********* ReadAList.cpp **************
---------------------------------------------------------------------
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
---------------------------------------------------------------------

********* maxentcrossprod.cpp **************
---------------------------------------------------------------------
Computes the maximum entropy of the cross-product permutations
a <cross> b of the form descibed below WHICH ARE CYCLES (Section 3):

On input a file that must be in a valid A-list format must be
specified. Optionally, a segment of the list can be specified
(the default segment is the whole list).

Then the program generates all products a <cross> b with:
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
This is a device to increase the speed of the program in critical cases when
we know in advance that the entropy must be larger than an a priori value.

Usage: filename entropy-filter-level [ini:end]
---------------------------------------------------------------------

********* maxentdotprod.cpp **************
---------------------------------------------------------------------
Computes the maximum entropy of the dot-product permutations
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
---------------------------------------------------------------------


********* singleprodinfo.cpp **************
---------------------------------------------------------------------
This program is necessary because maxentcrossprod.exe and
maxentdotprod.exe do not display enough information on the
cycles of maximum entropy. They only display the maximum
spectral radius and the position of a and b in the list.
They do not display the actual permutation of those elements
the permutation of the resulting cycle and the Markov matrix.

On input a file that must be in a valid A-list format (Section 3)
must be specified, two valid indices a and b in the file, and one
type of product:
x for cross
. for dot

The program tests that both products of these two elements are cycles
and maximodal. In that case it performs the chosen product and displays
the actual permutations corresponding to a and b, the product cycle in
map form, the topological entropy and the Markov matrix with all row
and colum sums, together with the total sum.

Usage: filename alphaindex producttype(x,.) betaindex
---------------------------------------------------------------------

