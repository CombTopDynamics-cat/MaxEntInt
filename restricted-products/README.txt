********* GenerateRestrictedAlphaList.cpp **************
---------------------------------------------------------------------
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
---------------------------------------------------------------------
