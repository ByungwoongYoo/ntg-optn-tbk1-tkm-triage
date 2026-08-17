# Reduction of the 14-element B_2[2] problem in Z_100 to the normalized gap-1 case

Let A be a 14-element subset of Z/100Z satisfying r_A(t) <= 2 for every nonzero t.

## Lemma: A contains a unit difference

Assume for contradiction that every difference a-b (a != b) is a nonunit modulo 100. Because 100=2^2*5^2, every nonunit difference is divisible by 2 or by 5.

Map every element of A to its residue modulo 10. Write a residue as the pair

    (parity, residue mod 5) in Z/2Z x Z/5Z.

For any two occupied residue classes, their difference is divisible by 2 or 5, so the two pairs agree in the first coordinate or in the second coordinate.

A clique in the rook graph K_2 square K_5 is contained in one row or one column: if two vertices of a clique are in distinct rows, they must share their column; every third clique vertex must then share that same column. Hence all residues of A have the same parity, or all have the same residue modulo 5.

If all elements have the same parity, all 14*13=182 ordered nonzero differences lie among the 49 nonzero even residues. Since each residue may occur at most twice, their total capacity is 2*49=98, contradiction.

If all elements have the same residue modulo 5, all 182 ordered nonzero differences lie among the 19 nonzero multiples of 5. Their total capacity is only 2*19=38, again a contradiction.

Therefore A contains a difference u with gcd(u,100)=1.

## Corollary: only normalized gap 1 needs to be solved

Translate A so that one endpoint of a unit-difference pair is 0, then multiply all residues by u^{-1} modulo 100. Translation and multiplication by a unit are automorphisms of Z/100Z and merely permute the nonzero difference residues, so the B_2[2] property is preserved.

The transformed set contains 0 and 1. Consequently its minimum positive cyclic gap is exactly 1.

Thus existence of a 14-element B_2[2] subset of Z_100 is equivalent to existence in the normalized minimum-gap-1 case. The g=2,...,7 computations are useful independent checks but are not needed for a global exclusion once the lemma is accepted.

## Claim boundary

This note gives the mathematical reduction only. It does not by itself prove that the gap-1 case is unsatisfiable. A positive gap-1 witness or a complete independently checkable gap-1 UNSAT certificate is still required for a decisive answer.
