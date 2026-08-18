# Sound redundant propagation cuts for the B2[2]/Z100 SAT encoding

The optimized proof-CNF generator adds only constraints that follow from the original ordered-difference condition. They improve propagation but do not change the set of valid 14-element candidates.

## Parity cut

Let `e` be the number of selected even residues and `14-e` the number of selected odd residues. An ordered difference is even exactly when its endpoints have the same parity. Therefore the number of ordered even differences is

```text
e(e-1) + (14-e)(13-e).
```

There are 49 nonzero even residues modulo 100, each with capacity 2, hence the expression is at most 98. Direct evaluation gives a violation for `e<=4` and `e>=10`, while `e=5,...,9` are the only possible counts. Thus every valid 14-set satisfies

```text
5 <= e <= 9.
```

## Quotient-capacity cuts

Fix a divisor `m` of 100 and reduce differences modulo `m`.

### Quotient difference zero

An unordered selected pair in the same residue class modulo `m` contributes two ordered differences whose values modulo 100 are nonzero multiples of `m`. There are `100/m - 1` such residues, each of ordered capacity 2. Since each unordered pair consumes two units of this total capacity, the number of same-class unordered pairs is at most

```text
100/m - 1.
```

### Nonzero quotient distance `q`, not self-inverse

For `1 <= q < m/2`, an unordered pair whose quotient difference is `+q` or `-q` contributes one ordered difference to each direction. There are `100/m` residue preimages of `q` modulo 100, each with capacity 2. Therefore the number of such unordered pairs is at most

```text
2*(100/m).
```

### Self-inverse quotient distance `m/2`

When `m` is even and `q=m/2`, both orientations of an unordered pair project to the same quotient difference. The ordered capacity is `2*(100/m)`, and each pair consumes two units. Hence at most

```text
100/m
```

such unordered pairs are possible.

The generator applies these cuts for

```text
m in {2,4,5,10,20,25,50}.
```

All cuts are necessary consequences of the original `r_A(d)<=2` constraints. Consequently an UNSAT certificate for the strengthened CNF remains a valid nonexistence certificate for the original canonical branch, provided the exact encoding and proof object are preserved and independently checked.
