# Ulrich's `u4` is not a single axiom for positive implication

## Statement

Let

\[
\Phi=((x\to y)\to z)\to((y\to(z\to u))\to(y\to u)).
\]

Equivalently, with

\[
C(A,B,R):=(A\to(B\to R))\to(A\to R),
\]

we have

\[
\Phi=((x\to y)\to z)\to C(y,z,u).
\]

A proof from substitution instances of \(\Phi\), using modus ponens as the only inference rule, is represented by an applicative term whose leaves are independently polymorphic copies of a constant of principal type \(\Phi\). The theorem below proves that no such term has a reflexive type \(P\to P\).

## Right-comb principal types

Let \(R_1=\Phi\) and \(R_{n+1}=\Phi\,R_n\), where juxtaposition denotes application.

Define principal types \(T_n\) as follows. First,

\[
T_1=\Phi.
\]

Choose pairwise distinct type variables \(a,b,c,r_2,r_3,\ldots\), and set

\[
A_2=a,\qquad B_2=C(b,a,c),\qquad T_2=C(A_2,B_2,r_2).
\]

For \(n\ge2\), recursively define

\[
A_{n+1}=B_n\to r_n,\qquad
B_{n+1}=A_n\to r_n,\qquad
T_{n+1}=C(A_{n+1},B_{n+1},r_{n+1}).
\]

### Lemma 1

For every \(n\ge1\), \(T_n\) is the principal type of \(R_n\).

### Proof

For \(n=1\) this is the definition. Applying a fresh copy of \(\Phi\) to a term of type \(T_n\) requires unifying the first argument type \((x\to y)\to z\) of \(\Phi\) with \(T_n\).

For \(n=1\), this yields

\[
y=A_2=a,\qquad z=B_2=C(b,a,c),
\]

up to renaming, so the result is \(T_2\).

For \(n\ge2\), write

\[
T_n=(A_n\to(B_n\to r_n))\to(A_n\to r_n).
\]

The most general unifier is

\[
x=A_n,\qquad y=B_n\to r_n,\qquad z=A_n\to r_n.
\]

The result type of the application is therefore

\[
C(B_n\to r_n,A_n\to r_n,r_{n+1})=T_{n+1}.
\]

This proves the claim. ∎

## The interlocking lemma

For finite simple types \(P,Q\), write \(E(P,Q)\) for the following property:

> There are no substitutions \(\sigma\) and no finite type \(W\) such that
> \[
> \sigma(P)=\sigma(Q)\to W.
> \]

### Lemma 2

For every \(n\ge2\), both \(E(A_n,B_n)\) and \(E(B_n,A_n)\) hold.

### Proof

At \(n=2\), we have \(A_2=a\) and

\[
B_2=C(b,a,c)=((b\to(a\to c))\to(b\to c)).
\]

The variable \(a\) occurs properly inside \(B_2\). Hence \(\sigma(a)\) cannot equal \(\sigma(B_2)\to W\), since that would make the finite tree \(\sigma(a)\) contain itself as a proper subterm. Thus \(E(A_2,B_2)\).

Conversely, if \(\sigma(B_2)=\sigma(a)\to W\), injectivity of the arrow constructor forces

\[
\sigma(b\to(a\to c))=\sigma(a).
\]

The right-hand side would again contain itself as a proper subterm. Thus \(E(B_2,A_2)\).

Now assume both directions for \((A_n,B_n)\). Since

\[
A_{n+1}=B_n\to r_n,\qquad B_{n+1}=A_n\to r_n,
\]

an equality

\[
\sigma(A_{n+1})=\sigma(B_{n+1})\to W
\]

would, by injectivity of \(\to\), imply

\[
\sigma(B_n)=\sigma(A_n)\to\sigma(r_n),
\]

contradicting \(E(B_n,A_n)\). The converse direction is symmetric and contradicts \(E(A_n,B_n)\). Induction completes the proof. ∎

## No derived right comb can be used as a major premise

The antecedent of \(T_m\), for \(m\ge2\), is

\[
D_m=A_m\to(B_m\to r_m).
\]

### Lemma 3

For every \(m\ge2\) and every \(n\ge1\), \(D_m\) and \(T_n\) have no common substitution instance.

### Proof: the case \(n\ge2\)

Suppose substitutions \(\sigma,\tau\) made \(\sigma(D_m)=\tau(T_n)\). Expanding both sides and using injectivity of \(\to\) gives

\[
\sigma(A_m)=\tau(A_n)\to(\tau(B_n)\to\tau(r_n)),
\]

and

\[
\sigma(B_m)=\tau(A_n).
\]

Consequently,

\[
\sigma(A_m)=\sigma(B_m)\to(\tau(B_n)\to\tau(r_n)),
\]

contradicting \(E(A_m,B_m)\) from Lemma 2.

### Proof: the case \(n=1, m\ge3\)

For \(m\ge3\), both \(A_m\) and \(B_m\) end in the same fresh codomain \(r_{m-1}\):

\[
A_m=B_{m-1}\to r_{m-1},\qquad B_m=A_{m-1}\to r_{m-1}.
\]

If \(D_m\) had a common instance with

\[
T_1=((x\to y)\to z)\to C(y,z,u),
\]

comparison of the first outer antecedents would give

\[
\sigma(r_{m-1})=\tau(z).
\]

Comparison of the antecedents inside the outer consequents would also give

\[
\sigma(r_{m-1})=\tau(z\to u).
\]

Thus \(\tau(z)=\tau(z)\to\tau(u)\), impossible for a finite type.

### Proof: the case \(n=1, m=2\)

Write \(A_2=a\) and \(B_2=C(b,a,c)\). A common instance of \(D_2\) and \(T_1\) would imply

\[
\sigma(a)=\tau((x\to y)\to z)
\]

and

\[
\sigma(C(b,a,c))=\tau(y\to(z\to u)).
\]

The second equality forces

\[
\tau(y)=\sigma(b\to(a\to c)),\qquad
\tau(z)=\sigma(b),\qquad
\tau(u)=\sigma(c).
\]

Substitution into the first equality makes \(\sigma(a)\) contain itself as a proper subterm, a contradiction. ∎

## Classification of all typable applicative terms

### Theorem 4

Every typable applicative term whose leaves are copies of \(\Phi\) is a right comb \(R_n\).

### Proof

Induct on the application tree. A leaf is \(R_1\). If a composite term is \(MN\), then both subterms are typable, so by induction \(M=R_m\) and \(N=R_n\). If \(m\ge2\), the application is impossible by Lemma 3. Therefore \(m=1\), and the term is

\[
\Phi R_n=R_{n+1}.
\]

Lemma 1 supplies its principal type. ∎

## Reflexivity is not among the right-comb types

### Lemma 5

No substitution instance of any \(T_n\) is of the form \(P\to P\).

### Proof

For \(n\ge2\), equality of an instance of

\[
T_n=(A_n\to(B_n\to r_n))\to(A_n\to r_n)
\]

with \(P\to P\) would force

\[
\sigma(B_n\to r_n)=\sigma(r_n),
\]

making the finite type \(\sigma(r_n)\) contain itself as a proper subterm.

For \(n=1\), equality of the antecedent and consequent of

\[
T_1=((x\to y)\to z)\to C(y,z,u)
\]

would force

\[
\sigma(y)=\sigma(z\to u),\qquad
\sigma(z)=\sigma(y\to u).
\]

The first equality makes \(\sigma(y)\) strictly larger than \(\sigma(z)\), while the second makes \(\sigma(z)\) strictly larger than \(\sigma(y)\), an impossibility. ∎

## Main theorem

### Theorem 6

Ulrich's `u4`

\[
((x\to y)\to z)\to((y\to(z\to u))\to(y\to u))
\]

is **not** a single axiom for positive implicational logic under substitution and modus ponens.

### Proof

Any derivation from substitution instances of `u4` by modus ponens unfolds to a typable applicative term whose leaves are copies of `u4`. By Theorem 4, every such term is a right comb and therefore has a substitution instance of one of the principal types \(T_n\). By Lemma 5, no such type is reflexive. Hence even \(P\to P\), a theorem of positive implication, is not derivable from `u4`. Therefore `u4` is not a single axiom. ∎

## Consequence

Fitelson and Peltier (2026) eliminated Ulrich's other three 15-symbol candidates and stated that the shortest single axioms for positive implication have length 17 if and only if `u4` is not a single axiom. Combining their classification with the theorem above yields the consequence that the minimum length is **17 symbols**. The Lean development directly verifies the `u4` separation; the minimum-length conclusion additionally depends on their published elimination and classification results.

## Claim boundary

This document gives a self-contained finite-type/unification proof. A public timestamp for this result has been established. Academic acceptance and any comparative priority claim remain subject to independent technical review, peer review, and a search for equivalent prior published or publicly timestamped work.
