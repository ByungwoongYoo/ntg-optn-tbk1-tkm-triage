# An explicit substitution-closed model refuting the single-axiomhood of `u4`

This note repackages the type argument in `U4_NONAXIOM_PROOF.md` as an explicit infinite model of the first-order encoding used by Fitelson and Peltier.

Let

\[
\Phi=((x\to y)\to z)\to C(y,z,u),
\qquad
C(A,B,R):=(A\to(B\to R))\to(A\to R).
\]

Let \(T_1,T_2,\ldots\) be the right-comb principal types defined in the main proof:

\[
T_1=\Phi,
\]

\[
A_2=a,\quad B_2=C(b,a,c),\quad T_2=C(A_2,B_2,r_2),
\]

and, for \(n\ge2\),

\[
A_{n+1}=B_n\to r_n,\qquad
B_{n+1}=A_n\to r_n,\qquad
T_{n+1}=C(A_{n+1},B_{n+1},r_{n+1}).
\]

For a type scheme \(S\), write \(\operatorname{Inst}(S)\) for its finite substitution instances, and define

\[
\mathcal P=\bigcup_{n\ge1}\operatorname{Inst}(T_n).
\]

Interpret the unary theorem predicate \(p\) by

\[
p(F)\quad\Longleftrightarrow\quad F\in\mathcal P.
\]

The binary implication symbol is interpreted as the free binary type constructor.

## 1. Every instance of `u4` is true

This is immediate because \(T_1=\Phi\). Thus every substitution instance of \(\Phi\) belongs to \(\mathcal P\).

## 2. \(\mathcal P\) is closed under modus ponens

Suppose

\[
X\in\mathcal P,\qquad X\to Y\in\mathcal P.
\]

Choose \(n,m\ge1\) and substitutions \(\alpha,\beta\) with

\[
X=\alpha(T_n),\qquad X\to Y=\beta(T_m).
\]

If \(m\ge2\), the antecedent of \(\beta(T_m)\) is \(\beta(D_m)\), where

\[
D_m=A_m\to(B_m\to r_m).
\]

Consequently \(X\) would be a common substitution instance of \(D_m\) and \(T_n\). Lemma 3 of the main proof shows that this is impossible. Therefore \(m=1\).

Hence

\[
X\to Y=\beta(T_1)
 =\beta\bigl(((x\to y)\to z)\to C(y,z,u)\bigr),
\]

so

\[
X=\beta((x\to y)\to z),\qquad
Y=\beta(C(y,z,u)).
\]

We now show that \(Y\in\mathcal P\).

### Case \(n=1\)

Rename the variables of the minor copy of \(T_1\) to \(x',y',z',u'\). Equality

\[
\alpha(T_1)=\beta((x\to y)\to z)
\]

and injectivity of \(\to\) yield

\[
\beta(y)=\alpha(z'),\qquad
\beta(z)=\alpha(C(y',z',u')).
\]

Therefore

\[
Y=C\bigl(\alpha(z'),\alpha(C(y',z',u')),\beta(u)\bigr),
\]

which is a substitution instance of

\[
T_2=C(a,C(b,a,c),r_2).
\]

Thus \(Y\in\operatorname{Inst}(T_2)\subseteq\mathcal P\).

### Case \(n\ge2\)

Write

\[
T_n=(A_n\to(B_n\to r_n))\to(A_n\to r_n).
\]

From

\[
\alpha(T_n)=\beta((x\to y)\to z)
\]

we obtain

\[
\beta(y)=\alpha(B_n\to r_n),\qquad
\beta(z)=\alpha(A_n\to r_n).
\]

Hence

\[
Y=C\bigl(\alpha(B_n\to r_n),\alpha(A_n\to r_n),\beta(u)\bigr).
\]

Because \(r_{n+1}\) is fresh, extend \(\alpha\) by sending \(r_{n+1}\) to \(\beta(u)\). Under this extended substitution the displayed formula is exactly an instance of

\[
T_{n+1}=C(B_n\to r_n,A_n\to r_n,r_{n+1}).
\]

Thus \(Y\in\operatorname{Inst}(T_{n+1})\subseteq\mathcal P\).

In either case, \(X,X\to Y\in\mathcal P\) implies \(Y\in\mathcal P\). Therefore the interpretation satisfies modus ponens.

## 3. Reflexivity is false in this model

Lemma 5 of the main proof shows that no finite substitution instance of any \(T_n\) has the form \(P\to P\). Therefore

\[
P\to P\notin\mathcal P
\]

for every finite implicational formula \(P\).

## 4. Conclusion

The free term algebra equipped with the predicate \(\mathcal P\) satisfies every substitution instance of `u4` and is closed under modus ponens, while it falsifies reflexivity. Hence

\[
\{\mathrm{u4},\mathrm{mp}\}\nvdash P\to P.
\]

Since reflexivity is a theorem of positive implicational logic, `u4` is not a single axiom.

This model is infinite: its true formulas are precisely the substitution instances of the right-comb schemes \(T_1,T_2,\ldots\). The construction explains why finite-matrix searches can fail even though `u4` is not a single axiom.
