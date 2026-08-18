% Ulrich u4, checked against Fitelson--Peltier (JAR 2026), Figure c.
fof(mp,axiom, ![X,Y]: ((p(i(X,Y)) & p(X)) => p(Y))).
fof(refl,conjecture, ![P]: p(i(P,P))).
fof(u4,axiom, ![X,Y,Z,U]: p(i(i(i(X,Y),Z), i(i(Y,i(Z,U)), i(Y,U))))).
