% Ulrich u4 in positive implicational logic.
% The candidate axiom is asserted; reflexivity is the target theorem.
% Modus ponens is encoded as a closure axiom over predicate p.

fof(mp,axiom,
    ! [X,Y] :
      ( (p(i(X,Y)) & p(X)) => p(Y) ) ).

fof(u4,axiom,
    ! [X,Y,Z,U] :
      p(i(i(i(X,Y),Z), i(i(Y,i(Z,U)),i(Y,U)))) ).

fof(refl,conjecture,
    ! [Q] : p(i(Q,Q)) ).
