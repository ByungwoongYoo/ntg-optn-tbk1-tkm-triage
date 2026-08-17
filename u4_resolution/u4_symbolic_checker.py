#!/usr/bin/env python3
"""Independent symbolic checks supporting the u4 non-single-axiom proof.

This is not a substitute for the general mathematical induction in
U4_NONAXIOM_PROOF.md. It independently implements finite simple-type
unification, constructs the right-comb principal types, checks every lemma on a
large finite prefix, and exhaustively enumerates all application-tree shapes up
to a configurable leaf count.
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
import argparse, json, hashlib
from typing import Dict, Iterable, Iterator, Optional, Tuple, Union

@dataclass(frozen=True)
class Var:
    n: int

@dataclass(frozen=True)
class Imp:
    a: "Ty"
    b: "Ty"

Ty = Union[Var, Imp]


def C(a: Ty, b: Ty, r: Ty) -> Ty:
    return Imp(Imp(a, Imp(b, r)), Imp(a, r))


def u4(offset: int = 0) -> Ty:
    x, y, z, u = [Var(offset + k) for k in range(4)]
    return Imp(Imp(Imp(x, y), z), C(y, z, u))


def vars_of(t: Ty) -> set[int]:
    if isinstance(t, Var):
        return {t.n}
    return vars_of(t.a) | vars_of(t.b)


def apply(t: Ty, s: Dict[int, Ty]) -> Ty:
    if isinstance(t, Var):
        trail = []
        while isinstance(t, Var) and t.n in s:
            trail.append(t.n)
            t = s[t.n]
        if not isinstance(t, Var):
            t = apply(t, s)
        for v in trail:
            s[v] = t
        return t
    return Imp(apply(t.a, s), apply(t.b, s))


def occurs(v: int, t: Ty, s: Dict[int, Ty]) -> bool:
    t = apply(t, s)
    if isinstance(t, Var):
        return t.n == v
    return occurs(v, t.a, s) or occurs(v, t.b, s)


def unify(a: Ty, b: Ty) -> Optional[Dict[int, Ty]]:
    s: Dict[int, Ty] = {}
    stack = [(a, b)]
    while stack:
        x, y = stack.pop()
        x, y = apply(x, s), apply(y, s)
        if x == y:
            continue
        if isinstance(x, Var):
            if occurs(x.n, y, s):
                return None
            s[x.n] = y
        elif isinstance(y, Var):
            if occurs(y.n, x, s):
                return None
            s[y.n] = x
        else:
            stack.append((x.a, y.a))
            stack.append((x.b, y.b))
    return s


def fresh(t: Ty, start: int) -> Tuple[Ty, int]:
    mp: Dict[int, int] = {}
    nxt = start
    def rec(x: Ty) -> Ty:
        nonlocal nxt
        if isinstance(x, Var):
            if x.n not in mp:
                mp[x.n] = nxt
                nxt += 1
            return Var(mp[x.n])
        return Imp(rec(x.a), rec(x.b))
    return rec(t), nxt


def canonical(t: Ty) -> Ty:
    mp: Dict[int, int] = {}
    nxt = 0
    def rec(x: Ty) -> Ty:
        nonlocal nxt
        if isinstance(x, Var):
            if x.n not in mp:
                mp[x.n] = nxt
                nxt += 1
            return Var(mp[x.n])
        return Imp(rec(x.a), rec(x.b))
    return rec(t)


def app_type(f: Ty, a: Ty) -> Optional[Ty]:
    ff, nxt = fresh(f, 0)
    aa, nxt = fresh(a, nxt)
    r = Var(nxt)
    s = unify(ff, Imp(aa, r))
    return None if s is None else canonical(apply(r, s))


def family(N: int) -> tuple[list[Optional[Ty]], list[Optional[Ty]], list[Optional[Ty]]]:
    """Return A[n], B[n], T[n] for 1<=n<=N using canonical global names.

    T1 uses v0..v3. For n>=2, A2=v2, B2=C(v1,v2,v3), r_n=v(n+2).
    """
    A: list[Optional[Ty]] = [None] * (N + 1)
    B: list[Optional[Ty]] = [None] * (N + 1)
    T: list[Optional[Ty]] = [None] * (N + 1)
    T[1] = u4()
    if N >= 2:
        A[2] = Var(2)
        B[2] = C(Var(1), Var(2), Var(3))
        T[2] = C(A[2], B[2], Var(4))
    for n in range(2, N):
        rn = Var(n + 2)
        A[n + 1] = Imp(B[n], rn)  # type: ignore[arg-type]
        B[n + 1] = Imp(A[n], rn)  # type: ignore[arg-type]
        T[n + 1] = C(A[n + 1], B[n + 1], Var(n + 3))  # type: ignore[arg-type]
    return A, B, T


def antecedent(t: Ty) -> Ty:
    assert isinstance(t, Imp)
    return t.a


def reflexive_unifiable(t: Ty) -> bool:
    tt, nxt = fresh(t, 0)
    q = Var(nxt)
    return unify(tt, Imp(q, q)) is not None


Tree = Optional[Tuple["Tree", "Tree"]]  # None is one u4 leaf

@lru_cache(maxsize=None)
def shapes(n: int) -> tuple[Tree, ...]:
    if n == 1:
        return (None,)
    out: list[Tree] = []
    for k in range(1, n):
        for l in shapes(k):
            for r in shapes(n-k):
                out.append((l, r))
    return tuple(out)


def infer_tree(sh: Tree) -> Optional[Ty]:
    if sh is None:
        return u4()
    lt = infer_tree(sh[0])
    if lt is None:
        return None
    rt = infer_tree(sh[1])
    if rt is None:
        return None
    return app_type(lt, rt)


def tree_text(sh: Tree) -> str:
    return "U" if sh is None else f"({tree_text(sh[0])} {tree_text(sh[1])})"


def right_comb(n: int) -> Tree:
    t: Tree = None
    for _ in range(1, n):
        t = (None, t)
    return t


def run(max_family: int, max_enum: int) -> dict:
    A, B, T = family(max_family + 1)
    checks = {
        "right_comb_recurrence": True,
        "interlocking_pairs": True,
        "major_application_exclusion": True,
        "reflexivity_exclusion": True,
    }
    failures = []

    for n in range(1, max_family):
        got = app_type(T[1], T[n])  # type: ignore[arg-type]
        want = canonical(T[n + 1])  # type: ignore[arg-type]
        if got != want:
            checks["right_comb_recurrence"] = False
            failures.append(["recurrence", n])

    for n in range(2, max_family + 1):
        aa, nxt = fresh(A[n], 0)  # type: ignore[arg-type]
        bb, nxt = fresh(B[n], nxt)  # independent here is stronger than needed
        w = Var(nxt)
        # The mathematical lemma uses a common substitution. To test it faithfully,
        # use the shared-variable originals as well.
        if unify(A[n], Imp(B[n], Var(max_family + 1000))) is not None:  # type: ignore[arg-type]
            checks["interlocking_pairs"] = False
            failures.append(["E(A,B)", n])
        if unify(B[n], Imp(A[n], Var(max_family + 1001))) is not None:  # type: ignore[arg-type]
            checks["interlocking_pairs"] = False
            failures.append(["E(B,A)", n])

    for m in range(2, max_family + 1):
        for n in range(1, max_family + 1):
            dm, nxt = fresh(antecedent(T[m]), 0)  # type: ignore[arg-type]
            tn, nxt = fresh(T[n], nxt)  # type: ignore[arg-type]
            if unify(dm, tn) is not None:
                checks["major_application_exclusion"] = False
                failures.append(["major", m, n])

    for n in range(1, max_family + 1):
        if reflexive_unifiable(T[n]):  # type: ignore[arg-type]
            checks["reflexivity_exclusion"] = False
            failures.append(["refl", n])

    enum_rows = []
    for n in range(1, max_enum + 1):
        typable = []
        for sh in shapes(n):
            ty = infer_tree(sh)
            if ty is not None:
                typable.append((tree_text(sh), canonical(ty)))
        expected = tree_text(right_comb(n))
        enum_rows.append({
            "leaves": n,
            "catalan_shapes": len(shapes(n)),
            "typable_shapes": len(typable),
            "unique_principal_types": len({repr(t) for _, t in typable}),
            "only_shape": typable[0][0] if len(typable) == 1 else None,
            "expected_right_comb": expected,
            "matches_theorem": len(typable) == 1 and typable[0][0] == expected,
        })
        if not enum_rows[-1]["matches_theorem"]:
            failures.append(["enumeration", n])

    payload = {
        "problem": "Ulrich u4 single-axiom question",
        "formula": "(((x->y)->z)->((y->(z->u))->(y->u)))",
        "max_family_index_checked": max_family,
        "max_application_tree_leaves_exhaustively_enumerated": max_enum,
        "checks": checks,
        "enumeration": enum_rows,
        "failures": failures,
        "all_checks_pass": all(checks.values()) and not failures,
        "claim_boundary": "Finite symbolic and exhaustive checks support, but do not replace, the all-n induction in U4_NONAXIOM_PROOF.md.",
    }
    return payload


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--family", type=int, default=250)
    ap.add_argument("--enumerate", type=int, default=11)
    ap.add_argument("--output", default="U4_SYMBOLIC_CHECK.json")
    args = ap.parse_args()
    data = run(args.family, args.enumerate)
    text = json.dumps(data, indent=2, sort_keys=True)
    print(text)
    open(args.output, "w", encoding="utf-8").write(text + "\n")
    if not data["all_checks_pass"]:
        raise SystemExit(1)

if __name__ == "__main__":
    main()
