#!/usr/bin/env python3
"""Search finite matrices for a countermodel to Ulrich u4.

A returned model is accepted only after the explicit independent verifier
checks every MP instance, every u4 valuation, and the failed reflexivity
instance.  UNSAT/UNKNOWN results are never treated as a mathematical
conclusion.
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import z3


def verify(n: int, k: int, witness: int, table: list[list[int]]) -> tuple[bool, dict]:
    def I(a: int, b: int) -> int:
        return table[a][b]
    if len(table) != n or any(len(r) != n for r in table):
        return False, {"failure": "shape"}
    if any(not (0 <= v < n) for r in table for v in r):
        return False, {"failure": "range"}
    # Modus ponens closure.
    for x in range(n):
        for y in range(n):
            if x < k and I(x, y) < k and y >= k:
                return False, {"failure": "mp", "x": x, "y": y, "Ixy": I(x, y)}
    # Every substitution instance of u4 must be designated.
    for x in range(n):
        for y in range(n):
            for z in range(n):
                for u in range(n):
                    left = I(I(x, y), z)
                    right = I(I(y, I(z, u)), I(y, u))
                    value = I(left, right)
                    if value >= k:
                        return False, {"failure": "u4", "valuation": [x, y, z, u], "value": value}
    diag = I(witness, witness)
    if diag < k:
        return False, {"failure": "refl_not_refuted", "witness": witness, "value": diag}
    return True, {"failed_reflexivity_witness": witness, "Iww": diag}


def solve_case(n: int, k: int, witness: int, timeout_ms: int, seed: int) -> dict:
    s = z3.Solver()
    s.set(timeout=timeout_ms)
    s.set(random_seed=seed)
    op = z3.Array(f"op_{n}_{k}_{witness}", z3.IntSort(), z3.IntSort())

    def I(a, b):
        return z3.Select(op, a * n + b)

    for a in range(n):
        for b in range(n):
            v = I(a, b)
            s.add(v >= 0, v < n)

    # If x is designated and y is not, I(x,y) must not be designated;
    # this is exactly the MP closure condition for the fixed designated set {0,...,k-1}.
    for x in range(k):
        for y in range(k, n):
            s.add(I(x, y) >= k)

    # u4: (((x->y)->z) -> ((y->(z->u))->(y->u))) is designated.
    for x in range(n):
        for y in range(n):
            for z in range(n):
                for u in range(n):
                    left = I(I(x, y), z)
                    right = I(I(y, I(z, u)), I(y, u))
                    s.add(I(left, right) < k)

    # Chosen witness falsifies reflexivity.
    s.add(I(witness, witness) >= k)

    t0 = time.time()
    ans = s.check()
    elapsed = time.time() - t0
    out = {
        "n": n,
        "designated_count": k,
        "designated_values": list(range(k)),
        "witness": witness,
        "timeout_ms": timeout_ms,
        "solver_result": str(ans),
        "elapsed_seconds": elapsed,
        "z3_version": z3.get_version_string(),
    }
    if ans == z3.sat:
        m = s.model()
        table = [[m.eval(I(a, b), model_completion=True).as_long() for b in range(n)] for a in range(n)]
        ok, detail = verify(n, k, witness, table)
        out.update({"table": table, "independent_verifier_ok": ok, "verification": detail})
    else:
        out["independent_verifier_ok"] = False
        out["note"] = "UNSAT or UNKNOWN is logged only; it is not accepted as a u4 conclusion."
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--case-seconds", type=int, default=120)
    ap.add_argument("--global-seconds", type=int, default=1500)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    n = args.n
    deadline = time.time() + args.global_seconds
    runs = []
    found = None
    for k in range(1, n):
        # Any witness can be moved by a designated-set-preserving permutation
        # to 0 if designated, or to k if non-designated.
        for witness in (0, k):
            if time.time() >= deadline:
                break
            remaining = max(1, int(deadline - time.time()))
            sec = min(args.case_seconds, remaining)
            r = solve_case(n, k, witness, sec * 1000, 20260818 + 100*n + 10*k + witness)
            runs.append(r)
            print(json.dumps({q: r.get(q) for q in ("n","designated_count","witness","solver_result","elapsed_seconds","independent_verifier_ok")}), flush=True)
            if r.get("independent_verifier_ok"):
                found = r
                break
        if found or time.time() >= deadline:
            break
    result = {
        "problem": "Ulrich u4 finite-matrix countermodel search",
        "n": n,
        "status": "VERIFIED_FINITE_COUNTERMODEL" if found else "NO_VERIFIED_COUNTERMODEL_IN_THIS_RUN",
        "verified_model": found,
        "runs": runs,
        "claim_boundary": "Only a SAT table passing verify() is a result. UNSAT, UNKNOWN, and timeout are not conclusions.",
    }
    p = Path(args.out)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return 0 if found is None else 10

if __name__ == "__main__":
    raise SystemExit(main())
