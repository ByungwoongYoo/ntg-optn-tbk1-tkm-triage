#!/usr/bin/env python3
"""Verify a CaDiCaL competition-format SAT witness against the exact CNF.

For the B2/Z100 independent encoder, variables 1..100 are x_0..x_99.
The checker does not trust solver status: it evaluates every DIMACS clause,
then independently checks the selected residue set against the ordered-
difference definition.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path


def read_cnf(path: Path):
    nvars = None
    clauses = []
    with path.open('r', encoding='ascii', errors='strict') as f:
        current = []
        for raw in f:
            line = raw.strip()
            if not line or line.startswith('c'):
                continue
            if line.startswith('p '):
                parts = line.split()
                if len(parts) != 4 or parts[1] != 'cnf':
                    raise ValueError(f'bad DIMACS header: {line}')
                nvars = int(parts[2]); expected = int(parts[3])
                continue
            for tok in line.split():
                lit = int(tok)
                if lit == 0:
                    clauses.append(tuple(current)); current = []
                else:
                    current.append(lit)
        if current:
            raise ValueError('unterminated DIMACS clause')
    if nvars is None:
        raise ValueError('missing DIMACS header')
    if len(clauses) != expected:
        raise ValueError(f'clause count mismatch {len(clauses)} != {expected}')
    return nvars, clauses


def read_model(path: Path):
    status = None
    assignment = {}
    with path.open('r', encoding='utf-8', errors='replace') as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith('c'):
                continue
            if line.startswith('s '):
                status = line[2:].strip()
            elif line.startswith('v '):
                for tok in line[2:].split():
                    lit = int(tok)
                    if lit == 0:
                        continue
                    var = abs(lit)
                    val = lit > 0
                    if var in assignment and assignment[var] != val:
                        raise ValueError(f'contradictory assignment for variable {var}')
                    assignment[var] = val
    return status, assignment


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cnf', type=Path, required=True)
    ap.add_argument('--model', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    nvars, clauses = read_cnf(args.cnf)
    status, a = read_model(args.model)
    missing = [i for i in range(1, nvars + 1) if i not in a]
    bad_clause = None
    if not missing:
        for idx, clause in enumerate(clauses, 1):
            if not any(a[abs(lit)] == (lit > 0) for lit in clause):
                bad_clause = idx
                break
    selected = [i for i in range(100) if a.get(i + 1, False)]
    counts = [0] * 100
    for x in selected:
        for y in selected:
            if x != y:
                counts[(x - y) % 100] += 1
    max_mult = max(counts[1:]) if len(selected) > 1 else 0
    result = {
        'schema': 'b2-z100-cadical-model-verification-v1',
        'solver_status_line': status,
        'cnf_variables': nvars,
        'cnf_clauses': len(clauses),
        'assigned_variables': len(a),
        'missing_assignment_count': len(missing),
        'first_missing_variables': missing[:20],
        'first_unsatisfied_clause': bad_clause,
        'selected_residues': selected,
        'selected_size': len(selected),
        'contains_0_1': 0 in selected and 1 in selected,
        'ordered_difference_max_multiplicity': max_mult,
        'all_cnf_clauses_satisfied': not missing and bad_clause is None,
        'valid_14_set_counterexample': (
            status == 'SATISFIABLE' and not missing and bad_clause is None and
            len(selected) == 14 and 0 in selected and 1 in selected and max_mult <= 2
        ),
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2, sort_keys=True))
    if status != 'SATISFIABLE' or missing or bad_clause is not None:
        raise SystemExit(2)

if __name__ == '__main__':
    main()
