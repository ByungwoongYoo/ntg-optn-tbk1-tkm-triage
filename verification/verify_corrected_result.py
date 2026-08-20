#!/usr/bin/env python3
"""Verify the corrected B2/Z100 synopsis without overstating its scope.

This checker independently validates the explicit 13-set and the elementary
mod-10 support classification.  It also checks the corrected synopsis for a
complete, internally consistent map of the *recorded* final search runs.  It
does not rerun the exhaustive searches and is not third-party validation.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
RESULT_PATH = PACKAGE_ROOT / "claim/B2_Z100_FINAL_RESULTS_CORRECTED_20260820.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def main() -> int:
    result = json.loads(RESULT_PATH.read_text(encoding="utf-8"))

    require(result["schema_version"] == "2.1", "unexpected schema version")
    require(result["problem"]["reference"] == "b2-two-set-z100", "wrong problem reference")
    require(result["answer"] == "NO", "corrected synopsis does not record answer NO")
    require(result["maximum_size"] == 13, "corrected synopsis does not record maximum 13")
    require(
        result["status"]
        == "CANDIDATE_COMPUTATIONAL_RESOLUTION_WITH_DISTINCT_THIRD_PARTY_SUPPORT_PENDING_PACKET_REVIEW",
        "status boundary changed",
    )
    require(result["evidence_grade"] == "executable", "evidence grade changed")
    require(
        result["independent_third_party_reproduction_of_this_raw_corpus"] is False,
        "unsupported reproduction claim for this raw corpus",
    )
    require(result["peer_reviewed"] is False, "unsupported peer-review claim")
    require(result["formally_verified"] is False, "unsupported formal-verification claim")

    public = result["official_public_status"]
    require(public["status"] == "Open", "official P2624 status changed")
    require(public["established_bound"] == "13 <= M <= 14", "official bound changed")
    require(
        public["work_items"] == public["packet_items"] + public["recent_work_items"] == 17,
        "live work-item counts disagree",
    )
    external = result["distinct_external_audit"]
    require(external["records"] == ["R6088", "R6089", "R6090"], "external record map changed")
    require(external["package_name"] == "P2624_certified_resolution.zip", "external package name changed")
    require(
        external["package_sha256"]
        == "facbdfec87f5bfb302f197af5c13cce6e984444bd8270adce250de5e8fab8a35",
        "external package SHA-256 changed",
    )
    require(
        [
            external["quotient_occupancy_types"],
            external["admissible_raw_occupancy_vectors"],
            external["affine_orbits"],
            external["proof_tree_nodes"],
        ]
        == [7, 204360, 1341, 299903736],
        "external audit counts changed",
    )
    require(external["corrupt_byte_rejection_recorded"] is True, "negative control not recorded")
    require(external["source_url"] is None and external["source_locator"] is None, "external locator boundary changed")
    require(external["badge"] == "Needs packet proposal", "external packet badge changed")
    require(external["same_artifact_as_this_raw_corpus"] == "UNESTABLISHED", "artifact-lineage boundary changed")

    lower = result["lower_bound"]
    witness = lower["witness"]
    require(lower["cardinality"] == 13, "lower-bound cardinality field is not 13")
    require(len(witness) == 13 and len(set(witness)) == 13, "witness is not a 13-element set")
    require(all(isinstance(value, int) and 0 <= value < 100 for value in witness), "witness value out of range")
    counts = [0] * 100
    for left in witness:
        for right in witness:
            if left != right:
                counts[(left - right) % 100] += 1
    observed_max = max(counts[1:])
    require(observed_max == 2, "witness maximum ordered-difference multiplicity is not 2")
    require(
        lower["max_nonzero_ordered_difference_multiplicity"] == observed_max,
        "witness multiplicity field disagrees with direct computation",
    )

    # Exhaustively audit the elementary mod-10 support classification used by
    # the unit-difference normalization argument.
    unit_residues = {1, 3, 7, 9}
    valid_supports = 0
    for mask in range(1, 1 << 10):
        support = [residue for residue in range(10) if (mask >> residue) & 1]
        if all((a - b) % 10 not in unit_residues for a in support for b in support):
            same_parity = len({residue % 2 for residue in support}) == 1
            opposite_pair = (
                len(support) == 2
                and support[0] % 2 != support[1] % 2
                and (support[0] - support[1]) % 10 == 5
            )
            require(same_parity or opposite_pair, f"unclassified support: {support}")
            valid_supports += 1
    require(valid_supports == 67, "mod-10 support count is not 67")
    require(14 * 13 == 182 and 49 * 2 < 182 and 19 * 2 < 182, "counting inequalities failed")

    upper = result["upper_bound_reduction"]
    require(upper["unit_difference_lemma"] is True, "unit-difference lemma flag is false")
    require(upper["normalized_third_range"] == [2, 45], "normalized third range changed")
    require(upper["normalized_third_branches"] == 44, "normalized branch count changed")
    require(
        upper["normalized_third_range"][1] - upper["normalized_third_range"][0] + 1
        == upper["normalized_third_branches"],
        "normalized range and branch count disagree",
    )
    require(upper["all_final_branches_completed"] is True, "final branches not recorded complete")
    require(upper["final_timeouts"] == 0, "final timeout count is nonzero")
    require(upper["final_incomplete_branches"] == 0, "final incomplete count is nonzero")
    require(upper["final_witnesses"] == 0, "a final 14-set witness is recorded")

    dependencies = result["final_exact_search_dependencies"]
    v10 = dependencies["g1_t2_t4_v10"]
    v8_head = dependencies["g1_t5_t10_v8"]
    v8_tail = dependencies["g1_t11_t45_v8"]
    require(v10["run_id"] == 32040699080 and v10["final_subbranches"] == 255, "v10 dependency changed")
    require(
        v10["regression_control"] == {"third": 7, "fourth": 8, "count": 1},
        "v10 regression-control record changed",
    )
    require(
        v8_head["run_id"] == 32038183046
        and v8_head["third_range"] == [5, 10]
        and v8_head["final_branches"] == 6,
        "v8 head dependency changed",
    )
    require(
        v8_tail["run_id"] == 32038657803
        and v8_tail["third_range"] == [11, 45]
        and v8_tail["final_branches"] == 35,
        "v8 tail dependency changed",
    )
    require(3 + v8_head["final_branches"] + v8_tail["final_branches"] == 44, "primary t coverage is not 44")

    redundant = result["redundant_checks"]
    require(
        redundant["g2_g3_v9"]["run_id"] == 32040180627
        and redundant["g2_g3_v9"]["branches"] == 67
        and redundant["g2_g3_v9"]["all_completed_without_witness"] is True,
        "g=2,3 cross-check record changed",
    )
    require(redundant["g4_g7_drat"]["cases"] == [4, 5, 6, 7], "DRAT case coverage changed")
    require(redundant["g4_g7_drat"]["cadical_unsat_exit_code"] == 20, "CaDiCaL exit code changed")
    require(redundant["g4_g7_drat"]["actual_cnf_and_drat_preserved"] is True, "proof preservation flag false")

    raw = result["canonical_raw_release"]
    require(raw["bytes"] == 219202341, "canonical raw byte count changed")
    require(
        raw["sha256"] == "254031fc2fab17027e389900ca63e704d170eba6ad6861e235db3fe9be46727a",
        "canonical raw SHA-256 changed",
    )
    require(raw["manifest_file_count"] == 4155, "raw manifest file count changed")
    require(raw["manifest_uncompressed_bytes"] == 518219226, "raw manifest byte count changed")

    print(
        json.dumps(
            {
                "schema": "b2-z100-corrected-result-check-v1",
                "status": "PASS",
                "claim_file": RESULT_PATH.relative_to(PACKAGE_ROOT).as_posix(),
                "witness_cardinality": len(witness),
                "max_ordered_difference_multiplicity": observed_max,
                "mod10_supports_checked": valid_supports,
                "normalized_third_branches_recorded": upper["normalized_third_branches"],
                "distinct_third_party_audit_records": len(external["records"]),
                "scope": "elementary checks plus internal consistency of recorded run map; no solver replay",
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        raise SystemExit(1)
