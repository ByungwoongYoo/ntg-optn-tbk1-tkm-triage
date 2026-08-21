import importlib.util
import contextlib
import hashlib
import io
import json
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "validate_panax_remote_archive.py"
SPEC = importlib.util.spec_from_file_location("panax_remote_validator", MODULE_PATH)
VALIDATOR = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(VALIDATOR)
TEST_CANDIDATES = ("PNX_Picorna_A1", "PNX_Picorna_A2", "PNX_Picorna_B")
TEST_CANDIDATE = TEST_CANDIDATES[0]
TEST_CONTROL = "PNX_Panax_L2_control"


def search_report(
    qid, qlen, *, hits=None, valid_stats=True,
    program=None, database=None, expect=None, version=None,
):
    hits = [] if hits is None else hits
    stat = {
        "db_num": 1000,
        "db_len": 500000,
        "kappa": 0.041,
        "lambda": 0.267,
        "entropy": 0.140,
    }
    if not valid_stats:
        stat.update({"kappa": 0, "lambda": 0, "entropy": 0})
    search = {
        "query_id": "Query_1",
        "query_title": qid,
        "query_len": qlen,
        "hits": hits,
        "stat": stat,
    }
    if not hits:
        search["message"] = "No hits found"
    report = {"results": {"search": search}}
    if version is not None:
        report["version"] = version
    report["reference"] = "mock BLAST reference"
    if program is not None:
        report["program"] = program
    if database is not None:
        report["search_target"] = {"db": database}
    if expect is not None:
        report["params"] = {"expect": expect}
    return {"report": report}


def hit_values(
    qid, qlen, *, saccver="REP00001.1", sallacc=None,
    sallseqid="ref|REP00001.1|",
):
    if sallacc is None:
        sallacc = saccver.split(".")[0]
    return {
        "qseqid": qid,
        "saccver": saccver,
        "sallacc": sallacc,
        "sallseqid": sallseqid,
        "pident": "100.0",
        "length": str(qlen),
        "qlen": str(qlen),
        "slen": str(qlen),
        "qstart": "1",
        "qend": str(qlen),
        "sstart": "1",
        "send": str(qlen),
        "evalue": "0.0",
        "bitscore": "500",
        "qcovs": "100",
        "staxids": "1",
        "sscinames": "N/A",
        "stitle": "mock protein",
        "qseq": "A" * qlen,
        "sseq": "A" * qlen,
    }


def hit_row(
    qid, qlen, *, saccver="REP00001.1", sallacc=None,
    sallseqid="ref|REP00001.1|",
):
    values = hit_values(
        qid, qlen, saccver=saccver, sallacc=sallacc,
        sallseqid=sallseqid,
    )
    return "\t".join(values[field] for field in VALIDATOR.HIT_FIELDS)


def json_hit(qlen, *, accession="REP00001", qseq=None, sseq=None, slen=None):
    qseq = "A" * qlen if qseq is None else qseq
    sseq = "A" * qlen if sseq is None else sseq
    slen = qlen if slen is None else slen
    return {
        "description": [{
            "id": f"ref|{accession}.1|", "accession": accession,
            "title": "mock protein", "taxid": 1,
        }],
        "len": slen,
        "hsps": [{
            "query_from": 1,
            "query_to": len(qseq.replace("-", "")),
            "hit_from": 1,
            "hit_to": len(sseq.replace("-", "")),
            "align_len": len(qseq),
            "identity": sum(
                q.upper() == s.upper() and q != "-"
                for q, s in zip(qseq, sseq)
            ),
            "evalue": 0.0,
            "bit_score": 500.0,
            "qseq": qseq,
            "hseq": sseq,
        }],
    }


class ValidatorTests(unittest.TestCase):
    def run_case(
        self, reports, hits_text, *, control=True, mode=None,
        omit_control_metadata=False, control_length=10,
    ):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            query_path = root / "queries.faa"
            query_records = [(candidate, "A" * 12) for candidate in TEST_CANDIDATES]
            if control:
                query_records.append((TEST_CONTROL, "A" * control_length))
            query_path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in query_records))
            query_specs = [
                {
                    "id": candidate, "length": 12,
                    "sequence_sha256": hashlib.sha256(b"A" * 12).hexdigest(),
                }
                for candidate in TEST_CANDIDATES
            ]
            if control:
                query_specs.append({
                    "id": TEST_CONTROL, "length": control_length,
                    "sequence_sha256": hashlib.sha256(
                        b"A" * control_length
                    ).hexdigest(),
                })
            expected = {
                "queries": query_specs,
                "query_file": str(query_path),
                "query_file_sha256": hashlib.sha256(query_path.read_bytes()).hexdigest(),
                "validation_control_ids": (
                    [] if omit_control_metadata else ([TEST_CONTROL] if control else [])
                ),
                "validation_controls": (
                    [
                        {
                            "id": TEST_CONTROL,
                            "expected_accession": "YP_009121238.1",
                            "min_query_coverage": 99.0,
                            "min_identity": 99.0,
                        }
                    ]
                    if control and not omit_control_metadata
                    else []
                ),
            }
            effective_mode = mode or ("protein_nonviral" if control else "protein_tsa")
            search_contract = VALIDATOR.MODE_SEARCH_CONTRACT[effective_mode]
            result_path = root / "RESULTS.json"
            expected_path = root / "EXPECTED_QUERIES.json"
            hits_path = root / "HITS.tsv"
            observed_report_ids = {
                item["report"]["results"]["search"]["query_title"].split()[0]
                for item in reports
            }
            reports = list(reports) + [
                search_report(candidate, 12)
                for candidate in TEST_CANDIDATES
                if candidate not in observed_report_ids
            ]
            for item in reports:
                report = item.get("report") if isinstance(item, dict) else None
                if not isinstance(report, dict):
                    continue
                report.setdefault("program", search_contract["program"])
                report.setdefault(
                    "version", f"{search_contract['program'].upper()} 2.17.0+"
                )
                report.setdefault("search_target", {"db": search_contract["database"]})
                report.setdefault("params", dict(search_contract["params"]))
            result_path.write_text(json.dumps({"BlastOutput2": reports}))
            expected_path.write_text(json.dumps(expected))
            hits_path.write_text(hits_text)
            return VALIDATOR.validate(
                result_path, expected_path,
                effective_mode, hits_path,
            )

    def test_zero_statistics_are_structural_before_control_failure(self):
        reports = [
            search_report(TEST_CANDIDATE, 12, valid_stats=False),
            search_report(TEST_CONTROL, 10, valid_stats=False),
        ]
        structural, control, observed = self.run_case(reports, "")
        self.assertEqual(observed, 4)
        self.assertTrue(any("invalid result statistics" in x for x in structural))
        self.assertEqual(control, [])

    def test_complete_archive_with_missing_control_is_deterministic(self):
        reports = [search_report(TEST_CANDIDATE, 12), search_report(TEST_CONTROL, 10)]
        structural, control, _ = self.run_case(reports, "")
        self.assertEqual(structural, [])
        self.assertEqual(control, [f"positive control has no hit: {TEST_CONTROL}"])

    def test_expected_accession_in_nr_aliases_is_accepted(self):
        candidate = search_report(TEST_CANDIDATE, 12)
        control_hit = json_hit(10)
        control_hit["description"].append({
            "id": "gb|YP_009121238.1|", "accession": "YP_009121238",
            "title": "mock protein", "taxid": 1,
        })
        control = search_report(TEST_CONTROL, 10, hits=[control_hit])
        row = hit_row(
            TEST_CONTROL,
            10,
            sallacc="REP00001;YP_009121238",
            sallseqid="ref|REP00001.1|;gb|YP_009121238.1|",
        )
        structural, control_errors, _ = self.run_case(
            [candidate, control], row + "\n"
        )
        self.assertEqual(structural, [])
        self.assertEqual(control_errors, [])

    def test_tsv_only_control_accession_injection_is_rejected(self):
        candidate = search_report(TEST_CANDIDATE, 12)
        control = search_report(TEST_CONTROL, 10, hits=[json_hit(10)])
        row = hit_row(
            TEST_CONTROL, 10,
            sallseqid="ref|REP00001.1|;gb|YP_009121238.1|",
        )
        structural, control_errors, _ = self.run_case(
            [candidate, control], row + "\n"
        )
        self.assertTrue(any("HSP identity mismatch" in x for x in structural))
        self.assertEqual(control_errors, [])

    def test_primary_tsv_accession_must_be_present_in_json(self):
        reports = [
            search_report(TEST_CANDIDATE, 12, hits=[json_hit(12)])
        ]
        row = hit_row(
            TEST_CANDIDATE, 12,
            saccver="FAKE99999.9",
            sallseqid="ref|FAKE99999.9|;ref|REP00001.1|",
        )
        structural, control, _ = self.run_case(
            reports, row + "\n", control=False
        )
        self.assertTrue(any("HSP identity mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_tsv_taxids_cannot_add_an_unbound_value(self):
        reports = [
            search_report(TEST_CANDIDATE, 12, hits=[json_hit(12)])
        ]
        values = hit_values(TEST_CANDIDATE, 12)
        values["staxids"] = "1;999999"
        row = "\t".join(values[field] for field in VALIDATOR.HIT_FIELDS)
        structural, control, _ = self.run_case(
            reports, row + "\n", control=False
        )
        self.assertTrue(any("HSP identity mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_hits_and_no_hit_message_are_mutually_exclusive(self):
        report = search_report(TEST_CANDIDATE, 12, hits=[json_hit(12)])
        report["report"]["results"]["search"]["message"] = "No hits found"
        structural, control, _ = self.run_case(
            [report], hit_row(TEST_CANDIDATE, 12) + "\n", control=False
        )
        self.assertTrue(any("contradictory no-hit" in x for x in structural))
        self.assertEqual(control, [])

    def test_hits_with_fatal_result_message_are_structural(self):
        report = search_report(TEST_CANDIDATE, 12, hits=[json_hit(12)])
        report["report"]["results"]["search"]["message"] = "BLAST Database error"
        structural, control, _ = self.run_case(
            [report], hit_row(TEST_CANDIDATE, 12) + "\n", control=False
        )
        self.assertTrue(any("error result message" in x for x in structural))
        self.assertEqual(control, [])

    def test_incomplete_json_alias_cannot_support_tsv_alias(self):
        candidate_hit = json_hit(12)
        candidate_hit["description"].append({
            "id": "gb|EVIL99999.1|", "accession": "EVIL99999",
        })
        reports = [search_report(TEST_CANDIDATE, 12, hits=[candidate_hit])]
        row = hit_row(
            TEST_CANDIDATE, 12,
            sallacc="REP00001;EVIL99999",
            sallseqid="ref|REP00001.1|;gb|EVIL99999.1|",
        )
        structural, control, _ = self.run_case(
            reports, row + "\n", control=False
        )
        self.assertTrue(any("valid subject identity/length" in x for x in structural))
        self.assertEqual(control, [])

    def test_json_tsv_count_mismatch_is_structural(self):
        candidate_hit = json_hit(12, accession="HIT00001")
        reports = [
            search_report(TEST_CANDIDATE, 12, hits=[candidate_hit]),
        ]
        structural, control, _ = self.run_case(reports, "", control=False)
        self.assertTrue(any("JSON/TSV HSP-count mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_multiple_hsps_are_rejected_under_max_hsps_one_contract(self):
        candidate_hit = json_hit(12, accession="HIT00001")
        candidate_hit["hsps"].append(dict(candidate_hit["hsps"][0]))
        reports = [
            search_report(TEST_CANDIDATE, 12, hits=[candidate_hit]),
        ]
        rows = "\n".join([hit_row(TEST_CANDIDATE, 12), hit_row(TEST_CANDIDATE, 12)])
        structural, control, _ = self.run_case(
            reports, rows + "\n", control=False
        )
        self.assertTrue(any("despite max_hsps=1" in x for x in structural))
        self.assertEqual(control, [])

    def test_empty_qseq_is_structural(self):
        candidate_hit = json_hit(12, accession="HIT00001")
        reports = [
            search_report(TEST_CANDIDATE, 12, hits=[candidate_hit]),
        ]
        values = hit_row(TEST_CANDIDATE, 12).split("\t")
        values[VALIDATOR.HIT_FIELDS.index("qseq")] = ""
        structural, control, _ = self.run_case(
            reports, "\t".join(values) + "\n", control=False
        )
        self.assertTrue(any("aligned query/subject" in x for x in structural))
        self.assertEqual(control, [])

    def test_json_subject_cannot_be_replaced_only_in_tsv(self):
        candidate_hit = json_hit(12, accession="OTHER0001")
        reports = [
            search_report(TEST_CANDIDATE, 12, hits=[candidate_hit]),
        ]
        structural, control, _ = self.run_case(
            reports, hit_row(TEST_CANDIDATE, 12) + "\n", control=False
        )
        self.assertTrue(any("HSP identity mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_json_tsv_accession_versions_must_match(self):
        candidate_hit = json_hit(12, accession="REP00001")
        reports = [search_report(TEST_CANDIDATE, 12, hits=[candidate_hit])]
        structural, control, _ = self.run_case(
            reports,
            hit_row(
                TEST_CANDIDATE, 12,
                saccver="REP00001.9", sallseqid="ref|REP00001.9|",
            ) + "\n",
            control=False,
        )
        self.assertTrue(any("HSP identity mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_saccver_cannot_drop_json_accession_version(self):
        reports = [
            search_report(TEST_CANDIDATE, 12, hits=[json_hit(12)])
        ]
        structural, control, _ = self.run_case(
            reports,
            hit_row(
                TEST_CANDIDATE, 12,
                saccver="REP00001", sallseqid="ref|REP00001.1|",
            ) + "\n",
            control=False,
        )
        self.assertTrue(any("HSP identity mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_sallacc_cannot_add_tsv_only_alias(self):
        reports = [
            search_report(TEST_CANDIDATE, 12, hits=[json_hit(12)])
        ]
        structural, control, _ = self.run_case(
            reports,
            hit_row(
                TEST_CANDIDATE, 12,
                sallacc="REP00001;EVIL99999",
            ) + "\n",
            control=False,
        )
        self.assertTrue(any("HSP identity mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_sallseqid_cannot_add_tsv_only_general_id(self):
        reports = [
            search_report(TEST_CANDIDATE, 12, hits=[json_hit(12)])
        ]
        structural, control, _ = self.run_case(
            reports,
            hit_row(
                TEST_CANDIDATE, 12,
                sallseqid="ref|REP00001.1|;gnl|namespace|EVILTAG",
            ) + "\n",
            control=False,
        )
        self.assertTrue(any("HSP identity mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_standard_pdb_chain_identifier_is_accepted(self):
        candidate_hit = json_hit(12, accession="3GRX_A")
        candidate_hit["description"][0]["id"] = "pdb|3GRX|A"
        reports = [search_report(TEST_CANDIDATE, 12, hits=[candidate_hit])]
        row = hit_row(
            TEST_CANDIDATE, 12,
            saccver="3GRX_A", sallacc="3GRX_A",
            sallseqid="pdb|3GRX|A",
        )
        structural, control, _ = self.run_case(
            reports, row + "\n", control=False
        )
        self.assertEqual(structural, [])
        self.assertEqual(control, [])

    def test_duplicate_primary_subject_is_structural(self):
        duplicate_hit = json_hit(12)
        reports = [
            search_report(
                TEST_CANDIDATE, 12,
                hits=[duplicate_hit, json.loads(json.dumps(duplicate_hit))],
            )
        ]
        row = hit_row(TEST_CANDIDATE, 12)
        structural, control, _ = self.run_case(
            reports, row + "\n" + row + "\n", control=False
        )
        self.assertTrue(any("duplicate primary subject" in x for x in structural))
        self.assertEqual(control, [])

    def test_reports_must_share_database_size_signature(self):
        reports = [search_report(TEST_CANDIDATE, 12)]
        control_report = search_report(TEST_CONTROL, 10)
        control_report["report"]["results"]["search"]["stat"].update({
            "db_num": 777, "db_len": 123456,
        })
        reports.append(control_report)
        structural, control, _ = self.run_case(reports, "")
        self.assertTrue(any("database-size signature" in x for x in structural))
        self.assertEqual(control, [])

    def test_fake_full_qcov_on_one_residue_control_is_structural(self):
        control_hit = json_hit(1, accession="REP00001", slen=10)
        reports = [
            search_report(TEST_CANDIDATE, 12),
            search_report(TEST_CONTROL, 10, hits=[control_hit]),
        ]
        values = hit_values(
            TEST_CONTROL, 10,
            sallacc="REP00001;YP_009121238",
            sallseqid="ref|REP00001.1|;gb|YP_009121238.1|",
        )
        values.update({
            "length": "1", "qend": "1", "send": "1",
            "qseq": "A", "sseq": "A", "qcovs": "100",
        })
        row = "\t".join(values[field] for field in VALIDATOR.HIT_FIELDS)
        structural, control, _ = self.run_case(reports, row + "\n")
        self.assertTrue(any("query-coverage/alignment mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_reported_identity_must_match_aligned_sequences(self):
        control_hit = json_hit(10, qseq="A" * 10, sseq="C" * 10)
        reports = [
            search_report(TEST_CANDIDATE, 12),
            search_report(TEST_CONTROL, 10, hits=[control_hit]),
        ]
        values = hit_values(
            TEST_CONTROL, 10,
            sallacc="REP00001;YP_009121238",
            sallseqid="ref|REP00001.1|;gb|YP_009121238.1|",
        )
        values["sseq"] = "C" * 10
        row = "\t".join(values[field] for field in VALIDATOR.HIT_FIELDS)
        structural, control, _ = self.run_case(reports, row + "\n")
        self.assertTrue(any("identity/alignment mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_reported_evalue_and_bitscore_must_bind_to_json_hsp(self):
        candidate_hit = json_hit(12)
        reports = [search_report(TEST_CANDIDATE, 12, hits=[candidate_hit])]
        values = hit_values(TEST_CANDIDATE, 12)
        values.update({"evalue": "999", "bitscore": "999999"})
        row = "\t".join(values[field] for field in VALIDATOR.HIT_FIELDS)
        structural, control, _ = self.run_case(
            reports, row + "\n", control=False
        )
        self.assertTrue(any("HSP identity mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_search_program_database_and_expect_are_mode_bound(self):
        reports = [
            search_report(
                TEST_CANDIDATE, 12,
                program="blastn", database="nt", expect=1.0,
            )
        ]
        structural, control, _ = self.run_case(
            reports, "", control=False, mode="protein_tsa"
        )
        self.assertTrue(any("program mismatch" in x for x in structural))
        self.assertTrue(any("database mismatch" in x for x in structural))
        self.assertTrue(any("expect mismatch" in x for x in structural))
        self.assertTrue(any("parameter mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_blast_version_prefix_must_match_program(self):
        reports = [search_report(
            TEST_CANDIDATE, 12, version="BLASTN 2.17.0+",
        )]
        structural, control, _ = self.run_case(
            reports, "", control=False, mode="protein_tsa"
        )
        self.assertTrue(any("version/program mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_search_target_cannot_carry_unasserted_extra_fields(self):
        reports = [search_report(TEST_CANDIDATE, 12)]
        reports[0]["report"]["search_target"] = {
            "db": "tsa_nr", "entrez_query": "txid10239[ORGN]",
        }
        structural, control, _ = self.run_case(
            reports, "", control=False, mode="protein_tsa"
        )
        self.assertTrue(any("database mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_unexpected_search_parameter_is_structural(self):
        reports = [search_report(TEST_CANDIDATE, 12)]
        reports[0]["report"]["params"] = {
            **VALIDATOR.MODE_SEARCH_CONTRACT["protein_tsa"]["params"],
            "word_size": 999,
        }
        structural, control, _ = self.run_case(
            reports, "", control=False, mode="protein_tsa"
        )
        self.assertTrue(any("parameter-key mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_control_threshold_uses_unrounded_query_coverage(self):
        control_hit = json_hit(985, slen=1000)
        control_hit["description"].append({
            "id": "gb|YP_009121238.1|", "accession": "YP_009121238",
            "title": "mock protein", "taxid": 1,
        })
        reports = [
            search_report(TEST_CANDIDATE, 12),
            search_report(TEST_CONTROL, 1000, hits=[control_hit]),
        ]
        values = hit_values(
            TEST_CONTROL, 1000,
            sallacc="REP00001;YP_009121238",
            sallseqid="ref|REP00001.1|;gb|YP_009121238.1|",
        )
        values.update({
            "length": "985", "qend": "985", "send": "985",
            "qseq": "A" * 985, "sseq": "A" * 985, "qcovs": "99",
        })
        row = "\t".join(values[field] for field in VALIDATOR.HIT_FIELDS)
        structural, control, _ = self.run_case(
            reports, row + "\n", control_length=1000
        )
        self.assertEqual(structural, [])
        self.assertEqual(
            control,
            [f"positive control did not recover a near-exact match: {TEST_CONTROL}"],
        )

    def test_standard_blastn_scoring_signature_cannot_be_megablast(self):
        reports = [search_report(TEST_CANDIDATE, 12)]
        structural, control, _ = self.run_case(
            reports, "", control=False, mode="nt_tsa"
        )
        self.assertEqual(structural, [])
        self.assertEqual(control, [])
        wrong = [search_report(TEST_CANDIDATE, 12)]
        wrong[0]["report"]["params"] = dict(
            VALIDATOR.MODE_SEARCH_CONTRACT["nt_megablast"]["params"]
        )
        structural, control, _ = self.run_case(
            wrong, "", control=False, mode="nt_tsa"
        )
        self.assertTrue(any("sc_match" in x for x in structural))
        self.assertTrue(any("sc_mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_hsp_evalue_above_command_threshold_is_structural(self):
        candidate_hit = json_hit(12)
        candidate_hit["hsps"][0]["evalue"] = 1.0
        reports = [search_report(TEST_CANDIDATE, 12, hits=[candidate_hit])]
        values = hit_values(TEST_CANDIDATE, 12)
        values["evalue"] = "1.0"
        row = "\t".join(values[field] for field in VALIDATOR.HIT_FIELDS)
        structural, control, _ = self.run_case(
            reports, row + "\n", control=False, mode="protein_tsa"
        )
        self.assertTrue(any("invalid JSON HSP numeric range" in x for x in structural))
        self.assertTrue(any("invalid numeric range" in x for x in structural))
        self.assertEqual(control, [])

    def test_external_json_integer_fields_are_strict(self):
        candidate_hit = json_hit(12)
        candidate_hit["len"] = 12.9
        candidate_hit["hsps"][0]["align_len"] = 12.9
        reports = [search_report(TEST_CANDIDATE, 12, hits=[candidate_hit])]
        structural, control, _ = self.run_case(
            reports, hit_row(TEST_CANDIDATE, 12) + "\n", control=False
        )
        self.assertTrue(any("valid subject identity/length" in x for x in structural))
        self.assertTrue(any("incomplete JSON HSP" in x for x in structural))
        self.assertEqual(control, [])

    def test_null_json_accession_is_not_coerced_to_text(self):
        candidate_hit = json_hit(12)
        candidate_hit["description"] = [{
            "id": None, "accession": None,
            "title": "mock protein", "taxid": 1,
        }]
        reports = [search_report(TEST_CANDIDATE, 12, hits=[candidate_hit])]
        structural, control, _ = self.run_case(
            reports,
            hit_row(
                TEST_CANDIDATE, 12,
                saccver="None", sallseqid="ref|None|",
            ) + "\n",
            control=False,
        )
        self.assertTrue(any("valid subject identity/length" in x for x in structural))
        self.assertEqual(control, [])

    def test_json_title_and_taxid_must_bind_to_same_description(self):
        candidate_hit = json_hit(12)
        candidate_hit["description"][0]["title"] = "different title"
        candidate_hit["description"][0]["taxid"] = 2
        reports = [search_report(TEST_CANDIDATE, 12, hits=[candidate_hit])]
        structural, control, _ = self.run_case(
            reports, hit_row(TEST_CANDIDATE, 12) + "\n", control=False
        )
        self.assertTrue(any("HSP identity mismatch" in x for x in structural))
        self.assertEqual(control, [])

    def test_json_taxid_must_be_an_integer_not_a_digit_string(self):
        candidate_hit = json_hit(12)
        candidate_hit["description"][0]["taxid"] = "1"
        reports = [search_report(TEST_CANDIDATE, 12, hits=[candidate_hit])]
        structural, control, _ = self.run_case(
            reports, hit_row(TEST_CANDIDATE, 12) + "\n", control=False
        )
        self.assertTrue(any("valid subject identity/length" in x for x in structural))
        self.assertEqual(control, [])

    def test_more_than_one_hundred_subjects_is_structural(self):
        hits = [json_hit(12, accession=f"REP{index:05d}") for index in range(101)]
        reports = [search_report(TEST_CANDIDATE, 12, hits=hits)]
        structural, control, _ = self.run_case(reports, "", control=False)
        self.assertTrue(any("max_target_seqs=100" in x for x in structural))
        self.assertEqual(control, [])

    def test_mode_contract_cannot_omit_required_control_metadata(self):
        reports = [
            search_report(TEST_CANDIDATE, 12),
            search_report(TEST_CONTROL, 10),
        ]
        with self.assertRaisesRegex(ValueError, "control ID/spec contract"):
            self.run_case(
                reports, "", control=True, omit_control_metadata=True
            )

    def test_nonviral_and_panax_nt_partitions_have_distinct_query_controls(self):
        nonviral = set(VALIDATOR.MODE_CONTROL_CONTRACT["nt_nonviral"])
        panax = set(VALIDATOR.MODE_CONTROL_CONTRACT["nt_panax"])
        self.assertEqual(
            nonviral - panax, {"PNX_NonPanax_mtDNA_control"}
        )
        self.assertEqual(panax - nonviral, set())

    def test_nucleotide_mode_rejects_non_iupac_subject_sequence(self):
        candidate_hit = json_hit(12, sseq="Z" * 12)
        reports = [search_report(TEST_CANDIDATE, 12, hits=[candidate_hit])]
        values = hit_values(TEST_CANDIDATE, 12)
        values.update({"sseq": "Z" * 12, "pident": "0.0"})
        row = "\t".join(values[field] for field in VALIDATOR.HIT_FIELDS)
        structural, control, _ = self.run_case(
            reports, row + "\n", control=False, mode="nt_tsa"
        )
        self.assertTrue(any("aligned query/subject" in x for x in structural))
        self.assertEqual(control, [])

    def test_json_hsp_sequences_must_be_strings(self):
        candidate_hit = json_hit(12)
        candidate_hit["hsps"][0]["hseq"] = None
        reports = [search_report(TEST_CANDIDATE, 12, hits=[candidate_hit])]
        structural, control, _ = self.run_case(
            reports, "", control=False
        )
        self.assertTrue(any("non-string JSON HSP sequence" in x for x in structural))
        self.assertEqual(control, [])

    def test_malformed_result_json_cli_returns_structural_exit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result_path = root / "RESULTS.json"
            expected_path = root / "EXPECTED_QUERIES.json"
            hits_path = root / "HITS.tsv"
            result_path.write_text("")
            expected_path.write_text(
                json.dumps({"queries": [{"id": "candidate", "length": 12}]})
            )
            hits_path.write_text("")
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                rc = VALIDATOR.main(
                    ["validator", str(result_path), str(expected_path), "mode", str(hits_path)]
                )
            self.assertEqual(rc, VALIDATOR.EXIT_STRUCTURAL)
            self.assertIn("malformed result JSON", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
