import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from virus_hunt.finalization import finalize_panax_sequence_gate as FINALIZER


ROOT = Path(__file__).parents[1]
RUNNER = ROOT / "run_panax_remote_search.sh"
FINALIZER_SCRIPT = ROOT / "finalize_panax_sequence_gate.py"
CANDIDATES = ("PNX_Picorna_A1", "PNX_Picorna_A2", "PNX_Picorna_B")


FAKE_BLAST = r'''#!/usr/bin/env python3
import json,pathlib,sys
program=pathlib.Path(sys.argv[0]).name
if "-version" in sys.argv:
    print(f"{program}: 2.17.0+")
    raise SystemExit(0)
query=sys.argv[sys.argv.index("-query")+1]
database=sys.argv[sys.argv.index("-db")+1]
out=pathlib.Path(sys.argv[sys.argv.index("-out")+1])
out.write_text(json.dumps({"query":query,"program":program,"database":database}))
'''


FAKE_FORMATTER = r'''#!/usr/bin/env python3
import json,pathlib,sys
if "-version" in sys.argv:
    print("blast_formatter: 2.17.0+")
    raise SystemExit(0)
archive=pathlib.Path(sys.argv[sys.argv.index("-archive")+1])
metadata=json.loads(archive.read_text())
records={}; name=None
for raw in pathlib.Path(metadata["query"]).read_text().splitlines():
    if raw.startswith(">"):
        name=raw[1:].split()[0]; records[name]=""
    elif name is not None:
        records[name]+=raw.strip().upper()
accessions={
    "PNX_Panax_L2_control":"YP_009121238.1",
    "PNX_Panax_cpDNA_control":"NC_026447.1",
    "PNX_NonPanax_mtDNA_control":"NC_012920.1",
}
out=pathlib.Path(sys.argv[sys.argv.index("-out")+1])
outfmt=sys.argv[sys.argv.index("-outfmt")+1]
if outfmt=="15":
    reports=[]
    for qid,sequence in records.items():
        hits=[]
        if qid in accessions:
            accession=accessions[qid]
            hits=[{
                "description":[{
                    "id":f"ref|{accession}|",
                    "accession":accession.rsplit(".",1)[0],
                    "title":"mock control", "taxid":1,
                }],
                "len":len(sequence),
                "hsps":[{
                    "query_from":1,"query_to":len(sequence),
                    "hit_from":1,"hit_to":len(sequence),
                    "align_len":len(sequence),"identity":len(sequence),
                    "evalue":0.0,"bit_score":500.0,
                    "qseq":sequence,"hseq":sequence,
                }],
            }]
        params=(
            {"matrix":"BLOSUM62","expect":1e-5,"gap_open":11,
             "gap_extend":1,"filter":"L;","cbs":2}
            if metadata["program"]=="blastp" else
            {"expect":1e-5,"sc_match":2,"sc_mismatch":-3,
             "gap_open":5,"gap_extend":2,"filter":"L;m;"}
        )
        search={
            "query_title":qid,"query_len":len(sequence),"hits":hits,
            "stat":{"db_num":1000,"db_len":500000,"kappa":0.041,
                    "lambda":0.267,"entropy":0.14},
        }
        if not hits: search["message"]="No hits found"
        reports.append({"report":{
            "program":metadata["program"],
            "version":f"{metadata['program'].upper()} 2.17.0+",
            "reference":"mock BLAST reference",
            "search_target":{"db":metadata["database"]},
            "params":params,"results":{"search":search},
        }})
    out.write_text(json.dumps({"BlastOutput2":reports}))
else:
    lines=[]
    for qid,sequence in records.items():
        if qid not in accessions: continue
        accession=accessions[qid]; length=str(len(sequence))
        lines.append("\t".join([
            qid,accession,accession.rsplit(".",1)[0],f"ref|{accession}|",
            "100",length,length,length,"1",length,"1",length,"0","500",
            "100","1","mock","mock control",sequence,sequence,
        ]))
    out.write_text("\n".join(lines)+( "\n" if lines else ""))
'''


class SplitRemoteRunnerTests(unittest.TestCase):
    def run_mode(self, mode):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tools = root / "tools"
            queries = root / "queries"
            out = root / f"remote-{mode}"
            tools.mkdir(); queries.mkdir()
            for name, payload in (
                ("blastp", FAKE_BLAST), ("blastn", FAKE_BLAST),
                ("blast_formatter", FAKE_FORMATTER),
            ):
                path = tools / name
                path.write_text(payload); path.chmod(0o755)
            suffix = "faa" if mode == "protein_nonviral" else "fna"
            filename = (
                "panax_three_partial_orfs.faa"
                if mode == "protein_nonviral" else "panax_three_contigs.fna"
            )
            sequence = "ACDEFGHIKLMN" if suffix == "faa" else "ACGTACGTACGT"
            (queries / filename).write_text(
                "".join(f">{candidate}\n{sequence}\n" for candidate in CANDIDATES)
            )
            env = dict(os.environ)
            env.update({
                "PATH": f"{tools}:{env['PATH']}",
                "PANAX_REMOTE_MAX_ATTEMPTS": "1",
                "PANAX_REMOTE_ATTEMPT_TIMEOUT_SECONDS": "60",
                "PANAX_PROTEIN_NONVIRAL_SPLIT_BUDGET_SECONDS": "300",
                "PANAX_NT_NONVIRAL_SPLIT_BUDGET_SECONDS": "300",
            })
            completed = subprocess.run(
                ["bash", str(RUNNER), mode, "queries", f"remote-{mode}"],
                cwd=root, env=env, text=True, capture_output=True, timeout=60,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            status = json.loads((out / "SEARCH_STATUS.json").read_text())
            self.assertTrue(status["technical_complete"])
            self.assertEqual(set(status["split_results"]), set(CANDIDATES))
            self.assertEqual(
                set(path.name for path in (out / "SPLITS").iterdir()),
                set(CANDIDATES),
            )
            verified = subprocess.run(
                ["sha256sum", "-c", "SHA256SUMS.txt"], cwd=out,
                text=True, capture_output=True,
            )
            self.assertEqual(verified.returncode, 0, verified.stderr)

            collected = root / "collected"
            preflight = collected / "panax-query-preflight"
            preflight.mkdir(parents=True)
            (preflight / filename).write_bytes((queries / filename).read_bytes())
            remote = collected / f"panax-remote-{mode}"
            out.rename(remote)
            failures = (
                FINALIZER.validate_protein_nonviral_split_contract(
                    collected, status
                )
                if mode == "protein_nonviral"
                else FINALIZER.validate_nt_nonviral_split_contract(
                    collected, status
                )
            )
            self.assertEqual(failures, [])

    def test_protein_nonviral_split_runner_end_to_end(self):
        self.run_mode("protein_nonviral")

    def test_nt_nonviral_split_runner_end_to_end(self):
        self.run_mode("nt_nonviral")

    def test_finalizer_direct_workflow_invocation_imports_validator(self):
        completed = subprocess.run(
            [sys.executable, str(FINALIZER_SCRIPT), "--help"],
            cwd=ROOT.parents[1], text=True, capture_output=True, timeout=10,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)


if __name__ == "__main__":
    unittest.main()
