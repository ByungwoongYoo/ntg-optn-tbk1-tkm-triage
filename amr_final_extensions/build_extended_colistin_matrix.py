#!/usr/bin/env python3
"""Extend the frozen colistin pathway scan before examining any new results.

The base audit already covers canonical PhoPQ/PmrAB/CrrAB, lipid-A modification,
and selected envelope genes. This wrapper adds genes and promoter regions reported
in prior experimental evolution, complementation, genomic, or transcriptomic studies
of K. pneumoniae colistin resistance/heteroresistance. Reference-relative calls remain
association-screen features; inclusion here is not evidence that a gene is causal.
"""
from __future__ import annotations

import sys
from pathlib import Path

# Running a file by path puts only this subdirectory on sys.path. Add the repository
# root explicitly so the sibling amr_discovery namespace is importable on Actions.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from amr_discovery import build_colistin_variant_matrix as base

# Literature-prespecified additions. Aliases are intentionally conservative and are
# matched against the pinned reference GenBank annotation by the base implementation.
EXTENDED_TARGETS: dict[str, tuple[str, ...]] = {
    "yciM": ("ycim", "lapb", "lipopolysaccharide assembly protein b", "negative regulator of lps biosynthesis"),
    "lptD": ("lptd", "lps-assembly protein lptd", "lipopolysaccharide assembly protein lptd"),
    "wcaJ": ("wcaj", "undecaprenyl-phosphate glucose phosphotransferase wcaj"),
    "wzc": ("wzc", "tyrosine-protein kinase wzc"),
    "rho": ("rho", "transcription termination factor rho"),
    "ecpR": ("ecpr", "e. coli common pilus transcriptional regulator ecpr", "common pilus transcriptional regulator"),
    "phnC": ("phnc", "phosphonate abc transporter atp-binding protein phnc"),
    "lpxL": ("lpxl", "htrb", "kdo2-lipid a lauroyl acyltransferase"),
    "eptB": ("eptb", "phosphoethanolamine transferase eptb"),
    "acrA": ("acra", "multidrug efflux pump subunit acra"),
    "tolC": ("tolc", "outer membrane channel protein tolc"),
    "ramR": ("ramr", "transcriptional repressor ramr"),
    "soxS": ("soxs", "superoxide response transcriptional activator soxs"),
    "ompA": ("ompa", "outer membrane protein a"),
    "lpxH": ("lpxh", "udp-2,3-diacylglucosamine hydrolase"),
    "lpxK": ("lpxk", "tetraacyldisaccharide 4'-kinase"),
    "waaC": ("waac", "heptosyltransferase i"),
    "waaF": ("waaf", "heptosyltransferase ii"),
}

base.TARGETS.update(EXTENDED_TARGETS)

# The base analysis previously queried only three promoter regions. The extension
# adds a uniform upstream screen for the principal regulatory and literature-defined
# loci. These are reference-relative promoter features, not expression measurements.
for gene in [
    "mgrB", "phoP", "phoQ", "pmrA", "pmrB", "pmrD", "crrA", "crrC",
    "eptA", "eptB", "ugd", "yciM", "lptD", "wcaJ", "wzc", "rho",
    "ecpR", "phnC", "lpxL", "ramA", "ramR", "soxS", "acrA", "acrB",
    "tolC", "ompA",
]:
    base.PROMOTER_UPSTREAM[gene] = 300
base.PROMOTER_UPSTREAM["mgrB"] = 500
base.PROMOTER_UPSTREAM["crrA"] = 500
base.PROMOTER_UPSTREAM["crrC"] = 500

if __name__ == "__main__":
    base.main()
