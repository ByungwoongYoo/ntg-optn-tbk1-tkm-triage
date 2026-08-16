#!/usr/bin/env python3
"""Tighten the frozen corpus to herbal/natural-product interventions only.

The v2 development smoke test exposed a synthetic-DMARD false inclusion caused by generic
words such as 'compound'. This preregistered patch removes those generic signals and requires
an explicit herbal, traditional-medicine, botanical, natural-product, plant-derived, or
phytochemical context before the held-out full analysis.
"""
from pathlib import Path
import hashlib

path = Path(__file__).resolve().parent / "audit.py"
source = path.read_text(encoding="utf-8")
old = '''INTERVENTION_PATTERNS = [
    r"\\bdecoction\\b", r"\\bformula\\b", r"\\bherbal\\b", r"\\bherb\\b",
    r"\\btraditional chinese medicine\\b", r"\\btcm\\b", r"\\btraditional medicine\\b",
    r"\\bphytochemical\\b", r"\\bmedicinal plant\\b", r"\\bnatural product\\b",
    r"\\bextract\\b", r"\\bcompound\\b", r"\\bflavonoid\\b", r"\\balkaloid\\b",
    r"\\bpolyphenol\\b", r"\\bsaponin\\b", r"\\bterpenoid\\b", r"\\bpolysaccharide\\b",
    r"\\bginseng\\b", r"\\bberberine\\b", r"\\bcurcumin\\b", r"\\bquercetin\\b",
    r"\\bresveratrol\\b", r"\\bastragalus\\b", r"\\bsalvia\\b", r"\\brhubarb\\b",
]'''
new = '''INTERVENTION_PATTERNS = [
    r"\\bdecoction\\b", r"\\bformula(?:tion)?\\b", r"\\bherbal\\b", r"\\bherb(?:al)? medicine\\b",
    r"\\btraditional chinese medicine\\b", r"\\btcm formula\\b", r"\\btraditional medicine\\b",
    r"\\bphytochemical\\b", r"\\bmedicinal plant\\b", r"\\bbotanical\\b",
    r"\\bnatural product\\b", r"\\bnatural compound\\b", r"\\bplant[- ]derived\\b",
    r"\\bherbal extract\\b", r"\\bplant extract\\b", r"\\bextract(?:ed)? from (?:the )?(?:herb|plant|root|rhizome|flower|fruit|seed|bark|leaf|leaves)\\b",
    r"\\bisolated from (?:the )?(?:herb|plant|root|rhizome|flower|fruit|seed|bark|leaf|leaves)\\b",
    r"\\bbioactive (?:component|constituent|ingredient)s?\\b", r"\\bactive ingredient[s]? of\\b",
    r"\\bflavonoid\\b", r"\\balkaloid\\b", r"\\bpolyphenol\\b", r"\\bsaponin\\b",
    r"\\bterpenoid\\b", r"\\bpolysaccharide\\b", r"\\bessential oil\\b",
    r"\\bginseng\\b", r"\\bberberine\\b", r"\\bcurcumin\\b", r"\\bquercetin\\b",
    r"\\bresveratrol\\b", r"\\bastragalus\\b", r"\\bsalvia\\b", r"\\brhubarb\\b",
]'''
if source.count(old) != 1:
    raise RuntimeError(f"Expected exactly one v2 intervention block, found {source.count(old)}")
source = source.replace(old, new)
path.write_text(source, encoding="utf-8")
print(f"patched={path} bytes={len(source.encode())} sha256={hashlib.sha256(source.encode()).hexdigest()}")
