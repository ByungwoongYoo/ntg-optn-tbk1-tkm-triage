# Development pilot — superseded

This directory contains provenance from the 500-record development pilot. It is **not** the final held-out analysis and must not be cited as such.

The authoritative result is now stored in:

- [`../final_v3/README.md`](../final_v3/README.md)
- [`../final_v3/FINAL_REPORT.md`](../final_v3/FINAL_REPORT.md)
- [`../final_v3/results_key.json`](../final_v3/results_key.json)

## Final held-out conclusion

The successful full frozen-corpus run (`31949335813`) produced the prespecified decision:

```text
no_reliable_positive_specificity_signal
```

In 114 strictly eligible articles across nine repeatedly represented diseases, mean Jaccard similarity was 0.1219115689 within diseases and 0.1214275645 between diseases (`Delta=0.0004840045`; one-sided label-permutation `p=0.4335`; bootstrap 95% CI `[-0.0124, 0.0140]`). The result did not become positive under list-size/year stratification, frequency-preserving randomization, alternative similarity metrics, independent symbol extraction, exact-set deduplication, recurrent-hub removal, or disease-label recoverability analysis.

The original development-pilot files remain only to document the analysis-development history.
