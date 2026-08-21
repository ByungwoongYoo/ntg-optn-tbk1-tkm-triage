# Panax public-source context audit

## Corrected verdict

`Zh` and `Bh` are retained only as archive group tokens. Their correspondence to the healthy and root-rot conditions is **unresolved**.

The source article describes primary-root RNA-seq from healthy and root-rot-affected *Panax notoginseng*, with three biological replicates per condition and five plants pooled per replicate. Its Figure 1 labels transcriptome samples `CK1`-`CK3` and `ROT1`-`ROT3`. No recovered article text, figure, supplementary file, run record, experiment record, BioSample record, or project record explicitly maps `CK/ROT` to `Zh/Bh`.

Consequently, abundance or detection differences between the two archive-token sets are not interpreted as health-status, root-rot, resistance, susceptibility, infection, or disease associations.

## Run and registry crosswalk

| Archive token | DDBJ/ENA run | NGDC run | DDBJ experiment | NGDC experiment | BioSample | ENA sample | NGDC sample |
|---|---|---|---|---|---|---|---|
| ZhA | DRR853907 | CRR2105996 | DRX831960 | CRX1962292 | SAMD01772848 | DRS610865 | SAMC5870565 |
| ZhB | DRR853908 | CRR2105997 | DRX831961 | CRX1962293 | SAMD01772849 | DRS610866 | SAMC5870566 |
| ZhC | DRR853909 | CRR2105998 | DRX831962 | CRX1962294 | SAMD01772850 | DRS610867 | SAMC5870567 |
| BhA | DRR853910 | CRR2105999 | DRX831963 | CRX1962295 | SAMD01772851 | DRS610868 | SAMC5870568 |
| BhB | DRR853911 | CRR2106000 | DRX831964 | CRX1962296 | SAMD01772852 | DRS610869 | SAMC5870569 |
| BhC | DRR853912 | CRR2106001 | DRX831965 | CRX1962297 | SAMD01772853 | DRS610870 | SAMC5870570 |

The six run, experiment, and BioSample accessions are separate archive objects. They do not establish six individual plants or six independent infections; the article says each biological replicate pooled five plants.

## Why the earlier explicit-mapping result was rejected

The earlier source-context package labeled the mapping `explicit_mapping_found`, but its displayed evidence did not contain an exact `Zh`/`Bh` condition mapping. The alleged `Zh` evidence consisted of substrings inside author and reference names such as `Zihan`, `Zhu`, `Zhang`, `Zheng`, and `Zhou`, captured near unrelated root-rot language. It found no corresponding `Bh` evidence. That substring procedure therefore produced false positives and contradicted the same package's interpretation boundary.

The corrected decision rule requires an exact archive token and an explicit condition assignment in the same source statement or table. No recovered source meets that rule.

## Metadata discrepancies retained for reporting

- The article reports RNA-seq on an Illumina HiSeq 4000. All six DDBJ/ENA experiment records report an Illumina NovaSeq 6000. Both records must be cited as an unresolved discrepancy.
- The article reports GSA accession `CRA029712`. Linked registry objects include NGDC BioProject `PRJCA045816`, DDBJ/ENA project `PRJDB39299`, and DDBJ study `DRP016351`. They are listed as cross-accessions rather than collapsed into one identifier.
- Archive sample titles contain only the `ZhA`-`ZhC` and `BhA`-`BhC` tokens; their deposited attributes do not include health status.
- Registry records were retrieved on 2026-08-21 UTC, but exact XML/JSON
  snapshots have not yet been frozen and hashed. The machine-readable
  crosswalk therefore marks snapshot provenance `not_frozen`; freezing those
  records is a submission blocker.

## Publication-safe source statement

> Six public RNA-seq libraries were deposited under the archive group tokens ZhA-ZhC and BhA-BhC. The source article describes primary-root RNA-seq from healthy and root-rot-affected *Panax notoginseng*, with three biological replicates per condition and five plants pooled per replicate, and labels the transcriptome samples CK1-CK3 and ROT1-ROT3. No recovered article, supplementary file, or archive record explicitly maps CK/ROT to Zh/Bh. Zh and Bh are therefore used only as archive group tokens, and candidate abundance differences between them are not interpreted as associations with health status or root rot.

## Official sources

- Source article: <https://doi.org/10.1186/s12870-026-08239-w>
- PubMed Central full text: <https://pmc.ncbi.nlm.nih.gov/articles/PMC12918092/>
- Genome Sequence Archive portal (`CRA029712`): <https://ngdc.cncb.ac.cn/gsa>
- DDBJ/ENA project (`PRJDB39299`): <https://www.ebi.ac.uk/ena/browser/view/PRJDB39299>
- DDBJ/ENA run records: <https://www.ebi.ac.uk/ena/browser/view/DRR853907> through <https://www.ebi.ac.uk/ena/browser/view/DRR853912>
- ENA XML API pattern: `https://www.ebi.ac.uk/ena/browser/api/xml/{accession}`
