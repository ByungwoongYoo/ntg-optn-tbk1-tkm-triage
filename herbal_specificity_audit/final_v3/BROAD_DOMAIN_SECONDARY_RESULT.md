# Secondary broad-domain analysis

## Why this matters

The primary prespecified analysis asked whether author-reported hub/core/key genes distinguish **exact disease labels**. It did not. A secondary analysis asked a lower-resolution question: do the lists distinguish broad domains such as cancer, cardiovascular, musculoskeletal, neurological, gastrointestinal, and renal/metabolic disease?

## Result

| Quantity | Result |
|---|---:|
| Articles | 229 |
| Broad domains | 6 |
| Within-domain mean Jaccard | 0.1243671152 |
| Between-domain mean Jaccard | 0.1129680619 |
| Delta | 0.0113990534 |
| One-sided label-permutation p | 0.0352 |
| Article-bootstrap 95% CI | [-0.0020, 0.0257] |
| Frequency-preserving null p | 0.0033 |
| Delta after top-10 hub removal | 0.0107213347 |
| Hub-removal permutation p | 0.0045 |
| Gene-only classifier balanced accuracy | 0.2326 |
| Balanced-accuracy chance | 0.1667 |
| Classifier permutation p | 0.0249 |

## Interpretation

The broad-domain result contains modest information: the genes distinguish coarse domains better than exact diseases. However, its article-bootstrap interval includes zero, so it did not satisfy the audit's strict robust-positive rule.

The most accurate overall interpretation is therefore a **disease-resolution ceiling**:

> Author-reported hub/core/key lists were not a reproducible high-resolution fingerprint of exact disease identity, although they retained some low-resolution information about broad disease domains.

This secondary signal prevents the overstatement that the lists contain no disease information whatsoever.
