# ChalMolDB v0.9 Performance Baseline

This document records the first v0.9 performance profiling cycle. Measurements are intended for relative comparison inside the same GitHub Actions environment, not as absolute guarantees for every deployment.

## Environment

- GitHub Actions hosted runner
- Ubuntu 24.04
- Python 3.12
- RDKit 2026.03.6
- pandas 3.0.6
- NumPy 2.5.3
- Benchmark analysis size: 200 molecular records

## Initial baseline

Measured before the v0.9 analysis-pipeline optimizations:

| Operation | Time |
|---|---:|
| Curated database load | 0.0617 s |
| Database search | 0.0117 s |
| Public table formatting | 0.0241 s |
| Descriptor processing | 0.2113 s |
| Fingerprint generation | 0.3193 s |
| Similarity search | 0.0312 s |

The main performance cost was fingerprint generation, followed by descriptor processing. Database loading and browsing operations were already small relative to molecular analysis.

## Optimizations

### Single RDKit parse

Descriptor and fingerprint generation previously parsed the same uploaded SMILES in separate passes. v0.9 adds a combined pipeline that reuses one RDKit molecule object per valid record.

Measured result:

- Separate descriptor + fingerprint total: 0.5198 s
- Combined pipeline total: 0.4621 s
- Improvement: 11.1%

### Vectorized fingerprint table construction

Fingerprint exports contain Morgan and MACCS bit columns. The original implementation built thousands of Python dictionary entries per molecule. v0.9 now collects fingerprint vectors as NumPy arrays and constructs the bit table in one matrix operation.

Measured after vectorization:

| Operation | Time |
|---|---:|
| Curated database load | 0.0444 s |
| Database search | 0.0140 s |
| Public table formatting | 0.0178 s |
| Descriptor processing | 0.1660 s |
| Fingerprint generation | 0.0541 s |
| Combined descriptor + fingerprint pipeline | 0.1978 s |
| Similarity search | 0.0138 s |

Fingerprint generation decreased from 0.3193 s to 0.0541 s in this benchmark, an approximately 83% reduction.

The end-to-end descriptor + fingerprint path decreased from the original separate total of about 0.5306 s to 0.1978 s, approximately 63% lower in the measured CI environment.

## Additional changes

- Removed an unnecessary full curated-database DataFrame copy before filtering.
- Replaced selected pandas `iterrows()` loops with lower-overhead tuple iteration.
- Added performance measurement to CI so future changes can be compared against this baseline.
- Added regression coverage to verify that the combined pipeline preserves descriptor, fingerprint, invalid-record, and summary outputs.

## Current interpretation

The curated database browser is not the main performance bottleneck at the current database size. v0.9 performance work should prioritize molecular analysis and rendering before attempting low-value database micro-optimizations.

Future performance changes should preserve the scientific definitions and exported values established by the Scientific Core.
