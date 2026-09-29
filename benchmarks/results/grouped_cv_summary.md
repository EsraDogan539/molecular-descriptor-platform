# Grouped cross-validation: donor/acceptor-aware splits

Development records: 3086 · unique InChIKeys: 2981 · donors: 60 · acceptors: 52

Merged identical fragments (alias -> group): donors [('D26', 'D12'), ('D50', 'D10')]; acceptors none

Excluded (no acceptor label in source name): EROL_1093 (D21), EROL_1249 (D24)

Values are mean ± SD over folds. Target: Eg_eV (calculated HOMO-LUMO gap).
Hyperparameters are fixed (no tuning); see make_models().

## Split: inchikey

| Representation | Model | MAE (eV) | RMSE (eV) | R² | test records |
|---|---|---|---|---|---|
| - | Mean baseline | 0.409 ± 0.015 | 0.529 ± 0.014 | -0.003 ± 0.004 | 3086 |
| Donor/acceptor one-hot | Ridge | 0.081 ± 0.003 | 0.124 ± 0.016 | 0.944 ± 0.016 | 3086 |
| General | RF | 0.223 ± 0.005 | 0.294 ± 0.006 | 0.689 ± 0.028 | 3086 |
| General | SVR | 0.193 ± 0.006 | 0.261 ± 0.011 | 0.755 ± 0.021 | 3086 |
| General | XGB | 0.196 ± 0.007 | 0.262 ± 0.003 | 0.754 ± 0.017 | 3086 |
| Chalcogen-aware | RF | 0.221 ± 0.011 | 0.296 ± 0.013 | 0.686 ± 0.023 | 3086 |
| Chalcogen-aware | SVR | 0.233 ± 0.012 | 0.312 ± 0.022 | 0.649 ± 0.046 | 3086 |
| Chalcogen-aware | XGB | 0.211 ± 0.008 | 0.285 ± 0.011 | 0.709 ± 0.022 | 3086 |
| Combined | RF | 0.187 ± 0.007 | 0.251 ± 0.008 | 0.773 ± 0.016 | 3086 |
| Combined | SVR | 0.161 ± 0.009 | 0.222 ± 0.016 | 0.823 ± 0.025 | 3086 |
| Combined | XGB | 0.161 ± 0.008 | 0.220 ± 0.011 | 0.826 ± 0.016 | 3086 |

## Split: donor

| Representation | Model | MAE (eV) | RMSE (eV) | R² | test records |
|---|---|---|---|---|---|
| - | Mean baseline | 0.414 ± 0.015 | 0.534 ± 0.020 | -0.056 ± 0.090 | 3086 |
| Donor/acceptor one-hot | Ridge | 0.228 ± 0.039 | 0.304 ± 0.042 | 0.654 ± 0.102 | 3086 |
| General | RF | 0.281 ± 0.027 | 0.370 ± 0.035 | 0.488 ± 0.126 | 3086 |
| General | SVR | 0.312 ± 0.084 | 0.404 ± 0.108 | 0.355 ± 0.391 | 3086 |
| General | XGB | 0.281 ± 0.033 | 0.366 ± 0.043 | 0.495 ± 0.147 | 3086 |
| Chalcogen-aware | RF | 0.247 ± 0.038 | 0.323 ± 0.048 | 0.605 ± 0.134 | 3086 |
| Chalcogen-aware | SVR | 0.270 ± 0.044 | 0.352 ± 0.054 | 0.530 ± 0.168 | 3086 |
| Chalcogen-aware | XGB | 0.249 ± 0.041 | 0.323 ± 0.048 | 0.605 ± 0.131 | 3086 |
| Combined | RF | 0.235 ± 0.031 | 0.315 ± 0.036 | 0.627 ± 0.103 | 3086 |
| Combined | SVR | 0.272 ± 0.058 | 0.350 ± 0.072 | 0.528 ± 0.221 | 3086 |
| Combined | XGB | 0.239 ± 0.038 | 0.312 ± 0.046 | 0.633 ± 0.122 | 3086 |

## Split: acceptor

| Representation | Model | MAE (eV) | RMSE (eV) | R² | test records |
|---|---|---|---|---|---|
| - | Mean baseline | 0.413 ± 0.082 | 0.530 ± 0.103 | -0.069 ± 0.098 | 3086 |
| Donor/acceptor one-hot | Ridge | 0.354 ± 0.104 | 0.452 ± 0.123 | 0.230 ± 0.163 | 3086 |
| General | RF | 0.317 ± 0.044 | 0.411 ± 0.050 | 0.339 ± 0.110 | 3086 |
| General | SVR | 0.292 ± 0.025 | 0.384 ± 0.038 | 0.387 ± 0.214 | 3086 |
| General | XGB | 0.302 ± 0.061 | 0.398 ± 0.076 | 0.382 ± 0.165 | 3086 |
| Chalcogen-aware | RF | 0.317 ± 0.086 | 0.413 ± 0.097 | 0.314 ± 0.319 | 3086 |
| Chalcogen-aware | SVR | 0.314 ± 0.074 | 0.407 ± 0.087 | 0.329 ± 0.320 | 3086 |
| Chalcogen-aware | XGB | 0.314 ± 0.090 | 0.406 ± 0.102 | 0.342 ± 0.288 | 3086 |
| Combined | RF | 0.304 ± 0.071 | 0.405 ± 0.079 | 0.366 ± 0.138 | 3086 |
| Combined | SVR | 0.297 ± 0.042 | 0.391 ± 0.040 | 0.382 ± 0.205 | 3086 |
| Combined | XGB | 0.288 ± 0.069 | 0.385 ± 0.076 | 0.431 ± 0.093 | 3086 |

## Split: donor+acceptor

| Representation | Model | MAE (eV) | RMSE (eV) | R² | test records |
|---|---|---|---|---|---|
| - | Mean baseline | 0.390 ± 0.065 | 0.509 ± 0.082 | -0.031 ± 0.037 | 622 |
| Donor/acceptor one-hot | Ridge | 0.390 ± 0.065 | 0.509 ± 0.082 | -0.031 ± 0.038 | 622 |
| General | RF | 0.360 ± 0.083 | 0.456 ± 0.097 | 0.146 ± 0.277 | 622 |
| General | SVR | 0.417 ± 0.055 | 0.518 ± 0.073 | -0.094 ± 0.248 | 622 |
| General | XGB | 0.378 ± 0.063 | 0.471 ± 0.078 | 0.083 ± 0.282 | 622 |
| Chalcogen-aware | RF | 0.352 ± 0.039 | 0.446 ± 0.044 | 0.193 ± 0.113 | 622 |
| Chalcogen-aware | SVR | 0.348 ± 0.028 | 0.445 ± 0.032 | 0.170 ± 0.218 | 622 |
| Chalcogen-aware | XGB | 0.349 ± 0.030 | 0.439 ± 0.036 | 0.201 ± 0.198 | 622 |
| Combined | RF | 0.352 ± 0.078 | 0.448 ± 0.097 | 0.203 ± 0.144 | 622 |
| Combined | SVR | 0.387 ± 0.069 | 0.488 ± 0.081 | 0.033 ± 0.212 | 622 |
| Combined | XGB | 0.347 ± 0.052 | 0.431 ± 0.063 | 0.242 ± 0.180 | 622 |

Runtime: 187 s
