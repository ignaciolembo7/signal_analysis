# Diagnostic report: apparent water-contrast collapse in the July 2026 acquisitions

Date: 2026-09-15
Scope: `20220610P` versus `20260706P`, the July 2026 brain acquisitions, local preprocessing and signal-analysis products, and the methodological literature stored under `research/`.
Constraint: this was a read-only diagnostic audit. No source code, processed data, configuration, or notebook output was changed.

> **Resolution adopted on 2026-09-15:** use `g` and `bvalue_g` as the default
> axes throughout signal analysis. These columns retain the non-equally-spaced
> sampling recorded in the b-values. Keep the Thorsten columns only for explicit
> legacy sensitivity checks.

## Executive conclusion

The low water-contrast maxima in `20260706P` are an **analysis-axis artifact**, not an effect introduced by topup and not evidence of restricted water.

The processed master table constructs `g_thorsten` as

```text
g_thorsten(point) = G_thorsten(maximum) * b_step / max(observed b_step)
```

This is correct only when the acquisition includes all ten nominal non-zero gradient steps. The July 2026 long-diffusion-time protocols deliberately acquire fewer steps (`nbvals=9`, `6`, or `8`) because of their gradient/b-value design. Nevertheless, `G thorsten` in the sequence-parameter workbook remains the value for the **nominal tenth step**. Dividing by the last observed step therefore stretches the gradient axis to a maximum that was not acquired.

The distortion is different for N=4 and N=8, so their signals are compared at mismatched physical gradients. This necessarily suppresses the maximum of the reconstructed N=8 minus N=4 contrast. The effect starts at TD=143.4 ms and becomes severe at TD=210 ms, exactly as seen in the figures.

Using the nominal denominator of 10, or using the gradient derived directly from the recorded b-values, restores the expected free-water maximum to approximately 0.46–0.48 at every TD. The same artifact is present in the July 2026 brain acquisitions and is especially clear in the lateral ventricles and the syringe reference.

## Evidence chain

### 1. The anomaly is global, not an ROI or direction problem

In `20260706P`, all four water ROIs and both derived directions show the same trajectory. The variation between water ROIs is only a few thousandths, while the TD-dependent loss is approximately two tenths. This excludes ROI placement, a local susceptibility distortion, and tensor direction as plausible primary causes.

The current fitted peak values are approximately:

| TD (ms) | Water peak, long | Water peak, transverse |
|---:|---:|---:|
| 90.0 | 0.476 | 0.474 |
| 120.0 | 0.478 | 0.476 |
| 143.4 | 0.413 | 0.411 |
| 210.0 | 0.287 | 0.280 |

In `20220610P`, the corresponding water maximum is approximately 0.47 for TD=97.1, 119.1, 142.5, and 209.1 ms. The lower value at TD=75.1 ms is boundary-limited: the maximum occurs at the end of the available gradient range.

### 2. The collapse exists without the signal or pseudo-Huber fits

The N=4 and N=8 measured points were linearly interpolated on their common corrected-gradient support and subtracted, without using the restricted-signal fit, correlation-time estimates, maximum-derived length, or pseudo-Huber model.

Mean direct-interpolation maxima across the four water ROIs were:

| Dataset | TD (ms) | Direct water peak |
|---|---:|---:|
| `20220610P` | 97.1 | 0.462–0.466 |
| `20220610P` | 119.1 | 0.467 |
| `20220610P` | 142.5 | 0.462 |
| `20220610P` | 209.1 | 0.462–0.469 |
| `20260706P` | 90.0 | 0.473–0.476 |
| `20260706P` | 120.0 | 0.478–0.479 |
| `20260706P` | 143.4 | 0.414–0.415 |
| `20260706P` | 210.0 | 0.286–0.292 |

The fitted curves reproduce an artifact already present after mapping the measured points onto the current gradient axis. Pseudo-Huber is downstream of this problem and cannot be its cause.

### 3. The supposedly missing points are deliberate protocol sampling, not truncated files

The new acquisitions contain fewer shells at the affected TDs, but the parameter workbooks explicitly record the same counts:

| TD (ms) | N | Workbook `nbvals` | Observed non-zero steps | Last observed `b_step` |
|---:|---:|---:|---:|---:|
| 90.0 | 4 | 10 | 10 | 10 |
| 90.0 | 8 | 10 | 10 | 10 |
| 120.0 | 4 | 10 | 10 | 10 |
| 120.0 | 8 | 10 | 10 | 10 |
| 143.4 | 4 | 9 | 9 | 9 |
| 143.4 | 8 | 10 | 10 | 10 |
| 210.0 | 4 | 6 | 6 | 6 |
| 210.0 | 8 | 8 | 8 | 8 |

The NIfTI volume counts agree with these designs (`2 + 6 * nbvals` volumes). Moreover, the same sampling pattern occurs in every July 2026 brain acquisition. The files were not accidentally cut off during preprocessing or ingestion.

### 4. Exact mechanism of the gradient-axis error

The `G thorsten` value is converted to an RMS-like magnitude by a factor of square root of two and then distributed across the observed b-step indices. For complete ten-step acquisitions, the present formula is harmless. For shorter acquisitions, the expected and currently assigned fractions are:

| TD (ms) | N | Physical fraction of nominal maximum | Current fraction assigned to final point | Gradient inflation | b-value inflation |
|---:|---:|---:|---:|---:|---:|
| 143.4 | 4 | 9/10 | 9/9 | 10/9 = 1.111 | 1.235 |
| 143.4 | 8 | 10/10 | 10/10 | 1.000 | 1.000 |
| 210.0 | 4 | 6/10 | 6/6 | 10/6 = 1.667 | 2.778 |
| 210.0 | 8 | 8/10 | 8/8 | 10/8 = 1.250 | 1.563 |

The squared factors follow because b is proportional to the square of gradient strength. They match the discrepancies observed in `bvalue_thorsten` numerically.

This also explains why the artifact has the observed time dependence:

- TD=90 and 120 ms: both branches have ten steps, so no error.
- TD=143.4 ms: only N=4 is stretched, producing a moderate loss.
- TD=210 ms: both branches are stretched, but by different factors, producing a large loss.

The relevant construction is in `src/data_processing/result_signals.py`, where `g_thorsten` is multiplied by `b_step / b_step_max`. The ESMRMB contrast workflow then uses `g_thorsten * grad_correction_factor`, so the incorrect scale enters both the signal fits and the common contrast grid. The empirical gradient-correction factors are close between N=4 and N=8 and do not undo this branch-dependent stretch.

### 5. The recorded signals themselves are consistent with free water

A simple monoexponential diagnostic fit was performed on the normalized water signals above the noise floor. On the recorded scanner b-values (`bvalue_orig`), the longitudinal apparent diffusion coefficients remain approximately constant across N and TD:

| TD (ms) | D from N=4 (`10^-3 mm2/s`) | D from N=8 (`10^-3 mm2/s`) |
|---:|---:|---:|
| 90.0 | 2.095 | 2.089 |
| 120.0 | 2.076 | 2.091 |
| 143.4 | 2.110 | 2.106 |
| 210.0 | 2.105 | 2.114 |

On the current `bvalue_thorsten` axis, the apparent coefficients become artificially inconsistent:

| TD (ms) | D from N=4 (`10^-3 mm2/s`) | D from N=8 (`10^-3 mm2/s`) |
|---:|---:|---:|
| 90.0 | 2.127 | 2.123 |
| 120.0 | 2.096 | 2.097 |
| 143.4 | 1.712 | 2.101 |
| 210.0 | 0.762 | 1.372 |

This is a strong internal consistency test: the measured attenuation behaves like water with D approximately `2.1e-3 mm2/s`; the nonphysical N- and TD-dependence appears only after the Thorsten-axis remapping.

### 6. Counterfactual axis reconstruction restores the water invariant

The measured signals were recombined three ways, without writing new derivatives:

1. current `g_thorsten` construction;
2. `g_thorsten` scaled by `max_observed_b_step / 10`, which restores the nominal ten-step denominator;
3. `g` derived directly from the recorded b-values.

Mean direct-interpolation peaks across all water ROIs and directions were:

| TD (ms) | Current `g_thorsten` | Nominal-ten-step `g_thorsten` | Recorded-b-value `g` |
|---:|---:|---:|---:|
| 90.0 | 0.475 | 0.475 | 0.475 |
| 120.0 | 0.478 | 0.478 | 0.479 |
| 143.4 | 0.415 | 0.477 | 0.478 |
| 210.0 | 0.289 | 0.463 | 0.470 |

The small residual difference between the last two reconstructions is compatible with rounded workbook parameters and gradient calibration. The decisive result is that both physically plausible axes recover the expected invariant height.

### 7. Topup is not responsible

The image preprocessing applies the topup field to all DWI volumes with `applytopup --method=jac`. A direct check of the `20260706P`, TD=143.4 ms, N=4 water ROI across raw, denoised, Gibbs-corrected, and topup-applied images found:

- raw water b=0 mean: approximately 1930.5;
- topup water b=0 mean: approximately 1913.7;
- maximum change in the shell-wise curve after normalization to each acquisition's own b=0: approximately `3e-5`.

Topup changes the absolute ROI intensity slightly through spatial resampling/Jacobian modulation, but that nearly common multiplicative effect cancels under the per-acquisition normalization. It cannot explain a decrease from approximately 0.47 to 0.41 or 0.29. More importantly, the exact TD- and N-dependent error factors above arise after signal extraction, in the gradient-coordinate construction.

The TD=210 ms phantom N=4/N=8 pair also has unmatched TE/TR (267/9900 ms versus 228/8800 ms), unlike the corresponding brain pair. This is a genuine protocol confound that should still be reported. However, normalization to each branch's own b=0 removes a branch-wide T2/T1 scale in free water, and the corrected-axis counterfactual restores the expected peak despite that mismatch. TE/TR mismatch is therefore not the primary cause of this particular water-height anomaly.

## Transfer to brains

The problem transfers directly to the July 2026 brains (`ADBN`, `ARVE`, and `SNVN`) because they use the same variable-step protocol and the same downstream axis construction.

### Internal water references

Median direct peaks across the three subjects are:

| ROI | TD (ms) | Current axis | Nominal-ten-step axis | Recorded-b-value axis |
|---|---:|---:|---:|---:|
| Left lateral ventricle | 90 | 0.465–0.472 | 0.465–0.472 | 0.464–0.471 |
| Left lateral ventricle | 143.4 | 0.412–0.416 | 0.469–0.472 | 0.467–0.470 |
| Left lateral ventricle | 210 | 0.273 | 0.460–0.462 | 0.467–0.470 |
| Right lateral ventricle | 143.4 | 0.407–0.415 | 0.464–0.473 | 0.462–0.471 |
| Right lateral ventricle | 210 | 0.269–0.272 | 0.461–0.463 | 0.469–0.470 |
| Syringe | 143.4 | 0.405 | 0.467–0.473 | 0.468–0.473 |
| Syringe | 210 | 0.275 | 0.466 | 0.473–0.474 |

The ventricle and syringe results reproduce the phantom anomaly and its correction almost exactly. This is the strongest evidence that the issue is protocol-coordinate handling rather than phantom-specific physics.

### Corpus callosum ROIs

Across the three July 2026 brains and five callosal ROIs, median model-free peaks are:

| TD (ms) | Direction | Current axis | Nominal-ten-step axis | Recorded-b-value axis |
|---:|---|---:|---:|---:|
| 90.0 | long | 0.459 | 0.459 | 0.458 |
| 90.0 | transverse | 0.299 | 0.299 | 0.300 |
| 120.0 | long | 0.447 | 0.447 | 0.448 |
| 120.0 | transverse | 0.316 | 0.316 | 0.317 |
| 143.4 | long | 0.363 | 0.409 | 0.407 |
| 143.4 | transverse | 0.281 | 0.320 | 0.318 |
| 210.0 | long | 0.253 | 0.390 | 0.397 |
| 210.0 | transverse | 0.230 | 0.337 | 0.341 |

Thus the current analysis underestimates the July 2026 callosal peak by roughly 14–16% at TD=143.4 ms and 48–56% at TD=210 ms in this direct interpolation audit. Any TD trajectory, curve area, centroid, inferred peak position, correlation time, filtered length, or pseudo-Huber result based on the current axis is contaminated for these acquisitions.

This does not imply that the corrected tissue peaks must be TD-invariant; restricted and heterogeneous tissue is expected to behave differently from free water. It means the present apparent late-TD attenuation cannot be interpreted biologically until the gradient coordinate is corrected and all downstream products are regenerated.

## Consistency with the local theory

The local theses provide a clear negative control:

- Saidman (2022) explains that reducing the free diffusion coefficient shifts and dilates the free/tortuosity curve horizontally while preserving its maximum height; adding a restricted component attenuates the maximum.
- Lembo Ferrari (2024), Section 4.1.1, describes a universal collapse for free diffusion and pure tortuosity across diffusion times, while restricted diffusion produces TD-dependent displacement and attenuation.

The observed water collapse would therefore imply an implausible new restriction component if taken literally. The axis counterfactual instead restores the theoretical free-water invariant, so the data and theory are mutually consistent once the coordinate error is removed.

## Consequences for the current ESMRMB analysis

Until the coordinate construction is corrected and the analysis rerun, the following July 2026 results should not be used for scientific interpretation:

- contrast peak and area at TD=143.4 and 210 ms;
- peak position and any transformation to `lcf`, `lc,min`, or correlation time;
- TD-dependent centroids and widths that use the affected gradient/length coordinate;
- pseudo-Huber fits to any of those derived trajectories;
- quantitative brain-versus-`20260706P` comparisons involving the affected TDs.

The raw/preprocessed images, extracted normalized signals, ROI definitions, and topup outputs do not need to be discarded on the basis of this audit. The fault is downstream and deterministic.

## Recommended correction and validation plan

The subsequent project decision was to use `g`/`bvalue_g`, rather than repair and retain the equally-spaced Thorsten axis as the default. The resulting validation plan is:

1. Make `g` the default gradient coordinate for plots, physical signal models, contrast construction, and transformed-length analyses.
2. Make `bvalue_g` the default b-value coordinate for monoexponential fits and gradient calibration.
3. Keep `g_thorsten` and `bvalue_thorsten` available only as explicit legacy alternatives.
4. Add a QC that compares `g` with `g_thorsten` and warns when the latter imposes a different spacing pattern.
5. Add free-water invariants: D should be stable across N/TD on `bvalue_g`, and the N=8 minus N=4 peak should be close to the theoretical free-water value when support contains the maximum.
6. Regenerate the ESMRMB notebook tables and figures for `20260706P` and the July 2026 brains.
7. Keep the TD=210 ms phantom TE/TR mismatch flagged separately after the axis migration.
8. Re-evaluate the general ESMRMB conclusions after regeneration; several previously reported late-TD differences are dominated by this artifact.

## Files audited

- `Data-BIDS/sequence_parameters_phantoms.xlsx`
- `Data-BIDS/sequence_parameters_brains.xlsx`
- `Data-BIDS/derivatives/signal_analysis/den_gr--manual/phantoms/ogse_experiments/master.long.parquet`
- `Data-BIDS/derivatives/signal_analysis/den_gr-topup--manual/phantoms/ogse_experiments/master.long.parquet`
- `Data-BIDS/derivatives/signal_analysis/den_gr-topup--plain/brains/ogse_experiments/master.long.parquet`
- `repos/signal_analysis/src/data_processing/result_signals.py`
- `repos/signal_analysis/src/presentation_analysis/abstract_figures.py`
- `repos/img_preprocessing/steps/step_applytopup.sh`
- `research/tesis_master_Ezequiel L. Saidman_2022.pdf`
- `research/tesis_master_Ignacio_Lembo_Ferrar_2024.pdf`
