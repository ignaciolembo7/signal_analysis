# ESMRMB 2026 phantom and brain analysis: audit and proposed revision

> **Correction dated 2026-09-15:** the later targeted water audit identified a
> deterministic `g_thorsten` scaling artifact in variable-`nbvals` acquisitions.
> It supersedes this report's earlier attribution of the July 2026 late-TD peak
> attenuation primarily to protocol/TE sensitivity. See
> `esmrmb_20260706_water_axis_diagnostic_report.md`. The TE mismatch remains a
> real protocol confound, but it is not the primary cause of the water-height
> collapse. All July 2026 TD=143.4/210 ms curve-coordinate results should be
> regenerated before the broader conclusions below are used.

Date: 2026-09-14
Scope: processed OGSE/NOGSE data in this project, the current ESMRMB comparison workflow, the submitted draft abstract, and the literature stored in `research/`.

## Executive conclusion

The project is ready for a reproducible **descriptive comparison** of the two phantom acquisitions and for a configurable comparison of brain data against either phantom, with `20260706P` as the default. It is not yet ready to support a strong biological interpretation of the pseudo-Huber transition parameter or of a correlation time inferred from the position of a fitted contrast maximum.

The most reproducible finding in the present phantom data is a displacement of the entire NOGSE contrast profile toward larger model-derived restriction length as diffusion time increases. The center of mass of the positive contrast curve in log-length space increases in every fiber/direction group in both phantom acquisitions. By contrast, peak amplitude behaves in opposite ways in the two acquisitions. The latter is therefore not a stable microstructural result and is particularly compromised by a 39 ms N=4/N=8 echo-time mismatch at TD=210 ms in `20260706P`.

The main recommendation is to make full-curve, minimally parametric outcomes primary and to demote pseudo-Huber `delta` to an exploratory sensitivity analysis. If a physical correlation time is required, it should be estimated by a joint generative fit of the measured signals, with one tissue-level restriction-time distribution and an explicit noise model, rather than by fitting two independent correlation times and transforming an interpolated maximum.

## Data readiness and corrections made

The following processed inputs are available:

| Cohort | Processing profile | Subject/session key | Diffusion times (ms) | Relevant N values | Rotated rows |
|---|---|---|---|---|---:|
| Older phantom | `den_gr--manual` | `20220610P` | 75.1, 97.1, 119.1, 142.5, 209.1 | 4, 8 | 17,280 |
| Newer phantom | `den_gr-topup--manual` | `20260706P` | 90, 120, 143.4, 210 | 4, 8 | 4,848 |
| Brains | `den_gr-topup--plain` by default | six subject keys, ten sheets | 76, 90, 120, 143.4, 210 | 4, 8 | 57,600 |

Two configuration errors or ambiguities had to be addressed before comparison:

1. The notebook's phantom profile pointed to raw DWI data and therefore did not discover the existing processed phantom derivative. It now uses `den_gr-topup` by default, while the cross-dataset section explicitly maps each phantom to its actual processing profile.
2. Both stored phantom alpha summaries had been computed with the brain free-diffusion reference, `D0 = 3.2e-3 mm2/s`, although the abstract specifies `2.3e-3 mm2/s` for phantoms. The analysis loader now harmonizes these values without modifying the source derivative, and the pipeline uses dataset-specific defaults for future runs.

Historical alpha workbooks contain both direct rotated and reconstructed aliases of the long/transverse directions. Those reconstructed rows used the brain axis convention and are not valid for phantoms. The loader now retains the direct rotated rows, and newly generated summaries do not reconstruct aliases unless explicitly requested for a legacy input.

The notebook exports a protocol-compatibility audit. Its most consequential result is:

| Cohort | TD (ms) | TE at N=4 (ms) | TE at N=8 (ms) | Difference |
|---|---:|---:|---:|---:|
| `20220610P` | 209.1 | 232 | 232 | 0 |
| `20260706P` | 210 | 267 | 228 | 39 |
| Comparable brains | 210 | 232 | 232 | 0 |

At the last time point of the newer phantom, the contrast therefore mixes modulation-number and T2-weighting effects. A shared Rician floor is also less defensible for that pair. This point should be excluded from protocol-matched inferential comparisons or reacquired with matched TE.

## What the current data show

### Phantom-to-phantom behavior

For each N=8 minus N=4 contrast curve, the revised analysis reports the peak, normalized-gradient area, equivalent width, gradient centroid/spread, log-length centroid/geometric width, and whether the maximum lies at a sampled boundary.

Across the endpoint diffusion times, the log-length centroid changes as follows:

| Phantom | ROI | Direction | Change in log-length centroid | Spearman rho with TD | Change in peak amplitude |
|---|---|---|---:|---:|---:|
| `20220610P` | fiber1 | long | +26.3% | 0.70 | +23.6% |
| `20220610P` | fiber1 | transverse | +22.4% | 0.70 | +38.1% |
| `20220610P` | fiber2 | long | +30.0% | 0.70 | +5.8% |
| `20220610P` | fiber2 | transverse | +29.2% | 0.70 | +14.0% |
| `20260706P` | fiber1 | long | +11.9% | 1.00 | -38.6% |
| `20260706P` | fiber1 | transverse | +6.0% | 0.80 | -41.8% |
| `20260706P` | fiber2 | long | +15.0% | 1.00 | -40.1% |
| `20260706P` | fiber2 | transverse | +9.2% | 0.80 | -17.9% |

This distinction matters. The spatial shift is directionally consistent in all eight groups; the amplitude trajectory is not. At approximately matched TD=120 ms, the fiber1 peak differs by less than 1% between acquisitions, whereas near 210 ms the newer phantom is approximately 22–44% lower, depending on fiber and direction. The late discrepancy is consistent with the TE confound and must not be presented as evidence of a different restriction regime.

The two phantom acquisitions are not independent biological replicates: they differ in date, protocol, processing branch, available diffusion times, and probably physical state. They can validate qualitative robustness and expose protocol sensitivity, but they do not by themselves supply a population-level variance estimate.

### Brain-to-phantom comparison

Using brains processed with the `plain` segmentation and `20260706P` as the selectable reference phantom, descriptive averages are:

| Cohort | Direction | Mean peak contrast | Mean log-length centroid (um) | Mean normalized-gradient area |
|---|---|---:|---:|---:|
| Brains | long | 0.409 | 8.46 | 0.221 |
| Brains | transverse | 0.321 | 8.36 | 0.173 |
| `20260706P` fibers | long | 0.366 | 7.40 | 0.213 |
| `20260706P` fibers | transverse | 0.351 | 7.22 | 0.203 |

These numbers are descriptive, not evidence that brain restrictions are physically larger. The length coordinate is derived from a model and uses different reference diffusivities for brain and phantom. Cross-material inference should therefore emphasize dimensionless curve coordinates and within-cohort trajectories, and treat micrometer-valued axes as model-dependent.

External-reference alpha values in `20260706P` are approximately 0.515/0.465 for fiber1 and 0.284/0.282 for fiber2 in the long/transverse directions. The main callosal brain regions are approximately 0.278–0.316 longitudinally and 0.162–0.183 transversely. The same-scan water-normalized phantom values are about 0.610/0.553 for fiber1 and 0.337/0.335 for fiber2. Internal normalization removes dependence on an assumed phantom free diffusivity but inherits any water-temperature, flow, partial-volume, or estimator bias.

An important quality-control observation is that externally normalized water alpha is only roughly 0.75–0.85 rather than one. This may reflect temperature/reference mismatch, residual noise-floor bias, or the high-b estimator. Both external and same-scan-normalized results should be shown until this discrepancy is resolved.

### Segmentation sensitivity in brains

The previous notebook results already show that ROI construction materially changes the proposed abstract outcomes. Relative to `plain`, mean alpha was lower by approximately 15% longitudinally and 37% transversely for `erode1`, and by 22% and 48% for `sket1`. Pseudo-Huber `delta` also changed by roughly 10–40%, with missing fits in the more restrictive variants. This is not a minor plotting detail: it means that segmentation uncertainty is currently comparable to or larger than several biological effects of interest.

The abstract should either predefine one segmentation as primary and treat the others as sensitivity analyses, or use a hierarchical model that carries ROI-definition uncertainty into the interval estimates. Choosing the segmentation with the most favorable curve would be circular.

## Why the present correlation-time inference is fragile

### Two independent branch fits do not define one tissue correlation time

The original contrast reconstruction fits the N=4 and N=8 signals independently, including a separate correlation time for each branch, and then finds the maximum of their difference. A single tissue measured at one diffusion time should not acquire two physical restriction times solely because N changes. Separate fits can be useful interpolators, but their correlation times should not be interpreted mechanistically.

The revised code provides a shared-correlation-time alternative and compares it with the separate model using BIC:

| Cohort | Signal pairs | Median R2 shared | Median R2 separate | Strong BIC preference for separate model |
|---|---:|---:|---:|---:|
| `20220610P` | 20 | 0.9983 | 0.9994 | 55.0% |
| `20260706P` | 16 | 0.9967 | 0.9989 | 81.3% |
| Brains | 250 | 0.9715 | 0.9945 | 89.6% |

The shared model's poorer performance is informative: the present single-time/single-compartment signal model cannot explain both branches for much of the data. It does not validate two physical correlation times. Plausible causes include heterogeneous restriction sizes, exchange, branch-specific TE/noise, imperfect gradient calibration, and model misspecification.

### Peak position is an unstable inverse estimator

A fitted contrast maximum is sensitive to the interpolation model, allowed gradient range, noise floor, and whether the maximum lies at a boundary. In the older phantom, 20% of the fitted maxima are at the resampled boundary. Transforming such a maximum into correlation time magnifies model assumptions while hiding identifiability.

This concern is also anticipated in the local methodological literature. Saidman's 2022 thesis explicitly notes that the capabilities and limitations of inferring correlation time from the maximum require deeper study. The original NOGSE work and subsequent size-distribution imaging papers motivate sensitivity to restrictions, but do not make a maximum-derived time automatically identifiable in sparse, heterogeneous in-vivo data.[1–3]

### The pseudo-Huber transition is outside the sampled experiment

The pseudo-Huber fit is useful as a smooth empirical curve,

`tc(TD) = c + alpha * delta * (sqrt(1 + (TD/delta)^2) - 1)`,

but `delta` has a transition interpretation only when the data span the crossover. They do not:

- `20260706P` fiber estimates are approximately 1.67–3.73 s, while maximum TD is 210 ms.
- `20220610P` fiber estimates are approximately 0.81–1.39 s, while maximum TD is 209.1 ms.
- Brain estimates are typically several hundred milliseconds (median approximately 0.65 s), again beyond the 210 ms acquisition.

On the sampled range the pseudo-Huber function is mainly on its quadratic limb. With alpha fixed, many large `delta` values give nearly indistinguishable curves. A high R2 therefore does not establish a measured transition. In the current default phantom execution, none of the fitted crossovers is sampled and the median `delta / TDmax` is about 25.6 when water ROIs are included.

Consequently, claims such as "earlier transition" or a quantitatively different tortuosity timescale are not supported by these measurements. Replacing pseudo-Huber with another two-parameter sigmoid on the same four or five TD points would rename the identifiability problem, not solve it.

## Recommended outcome hierarchy

### Primary outcomes now supportable

1. **Log-length centroid and geometric width of the positive contrast curve.** These use the full curve and separate a displacement from a broadening. Report the model dependence of the length mapping; also report the equivalent dimensionless gradient centroid.
2. **Normalized-gradient area and peak amplitude as separate outcomes.** Area captures total contrast; peak captures local amplitude. Neither should substitute for the other.
3. **Within-ROI anisotropy.** Use paired long-versus-transverse differences or log-ratios for alpha and curve-shape metrics. Pairing within subject/ROI reduces scanner- and normalization-level variation.
4. **External and internal alpha.** Report `Deff / D0` with the literature/reference D0 and `Deff / Deff_water` as a sensitivity analysis in phantoms. Do not silently mix the brain and phantom D0 values.
5. **Functional curve distance at common TD.** Compare complete normalized contrast curves on their common measured support using an L1/L2 distance or a functional mixed model, with bootstrap intervals. This avoids reducing every curve to one possibly unstable maximum.

For the brain cohort, inference should retain subject as the independent unit. A mixed-effects model or hierarchical bootstrap can include direction, TD and ROI as repeated factors, with subject-level resampling. Thousands of gradient samples or many ROIs from the same subject must not be treated as independent observations.

For the two phantom acquisitions, use matched-TD plots, paired differences, and Bland–Altman-style QC where matching is meaningful. With effectively two sessions, concordance coefficients or p-values would create false precision; protocol discrepancies should be displayed directly.

### Exploratory outcomes requiring qualification

- Maximum-derived correlation time: retain only as a sensitivity analysis, with boundary flags and bootstrap/profile-likelihood intervals.
- Pseudo-Huber `delta`: report only if the acquisition spans the estimated crossover and its profile-likelihood interval is bounded within a prespecified feasible domain.
- Direct size-distribution inversion: scientifically attractive and closer to the NOGSE literature, but it requires an explicit forward kernel, regularization choice, resolution/point-spread analysis, and simulation-based recovery tests before a recovered histogram can be called quantitative.[2,3]

## Proposed mechanistic analysis

If the goal remains a tissue correlation-time or restriction-size distribution, the next model should be fitted jointly to the original signals rather than to maxima extracted from separately fitted curves:

1. Use a shared tissue-level distribution of correlation times or restriction lengths for N=4 and N=8.
2. Permit direction-specific distribution parameters where anisotropy is hypothesized, while sharing parameters that are physically common.
3. Include the actual waveform, TD, gradient strength, TE, and any calibrated gradient correction in the forward model.
4. Fit branch-specific scale/noise terms when TE differs; otherwise constrain them only after testing the constraint.
5. Estimate uncertainty with subject/session bootstrap and profile likelihood or a Bayesian posterior. Validate parameter recovery on synthetic data at exactly the acquired gradient/TD grid.
6. Compare a single-size model, a log-normal distribution, and a minimal mixture using out-of-sample prediction or information criteria. Reject parameters whose recovery simulations are biased or whose intervals run to the bounds.

This follows the physical idea that time-dependent diffusion probes the structural correlation spectrum and that short- and long-time regimes carry different structural information.[4–6] It also prevents an empirical transition curve from being mistaken for a directly measured microstructural timescale.

## Noise model and acquisition recommendations

The current `sqrt(S^2 + C^2)` correction is a practical magnitude-noise approximation. MRI magnitude data are Rician in the single-coil idealization and can be noncentral-chi after coil combination; modern preprocessing can alter this distribution.[7,8] The fitted floor can therefore absorb model mismatch, partial volume, or TE differences.

Recommended next steps are:

- estimate spatial noise from repeated b=0 images or a validated MP-PCA method;
- inspect residuals versus gradient, N, TE, ROI and subject;
- use the appropriate Rician/noncentral-chi likelihood where raw magnitude statistics are retained, or document when preprocessing justifies an approximately Gaussian likelihood;
- acquire N=4 and N=8 at matched TE and preferably matched TR for every TD;
- add diffusion times that bracket the hypothesized transition. A transition above 0.6–1 s cannot be identified from measurements ending at 0.21 s;
- include repeated phantom scans within session and across days, with temperature recorded, to separate repeatability from protocol drift;
- acquire a common water reference and a stable calibration phantom in each session.

## Recommended revision of the abstract's scientific claim

The current data can support a claim of **time-dependent displacement and broadening of NOGSE contrast profiles, with direction-dependent attenuation and substantial protocol/segmentation sensitivity**. They do not yet support a robust numerical transition time or a definitive statement that an empirically fitted pseudo-Huber crossover measures microscopic tortuosity.

A defensible result paragraph would focus on:

- the consistently increasing contrast-curve centroid with TD in both phantom acquisitions;
- longitudinal/transverse differences in brain alpha and profile shape, with subject-level uncertainty;
- agreement and disagreement between phantom sessions at matched TD;
- the failure of the shared single-restriction model as evidence for heterogeneity/model mismatch, not as evidence for two N-dependent tissue times;
- sensitivity to segmentation, reference diffusivity, noise treatment and TE matching.

The strongest route to a genuinely improved abstract is to present pseudo-Huber as the analysis being superseded: the full-curve metrics show which part of the conclusion is reproducible, while the support audit explains why the former transition parameter was overinterpreted.

## Reproducible implementation

The expanded notebook now:

- selects either phantom acquisition independently;
- compares both phantoms without pretending they are segmentation variants;
- compares brains against a selectable phantom, defaulting to `20260706P`;
- harmonizes alpha reference diffusivity and provides same-scan water normalization;
- computes full-curve metrics and boundary flags;
- compares shared and separate correlation-time signal models with AIC/AICc/BIC;
- audits TE/TR and gradient support;
- labels pseudo-Huber as a legacy sensitivity analysis and tests whether its crossover was sampled;
- exports the plot-ready and audit CSV files under `notebooks/outputs/esmrmb_variant_comparison/cross_dataset/`.

## Sources

1. E. L. Saidman, *Master's thesis* (2022), especially Chapter 4.2 and Section 5.3, local file `research/tesis_master_Ezequiel L. Saidman_2022.pdf`.
2. N. Shemesh et al., “Nonuniform oscillating-gradient spin-echo (NOGSE) MRI,” *Journal of Magnetic Resonance* 237 (2013), [doi:10.1016/j.jmr.2013.09.009](https://doi.org/10.1016/j.jmr.2013.09.009).
3. N. Shemesh et al., “Size distribution imaging by non-uniform oscillating-gradient spin echo (NOGSE) MRI,” *NeuroImage* (2015), [PMC4509907](https://pmc.ncbi.nlm.nih.gov/articles/PMC4509907/).
4. L. M. Burcaw, E. Fieremans, and D. S. Novikov, “Mesoscopic structure of neuronal tracts from time-dependent diffusion,” *NeuroImage* 114 (2015), conceptual framework summarized in D. S. Novikov et al., [doi:10.1002/nbm.3998](https://doi.org/10.1002/nbm.3998).
5. D. S. Novikov et al., “Revealing mesoscopic structural universality with diffusion,” *PNAS* 111 (2014), [doi:10.1073/pnas.1316944111](https://doi.org/10.1073/pnas.1316944111).
6. P. P. Mitra et al., “Diffusion propagator as a probe of the structure of porous media,” *Physical Review Letters* 68 (1992), [doi:10.1103/PhysRevLett.68.3555](https://doi.org/10.1103/PhysRevLett.68.3555).
7. H. Gudbjartsson and S. Patz, “The Rician distribution of noisy MRI data,” *Magnetic Resonance in Medicine* 34 (1995), [PMC2254141](https://pmc.ncbi.nlm.nih.gov/articles/PMC2254141/).
8. J. Veraart et al., “Diffusion MRI noise mapping using random matrix theory,” *Magnetic Resonance in Medicine* 76 (2016), [doi:10.1002/mrm.26059](https://doi.org/10.1002/mrm.26059).
9. M. Capiglioni et al., “Probing axonal morphology using time-dependent diffusion MRI,” *Physical Review Applied* 15 (2021), [doi:10.1103/PhysRevApplied.15.014045](https://doi.org/10.1103/PhysRevApplied.15.014045).
10. Project-specific context and prior claims: `research/ESMRMB 2026 Abstract.docx`, `research/tesis_master_Ignacio_Lembo_Ferrar_2024.pdf`, `research/2023-ISMRM_Saidman et al.pdf`, `research/2024-ISMRM_Zwick et al.pdf`, and `research/Capiglioni et al, NOGSE @ 3T.pdf`.
