# Fourier-domain PTA Detection Statistic + Imhof p-value

This repository implements a **Fourier-domain PTA detection statistic** together with **right-tail p-value evaluation** under \(H_0\) using **Imhof’s method** (generalized \(\chi^2\)). It also includes a lightweight `get_phi` implementation (`PTA_Lite`) to avoid `enterprise`-level constraints, while keeping the same \(\phi\)-structure: per-pulsar red noise (RN) + a common GW process, optionally correlated through the Hellings–Downs ORF.

---

## Features

- **Portable parameter dump**: `param_list.txt` + `chain_1.txt` → `samples.json`
- **Portable PTA model I/O**: save/load nested `pta_model` dict via Feather (`pyarrow`)
- **Custom `get_phi` (PTA_Lite)**:
  - `mode="curn"`: per-pulsar diagonal \(\phi_N\) (no cross-pulsar correlations)
  - `mode="hd"`: full \(\phi_C\) with Hellings–Downs cross correlations
- **Fourier-domain detection statistic** with **dimension reduction** from \(\Delta\phi=\phi_C-\phi_N\):
  - `proj_method="mask"` (recommended): select active support of \(\Delta\phi\)
  - `proj_method="svd"`: thin-SVD basis (optional)
- **p-value** under \(H_0\) using:
  - PSD repair for `Sigma_y` (clip tiny negative eigenvalues)
  - whitening + diagonalization → generalized \(\chi^2\)
  - **Imhof integral** for CDF → right-tail p-value
- **(Optional) demo plot** of GX² PDF under \(H_0\)

---

## Dependencies

From the current implementation, you will need at least:

- `numpy`, `scipy`, `matplotlib`
- `pandas`, `pyarrow`
- `astropy`
- `dill`, `tqdm`

(Additional imports in your environment include `la_forge`.)

---

## Data structures

### `pta_model` (dict)

Top level: `pta_model[psr_name]` is a dict containing:

**For `PTA_Lite` (`get_phi`)**
- `phi`, `theta` : sky position angles in radians
- `Tspan` : timespan in Julian years (internally converted to seconds)
- `nfrequencies` : number of positive Fourier frequencies \(K\) (must be identical across pulsars)

**For `FourierDetectionStatistic`**
- `phiinv` : per-pulsar \(\Phi_0^{-1}\) (vector diag or square matrix)
- `Sigma` / `Sigma_inv` : per-pulsar \(\Sigma_0\), \(\Sigma_0^{-1}\) (vector diag or square matrix)
- `a_hat` : per-pulsar vector of length `fourier_num`

### `par_dict` (dict of floats)

Required keys:

**Global GW**
- `gw_log10_A`, `gw_gamma`

**Per-pulsar RN** (if `rn_name="red_noise"`)
- `{PSR}_red_noise_log10_A`
- `{PSR}_red_noise_gamma`
- `{PSR}_red_noise_log10_kappa`

---

## Quickstart

### 1) Export samples to JSON

Inputs:
- `param_list.txt` (parameter names)
- `chain_1.txt` (chain matrix)

Produces:
- `samples.json` (list of dicts)

### 2) Build \(\phi\) providers (H0 / H1)

```python
K = next(iter(pta_model.values()))["nfrequencies"]

pta_h0 = PTA_Lite(pta_model, mode="curn", components=K, gw_components=5, rn_name="red_noise")
pta_h1 = PTA_Lite(pta_model, mode="hd",   components=K, gw_components=5, rn_name="red_noise")
