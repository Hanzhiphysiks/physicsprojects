# Fourier-domain PTA Detection Statistic + p-value

This repository implements a **Fourier-domain PTA detection statistic** together with **right-tail p-value evaluation** under $H_0$ using **Imhof’s method** (generalized $\chi^2$). It also includes a lightweight `get_phi` implementation (`PTA_Lite`) to avoid `enterprise`-level constraints, while keeping the same $\phi$-structure: per-pulsar red noise (RN) + a common GW process with optional Hellings–Downs correlations.

> **Note (current repo layout):** all core implementations live inside Jupyter notebooks under `My Coding/`.
> To reproduce results, run the notebooks in the order listed in **Quickstart** below.

---

## Repository layout

- `My Coding/`
  - `Initial Import.ipynb` (imports + basic setup)
  - `get_phi function.ipynb` (defines `PTA_Lite`, `_rho_flat_tail`, `hellings_downs_orf`, and helpers)
  - `Fourier Detection Statistics.ipynb` (defines `FourierDetectionStatistic` and computes `os_val, y, Q, Sigma_y`)
  - `p value.ipynb` (defines `imhof`, `gx2cdf`, `spectral_and_pvalue_from_yQ`, and optional GX$^2$ PDF plotting)
- `some sample data that you can pl.../`
  - `param_list.txt`
  - `chain (just try the last few eleme...).txt`
- Root
  - `pta_model_gpta_2.feather`
  - `requirements.txt`
  - `README.md`

---

## What’s inside (conceptually)

### 1) Portable parameter dump (`param_list.txt` + chain → `samples.json`)
- Reads parameter names from `param_list.txt`
- Reads the chain matrix from the provided chain text file
- Aligns columns (default: `align="left"`) and writes a list of dict samples to `samples.json`  
  (parameter-name → float value, one dict per sample)

### 2) Portable PTA model I/O (Feather)
Helpers to write/read a nested per-pulsar dictionary to a Feather file:
- `write_fourier_dataframe(pta_model, filename)`
- `read_fourier_dataframe(filename)`

Arrays are converted into JSON-like objects that preserve shape (`{"data": ..., "shape": ...}` for ndim>1).

### 3) `PTA_Lite.get_phi`: custom $\phi$ builder (H0 / H1)
`PTA_Lite` constructs Fourier prior variances in a real (cos/sin) basis:

- Fourier grid:
  - $K =$ `nfrequencies` positive Fourier frequencies per pulsar
  - real basis duplication → vectors length `2K`
- Common GW process:
  - only first `gw_components = G` freqs active → length `2G`, then zero-pad to `2K`
  - uses `_rho_flat_tail(log10_A, gamma, f, df, log10_kappa=None)` (GW uses `log10_kappa=None`)
- Per-pulsar red noise (RN):
  - each pulsar has its own `log10_A`, `gamma`, `log10_kappa`
- Modes:
  - `mode="curn"` (H0): returns a list (len=P) of per-pulsar diagonals `rho_rn + rho_gw_padded`
  - `mode="hd"` (H1): returns the full correlated matrix
    $$\Phi_C = I\otimes \mathrm{diag}(\rho_\mathrm{RN}) + \Gamma\otimes \mathrm{diag}(\rho_\mathrm{GW})$$
    where $\Gamma$ is the Hellings–Downs ORF built from sky positions `(phi, theta)`.

**Required keys in `pta_model[pulsar]` for `PTA_Lite`:**
- `phi`, `theta` (radians)
- `Tspan` (in Julian years; converted internally to seconds)
- `nfrequencies` (same for all pulsars)

**Required keys in `par_dict` for `PTA_Lite.get_phi` (if `rn_name="red_noise"`):**
- Global GW: `gw_log10_A`, `gw_gamma`
- Per pulsar RN:
  - `{PSR}_red_noise_log10_A`
  - `{PSR}_red_noise_gamma`
  - `{PSR}_red_noise_log10_kappa`

### 4) Fourier-domain detection statistic (dimension-reduced)
`FourierDetectionStatistic` computes:
- `os_val`: detection statistic
- `y`: reduced data vector
- `Q`: reduced quadratic form matrix
- `Sigma_y`: covariance of `y` under $H_0$ (used for p-values)

**Expected keys in `pta_model[pulsar]` for the statistic:**
- `phiinv`  (vector-diagonal or square matrix)
- `Sigma` / `Sigma_inv` (vector-diagonal or square matrix)
- `a_hat` (shape `(fourier_num,)`)

**Core steps (as implemented):**
1. Build big block-diagonals: `BigPhi0Inv`, `Sigma0`, `Sigma0Inv`
2. Stack `ahat0` across pulsars
3. Build $\phi_N$ from H0 and $\phi_C$ from H1, then $\Delta\phi=\phi_C-\phi_N$  
   (slicing to the **last-K** coefficients per pulsar via `_idx_lastK`)
4. Transfer covariance and coefficients:
   $$\Sigma^{-1} = \Sigma_0^{-1} + \phi_N^{-1} - \Phi_0^{-1},\qquad
   \hat a = \Sigma\,\Sigma_0^{-1}\,\hat a_0$$
5. Dimension reduction from $\Delta\phi$:
   - `proj_method="mask"`: select active support (recommended/default)
   - `proj_method="svd"`: thin-SVD basis (optional)
6. Reduced statistic:
   - `y = G_A^T \hat a`
   - `Q_num = PhiN_red @ dphi_red @ PhiN_red.T`
   - normalize by `den` to get `Q = Q_num/den` and `os_val = (y^T Q_num y)/den`
7. Output `Sigma_y = G_A^T (phiN - Sigma) G_A`

### 5) p-value via Imhof (generalized $\chi^2$)
Given `Q`, `Sigma_y`, and observed `os_val`, the p-value routine:
- repairs `Sigma_y` to be PSD (clips tiny negative eigenvalues)
- whitens using Cholesky (adds a small diagonal jitter if needed)
- diagonalizes $S=L^T Q L$ so $D=\sum_i \lambda_i z_i^2$
- evaluates $F_D(\mathrm{os\_val})$ via Imhof and returns `p_right = 1 - F_D(os_val)`

Main function:
- `spectral_and_pvalue_from_yQ(Q, Sigma_y, os_val, ...) -> p_right`

### 6) (Optional) Demo: plot GX$^2$ PDF under $H_0$
`plot_pdf(Q, Sigma_y, os_val, ...)` evaluates the GX$^2$ PDF (via `gx2pdf`) and marks `os_val` and the mean $\mu=\sum\lambda$.

---

##
