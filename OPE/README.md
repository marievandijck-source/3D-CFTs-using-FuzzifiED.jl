# OPE Coefficients

The OPE coefficient $f_{\sigma\epsilon\sigma}$ is calculated using conformal perturbation theory, following the methods described in [arXiv:2510.09755](https://arxiv.org/abs/2510.09755).

The energy shift of the ground state due to a small perturbation $g_\varepsilon(h)$ of the coupling is given by:

$$\delta E_\sigma(R, h) = \frac{c}{R} \Delta_\sigma(R) + g_\varepsilon(h) \, f_{\sigma\varepsilon\sigma}(R)$$

The OPE coefficient is then extracted as:

$$f_{\sigma\varepsilon\sigma}(R) = \frac{\partial \, \delta E_\sigma(R)}{\partial \, g_\varepsilon(R)} \Bigg|_{g_\varepsilon(R)=0}$$

Since $g_\varepsilon(h) = 0$ at the critical point $h = h_c$, the function $g_\varepsilon(h)$ is evaluated at small deviations from $h_c$ for different states and system sizes.

## Code Overview

### `g_function.jl`
Calculates $g_\varepsilon(h)$ at small deviations from the critical coupling $h_c$, for different states and system sizes.

### `OPE.jl`
Calculates the OPE coefficients $f_{\sigma\varepsilon\sigma}$ using the values of $g_\varepsilon(h)$ obtained from `g_function.jl`.

### `Functions.jl`
Contains helper functions used across the other scripts.

### `combine.jl`
Combines the `.jld2` output files from `OPE.jl` into a single file for use in figure generation.

### `plot.jl`
Generates figures of the OPE coefficients for representations with angular momentum $L = 0, 1, 2$ and spin $S = 0, 1, 2$. The resulting figures are saved as `.pdf` files.
