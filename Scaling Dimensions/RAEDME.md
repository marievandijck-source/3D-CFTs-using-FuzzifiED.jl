## Code Overview

### `O3_ED.jl`
Calculates the scaling dimensions of the O(3) Wilson-Fisher CFT using Exact Diagonalization (ED) on the fuzzy sphere.

### `O3_DMRG.jl`
Calculates the scaling dimensions using the Density Matrix Renormalization Group (DMRG) method.

### `O3_combine.jl`
Collects the `.jld2` output files from both `O3_ED.jl` and `O3_DMRG.jl` into a single combined data file, which serves as input for the figure scripts.

### `Figures.jl`
Generates figures showing the scaling dimensions as a function of system size $R^{-\omega}$ for representations with angular momentum $L = 0, 1, 2$ and spin $S = 0, 1, 2$. The resulting figures are saved as `.pdf` files.
