# Automatic differentiation for Lax-Wendroff-type discretizations

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15607814.svg)](https://zenodo.org/doi/10.5281/zenodo.15607814)

This repository contains information and code to reproduce the results
presented in the article
```bibtex
@online{babbar2025automatic,
  title={Automatic differentiation for {L}ax-{W}endroff-type discretizations},
  author={Babbar, Arpit and Churavy, Valentin and Schlottke-Lakemper, Michael and Ranocha,  Hendrik},
  year={2025},
  eprint={TODO},
  eprinttype={arxiv},
  eprintclass={math.NA}
}
```

If you find these results useful, please cite the article mentioned above.
If you use the implementations provided here, please **also** cite this
repository as
```bibtex
@misc{babbar2025automaticRepro,
  title={Reproducibility repository for
         "Automatic differentiation for {L}ax-{W}endroff-type
         discretizations"},
  author={Babbar, Arpit and Churavy, Valentin and Schlottke-Lakemper, Michael and Ranocha,  Hendrik},
  year={2025},
  howpublished={\url{https://github.com/Arpit-Babbar/2025_paper_lw_ad}},
  doi={10.5281/zenodo.15607814}
}
```

## Abstract

Lax-Wendroff methods combined with discontinuous Galerkin/flux reconstruction spatial discretization
provide a high-order, single-stage, quadrature-free method for solving hyperbolic conservation laws.
In this work, we introduce automatic differentiation (AD) in the element-local time average flux
computation step (the predictor step) of Lax-Wendroff methods. The application of AD is similar
for methods of any order and does not need positivity corrections during the predictor step.
This contrasts with the approximate Lax-Wendroff procedure, which requires different finite difference
formulas for different orders of the method and positivity corrections in the predictor step for
fluxes that can only be computed on admissible states. The method is Jacobian-free and problem-independent,
allowing direct application to any physical flux function. Numerical experiments demonstrate the order
and positivity preservation of the method. Additionally, performance comparisons indicate that the
wall-clock time of automatic differentiation is always on par with the approximate Lax-Wendroff method.


## Numerical experiments

In order to generate the results from this repository, you need to install [Julia](https://julialang.org).
We recommend using `juliaup`, as detailed in the official website [https://julialang.org](https://julialang.org).

The results have been generated using Julia version 1.10.8, and we recommend installing the same.
Once you have installed Julia, you can clone this repository, enter this directory and start the executable
`julia` with the following steps

```shell
git clone https://github.com/Arpit-Babbar/2025_paper_lw_ad.git
cd 2025_paper_lw_ad
julia --project=.
```

Then enter the following commands to generate all the results

```julia
julia> import Pkg; Pkg.instantiate() # Does not need to be re-run the next time you enter the REPL
julia> include("generate_all.jl") # Generate all data, postprocess 1D profiles and convergence plots
julia> include("plotting/plot_wct.jl") # See wall clock time performance comparing ALW and AD on screen
```

You should start Julia with a single thread to get reliable timings.
You can do this by running `julia --project=. --threads=1` in the shell.
However, the last command in `generate_all.jl` (the Mach 2000 jet flow)
can take quite a while if using a single thread. Thus, you can also execute
the first few commands in `generate_all.jl` in serial, start Julia again
with multiple threads (e.g., `julia --project=. --threads=auto`), and then
run the last command in `generate_all.jl` to generate the Mach 2000 jet flow.

If you wish to visualize the 2D figure, you need [ParaView](https://www.paraview.org)
and its command line version `pvpython`. Then, in your shell, you can run

```shell
pvpython m2000.py
```

All the figures are now ready and available in the following locations:
1. Double rarefaction test: `paper_figures/isentropic/density.pdf`
2. RHD first Riemann problem: `paper_figures/rhd/density.pdf`
3. Convergence analysis of isentropic vortex test for 2-D compressible Euler's equations: `paper_figures/isentropic_error.pdf`
4. Density profile of Mach 2000 astrophysical jet flow: `paper_figures/m2000.png`.


## Authors

- [Arpit Babbar](https://arpit-babbar.github.io) (Johannes Gutenberg University Mainz, Germany)
- [Valentin Churavy](https://vchuravy.dev) (Johannes Gutenberg University Mainz, Germany and University of Augsburg, Germany)
- [Michael Schlottke-Lakemper](https://lakemper.eu) (University of Augsburg, Germany)
- [Hendrik Ranocha](https://ranocha.de) (Johannes Gutenberg University Mainz, Germany)


## License

The code in this repository is published under the MIT license, see the
`LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
