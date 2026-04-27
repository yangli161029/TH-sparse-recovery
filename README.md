# Truncated Huber Penalty for Sparse Signal Recovery

MATLAB code accompanying the paper

> **L. Yang, S. Morigi, M. K. Ng, Y.-W. Wen.**
> *Truncated Huber Penalty for Sparse Signal Recovery with Convergence Analysis.*
> SIAM Journal on Scientific Computing, **48** (2026), pp. A929–A957.
> [doi:10.1137/25M1748184](https://doi.org/10.1137/25M1748184)

The paper introduces a truncated Huber (TH) penalty that bridges
nonconvex unbiased sparse recovery and differentiable optimization,
and develops a block coordinate descent algorithm with finite-step
convergence guarantees under spark conditions.

This repository contains a **minimal, self-contained subset** of the
authors' code: just enough to reproduce the TH curves of Figure 3
of the paper (decaying-signal recovery, Gaussian and oversampled DCT
sensing matrices, noise-free).

## What is *not* included

For copyright reasons, the implementations of the eight baseline
methods compared against in the paper (L1, TL1, L0/AIHT, MCP, Lp/IRLS,
L1-L2, L1/L2, sorted L1/L2) are **not** redistributed here. Pointers:

| Method | Reference |
|---|---|
| L1 (BPDN) | Glowinski & Marrocco, RAIRO, 1975 — [14] |
| TL1 | Zhang & Xin, *Math. Program.*, 2018 — [58] |
| L0 / AIHT | Blumensath & Davies, *ACHA*, 2009 — [2] |
| MCP | Sun, Chen, Tao, *IET Signal Process.*, 2018 — [41] |
| Lp / IRLS | Lai, Xu, Yin, *SINUM*, 2013 — [25] |
| L1−L2 | Yin, Lou, He, Xin, *SISC*, 2015 — [53] |
| L1/L2 (ADMM) | Tao, *SISC*, 2022 — [43] |
| Sorted L1/L2 | Wang, Yan, Yu, *J. Sci. Comput.*, 2024 — [49] |

## Repository layout

```
TH-sparse-recovery/
├── README.md
├── LICENSE
├── reproduce_figure3_Gaussian.m   driver: Gaussian sensing matrix
├── reproduce_figure3_DCT.m        driver: oversampled DCT sensing matrix
├── outputs/                       created on first run; PNG figures land here
└── functions/
    ├── data_generator_A.m         build Gaussian / oversampled-DCT matrix
    ├── ADMM_L1.m                  ADMM for the L1 problem (warm start)
    ├── shrink.m                   soft-thresholding helper
    └── THuberLagrange.m           the TH method (Algorithm 3.1 / 3.2)
```

## Requirements

- MATLAB R2019b or newer.
- Statistics and Machine Learning Toolbox (uses `mvnrnd` to build
  the correlated Gaussian sensing matrix).

No other toolboxes or external packages are needed.

## How to run

Clone the repository, open MATLAB, change the current folder to the
repository root, and run either driver:

```matlab
reproduce_figure3_Gaussian
reproduce_figure3_DCT
```

Each driver prints the relative recovery error (RRE), runtime, and
iteration count, and saves a 1×2 PNG figure (ground truth + TH
recovery) under `outputs/`.

### Expected output

| Driver | RRE | Iterations | Runtime |
|---|---|---|---|
| `reproduce_figure3_Gaussian` | ≈ 9.6 × 10⁻¹⁶ | ≈ 25 | < 0.1 s |
| `reproduce_figure3_DCT`      | ≈ 6.6 × 10⁻¹⁶ | ≈ 34 | < 0.1 s |

These match the values reported in the TH panels of Figure 3 in the
paper. Tiny variations across MATLAB versions / hardware are normal.

## Citing

If you use this code, please cite the paper:

```bibtex
@article{Yang2026TruncatedHuber,
  author  = {Yang, Li and Morigi, Serena and Ng, Michael K. and Wen, You-Wei},
  title   = {Truncated {H}uber Penalty for Sparse Signal Recovery
             with Convergence Analysis},
  journal = {SIAM Journal on Scientific Computing},
  volume  = {48},
  number  = {2},
  pages   = {A929--A957},
  year    = {2026},
  doi     = {10.1137/25M1748184}
}
```

## License

The code in this repository is released under the MIT License — see
[`LICENSE`](LICENSE).

## Contact

Li Yang — `liyang161029@gmail.com`
School of Mathematics and Statistics, Hunan Normal University.
