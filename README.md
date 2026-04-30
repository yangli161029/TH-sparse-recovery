# TH-sparse-recovery

Official MATLAB implementation of the paper:

> **L. Yang, S. Morigi, M. K. Ng, Y.-W. Wen.**  
> *Truncated Huber Penalty for Sparse Signal Recovery with Convergence Analysis.*  
> SIAM J. Sci. Comput., **48** (2026), pp. A929–A957.  
> DOI: [10.1137/25M1748184](https://doi.org/10.1137/25M1748184)

**Abstract.** Sparse signal recovery from underdetermined systems presents significant challenges when using conventional L0 and L1 penalties, primarily due to computational complexity and estimation bias. This paper introduces a truncated Huber penalty, a nonconvex metric that effectively bridges the gap between unbiased sparse recovery and differentiable optimization. The proposed penalty applies quadratic regularization to small entries while truncating large magnitudes, avoiding nondifferentiable points at optimal solutions. Theoretical analysis demonstrates that, for an appropriately chosen threshold, any (s,μ)-sparse solution recoverable via conventional penalties remains a local optimum under the truncated Huber function. To solve the optimization problem, we develop a block coordinate descent (BCD) algorithm with finite-step convergence guarantees under spark conditions. Numerical experiments validate the effectiveness and robustness of the proposed method in both synthetic and real scenarios.

---

## Repository layout

```
TH-sparse-recovery/
├── data/
│   ├── blocks.txt               clean 1D piecewise-constant signal (Figure 8)
│   └── blocks_noisy.txt         noisy observation (sigma = 0.5)
├── functions/
│   ├── THuberLagrange.m         TH algorithm (Algorithm 3.1 / 3.2)
│   ├── ADMM_L1.m                ADMM solver for L1 warm start
│   ├── algvetgradiffbdy.m       TH in the gradient domain (Figure 8)
│   ├── Diffmatrixboundarycond.m difference matrix with boundary conditions
│   ├── vecsubprb.m              linear subproblem solver
│   ├── data_generator_A.m       sensing matrix generator (Gaussian / DCT)
│   ├── data_generator_xg.m      sparse signal generator
│   └── shrink.m                 soft-thresholding helper
├── reproduce_figure3_Gaussian.m
├── reproduce_figure3_DCT.m
├── reproduce_figure6_Gaussian.m
├── reproduce_figure6_oDCT.m
├── reproduce_figure8.m
├── LICENSE
└── README.md
```

Only the TH method is included; baseline implementations are not redistributed due to licensing. Output figures are saved to `output/` (created automatically).

---

## How to run

Requires MATLAB R2019b or newer. Clone the repository, set the root as the working directory, and run any driver script:

```matlab
reproduce_figure3_Gaussian
reproduce_figure3_DCT
reproduce_figure6_Gaussian
reproduce_figure6_oDCT
reproduce_figure8
```

---

## Citation

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

---

## Contact

Questions and comments are welcome — feel free to open an issue or reach out directly at `liyang161029@gmail.com`.
