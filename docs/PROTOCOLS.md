# Experimental protocols

All practical solves use float64, maximum 20 outer cycles/macrocyles, and the sole
stopping test $\|F(x_k)\|_2\le\mathrm{TOL}$. The known root is used only for a
reported error. No inner correction terminates a block early. Terminal residuals
are cached for the next block. All counters include the final F evaluation,
except stationary one-macrocycle measurements, which do not test convergence.

## Suites

| CLI suite | Fields and dimensions | Parameters | Warmups / repetitions |
|---|---|---|---|
| `initialization` | P3/H1, P9/H5, P5 and P7 on H1--H5 | delta=0.10, TOL=1e-12; I versus D | 3 / 51 |
| `sweep` | H1 n=20; H5 n=200 | P2--P20, S1--S20; three deltas; TOL=1e-12 | 2 / 31 |
| `p5-m6` | H1--H5 | P5/M6; three deltas and three tolerances | 2 / 31 |
| `p7-m8` | H1--H5 | P7/M8; three deltas and three tolerances | 2 / 31 |
| `hammerstein` | n=8 | x0=0.5 times the all-ones vector; P5/M6 and P7/M8; three tolerances | 2 / 31 |
| `fusion` | G1,G5; n=40,55,100,200 | Phat2/P2composeP2; three deltas and three tolerances | 2 / 31 |
| `stationary` | G1,G5; n=40,55,100,200 | Phat2/P2composeP2; inherited state prepared outside timing | 2 / 31 |
| `phat2-p10` | H3 n=100; H5 n=200 | Phat2/P10; three deltas and three tolerances | 2 / 31 |

The three deltas are 0.10, 0.30 and 0.60; the three tolerances are 1e-8, 1e-10
and 1e-12. Phat2 always means $\widehat P_{2,4}$. All performance comparisons use
delayed initialization D, except the explicit initialization sensitivity study.
G1 and G5 use the H1 and H5 field definitions at the indicated varying dimensions.

The fresh suite contains 646 method configurations in total: 24 initialization,
234 family-sweep, 90 P5/M6, 90 P7/M8, 12 Hammerstein, 144 fusion end-to-end,
16 stationary, and 36 Phat2/P10. These are not 646 distinct roots/problems.

## Randomness and timing

Each dimension uses one PCG64 stream with seed `20260923+n`; matrices RA, RB and
RC are sampled sequentially. Matrix preparation and fixed scalar matrix multiples
are outside the timed region. The vector initial guess follows the manuscript
exactly, independently of later timings.

The English timing runner interleaves methods within a group, with shuffle seed
`20260923+n+round(1000*delta)+round(-log10(TOL))`. Initialization alternates I/D
and D/I. The order and raw elapsed time of each repetition are recorded. This is
a documented fresh schedule; identity with missing historical schedules is not
asserted. A result that fails the stopping test is marked nonconvergent and never
receives a resolved temporal-win flag; ratios for such pairs are only diagnostics.

No BLAS thread count is forced unless `--threads` is supplied. `--quick` uses only
the first group of each requested suite, at most three repetitions and one warmup;
its output is an execution check, not evidence for article timing conclusions.

The descriptive separation filter is

$$|\widetilde t_A-\widetilde t_B|>2(\operatorname{IQR}_A+\operatorname{IQR}_B).$$

It is not a statistical significance test. Relative times use A/B as denominator
convention: $100(\widetilde t_A/\widetilde t_B-1)$. Negative values favor A.
Never turn a reported 37% excess of A over B into a 37% saving of B over A.

## High-precision verification

The stress field is $F(u,v)=(u+u^2+uv,\ v+u^2+uv+v^2)^T$, with root (0,0).
The main protocol starts at (0.02,0.04), uses 6,000 decimal digits for P_N and
30,000 for Phat_N, and reports the last non-saturated COC. The q-transition program
keeps its distinct precisions and fixed macrocycle counts. Read the two programs
for their explicit loops. Neither program contributes performance timing data.

The clock and computational scope differ deliberately: the untimed order driver
may solve its 2-by-2 linear systems through mpmath's own solve, whereas float64
performance methods explicitly reuse SciPy's LU factors and pivots.
