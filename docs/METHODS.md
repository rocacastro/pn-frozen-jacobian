# Methods and cost model

This document specifies the implemented operations; it is not a replacement for
the convergence proofs in the accompanying article. All inverses below denote linear solves.
The actual float64 code uses `scipy.linalg.lu_factor` and `lu_solve`, preserving
both factors and pivot indices for reuse.

## Ordinary family and initialization

For $N\ge2$, the stationary P_N cycle is

$$u_k=x_k-A_{k-1}^{-1}F(x_k),\quad A_k=J(u_k),\quad
z_{k,0}=x_k,\quad z_{k,j+1}=z_{k,j}-A_k^{-1}F(z_{k,j}),\quad x_{k+1}=z_{k,N}.$$

In delayed initialization D, $u_0=x_0$, so the first cycle is S_N. In immediate
initialization I, J(x0) is factored before constructing the first predictor.
The reference R-order is $\rho_N=(N+2+\sqrt{(N+2)^2-8})/2$, not an assertion of
an exact order for every trajectory. Shamanskii S_m is the same frozen correction
block with $A_k=J(x_k)$ and no inherited predictor.

## Fused family

Define $E=I-A^{-1}K$ and $\mathcal H_q=(I+E+\cdots+E^{q-1})A^{-1}$.
To apply H_q to b, solve $Av_0=b$, set $h=v_0$, and repeat
$v_l=v_{l-1}-A^{-1}Kv_{l-1}$, $h\leftarrow h+v_l$ for $l=1,\ldots,q-1$.
This uses q solves and q-1 matvecs, without forming an inverse or matrix powers.

A fused macrocycle first predicts with the inherited H_q, factors the Jacobian
there, performs N frozen corrections to y, predicts $v=y-A^{-1}F(y)$,
evaluates $K=J(v)$ without factoring K, and performs N corrections with the new
H_q starting at y. D omits only the inherited predictor of the first macrocycle.
Phat_N denotes q=N+2. The reference order is

$$p_{N,q}=\begin{cases}
[(N+1)(q+1)+\sqrt{(N+1)^2(q+1)^2-8q}]/2,&2\le q\le N,\\
[(N+1)(N+2)+\sqrt{(N+1)^2(N+2)^2-16(N+1)}]/2,&q=N+1,\\
\rho_N^2,&q\ge N+2.
\end{cases}$$

Pcompose executes two ordinary P_N cycles and tests the residual only after the
complete pair. It is the structural reference for fusion, not a single P_N cycle.

## External comparators

M6 (Wang--Li) uses y=x-J(x)^{-1}F(x) and the action
$H=(2I-J(x)^{-1}J(y))J(x)^{-1}$, then z=y-HF(y) and xnew=z-HF(z).
M8 (Cordero--Torregrosa--Vassileva) uses s=J(x)^{-1}F(x), y=x-s/2,
z=x-2s/3, B=J(x)-3J(z), u=y+B^{-1}F(x), v=u+2B^{-1}F(u), and
xnew=v+2B^{-1}F(v). B is factored once and reused. Constant vector factors are
precomputed. The code counts the n-by-n scaling by 3 and vector scalings.
The bibliographic identities are those of the manuscript, not new priority claims.

## Equivalent products

$\mathcal L_\kappa=n(2n-1)(n-1)/6+\kappa n(n-1)/2$ and
$\mathcal T_\kappa=n(n-1)+\kappa n$. A full F evaluation costs n mu0,
a full J costs n^2 mu1, and a dense matvec costs n^2. Sums, pivot searches,
permutations, memory traffic and cache behavior are outside this algebraic model;
they are not outside measured wall-clock time.

$$\begin{aligned}
C_{P_N}&=Nn\mu_0+n^2\mu_1+\mathcal L_\kappa+(N+1)\mathcal T_\kappa,\\
C_{S_m}&=mn\mu_0+n^2\mu_1+\mathcal L_\kappa+m\mathcal T_\kappa,\\
C_{\widehat P_{N,q}}&=2Nn\mu_0+2n^2\mu_1+\mathcal L_\kappa+(N+1)(q+1)\mathcal T_\kappa+(N+1)(q-1)n^2,\\
C_{M6}&=3n\mu_0+2n^2\mu_1+\mathcal L_\kappa+5\mathcal T_\kappa+2n^2+2n,\\
C_{M8}&=3n\mu_0+2n^2\mu_1+2\mathcal L_\kappa+4\mathcal T_\kappa+n^2+4n.
\end{aligned}$$

The log index is eta=log(p_ref)/C. Continuous unimodality selects adjacent
integer release preparations for the optimum of P_N and S_m. The optimum q for fixed N
is searched only over 2,...,N+2; larger q cannot improve this reference index.

For r delayed ordinary cycles, resources are `(rN+1,r,r,r(N+1)-1,0)` in the
order `(F,J,LU,solves,matvecs)`, including one terminal F. For r delayed fused
macrocycles they are
`(2rN+1,2r,r,r(N+1)(q+1)-q,r(N+1)(q-1)-(q-1))`.
For r delayed composed macrocycles they are `(2rN+1,2r,2r,2r(N+1)-1,0)`.
These formulas count a complete execution, not each iteration plus an extra F.

## Canonical fields

Use alpha_j=j/10, epsilon=x-alpha, and matrices A=4I+RA/sqrt(n), B=RB/sqrt(n),
C=RC/sqrt(n). The initial point is $x^{(0)}=\alpha+\delta d/\|d\|_2$, where
$d_j=\sin(\sqrt{2}j)+\tfrac12\cos(\sqrt{3}j)$.
The random construction is specified in PROTOCOLS.md. Functions
and powers below act componentwise.

$$\begin{aligned}
H_1&=A\epsilon+0.8B\epsilon^2+0.3C\epsilon^3,\\
H_2&=A\epsilon+0.5B(\exp\epsilon-1-\epsilon)+0.4C(\sin\epsilon-\epsilon),\\
H_3&=A\epsilon+0.6B(\arctan\epsilon-\epsilon)+0.25C\epsilon^2,\\
H_4&=A\epsilon+0.5B\log(1+\epsilon^2)+0.3C(\cos\epsilon-1+\epsilon^2/2),\\
H_5&=A\epsilon+0.4B(\exp\epsilon-1-\epsilon)+0.3C\log(1+\epsilon^2).
\end{aligned}$$

The pairs (mu0,mu1) are `(3n+2,2+1/n)`, `(3n+30,2+30/n)`,
`(3n+15,2+(1+kappa)/n)`, `(3n+30,2+(15+kappa)/n)`, and
`(3n+31,2+(17+kappa)/n)`. These are declared equivalent-product weights, not
machine-measured function latencies. `problems.py` precomputes fixed matrix
multiples and evaluates analytic Jacobians; `structure` checks them independently
by complex-step differentiation.

Hammerstein uses eight Gauss--Legendre nodes on [0,1],
`a_ij=w_j*min(t_i,t_j)*(1-max(t_i,t_j))`, F_i=5x_i-5-sum(a_ij*x_j^3),
and J_ij=5delta_ij-3a_ij*x_j^2. The initial vector is 0.5 times the all-ones vector.
A and 3A are precomputed, giving mu0=n+3 and mu1=1+1/n. The fresh solver builds
its reference root outside timing; it never uses that root to stop. The archived
high-precision reference text is retained separately.
