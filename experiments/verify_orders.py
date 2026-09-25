"""High-precision verification of P_N and Phat_N reference orders.
Faithful English rewrite of the supplied order driver, not a CPU benchmark.
The archived driver constructs X0 at mpmath's import precision (53 bits).
That initialization is retained for numerical traceability. Working precision
is then set to 6000 digits for P_N and 30000 for Phat_N. Scalar errors may
have large negative exponents even when the working significand is shorter.
The small linear systems use mp.lu_solve, as in the original untimed driver.
This is intentionally not the float64 LU-reuse performance implementation."""
from __future__ import annotations
import csv
import sys
import mpmath as mp
try:
    sys.set_int_max_str_digits(0)
except AttributeError:
    pass
ALPHA = mp.matrix([mp.mpf('0'), mp.mpf('0')])
X0 = mp.matrix([mp.mpf('0.02'), mp.mpf('0.04')])

def F(x):
    u, v = (x[0], x[1])
    return mp.matrix([u + u * u + u * v, v + u * u + u * v + v * v])

def J(x):
    u, v = (x[0], x[1])
    return mp.matrix([[1 + 2 * u + v, u], [2 * u + v, 1 + u + 2 * v]])

def nrm(x):
    return mp.sqrt(sum((z * z for z in x)))

def rho(N):
    return (N + 2 + mp.sqrt((N + 2) ** 2 - 8)) / 2

def solve(A, b):
    return mp.lu_solve(A, b)

def poly_action(A, K, b, q):
    v = solve(A, b)
    h = mp.matrix(v)
    for _ in range(1, q):
        t = solve(A, K * v)
        v = v - t
        h = h + v
    return h

def coc(xs):
    errors = [nrm(x - ALPHA) for x in xs]
    values = []
    for k in range(1, len(errors) - 1):
        values.append(mp.log(errors[k + 1] / errors[k]) / mp.log(errors[k] / errors[k - 1]))
    return (values, errors)

def run_PN(N, cycles):
    x = mp.matrix(X0)
    Aprev = J(x)
    z = mp.matrix(x)
    for _ in range(N):
        z = z - solve(Aprev, F(z))
    x = z
    xs = [mp.matrix(X0), mp.matrix(x)]
    while len(xs) <= cycles:
        Fx = F(x)
        u = x - solve(Aprev, Fx)
        A = J(u)
        z = x - solve(A, Fx)
        for _ in range(1, N):
            z = z - solve(A, F(z))
        x = z
        Aprev = A
        xs.append(mp.matrix(x))
    return xs

def run_Phat(N, macrocycles):
    q = N + 2
    x = mp.matrix(X0)
    A = J(x)
    z = mp.matrix(x)
    for _ in range(N):
        z = z - solve(A, F(z))
    y = z
    Fy = F(y)
    v = y - solve(A, Fy)
    K = J(v)
    w = y - poly_action(A, K, Fy, q)
    for _ in range(1, N):
        w = w - poly_action(A, K, F(w), q)
    x = w
    Aprev, Kprev = (A, K)
    xs = [mp.matrix(X0), mp.matrix(x)]
    while len(xs) <= macrocycles:
        Fx = F(x)
        u = x - poly_action(Aprev, Kprev, Fx, q)
        A = J(u)
        z = x - solve(A, Fx)
        for _ in range(1, N):
            z = z - solve(A, F(z))
        y = z
        Fy = F(y)
        v = y - solve(A, Fy)
        K = J(v)
        w = y - poly_action(A, K, Fy, q)
        for _ in range(1, N):
            w = w - poly_action(A, K, F(w), q)
        x = w
        Aprev, Kprev = (A, K)
        xs.append(mp.matrix(x))
    return xs

def main():
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', default='build/high_precision')
    args = parser.parse_args()
    output_directory = Path(args.outdir)
    output_directory.mkdir(parents=True, exist_ok=True)
    rows = []
    mp.mp.dps = 6000
    for N in range(2, 9):
        cycles = 5 if N <= 6 else 4
        xs = run_PN(N, cycles)
        c, e = coc(xs)
        ref = rho(N)
        obs = c[-1]
        rows.append(['P_N', N, mp.nstr(ref, 20), mp.nstr(obs, 20), mp.nstr(abs(obs - ref), 8), cycles, mp.nstr(-mp.log10(e[-1]), 12), mp.nstr(-mp.log10(nrm(F(xs[-1]))), 12)])
        print(f'P_{N}: ref={mp.nstr(ref, 12)}, COC={mp.nstr(obs, 12)}')
    mp.mp.dps = 30000
    for N in range(2, 5):
        macrocycles = 4 if N == 2 else 3
        xs = run_Phat(N, macrocycles)
        c, e = coc(xs)
        ref = rho(N) ** 2
        obs = c[-1]
        rows.append(['Phat_N', N, mp.nstr(ref, 20), mp.nstr(obs, 20), mp.nstr(abs(obs - ref), 8), macrocycles, mp.nstr(-mp.log10(e[-1]), 12), mp.nstr(-mp.log10(nrm(F(xs[-1]))), 12)])
        print(f'Phat_{N}: ref={mp.nstr(ref, 12)}, COC={mp.nstr(obs, 12)}')
    with open(str(output_directory / 'high_precision_orders.csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['method', 'N', 'reference_order', 'final_COC', 'abs_COC_minus_ref', 'cycles_or_macrocycles', 'minus_log10_error_final', 'minus_log10_residual_final'])
        w.writerows(rows)
if __name__ == '__main__':
    main()
