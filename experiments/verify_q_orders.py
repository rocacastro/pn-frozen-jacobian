"""Independent q-transition COC experiment, with three delayed macrocycles.
Faithful English rewrite of the supplied driver. This untimed 2D verification
uses an explicit 2-by-2 solve, not the timed float64 implementation.
X0 is created from decimal strings after setting each working precision.
N=2,3,4 use 5000,10000,25000 digits, respectively."""
from __future__ import annotations
import csv
import sys
import mpmath as mp
try:
    sys.set_int_max_str_digits(0)
except AttributeError:
    pass

def F(x):
    u, v = x
    return [u + u * u + u * v, v + u * u + u * v + v * v]

def J(x):
    u, v = x
    return [[1 + 2 * u + v, u], [2 * u + v, 1 + u + 2 * v]]

def solve(A, b):
    a, b0 = A[0]
    c, d = A[1]
    det = a * d - b0 * c
    return [(d * b[0] - b0 * b[1]) / det, (-c * b[0] + a * b[1]) / det]

def mv(A, x):
    return [A[0][0] * x[0] + A[0][1] * x[1], A[1][0] * x[0] + A[1][1] * x[1]]

def sub(x, y):
    return [x[0] - y[0], x[1] - y[1]]

def norm(x):
    return mp.sqrt(x[0] * x[0] + x[1] * x[1])

def poly(A, K, b, q):
    v = solve(A, b)
    h = v[:]
    for _ in range(1, q):
        t = solve(A, mv(K, v))
        v = sub(v, t)
        h = [h[0] + v[0], h[1] + v[1]]
    return h

def rho(N):
    return (N + 2 + mp.sqrt((N + 2) ** 2 - 8)) / 2

def pref(N, q):
    if q <= N:
        tr = (N + 1) * (q + 1)
        return (tr + mp.sqrt(tr ** 2 - 8 * q)) / 2
    if q == N + 1:
        tr = (N + 1) * (N + 2)
        return (tr + mp.sqrt(tr ** 2 - 16 * (N + 1))) / 2
    return rho(N) ** 2

def run(N, q, macrocycles=3):
    x0 = [mp.mpf('0.02'), mp.mpf('0.04')]
    x = x0[:]
    A = J(x)
    z = x[:]
    for _ in range(N):
        z = sub(z, solve(A, F(z)))
    y = z
    fy = F(y)
    v = sub(y, solve(A, fy))
    K = J(v)
    w = sub(y, poly(A, K, fy, q))
    for _ in range(1, N):
        w = sub(w, poly(A, K, F(w), q))
    x = w
    Aprev = A
    Kprev = K
    xs = [x0[:], x[:]]
    while len(xs) <= macrocycles:
        fx = F(x)
        u = sub(x, poly(Aprev, Kprev, fx, q))
        A = J(u)
        z = sub(x, solve(A, fx))
        for _ in range(1, N):
            z = sub(z, solve(A, F(z)))
        y = z
        fy = F(y)
        v = sub(y, solve(A, fy))
        K = J(v)
        w = sub(y, poly(A, K, fy, q))
        for _ in range(1, N):
            w = sub(w, poly(A, K, F(w), q))
        x = w
        Aprev = A
        Kprev = K
        xs.append(x[:])
    e = [norm(x) for x in xs]
    c = []
    for k in range(1, len(e) - 1):
        c.append(mp.log(e[k + 1] / e[k]) / mp.log(e[k] / e[k - 1]))
    return c[-1]

def main():
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', default='build/high_precision')
    args = parser.parse_args()
    output_directory = Path(args.outdir)
    output_directory.mkdir(parents=True, exist_ok=True)
    rows = []
    for N, dps in [(2, 5000), (3, 10000), (4, 25000)]:
        mp.mp.dps = dps
        for q in [N + 1, N + 2]:
            obs = run(N, q, 3)
            ref = pref(N, q)
            rows.append([N, q, mp.nstr(ref, 20), mp.nstr(obs, 20)])
            print(N, q, 'ref=', mp.nstr(ref, 14), 'COC=', mp.nstr(obs, 14))
    with open(str(output_directory / 'q_order_verification.csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['N', 'q', 'p_reference', 'COC_observed'])
        w.writerows(rows)
if __name__ == '__main__':
    main()
