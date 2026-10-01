"""Canonical reference indices; these weights are not measured MP timings."""
import math
from .core import method, expected


def mu(name,n,k=1.):
    return {'H1':(3*n+2,2+1/n),'H2':(3*n+30,2+30/n),'H3':(3*n+15,2+(1+k)/n),
            'H4':(3*n+30,2+(15+k)/n),'H5':(3*n+31,2+(17+k)/n)}[name]


def rho(N): return (N+2+math.sqrt((N+2)**2-8))/2


def order(key):
    m=method(key)
    if m.family=='P': return rho(m.N)
    if m.family=='S': return m.N+1
    if m.family=='M6': return 6
    if m.family=='M8': return 8
    if m.family=='COMP' or m.q>=m.N+2: return rho(m.N)**2
    t,d=((m.N+1)*(m.q+1),2*m.q) if m.q<=m.N else ((m.N+1)*(m.N+2),4*(m.N+1))
    return (t+math.sqrt(t*t-4*d))/2


def work(c,name,n,k=1.):
    f,j=mu(name,n,k); L=n*(2*n-1)*(n-1)/6+k*n*(n-1)/2; T=n*(n-1)+k*n
    return c['F']*n*f+c['J']*n*n*j+c['LU']*L+c['solves']*T+c['matvecs']*n*n+c['scales']


def eta(key,name,n,k=1.):
    a,b=expected(method(key),1,n),expected(method(key),2,n)
    return math.log(order(key))/work({s:b[s]-a[s] for s in a},name,n,k)


def optimum(fam,name,n,k=1.):
    f,j=mu(name,n,k); L=n*(2*n-1)*(n-1)/6+k*n*(n-1)/2; T=n*(n-1)+k*n
    B=n*f+T; A=n*n*j+L+(T if fam=='P' else 0); lo=2. if fam=='P' else 1.
    fun=(lambda x: math.sqrt((x+2)**2-8)*math.log(rho(x))-x) if fam=='P' else (lambda x:(x+1)*math.log(x+1)-x)
    if A/B<=fun(lo): candidates=[int(lo)]
    else:
        hi=2*lo
        while fun(hi)<A/B: hi*=2
        for _ in range(100):
            mid=(lo+hi)/2
            if fun(mid)<A/B: lo=mid
            else: hi=mid
        candidates=sorted(set([math.floor((lo+hi)/2),math.ceil((lo+hi)/2)]))
    return max(candidates,key=lambda i:eta(f'{fam}{i}',name,n,k))
