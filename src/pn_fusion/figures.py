"""English figures from archived measurements. Each panel is a separate plot."""
from __future__ import annotations
from pathlib import Path
import math
import statistics
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .models import CostModel,DIMENSIONS,boundary_p5_m6,rho
from .io import read_csv,write_json


def create_figures(output: Path) -> list[dict]:
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    inventory=[]
    def save(name,description,sources):
        fig=plt.gcf();fig.tight_layout()
        fig.savefig(output/(name+'.pdf'),bbox_inches='tight')
        fig.savefig(output/(name+'.png'),dpi=180,bbox_inches='tight')
        plt.close(fig)
        inventory.append(dict(name=name,description=description,sources=sources))
    def figure():
        plt.figure(figsize=(7.2,4.6))
    figure()
    for field,marker in [('H1','o'),('H5','s')]:
        model=CostModel.for_field(field);Ns=np.arange(2,21)
        values=np.array([model.eta_p(int(N)) for N in Ns]);values/=values.max()
        plt.plot(Ns,values,marker=marker,label=field)
    plt.axhline(.995,linestyle='--',linewidth=1,label='99.5% of maximum')
    plt.xlabel('Number of frozen corrections N');plt.ylabel(r'$\eta_P(N)/\eta_P^{\max}$')
    plt.title(r'Flatness of the stationary optimum ($\kappa=1$)');plt.grid(alpha=.25);plt.legend()
    save('fig01_normalized_index','Manuscript Figure 1',['models.CostModel'])
    rows=read_csv('pn_shamanskii_results.csv')
    for delta in (.1,.3,.6):
        selection=[r for r in rows if r['system']=='H5' and math.isclose(float(r['delta']),delta)]
        for what,ylabel,key,factor in [('time','Median elapsed time [ms]','time_median_s',1000),('cycles','Outer cycles','cycles',1)]:
            figure()
            for family,marker in [('P','o'),('S','s')]:
                rr=sorted([r for r in selection if r['family']==family],key=lambda r:int(r['index']))
                plt.plot([int(r['index']) for r in rr],[factor*float(r[key]) for r in rr],marker=marker,label=r'$P_N$' if family=='P' else r'$S_m$')
            plt.xlabel('Family index N or m');plt.ylabel(ylabel);plt.title(rf'$H_5$, $n=200$, $\delta={delta:.2f}$, TOL $=10^{{-12}}$')
            plt.grid(alpha=.25);plt.legend()
            save('fig02_h5_'+what+'_delta_'+f'{delta:.2f}'.replace('.','p'),'Manuscript Figure 2 panel',['pn_shamanskii_results.csv'])
    figure();x=np.linspace(0,700,600)
    for (field,n),ls in zip(DIMENSIONS.items(),['-','--','-.',':',(0,(5,2,1,2))]):
        a0,a1,a2=boundary_p5_m6(n)
        plt.plot(x,-(a0*x+a2)/a1,linestyle=ls,label=f'n={n} ({field})')
    offsets={'H1':(7,-18),'H2':(7,8),'H3':(7,-18),'H4':(7,8),'H5':(-50,-28)}
    for field,n in DIMENSIONS.items():
        m=CostModel.for_field(field);adv=100*(m.eta_p(5)/(math.log(6)/m.m6())-1)
        plt.scatter([m.mu0],[m.mu1],s=38)
        plt.annotate(f'{field}: {adv:+.3f}%',(m.mu0,m.mu1),textcoords='offset points',xytext=offsets[field],fontsize=8)
    plt.xlim(0,700);plt.ylim(0,5);plt.xlabel(r'$\mu_0$');plt.ylabel(r'$\mu_1$')
    plt.title(r'Equal-reference-index boundaries: $P_5/M6$, $\kappa=1$')
    plt.grid(alpha=.25);plt.legend(fontsize=8,loc='upper left')
    save('fig03_efficiency_region','Manuscript Figure 3; compare each field only with its own dimension',['models.boundary_p5_m6','p5_m6_efficiency_region.csv'])
    station=read_csv('fusion_stationary_pairs.csv')
    for what,key in [('work','work_relative_hat_vs_comp'),('time','CPU_relative_hat_vs_comp')]:
        figure()
        for family,marker in [('G1','o'),('G5','s')]:
            rr=sorted([r for r in station if r['family']==family],key=lambda r:int(r['n']))
            plt.plot([int(r['n']) for r in rr],[100*float(r[key]) for r in rr],marker=marker,label=family)
        plt.axhline(0,linewidth=1);plt.axvline(55,linestyle='--',linewidth=1,label='Algebraic threshold n=55')
        plt.xlabel('Dimension n');plt.ylabel('Relative '+('work' if what=='work' else 'elapsed time')+' [%]')
        plt.title('One stationary macrocycle: fused / composed');plt.grid(alpha=.25);plt.legend(fontsize=8)
        save('fig04_stationary_'+what,'Manuscript Figure 4 panel; negative favors fusion',['fusion_stationary_pairs.csv'])
    figure()
    rr=[r for r in read_csv('p7_m8_indices.csv') if float(r['kappa'])==1]
    plt.plot([int(r['n']) for r in rr],[100*float(r['relative_eta_P7_vs_M8']) for r in rr],marker='o')
    plt.axhline(0,linewidth=1);plt.xlabel('Dimension n');plt.ylabel(r'$100(\eta_{P_7}/\eta_{M8}-1)$ [%]')
    plt.title('Stationary index comparison: P7 / M8');plt.grid(alpha=.25)
    save('supp_figS1_index','Supplement Figure S1, index panel',['p7_m8_indices.csv'])
    figure()
    vals=[];low=[];high=[]
    for field in DIMENSIONS:
        times=[100*float(r['CPU_relative_P7_vs_M8']) for r in read_csv('p7_m8_pairs.csv') if r['system']==field]
        med=statistics.median(times);vals.append(med);low.append(med-min(times));high.append(max(times)-med)
    plt.errorbar(list(DIMENSIONS.values()),vals,yerr=[low,high],fmt='o',capsize=4,label='Median and configuration range')
    plt.axhline(0,linewidth=1);plt.xlabel('Dimension n');plt.ylabel('Relative elapsed time [%]')
    plt.title('P7 / M8: nine configurations per field');plt.grid(alpha=.25);plt.legend(fontsize=8)
    save('supp_figS1_time','Supplement Figure S1, time panel; ranges are not confidence intervals',['p7_m8_pairs.csv'])
    figure()
    for family,marker in [('G1','o'),('G5','s')]:
        means=[];low=[];high=[];dims=(40,55,100,200)
        for n in dims:
            times=[100*float(r['CPU_relative_hat_vs_comp']) for r in read_csv('fusion_end_to_end_pairs.csv') if r['family']==family and int(r['n'])==n]
            med=statistics.median(times);means.append(med);low.append(med-min(times));high.append(max(times)-med)
        plt.errorbar(dims,means,yerr=[low,high],marker=marker,capsize=3,label=family)
    plt.axhline(0,linewidth=1);plt.axvline(55,linestyle='--',linewidth=1,label='Stationary algebraic threshold')
    plt.xlabel('Dimension n');plt.ylabel('Relative elapsed time [%]');plt.title('Fused / composed: end-to-end median and range')
    plt.grid(alpha=.25);plt.legend(fontsize=8)
    save('supp_figS2_fusion','Supplement Figure S2; lines do not establish monotonicity',['fusion_end_to_end_pairs.csv'])
    write_json(output/'figure_manifest.json',inventory)
    return inventory
