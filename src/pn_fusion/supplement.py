"""Build the English translation of the computational supplement (not the paper)."""
from __future__ import annotations
from pathlib import Path
import subprocess
import shutil
from .io import read_csv
from .models import CostModel,KAPPA_POINTS


def generate_supplement(output: Path, figure_directory: Path, compile_pdf: bool=False) -> Path:
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    figure_directory=Path(figure_directory)
    # A self-contained supplement folder: copy only its three English plot panels.
    panels=('supp_figS1_index','supp_figS1_time','supp_figS2_fusion')
    for name in panels:
        source=figure_directory/(name+'.pdf')
        target=output/(name+'.pdf')
        if source.resolve()!=target.resolve(): shutil.copy2(source,target)
    text=r'''\documentclass[11pt]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage[english]{babel}
\usepackage{lmodern,microtype,amsmath,amssymb,booktabs,graphicx,geometry,hyperref}
\geometry{a4paper,margin=2.4cm}
\hypersetup{hidelinks,pdftitle={Supplementary Computational Material: Manuscript v0.41},pdfauthor={Rodrigo Castro Marin}}
\renewcommand{\thetable}{S\arabic{table}}
\renewcommand{\thefigure}{S\arabic{figure}}
\title{Supplementary Computational Material\\\large Predictor-memory frozen-Jacobian methods, manuscript v0.41}
\author{Rodrigo Castro Mar\'in}
\date{}
\begin{document}\maketitle
This document is the English counterpart of the computational supplement accompanying manuscript v0.41. It preserves the numerical breakdowns removed from the main article during editorial compression; no essential proof has been moved here. All times below are archived measurements, not fresh measurements from the repository build. The manuscript itself remains a separate document.

\section{Detailed external-comparator timings}
Tables~\ref{tab:p5} and~\ref{tab:p7} retain the results at $\mathrm{TOL}=10^{-12}$. The main article summarizes all nine configurations per field. Full-precision CSV files, including the other tolerances, are supplied in \texttt{data/reference/p5\_m6\_results.csv} and \texttt{p7\_m8\_results.csv}.
The relative quantities are $\Delta t=100(t_A/t_B-1)$ and $\Delta W=100(W_A/W_B-1)$ for $\kappa=1$; negative values favor method $A$. Parentheses contain the IQR, in the same units as the median. The time separation rule is $|\widetilde t_A-\widetilde t_B|>2(\mathrm{IQR}_A+\mathrm{IQR}_B)$; it is descriptive, not a formal statistical test.
'''
    for N,M,file,label in [(5,6,'p5_m6_pairs.csv','p5'),(7,8,'p7_m8_pairs.csv','p7')]:
        rows=[r for r in read_csv(file) if float(r['tol'])==1e-12]
        rows.sort(key=lambda r:(r['system'],float(r['delta'])))
        text+=rf'''\begin{{table}}[htbp]\centering\small
\caption{{Comparison $P_{N}/M{M}$ at $\mathrm{{TOL}}=10^{{-12}}$. Times and IQRs are in milliseconds.}}\label{{tab:{label}}}
\setlength{{\tabcolsep}}{{4pt}}
\begin{{tabular}}{{cccrrrr}}\toprule
Field & $\delta$ & Cycles & $P_{N}$ time (IQR) & $M{M}$ time (IQR) & $\Delta t$ [\%] & $\Delta W$ [\%]\\\midrule
'''
        for r in rows:
            text+=(f"${r['system'][0]}_{r['system'][1:]}$ & {float(r['delta']):.2f} & {r[f'cycles_P{N}']}/{r[f'cycles_M{M}']} & "
                   f"{float(r[f'CPU_P{N}_ms']):.4f} ({float(r[f'IQR_P{N}_ms']):.4f}) & "
                   f"{float(r[f'CPU_M{M}_ms']):.4f} ({float(r[f'IQR_M{M}_ms']):.4f}) & "
                   f"{100*float(r[f'CPU_relative_P{N}_vs_M{M}']):+.2f} & {100*float(r[f'work_relative_P{N}_vs_M{M}']):+.2f}\\\\\n")
        text+='\\bottomrule\n\\end{tabular}\n\\end{table}\n'
        if N==5: text+='\\clearpage\n'
    text+=r'''\section{Full sensitivity of the optimal members}
Table~\ref{tab:optima} retains all six control points. The main article displays two representative endpoints. Stability throughout $[1,4.5]$ follows from affine normalized-cost differences and the unimodality argument, not from this finite list alone.
\begin{table}[htbp]\centering\scriptsize
\caption{Stationary family optima and 0.5\% near-optimal bands.}\label{tab:optima}
\setlength{\tabcolsep}{4pt}
\begin{tabular}{ccrrrrrrl}\toprule
Field & $\kappa$ & $N_P^*$ & $m_S^*$ & $10^6\eta_P^*$ & $10^6\eta_S^*$ & Advantage [\%] & & P/S bands\\\midrule
'''
    for field in ('H1','H5'):
        for k in KAPPA_POINTS:
            m=CostModel.for_field(field,kappa=k);N,_=m.optimal('P');s,_=m.optimal('S')
            bp=','.join(map(str,m.near_optimal('P')));bs=','.join(map(str,m.near_optimal('S')))
            text+=(f"$H_{field[1]}$ & {k:g} & {N} & {s} & {1e6*m.eta_p(N):.6f} & {1e6*m.eta_s(s):.6f} & "
                   f"{100*(m.eta_p(N)/m.eta_s(s)-1):.2f} & & $\\{{{bp}\\}}/\\{{{bs}\\}}$\\\\\n")
    text+=r'''\bottomrule\end{tabular}\end{table}
\clearpage
\section{Secondary graphical representations}
Figures~\ref{fig:external} and~\ref{fig:fusion} preserve the secondary comparisons. All labels are in English; the underlying measurements are unchanged. A configuration range is not an IQR or confidence interval. No interpolation establishes a universal time crossover or monotonicity.
\begin{figure}[htbp]\centering
\includegraphics[width=.49\textwidth]{supp_figS1_index.pdf}\hfill
\includegraphics[width=.49\textwidth]{supp_figS1_time.pdf}
\caption{$P_7/M8$: stationary index advantage (left) and median with range of relative elapsed time over nine configurations per field (right). Positive index differences and negative time differences favor $P_7$. Time-separation decisions are made separately for each configuration.}\label{fig:external}
\end{figure}
\begin{figure}[htbp]\centering
\includegraphics[width=.79\textwidth]{supp_figS2_fusion.pdf}
\caption{End-to-end comparison of $\widehat P_2$ with $P_2\circ P_2$: median and range of relative elapsed time. The line $n=55$ is the stationary algebraic threshold, not an estimated temporal crossover. Negative values favor fusion.}\label{fig:fusion}
\end{figure}
\paragraph{Provenance note.} The archived fusion work columns for $G_5$ retain a historical evaluation-cost convention described in \texttt{docs/KNOWN\_ISSUES.md}. The elapsed-time values in Figure~\ref{fig:fusion} are unaffected. The repository preserves those inputs and provides the separately recomputed canonical costs; no measurement is silently replaced.
\end{document}
'''
    path=output/'supplementary_computational_material.tex';path.write_text(text,encoding='utf-8')
    if compile_pdf:
        executable=shutil.which('pdflatex')
        if executable is None: raise RuntimeError('pdflatex is required to rebuild the supplement PDF; a compiled PDF is included')
        for _ in range(2):
            run=subprocess.run([executable,'-interaction=nonstopmode','-halt-on-error',path.name],cwd=output,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
            if run.returncode: raise RuntimeError(run.stdout[-5000:])
    return path
