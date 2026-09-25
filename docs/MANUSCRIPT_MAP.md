# Manuscript v0.41: table and figure map

Table numbers and LaTeX labels below refer to the current 34-page manuscript, not
to v0.40 or an earlier draft. `run.py reproduce` regenerates every export. All
measurements come from the archived data; the program does not fabricate new
timings to fill a manuscript table. CSV exports preserve more digits and may use
one method per row instead of a typeset A/B cell.

| Table | v0.41 label | English content | Export | Sources |
|---|---|---|---|---|
| 1 | `tab:kappa-threshold-sensitivity` | Stationary fusion threshold sensitivity | [`table_01.csv`](../tables/table_01.csv) | `models.fusion_threshold` |
| 2 | `tab:q-transition` | Reference orders and the q-transition trade-off | [`table_02.csv`](../tables/table_02.csv) | `models.rho`, `models.fused_order` |
| 3 | `tab:qstar-fields` | Optimal number of Neumann terms | [`table_03.csv`](../tables/table_03.csv) | `models.CostModel.optimal_q` |
| 4 | `tab:precision-regimes` | Experimental protocols | [`table_04.csv`](../tables/table_04.csv) | `docs/PROTOCOLS.md` |
| 5 | `tab:H-fields` | Canonical dense fields and evaluation costs | [`table_05.csv`](../tables/table_05.csv) | `models.field_costs` |
| 6 | `tab:elementary-products` | Elementary equivalent-product weights | [`table_06.csv`](../tables/table_06.csv) | `Manuscript Table 6; declared weights, not new timings` |
| 7 | `tab:H-structural-audit` | Structural audit of the dense fields | [`table_07.csv`](../tables/table_07.csv) | `structural_summary.csv` |
| 8 | `tab:PN-initialization-results` | Immediate and delayed initialization | [`table_08.csv`](../tables/table_08.csv) | `data/manuscript/initialization_summary.csv; transcribed from Table 8, not raw samples` |
| 9 | `tab:hammerstein-stationary` | Hammerstein stationary model | [`table_09.csv`](../tables/table_09.csv) | `hammerstein_stationary_model.csv` |
| 10 | `tab:hammerstein-practical` | Hammerstein end-to-end comparison | [`table_10.csv`](../tables/table_10.csv) | `hammerstein_pairs.csv` |
| 11 | `tab:optimal-P-S-kappa` | Independently optimized families | [`table_11.csv`](../tables/table_11.csv) | `kappa_pn_shamanskii.csv`, `models.CostModel.optimal` |
| 12 | `tab:obj12-flatness` | Flatness of the stationary optimum | [`table_12.csv`](../tables/table_12.csv) | `continuous_curvature.csv` |
| 13 | `tab:asymptotic-members-finite` | Stationarily optimal members at finite tolerance | [`table_13.csv`](../tables/table_13.csv) | `pn_shamanskii_results.csv` |
| 14 | `tab:practical-minima-P-S` | Observed practical minima of the family sweep | [`table_14.csv`](../tables/table_14.csv) | `pn_shamanskii_results.csv` |
| 15 | `tab:external-indices` | Stationary indices of the external comparators | [`table_15.csv`](../tables/table_15.csv) | `p5_m6_indices.csv`, `p7_m8_indices.csv` |
| 16 | `tab:external-cpu` | External comparisons over all configurations | [`table_16.csv`](../tables/table_16.csv) | `p5_m6_pairs.csv`, `p7_m8_pairs.csv` |
| 17 | `tab:fusion-resource-accounting` | Stationary and transient fusion resources | [`table_17.csv`](../tables/table_17.csv) | `fusion_resources.csv` |
| 18 | `tab:fusion-stationary` | Stationary fusion comparison | [`table_18.csv`](../tables/table_18.csv) | `fusion_stationary_pairs.csv; see KNOWN_ISSUES.md for G5 work coefficients` |
| 19 | `tab:fusion-practical` | End-to-end fusion comparison | [`table_19.csv`](../tables/table_19.csv) | `fusion_end_to_end_pairs.csv; see KNOWN_ISSUES.md for G5 work coefficients` |
| 20 | `tab:Phat2-P10` | Secondary comparison: Phat2 versus P10 | [`table_20.csv`](../tables/table_20.csv) | `phat2_p10_field_summary.csv` |
| 21 | `tab:order-verification` | Independent high-precision order verification | [`table_21.csv`](../tables/table_21.csv) | `high_precision_orders.csv` |
| 22 | `tab:q-coc` | Order verification across the q transition | [`table_22.csv`](../tables/table_22.csv) | `q_order_verification.csv` |

## Figures

| Manuscript item | English files in `figures/` | Basis |
|---|---|---|
| Figure 1 | `fig01_normalized_index` | Normalized P_N stationary index; H1/H5 |
| Figure 2 | `fig02_h5_time_delta_*`, `fig02_h5_cycles_delta_*` | Six independent panels, all three deltas |
| Figure 3 | `fig03_efficiency_region` | P5/M6 dimension-matched boundaries and all five actual field points |
| Figure 4 | `fig04_stationary_work`, `fig04_stationary_time` | Work and time shown separately |
| Supplement Figure S1 | `supp_figS1_index`, `supp_figS1_time` | P7/M8 stationary index and temporal ranges |
| Supplement Figure S2 | `supp_figS2_fusion` | Fusion end-to-end temporal ranges |

Each of the 13 independent panels has a PDF and a PNG. The four logical main
figures account for ten panels; the two supplementary figures account for three.
The supplement assembles S1 into its original two-panel layout. Plot labels are
in English. Exact file-level source mappings are in `figures/figure_manifest.json`.

## Supplementary tables

S1 and S2 retain all fifteen rows of the P5/M6 and P7/M8 TOL=1e-12 breakdowns.
S3 retains all twelve field/kappa rows of stationary family optima. They are
rebuilt from `p5_m6_pairs.csv`, `p7_m8_pairs.csv`, their full method outputs, and
`kappa_pn_shamanskii.csv`. The three-page English PDF does not replace any proof.
The archived G5 work convention is explicitly disclosed in its provenance note.

## Material intentionally not copied into the public repository

The Spanish manuscript, superseded drafts, private referee responses, editorial
reduction reports, and images containing Spanish text are excluded. Their
absence does not remove scientific numerical datasets: the data provenance
manifest records the 46 retained CSVs and the separate initialization transcription.
The original source ZIP is identified by hash but is not duplicated in this
all-English public package.
