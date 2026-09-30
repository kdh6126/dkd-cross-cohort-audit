# Reproducing the analysis

Everything runs through one entry point. There is no notebook and no manual step.

```bash
pip install -r requirements.txt
# stage `wgcna_ref` additionally needs R (4.6) with the WGCNA package (1.74);
# set RSCRIPT to the Rscript path if it is not the default Windows location.
python run_all.py --list          # 133 stages, their status and estimated cost
python run_all.py --all           # the 42 stages that produce the reported results
```

Cold-start cost is roughly 9 hours on one CPU, dominated by three stages: `baselines`
(120 min), `bootstrap_ci` (90 min) and `deconfound` (45 min).

## Options

| flag | effect |
|---|---|
| `--list` | print the stage table with status and estimated minutes (also the default) |
| `--stage NAME` | run one stage; repeatable |
| `--from NAME` | run that stage and everything after it |
| `--all` | run every reproduction-path stage, in order (`--paper` runs the 104 the paper rests on) |
| `--paper` | only the 27 stages the manuscript's figures and tables depend on |
| `--force` | re-run even when the declared outputs already exist |
| `--dry-run` | print the plan and the exact commands without executing |
| `--skip-network` | skip the 16 stages that need the internet (downloads and API queries) |

A stage is "done" when every file it declares as an output exists, so an interrupted run
resumes where it stopped and `--all` is safe to repeat. `--force` is the only way to overwrite.
A stage that declares no outputs is never treated as done.

Each stage's combined stdout and stderr is copied to `logs/<stage>.log`, overwritten per run.

## Stage groups

`acquire` downloads. `analysis` is everything the paper rests on. `extra` (prefixed `x_`) are
five side investigations that were run during the work and are kept reproducible, but no
manuscript claim depends on them: an RBS component ablation, GO/KEGG over-representation, a
first simulation superseded by `sim_blocks`, an early candidate validation superseded by
`master_table`, and a v1-vs-v2 selector figure. They are excluded from `--all` and `--paper`
and are reachable with `--stage`.

## Network

Four stages reach outside: `download` (GEO series, KPMP atlas, annotations, GWAS catalog),
`literature` (PubMed E-utilities), `kpmp_evidence` (KPMP API) and `x_enrich` (Enrichr).
`--skip-network` runs the rest against what is already in `data/`. No credentials are needed;
every source is public.

## Order

Stages are declared in dependency order and `--all` respects it. The spine:

```
download -> genespace -> harmonise -> {null, baselines, proposed} -> compare
harmonise -> artifact* -> module_size -> sim_blocks -> sim_real
compare   -> bootstrap_ci -> cohort_sensitivity            (the paper's lead result)
compare   -> deconfound -> ruv -> candidates_v2 -> spec_null -> gate_fdr -> master_table
all above -> paper_figures, database, deck
```

`gate_fdr` must run after `spec_null` and before `paper_figures`: it writes the background
pass-rate and the gate's empirical FDR into every table that quotes a pass-rate, so that no
figure or table reports a bare percentage.

## What lands where

| path | contents |
|---|---|
| `data/` | raw downloads, one directory per accession |
| `results/` | every table the paper quotes, as TSV |
| `results/figures/` | `P1`-`P6` are the manuscript figures; `F1`-`F12` are the exploratory set |
| `db/dkd.sqlite` | relational store, 14 tables and 3 audit views |
| `db/columnar/` | 22 Parquet files (zstd) holding the expression matrices |
| `logs/` | per-stage combined stdout and stderr |

## Verifying a specific claim

Each numbered claim in the manuscript is backed by a TSV under `results/`; the mapping is in
[docs/worklog/METHOD_AND_RESULTS.md](docs/worklog/METHOD_AND_RESULTS.md). The three a reviewer
is most likely to re-check:

```bash
python run_all.py --stage cohort_sensitivity   # results/cohort_sensitivity/  -> Fig P1, S3.2
python run_all.py --stage null                 # results/null_control.tsv     -> Fig P2, S3.1
python run_all.py --stage gate_fdr             # results/gate_characteristics.tsv -> S4.1
```

## Determinism

Every stage that resamples takes a fixed seed, so a rerun reproduces the reported numbers
exactly. The one exception is `generative` (`torch`), reported in the manuscript as a negative
result; no claim depends on it.
