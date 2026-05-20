# Medicaid Expansion and Pharmaceutical Innovation

This repository contains the full research pipeline my MPhil thesis. The project investigates how the Affordable Care Act’s Medicaid expansion affected incentives to innovate in the pharmaceutical industry, using FDA approval records, patent filings, Medicaid drug utilisation data, and state-level preferred drug lists (PDLs).

## Data availability

Raw data files are not included in this repository due to size constraints.
Contact `mohammad.karim@economics.ox.ac.uk` for access.

## Research questions

1. **Aggregate innovation response** — Did the Medicaid expansion change the volume of FDA approvals and patent filings in therapeutic areas with greater Medicaid exposure?
2. **Medicaid share vs. PDL inclusion** — How do innovation responses differ when separating Medicaid market exposure and a firm's exposure to the state's PDL exposure.
3. **Heterogeneity** — How do these effects vary across drug classes or firm types?

## Repository structure

```
mphil-thesis/
├── raw_data/               # Source datasets (gitignored)
│   ├── Drugs@FDA/          #   FDA approval & product files
│   ├── patents/            #   USPTO patent grant data
│   ├── state_drug_utilisation_data/  #   CMS SDUD extracts
│   ├── PDL/                #   State preferred drug lists
│   ├── clinical_trials/    #   ClinicalTrials.gov data
│   ├── orange_book/        #   FDA Orange Book
│   ├── RxNorm/             #   RxNorm concept mappings
│   └── UMLS/               #   UMLS Metathesaurus files
│
├── aux_data/               # Crosswalks & lookup tables
│   ├── ndc_atc_mapping.csv
│   ├── cpc_atc_mapped.csv
│   ├── ndc_firm_mapping.csv
│   └── medicaid_expansion_proportion_by_*.csv
│
├── scripts/                # Data-processing pipeline (R)
│   ├── fda_approval_data.r
│   ├── patent_data.r
│   ├── state_drug_utilisatation_data.r
│   ├── pdl_data.r
│   ├── clinical_trials_data.r
│   ├── medicaid_expansion_dates.r
│   ├── rollout_map.r
│   └── summary_statistics.r
│
├── regressions/            # Estimation scripts
│   ├── RQ1.r               #   Aggregate event-study (fixest)
│   ├── RQ2.r               #   Medicaid share vs. PDL channel
│   ├── RQ3.r               #   Heterogeneity analysis
│   └── pdl_counterfactual_analysis.r
│
├── robustness/             # Robustness & sensitivity checks
│   ├── contdid.r
│   ├── binned_treatment.r
│   └── patent_matching.r
│
├── theoretical_model/      # Game-theoretic model (Python + LaTeX)
│
├── output/                 # Tables (.tex) and figures (.png)
│   ├── RQ1/ … RQ3/
│   ├── robustness/
│   ├── summary_statistics/
│   └── PDL_counterfactual/
│
├── exploration/            # Ad-hoc notebooks & exploratory work
│
├── drafting/               # Thesis chapters (LaTeX, modular)
│   ├── draft.tex           #   Master document
│   ├── introduction.tex
│   ├── literature_review.tex
│   ├── theoretical_model.tex
│   ├── methodology.tex
│   ├── data.tex
│   ├── empirical_results.tex
│   ├── robustness.tex
│   ├── conclusion.tex
│   └── references.bib
│
├── submission/             # Final compiled thesis
│
├── presentations/          # Quarto slide decks (.qmd → .html)
│
├── figures/                # Standalone visualisations
│
├── literature/             # Reading tracker (gitignored)
│
├── references/             # Term-by-term .bib files
│
└── tex_resources/          # Oxford thesis class & assets
```

## Reproducing the analysis

### Prerequisites

- **R ≥ 4.3** — all R package dependencies are managed by [`renv`](https://rstudio.github.io/renv/) and pinned in `renv.lock`
- **Python 3** (for the theoretical model simulations)
- **LaTeX** distribution with `latexmk` (for compiling the thesis)

### Steps

1. Clone the repo and open the project in R or VS Code.
2. Restore the R environment:
   ```r
   # install.packages("renv")
   renv::restore()
   ```
3. Place source data in `raw_data/` following the sub-folder layout above.
4. Run the data-processing scripts in `scripts/` — each script is self-contained and writes to `processed_data/` (gitignored).
5. Run the regression scripts in `regressions/` to produce tables and figures in `output/`.
6. Compile the thesis: `cd drafting && latexmk -pdf draft.tex`.

## Licence

MIT — see [LICENSE](LICENSE).
