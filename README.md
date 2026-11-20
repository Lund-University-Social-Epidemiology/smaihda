# smaihda

smaihda implements Simple-Means MAIHDA (S-MAIHDA), the simple-means form of the Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy, for binary outcomes. It reads the data as individual records or as a stratum table with the number of individuals and the number of cases in each stratum, in a two-level design (individuals in contexts) or a three-level design (individuals in sociodemographic-by-area strata in areas). All measures are computed from observed counts and proportions, without a random-effects model: Wilson 95% intervals for the proportions, the composition-standardised area proportions (EAR), the deviation of each stratum from its area, the variance partition coefficients (VPC, crude and standardised, uncorrected and corrected for binomial noise), the contextual clustering of inequalities (CCI) and its split into CCI, G and CMI, the contextual AUC (coAUC, coAUC* and coAUC*EAR), and the classification at a threshold. Intervals of the summary measures come from a within-context bootstrap that keeps the contexts fixed and resamples individuals.

The method and the measures are described in the S-MAIHDA preprint. The ESS preprint applies them to the European Social Survey and introduced the option for survey weights (see How to cite).

## Versions

Each published version is one commit with its tag. The files of a version are those published with it, unchanged.

| Tag | Date | Published with | Change |
|---|---|---|---|
| `v1.9.0` | 6 October 2026 | S-MAIHDA preprint, version 2: S1 File (Stata package), S2 File (Python), S3 File (R) | |
| `v1.10.0` | 16 October 2026 | ESS preprint, version 2: S1 File, folders SI-1 (Stata), SI-2 (Python) and SI-4 (R) | Survey weights as an option. Without weights every number, interval limits included, is the v1.9.0 number |

v1.10.0 was published without the Stata installation manifest (`smaihda.pkg`) and the test do-file of v1.9.0, so its `stata/` folder holds the command and its help file only.

## This version

Version 1.10.0.

| Folder | File | Language | Content |
|---|---|---|---|
| `stata/` | `smaihda.ado` | Stata 14 or later | The command `smaihda` |
| `stata/` | `smaihda.sthlp` | Stata | Help file |
| `python/` | `S2_File_smaihda_python.py` | Python 3 with numpy | The same measures |
| `r/` | `S3_File_smaihda_R.R` | R, base R only | The same measures |

The three programs give the same point estimates. The bootstrap limits agree within Monte Carlo error, since each language draws its own random numbers.

## Installation and use

Get a version with `git clone --branch <tag> https://github.com/Lund-University-Social-Epidemiology/smaihda`, for example `--branch v1.10.0`, or download the files of the tag.

**Stata.** For v1.10.0, copy `stata/smaihda.ado` and `stata/smaihda.sthlp` to the PERSONAL ado folder (shown by `sysdir` in Stata), or add the folder to the ado path with `adopath + "<path to the clone>/stata"`. Then `help smaihda`. For v1.9.0, `net install smaihda, from("<path to the clone>/stata") replace`.

**Python.** Requires numpy. `python python/S2_File_smaihda_python.py <stratum table or records.csv> <outcome> [<outcome> ...] [minn=10] [B=2000] [seed=2026] [weights=...]`. For the Malmö stratum table (S5 File of the S-MAIHDA preprint): `python python/S2_File_smaihda_python.py dataset-smaihda_stratum_table.csv psycmed private`.

**R.** Requires base R only. `Rscript r/S3_File_smaihda_R.R <stratum table or records.csv> <outcome> [<outcome> ...] [minn=10] [B=2000] [seed=2026] [weights=...]`.

The stratum table has one row per stratum and outcome, with the columns `area`, `strata_sociodem_id`, `strata_global_id`, `outcome`, `n` and `cases`. Individual records have the columns `area`, `strata_sociodem_id`, `strata_global_id` and one 0/1 column per outcome. With `weights=<prefix>` the stratum table also holds `<prefix>_sum`, `<prefix>_cases`, `<prefix>_sq_cases` and `<prefix>_sq_noncases` (sum of the weights, weighted cases, and the sums of the squared weights over the cases and over the non-cases). With individual records, `weights=<column>` names the weight variable.

## Authors

Juan Merlo, Nasir Z Bashir, Merida Rodriguez-Lopez, Kani Khalaf, Johan Öberg and Raquel Perez-Vicente, the authors of the S-MAIHDA preprint. Contact: Juan Merlo, Unit for Social Epidemiology, Faculty of Medicine, Lund University (juan.merlo@med.lu.se).

## How to cite

Merlo J, Bashir NZ, Rodriguez-Lopez M, Khalaf K, Öberg J, Perez-Vicente R. Describing health inequalities without distortion: Simple-Means MAIHDA vs Random-Effects MAIHDA. Version 2. medRxiv. 2026. doi:10.64898/2026.08.17.26360592.

The option for survey weights was published with: Merlo J, Öberg J, Rodriguez-Lopez M, Perez-Vicente R. Inequalities in self-rated health across European Union countries: a Simple-Means MAIHDA demonstration using the European Social Survey. Version 2. medRxiv. 2026. doi:10.64898/2026.10.20.26366006.

`CITATION.cff` gives both references in machine-readable form.

## Licence

smaihda is licensed under the Creative Commons Attribution 4.0 International licence (CC BY 4.0), see `LICENSE`, the licence under which its authors published it with the S-MAIHDA preprint and the ESS preprint. Reuse requires attribution to its authors and to the two preprints.
