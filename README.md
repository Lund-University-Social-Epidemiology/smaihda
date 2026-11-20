# smaihda

smaihda implements Simple-Means MAIHDA (S-MAIHDA), the simple-means form of the Multilevel Analysis of Individual Heterogeneity and Discriminatory Accuracy, for binary outcomes. It reads the data as individual records or as a stratum table with the number of individuals and the number of cases in each stratum, in a two-level design (individuals in contexts) or a three-level design (individuals in sociodemographic-by-area strata in areas). All measures are computed from observed counts and proportions, without a random-effects model: Wilson 95% intervals for the proportions, the composition-standardised area proportions (EAR), the deviation of each stratum from its area, the variance partition coefficients (VPC, crude and standardised, uncorrected and corrected for binomial noise), the contextual clustering of inequalities (CCI) and its split into CCI, G and CMI, the contextual AUC (coAUC, coAUC* and coAUC*EAR), and the classification at a threshold. Intervals of the summary measures come from a within-context bootstrap that keeps the contexts fixed and resamples individuals.

The method and the measures are described in the S-MAIHDA preprint (see How to cite).

## This version

Version 1.9.0 (6 October 2026), published as S1 File, S2 File and S3 File of the S-MAIHDA preprint, version 2. The files are unchanged.

| Folder | File | Language | Content |
|---|---|---|---|
| `stata/` | `smaihda.ado` | Stata 14 or later | The command `smaihda` |
| `stata/` | `smaihda.sthlp` | Stata | Help file |
| `stata/` | `smaihda.pkg` | Stata | Installation manifest |
| `stata/` | `smaihda_v1.9.0_test.do` | Stata | Test on the Malmö stratum table |
| `python/` | `S2_File_smaihda_python.py` | Python 3 with numpy | The same measures, the listing of S5 Text |
| `r/` | `S3_File_smaihda_R.R` | R, base R only | The same measures, the listing of S6 Text |

The three programs give the same point estimates. The bootstrap limits agree within Monte Carlo error, since each language draws its own random numbers.

## Installation and use

Get this version with `git clone --branch v1.9.0 https://github.com/Lund-University-Social-Epidemiology/smaihda`, or download the files of the tag `v1.9.0`.

**Stata.** In Stata, `net install smaihda, from("<path to the clone>/stata") replace`, then `help smaihda`. The test do-file expects the Malmö stratum table (S5 File of the S-MAIHDA preprint) in a folder `S5_File_stratum_table` beside the folder it runs from.

**Python.** Requires numpy. `python python/S2_File_smaihda_python.py <stratum table or records.csv> <outcome> [<outcome> ...] [minn=10] [B=2000] [seed=2026]`. For the Malmö stratum table: `python python/S2_File_smaihda_python.py dataset-smaihda_stratum_table.csv psycmed private`.

**R.** Requires base R only. `Rscript r/S3_File_smaihda_R.R <stratum table or records.csv> <outcome> [<outcome> ...] [minn=10] [B=2000] [seed=2026]`.

The stratum table has one row per stratum and outcome, with the columns `area`, `strata_sociodem_id`, `strata_global_id`, `outcome`, `n` and `cases`. Individual records have the columns `area`, `strata_sociodem_id`, `strata_global_id` and one 0/1 column per outcome.

## Authors

Juan Merlo, Nasir Z Bashir, Merida Rodriguez-Lopez, Kani Khalaf, Johan Öberg and Raquel Perez-Vicente, the authors of the S-MAIHDA preprint. Contact: Juan Merlo, Unit for Social Epidemiology, Faculty of Medicine, Lund University (juan.merlo@med.lu.se).

## How to cite

Merlo J, Bashir NZ, Rodriguez-Lopez M, Khalaf K, Öberg J, Perez-Vicente R. Describing health inequalities without distortion: Simple-Means MAIHDA vs Random-Effects MAIHDA. Version 2. medRxiv. 2026. doi:10.64898/2026.08.17.26360592.

`CITATION.cff` gives the same reference in machine-readable form.

## Licence

smaihda is licensed under the Creative Commons Attribution 4.0 International licence (CC BY 4.0), see `LICENSE`, the licence under which its authors published it with the S-MAIHDA preprint. Reuse requires attribution to its authors and to the S-MAIHDA preprint.
