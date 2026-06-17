# Monte Carlo Methods for Statistical Hypothesis Testing:
## Applications in Bioinformatics

This repository contains the code and data used in the course project:

"Monte Karlo metodai statistinių hipotezių tikrinimui. Taikymai bioinformatikoje"

## Overview

The project investigates statistically significant DNA sequence motifs in complete genomes using Monte Carlo hypothesis testing.

The workflow consists of:

1. k-mer enumeration
2. motif identification
3. first-order Markov null model construction
4. Monte Carlo simulation
5. empirical p-value estimation
6. biological motif localization

## Organisms

- Escherichia coli K12
- Escherichia coli Sakai
- Mycoplasmoides genitalium
- Streptomyces coelicolor
- Saccharomyces cerevisiae

## Main Results

Monte Carlo simulations revealed multiple motifs whose frequencies significantly exceeded expectations under first-order Markov null models.

## Reproducibility

Run:

python src/multi_motif_montecarlo.py

to generate null distributions.

Run:

python src/generate_five_genome_histograms.py

to reproduce Figure X.

Run:

python src/distribution_table.py

to reproduce Table X.
