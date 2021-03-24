#!/bin/bash

set -e

echo "Plot density distribution ecoli"
Rscript ./figure_plotting/density.R

echo "Plot density distribution for edge cases"
Rscript ./figure_plotting/density_edge_cases.R

echo "Plot detected proteins over k"
Rscript ./figure_plotting/detected_protein_k_psm.R

echo "Plot original length retrieval"
Rscript ./figure_plotting/length_retrieval.R

echo "Plot s-hat over e-value"
Rscript ./figure_plotting/scatterplot_s-hat.R

echo "Barplot species composition"
Rscript ./figure_plotting/species_composition.R

echo "Venn diagramms"
Rscript ./figure_plotting/venn.R

echo "Plot spectra"
python ./figure_plotting/plot_spectra.py