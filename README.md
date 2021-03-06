#  Analysis for: [Native Defects and their Doping Response in the Lithium Solid Electrolyte Li<sub>7</sub>La<sub>3</sub>Zr<sub>2</sub>O<sub>12</sub>](https://pubs.acs.org/doi/abs/10.1021/acs.chemmater.9b04319)

[Alexander G. Squires](https://orcid.org/0000-0001-6967-3690)  
[David O. Scanlon](https://orcid.org/0000-0001-9174-8601)  
[Benjamin J. Morgan.](http://orcid.org/0000-0002-3056-8233)

This repository contains supporting code for the paper [Native Defects and their Doping Response in the Lithium Solid Electrolyte Li<sub>7</sub>La<sub>3</sub>Zr<sub>2</sub>O<sub>12</sub>](https://pubs.acs.org/doi/abs/10.1021/acs.chemmater.9b04319) [1].

The repository contains:
1. a [Jupyter notebook](https://github.com/alexsquires/native-defects-in-llzo/blob/master/analysis/defects_in_llzo.ipynb) containing code used to generate Figures 1,2,4,5,6 in the manuscript, and the accompnaying analysis.
2. calculation metadata needed to run the analysis.

## Note

To run the complete analysis, you will need a compiled version of the Fortran code: [sc-fermi](https://github.com/jbuckeridge/sc-fermi), and obtain the relevant `DOSCAR` files available as part of the accompanying [dataset](https://researchportal.bath.ac.uk/en/datasets/dataset-for-native-defects-and-their-doping-response-in-the-lithi).

## Contents

`data/`: This folder contains a series of .yaml files, containing data extracted from VASP calculations. The inputs and outputs for the source VASP calculations, along with instructions for extracting the relevant data and generating these files, are available at the [University of Bath Data Archive](https://researchportal.bath.ac.uk/en/datasets/dataset-for-native-defects-and-their-doping-response-in-the-lithi).

`analysis/`: This folder contains Jupyter notebooks that perform the analysis of the DFT data, using the input data in the `metadata` folder. These notebooks also generate relevant figures for publication.

`figures/`: This folder contains figures for publication, produced by the analysis scripts.

`doscars/`: Should you want to replicate the full analysis, this directory should contain the relevant oxygen vacancy doscars obtained from the accompanying [dataset](https://researchportal.bath.ac.uk/en/datasets/dataset-for-native-defects-and-their-doping-response-in-the-lithi)

`cplap_out/`: contains the output from thermodynamic stability analysis performed using the code [`CPLAP`](https://github.com/jbuckeridge/cplap)

## Dataset

The dataset in full is available at: https://researchportal.bath.ac.uk/en/datasets/dataset-for-native-defects-and-their-doping-response-in-the-lithi

## References

1. A. G. Squires, D. O. Scanlon, B. J. Morgan *Native Defects and their Doping Response in the Lithium Solid Electrolyte Li<sub>7</sub>La<sub>3</sub>Zr<sub>2</sub>O<sub>12</sub>*
