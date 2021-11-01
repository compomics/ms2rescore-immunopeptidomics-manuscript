# ms2rescore-immunopeptidomics-manuscript

This github Repository contains code for the ms2rescore immunopeptidomics manuscript (DOI)\
All data used can be found and downloaded from Zenodo:

## Abstract

"Immunopeptidomics aims to identify immunopeptides, which are presented on Major Histocompatibility Complexes (MHC) on every cell and can be used to develop vaccines against pathogens and cancer. However, existing immunopeptidomics data analysis pipelines have some major hurdles to overcome, mostly due to the non-tryptic nature of immunopeptides, which complicates their identification. Previously, the machine and deep learning tools MS²PIP and DeepLC have shown to improve tryptic peptide identifications by accurately predicting tandem mass spectrometry (MS2)  peak intensities and retention times, respectively, and by using these predictions to rescore peptide-spectrum matches (PSMs) with post-processing tool Percolator. However, MS²PIP was still tailored towards tryptic peptides and fragmentation patterns of immunopeptides are drastically different. To enable MS²PIP-based rescoring of immunopeptide PSMs, we have retrained MS²PIP to include non-tryptic peptides  . These newly trained MS²PIP models greatly improve the predictions for immunopeptides and, surprisingly, also for tryptic peptides. Next, the new MS²PIP models, DeepLC, and Percolator were integrated into one software package, called MS²Rescore. Using MS²Rescore, 46% more spectra and 36% more unique peptides were identified at 1% false discovery rate (FDR), with even more extreme differences at 0.1% FDR, in comparison with standard Percolator rescoring. Due to the innovative extraction of MS²PIP-, DeepLC and search engine-based features, MS²Rescore even outperforms current state-of-the-art immunopeptide rescoring efforts. Thus, the integration of the new immunopeptide MS²PIP models, DeepLC, and Percolator into MS²Rescore shows great promise to substantially improve the identification of novel epitopes from immunopeptidomics workflows."

## images
ALl manuscript images can be found in data/Figures.\
![MS²PIP model performances](https://github.com/compomics/ms2rescore-immunopeptidomics-manuscript/blob/main/notebooks/data/Figures/Figure1A.svg)
![MS²Rescore feature contribution](https://github.com/compomics/ms2rescore-immunopeptidomics-manuscript/blob/main/notebooks/data/Figures/Figure3.svg)

## Notebooks
### Model training notebooks
These notebooks contain code used to retrain new non-tryptic models for MS²PIP
* Immunopeptide_chymotrypsin_ms2pipB
* Immunopeptide_chymotrypsin_ms2pipY
* Immunopeptide_ms2pipB
* Immunopeptide_ms2pipY
* No_tryptic_ms2pip

### Analyses notebooks
These notebooks contain code used to test the new models and the implementation into MS²Rescore
* Data_analysis
* MS2PIP_model_comparison
* MS2rescore_output_analysis

## Pipeline
Nextflow pipeline to create spectral libraries required to train MS²PIP.\
This pipeline requires following parameters:
1. MassIVE identifier | Pride identifier | local folder containing raw files
1. Search engine identification file
1. Search engine used :
    * Maxquant
    * Comet
    * Peaks
    * Spectrum mill

Optionally:
* Path to store output can be given
* Path to store downloaded RAW files (if MassIVE|PRIDE identifier is used)
* path to store converter MGF files

## Tools
This folder contains usefull tools for filehandeling, id-file parsing, downloading RAW files and converting xgboost models to C code

