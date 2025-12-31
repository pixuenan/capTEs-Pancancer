## Oncogenic transcription drives tumor-intrinsic dsRNA formation through TE activation and TE-anchored paired splicing

## Usage
### Pipeline code
To use the following code, python package `snakemake` need to be installed.
+ Code under `snakemake_fq` could generate gff files by StringTie and flair from fastq files.
+ Code under `snakemake_trans` could assemble transcripts and calculate TE coverage.

### TaPs Identification
The file `TaPs_Identification/TaPs_Identification.py` is used to identify and determine all TaPs in the samples.

### XGBoost model
The training and evaluation of the XGBoost model for the prediciton of the reponse to PD-1 inhibition in STAD patients was conducted by `model.py`.
