# How error correction affects PCR deduplication: a survey on UMI datasets of short reads

This repository contains the scripts, the state-of-the-art methods and datasets (links) used to conduct the analysis for the study "How error correction affects PCR deduplication: a survey on UMI datasets of short reads"

## Table of Contents

- [Reproducing results](#reproducing-results)
    - [Datasets](#datasets)
    - [Comparative Analysis](#comparative-analysis)
        - [Deduplication](#deduplication)
        - [Error-Correction](#error-correction)
        - [Deduplication after Error-Correction](#deduplication-after-error-correction)
    - [Reference](#reference)
<!-- - License -->
<!-- - How to cite -->
- [Contact](#contact)

## Reproducing results

### Datasets
We take 8 UMI-based single-end sequencing datasets of the accession IDs [SRR1543964 - SRR1543971](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP045430&o=acc_s%3Aa) together with a pair-end sequencing dataset of the accession ID [SRR11207257](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR11207257&display=metadata) for the comparative analysis.
#### Data pre-processing
For processing datasets SRR1543964 - SRR1543971,
```bash
$ python ./src/preprocessing_se.py
```
For processing pair-end datasets SRR11207257,
```bash
$ python ./src/preprocessing_pe.py
```

### Comparative Analysis
Clone this repository first
```bash
git clone https://github.com/Jappy0/Deduplication_ErrorCorrection
cd Deduplication_ErrorCorrection
```

#### Deduplication
Install some deduplication tools via Anaconda
```bash
conda 
```
Then, run
```bash
python ./src/deduplication.py
```

#### Error-Correction
```bash
python ./src/error_correction.py
```
#### Deduplication after Error-Correction
```bash
python ./src/deduplication_after_error_correction.py
```

### Reference


<!-- ## License -->

<!-- ## How to cite -->

## Contact
Feel free to contact me (ping.pengyao@gmail.com) if you have any questions, comments, suggestions etc regarding this study.

