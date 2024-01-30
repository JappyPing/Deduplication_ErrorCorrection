# How error correction affects PCR deduplication: a survey based on UMI datasets of short reads

This repository contains the scripts, the state-of-the-art methods and datasets (links) used to conduct the analysis for the study "How error correction affects PCR deduplication: a survey based on UMI datasets of short reads".

**Notes:** This repository contains other researchers' source codes and executable software of PCR-deduplication and error-correction on short reads. The copyright of these source codes and programs belongs to the original authors. If you use all or part of them, please remember to cite these works list in [Reference](#reference).

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
All the experiments were running on the Red Hat Enterprise Linux release 8.6 with the CPU model of Intel(R) Xeon(R) Gold 6238R CPU @ 2.20GHz.
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
Install some deduplication and error-correction tools via Anaconda and pip
```bash
pip install git+https://github.com/pinellolab/AmpUMI.git
```
```
conda install bioconda::fastp
```
```
conda install bioconda::minirmd
```
```bash
conda install bioconda::cd-hit
```
```
conda install bioconda::calib
```
```
conda install bioconda::fastuniq
```
```
conda install bioconda::bcool
```
Then, run the following scripts

**Notes**: please remember to replace the input directories to your owns before running the following scripts.

```bash
python ./src/deduplication.py
```

#### Error-Correction
```bash
python ./src/error_correction.py
```
for error correction by bcool, you need to run the separate script 
```bash
python ./src/bcool.py
```

#### Deduplication after Error-Correction
```bash
python ./src/deduplication_after_error_correction.py
```

### Reference
Kendell Clement, Rick Farouni, Daniel E. Bauer, and Luca Pinello. AmpUMI: Design and analysis of unique molecular identifiers for deep amplicon sequencing. Bioinformatics (Oxford, England), 34(13):i202–i210, 2018.

Tom Smith, Andreas Heger, and Ian Sudbery. UMI-tools: Modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy. Genome Research, 27(3):491–499, 2017.

Maria Tsagiopoulou, Maria Christina Maniou, Nikolaos Pechlivanis, Anastasis Togkousidis, Michaela Kotrová, Tobias Hutzenlaub, Ilias Kappas, Anastasia Chatzidimitriou, and Fotis Psomopoulos. UMIc: A Preprocessing Method for UMI Deduplication and Reads Correction. Frontiers in Genetics, 12:660366, 2021.

Antonio Sérgio Cruz Gaia, Pablo Henrique Caracciolo Gomes de Sá, Mônica Silva de Oliveira, and Adonney Allan de Oliveira Veras. NGSReadsTreatment – A Cuckoo Filter-based Tool for Removing Duplicate Reads in NGS Data. Scientific Reports, 9(1):11681, 2019.

Hang Dai and Yongtao Guan. Nubeam-dedup: A fast and RAMefficient tool to de-duplicate sequencing reads without mapping. Bioinformatics, 36(10):3254–3256, 2020.

Gianvito Urgese, Emanuele Parisi, Orazio Scicolone, Santa Di Cataldo, and Elisa Ficarra. BioSeqZip: A collapser of NGS redundant reads for the optimization of sequence analysis. Bioinformatics, 36(9):2705–2711, 2020.

Shifu Chen, Yanqing Zhou, Yaru Chen, and Jia Gu. Fastp: An ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17):i884–i890, 2018.

Jason A. Vander Heiden, Gur Yaari, Mohamed Uduman, Joel N.H. Stern, Kevin C. O’Connor, David A. Hafler, Francois Vigneault, and Steven H. Kleinstein. pRESTO: A toolkit for processing high-throughput sequencing raw reads of lymphocyte receptor repertoires. Bioinformatics, 30(13):1930–1932, 2014.

Ying Huang, Beifang Niu, Ying Gao, Limin Fu, and Weizhong Li. CD-HIT Suite: A web server for clustering and comparing biological sequences. Bioinformatics (Oxford, England), 26(5):680–682, 2010.

Jorge González-Domínguez and Bertil Schmidt. ParDRe: Faster parallel duplicated reads removal tool for sequencing studies. Bioinformatics (Oxford, England), 32(10):1562–1564, 2016.

Yuansheng Liu, Xiaocai Zhang, Quan Zou, and Xiangxiang Zeng. Minirmd: accurate and fast duplicate removal tool for short reads via multiple minimizers. Bioinformatics, 37(11):1604– 1606, 2021.

Baraa Orabi, Emre Erhan, Brian McConeghy, Stanislav V Volik, Stephane Le Bihan, Robert Bell, Colin C Collins, Cedric Chauve, and Faraz Hach. Alignment-free clustering of umi tagged dna molecules. Bioinformatics, 35(11):1829–1836, 2019.

Haibin Xu, Xiang Luo, Jun Qian, Xiaohui Pang, Jingyuan Song, Guangrui Qian, Jinhui Chen, and Shilin Chen. FastUniq: A fast de novo duplicates removal tool for paired short reads. PloS One, 7(12):e52249, 2012.

Heng Li. BFC: Correcting Illumina sequencing errors. Bioinformatics (Oxford, England), 31(17):2885–2887, 2015.

Felix Kallenborn, Andreas Hildebrandt, and Bertil Schmidt. CARE: Context-aware sequencing read error correction. Bioinformatics, 37(7):889–895, 2021.

Felix Kallenborn, Julian Cascitti, and Bertil Schmidt. CARE 2.0: Reducing false-positive sequencing error corrections using machine learning. BMC Bioinformatics, 23(1):227, 2022.

Leena Salmela and Jan Schröder. Correcting errors in short reads by multiple alignments. Bioinformatics, 27(11):1455–1461, 2011.

Marcel H. Schulz, David Weese, Manuel Holtgrewe, Viktoria Dimitrova, Sijia Niu, Knut Reinert, and Hugues Richard. Fiona: A parallel and automatic strategy for read error correction. Bioinformatics, 30(17):i356–i363, 2014.

Li Song, Liliana Florea, and Ben Langmead. Lighter: Fast and memory-efficient sequencing error correction without counting. Genome Biology, 15(11):509, 2014.

Eric Marinier, Daniel G. Brown, and Brendan J. McConkey. Pollux: Platform independent error correction of single and mixed genomes. BMC Bioinformatics, 16(1):10, 2015.

Lucian Ilie and Michael Molnar. RACER: Rapid and accurate correction of errors in reads. Bioinformatics (Oxford, England), 29(19):2490–2493, 2013.

Antoine Limasset, Jean-François Flot, and Pierre Peterlongo.Toward perfect reads: Self-correction of short reads via mapping on de Bruijn graphs. Bioinformatics, 36(2):651, 2020.

<!-- ## License -->

<!-- ## How to cite -->

## Contact
Feel free to contact me (ping.pengyao@gmail.com) if you have any questions, comments, suggestions etc regarding this study.

