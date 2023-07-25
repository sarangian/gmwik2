# GMWI2: Gut Microbiome Wellness Index 2
![poop on a chip](./poop.png)

[![Anaconda-Server Badge](https://anaconda.org/bioconda/gmwi2/badges/version.svg)](https://anaconda.org/bioconda/gmwi2)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gmwi2/badges/platforms.svg)](https://anaconda.org/bioconda/gmwi2)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gmwi2/badges/license.svg)](https://anaconda.org/bioconda/gmwi2)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gmwi2/badges/downloads.svg)](https://anaconda.org/bioconda/gmwi2)

### Description

GMWI2 (Gut Microbiome Wellness Index 2) is a robust and biologically interpretable predictor of health status based on gut microbiome taxonomic profiles.

On a stool metagenome sample, this command-line tool performs four major steps:
1. Quality control
   1. Removal of overrepresented sequences (probable adapter sequences) using [fastqc](https://github.com/s-andrews/FastQC)
   2. Removal of human DNA contaminants (reads that map to GRCh38/hg38) using [Bowtie 2](https://github.com/BenLangmead/bowtie2)
   3. Removal of adapter sequences and low quality reads using [Trimmomatic](https://github.com/timflutre/trimmomatic)
2. Taxonomic profiling using MetaPhlAn3 (v3.0.13) with the mpa_v30_CHOCOPhlAn_201901 marker database
3. Transformation of taxonomic relative abundances into a binary presence/absence profile
4. Computation of the GMWI2 score using a Lasso-penalized logistic regression model trained on a meta-dataset of 8,069 health status labeled stool shotgun metagenomes

If you use GMWI2, please cite:

TODO

### Installation

To avoid dependency conflicts, please create an isolated conda environment and 
install the GMWI2 package. Installation via conda/mamba automatically installs GMWI2 and 
its dependencies (MetaPhlAn3).

1. Create new conda environment and install mamba
```bash
conda create --name gmwi2_env -c conda-forge mamba python=3.8
```

2. Activate environment
```bash
conda activate gmwi2_env
```

3. Install GMWI2 package with mamba
```bash
mamba install -c danielchang2002 -c bioconda -c conda-forge gmwi2=1.4
```

### Usage

Try downloading and running GMWI2 on an [example stool metagenome](./example).


```bash
Input: Two (forward/reverse) raw fastq (or fastq.gz) files generated from paired-end stool metagenome reads
Output: The GMWI2 (Gut Microbiome Wellness Index 2) score

usage: gmwi2 [-h] -n NUM_THREADS -i INPUT -o OUTPUT [-v]

* Example usage:

$ ls
.
├── forward.fastq
└── reverse.fastq

$ gmwi2 -i forward.fastq,reverse.fastq -n 8 -o output_prefix

$ ls
.
├── forward.fastq
├── reverse.fastq
├── output_prefix_GMWI2.txt
├── output_prefix_GMWI2_taxa.txt
└── output_prefix_metaphlan.txt

The three output files are: 
(i) output_prefix_GMWI2.txt: GMWI2 score
(ii) output_prefix_GMWI2_taxa.txt: A list of the taxa present in the sample used to compute GMWI2
(iii) output_prefix_metaphlan.txt: Raw MetaPhlAn3 taxonomic profiling output

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

required named arguments:
  -n NUM_THREADS, --num_threads NUM_THREADS
                        number of threads
  -i INPUT, --input INPUT
                        metagenome (.fastq) file(s)
  -o OUTPUT, --output OUTPUT
                        prefix to designate output file names
```

### Reproducing manuscript results

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/danielchang2002/GMWI2/blob/main/manuscript/GMWI2_manuscript.ipynb)

Please use the colab notebook linked above to reproduce all downstream analyses on the pooled dataset. 
See the [manuscript directory](./manuscript) for more details.

### Poop on a chip??

The top image was generated via [OpenAI DALL·E 2](https://openai.com/dall-e-2) using the prompt: "3D render of GPU chip in the form of a poop emoji, digital art".
The image was then widened using the [Runway Infinite Image tool](https://runwayml.com/ai-magic-tools/infinite-image/).
