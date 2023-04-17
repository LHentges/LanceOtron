<p align="center">
    <img src="assets/LanceOtron_logo_shadow_dark.png" alt="LanceOtron logo">
</p>

-----------------------

[![PyPI version](https://badge.fury.io/py/lanceotron.svg)](https://badge.fury.io/py/lanceotron) [![Downloads](https://pepy.tech/badge/lanceotron)](https://pepy.tech/project/lanceotron) [![CircleCI](https://circleci.com/gh/Chris1221/lanceotron/tree/main.svg?style=svg)](https://circleci.com/gh/Chris1221/lanceotron/tree/main) [![codecov](https://codecov.io/gh/Chris1221/lanceotron/branch/main/graph/badge.svg?token=yhL3YI00UP)](https://codecov.io/gh/Chris1221/lanceotron)

[![Conda version](http://anaconda.org/sgriva/lanceotron/badges/version.svg)](https://anaconda.org/sgriva/lanceotron) [![Anaconda-Server Badge](https://anaconda.org/sgriva/lanceotron/badges/downloads.svg)](https://anaconda.org/sgriva/lanceotron)

**LanceOtron** is a machine learning, genomic data extraction and analysis tool trained for ATAC-seq, ChIP-seq, and DNase-seq peak calling. A freely available and fully-featured webtool version, utilising the graphical user interface [MLV](https://mlv.molbiol.ox.ac.uk) and hosted at the [MRC WIMM Centre of Computational Biology, University of Oxford](https://www.imm.ox.ac.uk/research/units-and-centres/mrc-wimm-centre-for-computational-biology), can be found at [LanceOtron.molbiol.ox.ac.uk](https://lanceotron.molbiol.ox.ac.uk).

## Python Package and Requirements

LanceOtron was built using Python 3.8.3 and TensorFlow 2. The models have been saved such that a TensorFlow 1 setup could be used making only minor amendments to the scripts (see note in modules folder). Additional packages were used for benchmarking LanceOtron - see [requirements.txt](lanceotron/requirements.txt) for specific version numbers used. 

LanceOTron is currently tested on Linux and MacOSX, please see the CircleCI implementation for more information if interested. We do not currently support native Windows installations, nor do we have any immediate plans due to the amount of work involved to maintain it. Windows subsystem for Linux should function as expected, though if you have issues using LanceOTron on your specific installation, please [raise an issue](https://github.com/LHentges/LanceOtron/issues/new/choose) and we will help you resolve them. 

**Additional Python Packages for Benchmarking:**

N.B. [bedtools](https://github.com/arq5x/bedtools2) needs to be installed to use *pybedtools*.
> * pandas
> * matplotlib
> * pybedtools
> * seaborn

## Command Line Installation

There are three ways to install LanceOTron. The first and second methods are recommended for general users while the second is recommended for developers interested in extending the source code.

### Method 1: Installing from Pypi 

LanceOTron is [hosted on Pypi](https://pypi.org/project/lanceotron/) and can be easily installed using `pip`. 

```sh
pip install lanceotron
```

### Method 2: Conda installation

Conda installation is available via [sgriva's](https://anaconda.org/sgriva) channel.

```sh
conda install -c sgriva lanceotron
```

### Method 3: Local installation

We recommend using a fresh virtual environment with Python 3.7+. 

1. Clone the repository and navigate to the CLI package
2. Install dependencies with pip.
3. Install the package.
4. Run tests to ensure that everything is working.

```{sh}
git clone git@github.com:LHentges/LanceOtron.git # Step 1
cd LanceOTron/lanceotron

pip install -r requirements.txt # Step 2

pip install -e . # Step 3

python -m unittest # Step 4
```

## Usage

Currently there are 3 LanceOtron modules available. By default the modules return a bed file of **all candidate peaks, good, bad and otherwise**, along with their associated scores. For more details regarding candidate peak selection and the deep neural network scoring please see the citation below. 

Detailed usage instructions for each of the three modules along with a tutorial are available in the CLI directory [here](lanceotron/).

### Preparing bigwig files

All LanceOtron modules require a bigwig file to supply the model with coverage data. We recommend directly converting BAM files to bigwigs with [deepTools](https://github.com/deeptools/deepTools/tree/develop) using the following command:

> `bamCoverage --bam filename.bam.sorted -o filename.bw --extendReads -bs 1 --normalizeUsing RPKM`

The options used in this command are important, as they affect the shape of peaks and therefore the neural network's assessment. Extending the reads out to the fragment length represents a more accurate picture of coverage (N.B for paired end sequencing the extension length is automatically determined, single end tracks will require the user to specify the `--extendReads` length), as does using a bin size of 1 (the `--bs` flag). We recommend RPKM normalisation, as this was also used for the training data.


## Citation

Please see our [Bioinformatics article](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac525/6648462) for further details on this project or if you use it in your own work.

```bibtex
@article{10.1093/bioinformatics/btac525,
    author = {Hentges, Lance D and Sergeant, Martin J and Cole, Christopher B and Downes, Damien J and Hughes, Jim R and Taylor, Stephen},
    title = "{LanceOtron: a deep learning peak caller for genome sequencing experiments}",
    journal = {Bioinformatics},
    year = {2022},
    month = {07},
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btac525},
    url = {https://doi.org/10.1093/bioinformatics/btac525},
    note = {btac525},
    eprint = {https://academic.oup.com/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btac525/45048211/btac525.pdf},
}
```
