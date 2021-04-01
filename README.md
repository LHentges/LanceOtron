<p align="center">
    <img src="LanceOtron_logo_shadow_dark.png" alt="LanceOtron logo">
</p>

-----------------------

**LanceOtron** is a machine learning, genomic data extraction and analysis tool trained for ATAC-seq, ChIP-seq, and DNase-seq peak calling. A freely available and fully-featured webtool version, utilising the graphical user interface [MLV](https://mlv.molbiol.ox.ac.uk) and hosted at the [MRC WIMM Centre of Computational Biology, University of Oxford](https://www.imm.ox.ac.uk/research/units-and-centres/mrc-wimm-centre-for-computational-biology), can be found at [LanceOtron.molbiol.ox.ac.uk](https://lanceotron.molbiol.ox.ac.uk).

## Python Packages

LanceOtron was built using Python 3.8.3 and TensorFlow 2. The models have been saved such that a TensorFlow 1 setup could be used making only minor amendments to the scripts (see note in modules folder). Additional packages were used for benchmarking LanceOtron - N.B. [bedtools](https://github.com/arq5x/bedtools2) needs to be installed to use *pybedtools*. See [requirements.txt](requirements.txt) for specific version numbers used. 

**Required Python Packages for LanceOtron:**
> * scipy
> * numpy
> * pyBigWig
> * scikit\_learn==0.23.1
> * tensorflow==2.X.X

**Additional Python Packages for Benchmarking:**
> * pandas
> * matplotlib
> * pybedtools
> * seaborn

## Command Line Installation

We recommend using a fresh virtual environment with Python 3.7+ (older versions of Python 3 may work, but are untested), and installing the required packages with pip.

1. Clone/download repository
> `cd directory_you_want_LanceOtron_installed/`

> `git clone https://github.com/LHentges/LanceOtron.git`

2. pip install requirements
> `pip install --user -r directory_you_want_LanceOtron_installed/LanceOtron/requirements.txt`

3. Run commands using your computer's Python interpreter (see next section for common use cases)
> `python directory_you_want_LanceOtron_installed/LanceOtron/modules/find_and_score_peaks.py path/to/bigwig/my_experiment.bw`

## Usage

Currently there are 3 LanceOtron modules available. By default the modules return a bed file of **all candidate peaks, good, bad and otherwise**, along with their associated scores. For more details regarding candidate peak selection and the deep neural network scoring please see the citation below. 

### Preparing bigwig files

All LanceOtron modules require a bigwig file to supply the model with coverage data. We recommend directly converting BAM files to bigwigs with [deepTools](https://github.com/deeptools/deepTools/tree/develop) using the following command:

> `bamCoverage --bam filename.bam.sorted -o filename.bw --extendReads -bs 1 --normalizeUsing RPKM`

The options used in this command are important, as they affect the shape of peaks and therefore the neural network's assessment. Extending the reads out to the fragment length represents a more accurate picture of coverage (N.B for paired end sequencing the extension length is automatically determined, single end tracks will require the user to specify the `--extendReads` length), as does using a bin size of 1 (the `--bs` flag). We recommend RPKM normalisation, as this was also used for the training data.

### Modules

Module | Operation | Files Used
------ | --------- | ----------
Find and Score Peaks | Find enriched regions from coverage track, score regions with neural network | bigwig file
Find and Score Peaks with Input | Find enriched regions from coverage track, score regions with neural network, calculate pvalues for enrichment over control | experimental bigwig file, input (control) bigwig file
Score Peaks | Score user-supplied regions with neural network | bed file, bigwig file

#### Find and Score Peaks

This module first finds candidate peaks using an algorithm taking 2 parameters: 1) threshold and 2) window; default parameters are recommended. Signal is extracted from the bigwig file for each candidate peak, then passed to LanceOtron's deep neural network.

##### Basic Command

>  `python find_and_score_peaks.py my_experiment.bw`

##### Options

Flag | Description | Default
---- | ----------- | -------
-h, --help | Display arguments | *none*
-t, --threshold | Initial threshold used for candidate peak selection algorithm | 4
-w, --window | Window size for rolling mean used for candidate peak selection algorithm | 400
-f, --folder | Folder to write results to | current directory
--skipheader | Skip writing out header for results | *none*

#### Find and Score Peaks with Input

This module build on the **Find and Score Peaks** module, but additionally calculates the enrichment p-value of the experimental track above the control using the Poisson distribution.

##### Basic Command

>  `python find_and_score_peaks_with_input.py my_experiment.bw -i my_input_control.bw`

##### Options

Flag | Description | Default
---- | ----------- | -------
-i, --input | bigwig control input track | *none*
-h, --help | Display arguments | *none*
-t, --threshold | Initial threshold used for candidate peak selection algorithm | 4
-w, --window | Window size for rolling mean used for candidate peak selection algorithm | 400
-f, --folder | Folder to write results to | current directory
--skipheader | Skip writing out header for results | *none*

#### Score Peaks

Rather than finding candidate peaks based on enrichment from the bigwig coverage track, the user provides a bed file of regions to be scored by LanceOtron's deep neural network. 

##### Basic Command

>  `python score_peaks.py my_experiment.bw -b my_regions.bed`

##### Options

Flag | Description | Default
---- | ----------- | -------
-b, --bed |  bed file of regions to be scored using LanceOtron's neural network | *none* 
-h, --help | Display arguments | *none*
-f, --folder | Folder to write results to | current directory
--skipheader | Skip writing out header for results | *none*


## Citation

Please see our [bioRxiv article](https://www.biorxiv.org/content/10.1101/2021.01.25.428108v2) for further details on this project.

Hentges, L. D., Sergeant, M. J., Downes, D. J., Hughes, J. R., & Taylor, S. (2021). LanceOtron : a deep learning peak caller for ATAC-seq , ChIP-seq , and DNase-seq. *bioRxiv*.
